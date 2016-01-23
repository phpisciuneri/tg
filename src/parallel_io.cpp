#include "parallel_io.hpp"

#include <iostream>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/multi_array.hpp>
#include <mytools/mpi/opaque_objects.hpp>
#include <boost/iterator/iterator_traits.hpp>

#include "diagnostic/tracer.hpp"

namespace iplmcfd {

  // Use lower precision for output to save space
  typedef float oreal;

  //  +--------------------+
  //  | FILL OUTPUT BUFFER |
  //  +--------------------+
  // Fills up in an order consistent with the names inserted above, 
  // the post-processing data buffer corresponding to a preset 
  // home grid index (set at the caller site) 
  template< class T >
  inline T op_fill_buffer( lid_t I, const Field& fields, T pd )
  {

    for ( std::size_t idx=0; idx<fields.size(); idx++ )
      *pd++ = fields.data( I, idx );

    return pd; // for checks in the caller site
  }

  //  +-------------------------+
  //  |  IMPLEMENTATION STRUCT  |
  //  +-------------------------+
  /**
  \par Implementation notes:
  - Data is not output in place (i.e. not from simulation class members), instead copied 
  to a contiguous buffer of home grids. This has the disadvantage of requiring extra memory, 
  One big advantage is that data can be output asyncroniously (with non-blocking collective write), 
  allowing simulation to proceed as the i/o operation continues. Moreover, mpi data 
  structures to represent data in memory are much less complicated.
  - VTI appended data is comprised of a 4-byte field byte-count (referred to as magicbytes from hereforth) 
  and the field data for each field output. 
  - 4byte data is output by ioroot only as a blocking write
  - The ioroot is the rank with the 0-th home grid
  */
  struct Parallel_io::Impl
  {
    //////////////////////////////////////////////////////////////////////////
    // TYPES 

    /// multi-dimensional array type used for data buffer
    typedef boost::multi_array<oreal,2> MArray;

    /// scoped mpi datatype (for RAII call MPI_Type_free)
    typedef ScopedMpiDataType MpiType;

    //////////////////////////////////////////////////////////////////////////
    // DATA

    std::string header;
    bool ioroot;
    int rank;
    mytools::mpi::Comm comm;
    boost::filesystem::path output_dir;
    MPI_File mpifile;
    size_t num_home_grids;
    size_t num_all_grids;
    bool closed; ///< flag to track if MPI_File is closed 

    // i/o data buffers and file view types
    MArray postdata_buffer; ///< the post processing data is copied onto this buffer for asynchronous output
    std::vector< int >   magicbytes_buffer; ///< the bytecounts field 
    MpiType postdata_view; 
    MpiType magicbytes_view;
  };

  //  +--------------+
  //  |  DESTRUCTOR  |
  //  +--------------+
  // Close any open non-blocking collective I/O
  // dtor also cannot be inline or default because 
  // of scoped ptr (needs to see struct Impl defined)
  Parallel_io::~Parallel_io() { close(); }

  //  +---------------+
  //  |  CONSTRUCTOR  |
  //  +---------------+
  Parallel_io::Parallel_io( 
    const Domain& d, 
    const Subdomain& sd,
    const Field& fields ) : p(new Impl)
  {
    tracer::scope log("construct Parallel_io");

    //  +---------------------------+
    //  |  Initialize some members  |
    //  +---------------------------+
    p->ioroot = sd.homes().find(0) != sd.homes().end();
    MPI_Comm_dup(MPI_COMM_WORLD, p->comm.p());
    p->output_dir = simparam::instance().initial_path() / "out";
    p->num_home_grids = sd.nhomes();
    p->num_all_grids  = d.npoints();
    p->closed = true;
    MPI_Comm_rank(p->comm, &p->rank);

    // create output directory
    //   not an error if output_dir already exists
    boost::filesystem::create_directory(p->output_dir);

    //  +------------------+
    //  |  Compose Header  |
    //  +------------------+
    using std::string;
    string space_joined_op_fields;
    std::vector<string> op_names = fields.names();
    const int OPnum = op_names.size();
    BOOST_FOREACH( const std::string& name , op_names )
      space_joined_op_fields += name + " ";

    p->header = (
      boost::format(
      "<?xml version=\"1.0\"?>"
      /// \todo fix endianness
      "<VTKFile\n"
      "   type    = \"ImageData\"   \n"
      "   version = \"0.1\"         \n"
      "   byte_order = \"LittleEndian\"  \n"
      ">\n"
      "<ImageData\n"
      "   WholeExtent = \"0 %1% 0 %2% 0 %3%\"  \n"
      "   Origin      = \"%4% %5% %6%\"        \n"
      "   Spacing     = \"%7% %8% %9%\"        \n"
      ">\n"
      "<Piece Extent=\"0 %1% 0 %2% 0 %3%\">\n"
      "<PointData Scalars = \"%10%\" >\n"
      ) %
      (d.shape()[XDIM]-1) % (d.shape()[YDIM]-1) % (d.shape()[ZDIM]-1) %
      d.xmin()[XDIM]     %  d.xmin()[YDIM]     %  d.xmin()[ZDIM]     %
      d.dx()[XDIM]       %  d.dx()[YDIM]       %  d.dx()[ZDIM]       %
      space_joined_op_fields ).str();
    ;

    // write dataarray tags with offsets
    boost::format data_array_tag(
      "<DataArray NumberOfComponents=\"1\" format=\"appended\" "
      "   type=\"Float%1%\" Name=\"%2%\" offset=\"%3%\" />\n" );

    // compute offsets by accumulating
    const int vti_data_count_sz = 4; // 4 is the number VTK expects
    assert( vti_data_count_sz == sizeof(int)); // assumed in other locations in the implementation
    size_t offset = 0;
    for(int i=0; i<OPnum; i++ ) 
    {
      p->header += (data_array_tag%(sizeof(oreal)*8)% op_names[i] % offset).str();
      offset += sizeof(oreal)*p->num_all_grids + vti_data_count_sz;
    }

    // wrap up
    p->header += "</PointData></Piece></ImageData><AppendedData encoding=\"raw\">\n_";

    //  +------------------------------------+
    //  |  CREATE MESSAGE/FILEBUFFERS/TYPES  |
    //  +------------------------------------+
    // Allocate postprocess data buffer
    boost::array<int,2> shape = { OPnum, sd.nhomes() };
    p->postdata_buffer.resize(shape);

    std::fill_n(p->postdata_buffer.data(), p->postdata_buffer.num_elements(), oreal(0));

    // Create FileDataRef for the postprocessing data
    //   Copy home gids to a int vector
    std::vector<int> gid_vector( sd.nhomes() );
    std::vector<int>::iterator pgid = gid_vector.begin();
    BOOST_FOREACH( Subdomain::gidlid_mp::const_reference gl, sd.homes() )
      *pgid++ = gl.first;

    //    Create datatype for one field
    Impl::MpiType one_field_type;
    MPI_Type_create_indexed_block(
      sd.nhomes(),
      1,
      &gid_vector[0],
      mytools::mpi::datatype<oreal>(),
      one_field_type.ref()
      );
    one_field_type.commit();

    //    Create datatype for all fields
    std::vector<int> blen(OPnum,1);
    std::vector<MPI_Aint> displ(OPnum);
    MPI_Aint one_size = sizeof(oreal)*d.npoints() + vti_data_count_sz;
    displ[0] = vti_data_count_sz;
    for(int i=1; i<OPnum; ++i)
      displ[i] = displ[i-1] + one_size;
    MPI_Type_hindexed(
      OPnum,
      &blen[0],
      &displ[0],
      one_field_type,
      p->postdata_view.ref()
      );
    p->postdata_view.commit();

    // Create FileDataRef for magicbytes
    // there are OPnum many magicbytes, each with same value (same type, same count for each OP field)
    if (p->ioroot) {
      p->magicbytes_buffer.resize( OPnum, sizeof(oreal)*d.npoints() );
      MPI_Type_hvector(
        OPnum,
        1,
        one_size,
        MPI_INT,
        p->magicbytes_view.ref()
        );
      p->magicbytes_view.commit();

    }
  }

  //  +-------------------+
  //  |  OUTPUT_PARAVIEW  |
  //  +-------------------+
  void Parallel_io::output_paraview( 
    const Field& fields,
    int iter,
    int iter_stats )
  {
    tracer::scope _("Parallel_io::output_paraview");

    //  +--------+
    //  |  OPEN  |
    //  +--------+
    // Complete last output
    close();

    // Set filename, eg. for iter=123 we get "f.000123.vti"
    // prefix (f.) is just an arbitrary letter. If no such prefix, Paraview doesn't recognize
    // counted output. 
    std::string filename = 
      (p->output_dir / ((boost::format("f.%06d")%iter).str() + ".vti")).string();

    // Open file
    CHECK_MPI_CALL(MPI_File_open( 
      p->comm,
      (char*) filename.c_str(),
      MPI_MODE_WRONLY | MPI_MODE_CREATE, // open, create, overwrite if exists
      MPI_INFO_NULL,
      &p->mpifile
      ));

    //  +-------------------------------------+
    //  |  OUTPUT HEADER AND ROOT-ONLY STUFF  |
    //  +-------------------------------------+

    // Output header and magic bytes
    if ( p->ioroot )
    {
      CHECK_MPI_CALL(MPI_File_write(
        p->mpifile,
        (void*) p->header.c_str(),
        p->header.size(),
        MPI_CHAR,
        MPI_STATUS_IGNORE
        ));
    }

    // set_view is a collective operation, 
    // needs to be called by all even if only ioroot will do the output
    CHECK_MPI_CALL(MPI_File_set_view(
      p->mpifile,
      p->header.size(),
      MPI_INT,
      p->ioroot? p->magicbytes_view : MPI_INT, // just a dummy view type for non-root
      "native",
      MPI_INFO_NULL
      ));

    if (p->ioroot)
    {
      CHECK_MPI_CALL(MPI_File_write(
        p->mpifile,
        &p->magicbytes_buffer[0],
        p->magicbytes_buffer.size(),
        MPI_INT,
        MPI_STATUS_IGNORE
        ));
    };

    //  +------------------+
    //  |  UPDATE BUFFERS  |
    //  +------------------+
    for( size_t i=0; i< p->num_home_grids; i++ ) 
    {
      typedef Impl::MArray::index_range range;
      typedef Impl::MArray::array_view<1>::type slice_type;
      slice_type ith_slice = p->postdata_buffer[ boost::indices[range()][i] ];

      slice_type::iterator pd = op_fill_buffer( i, fields, ith_slice.begin() );
      BOOST_ASSERT( pd == ith_slice.end() ) ; // make sure all of read buffer is visited
    }

    //  +---------------+
    //  |  OUTPUT DATA  |
    //  +---------------+
    CHECK_MPI_CALL( MPI_File_set_view(
      p->mpifile,
      p->header.size(),
      mytools::mpi::datatype<oreal>(),
      p->postdata_view,
      "native",
      MPI_INFO_NULL ) );

    CHECK_MPI_CALL( MPI_File_write_all_begin(
      p->mpifile,
      p->postdata_buffer.data(),
      p->postdata_buffer.num_elements(),
      mytools::mpi::datatype<oreal>() ) );

    p->closed = false;

  }

  //  +---------+
  //  |  CLOSE  |
  //  +---------+
  void Parallel_io::close()
  {
    if ( p->closed ) return;
    tracer::scope log("Parallel_io::close()");

    CHECK_MPI_CALL(MPI_File_write_at_all_end(
      p->mpifile,
      p->postdata_buffer.data(),
      MPI_STATUS_IGNORE
      ));

    CHECK_MPI_CALL(MPI_File_close(
      &p->mpifile
      ));

    p->closed = true;  

  }

} // namespace iplmcfd