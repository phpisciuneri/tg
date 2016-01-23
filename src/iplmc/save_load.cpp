#include "../iplmc.hpp"

#include <sstream>

#include <mytools/mpi/datatype.hpp>
#include <mytools/mpi/opaque_objects.hpp>
#include <boost/assert.hpp>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <mpi.h>

#define CHECK CHECK_MPI_CALL



namespace iplmcfd {

// Auxiliary structs are defined in an unnamed namespace to 
// avoid potential linkage clash
namespace{ 

//  +----------------------+
//  |  FixedParticleStuff  |
//  +----------------------+

/// Contains information that is fixed for all particles
struct FixedParticleStuff
{
   size_t fixed_size;
   size_t phi_size;
   int    size_bytes;
   ScopedMpiDataType etype;

   // ---
public:
   static const FixedParticleStuff& get() {
      static FixedParticleStuff c;
      return c;
   }
private:
   FixedParticleStuff() {
      particle sample_particle;
      fixed_size = particle::fixed_num_reals();
      phi_size   = sample_particle.phi.size();
      MPI_Type_contiguous( fixed_size+phi_size, mpi_real, etype.ref());
      etype.commit();
      MPI_Type_size(etype, &size_bytes);
   }
};


//  +----------------+
//  |  Fixed Header  |
//  +----------------+
struct CheckPointHeader
{
   char tagb[8];
   int num_cells;   // global number of cells
   int num_particles; // global number of particles 
   int max_particles_per_cell; // maximum number of particles per cell, globally
   int particle_size_bytes;  // number of bytes per particle
   char tage[8];

   //---
   CheckPointHeader(int num_global_cells = 0)
   : num_cells(num_global_cells)
   {
      //           0123456
      strcpy(tagb,"[ipmchd");
      strcpy(tage,"ipmchd]");
      particle_size_bytes = FixedParticleStuff::get().size_bytes;
      num_cells = num_global_cells;
   }
   friend bool operator==(const CheckPointHeader& left, const CheckPointHeader& right) {
      return std::equal(left.tagb, left.tagb+sizeof(left.tagb), right.tagb)
      &&  std::equal(left.tage, left.tage+sizeof(left.tage), right.tage)
      &&  left.particle_size_bytes == right.particle_size_bytes
      &&  left.num_cells == right.num_cells
      ;
   }

   static MPI_Datatype mpitype() {
      static MPI_Datatype type = MPI_DATATYPE_NULL;
      if ( type == MPI_DATATYPE_NULL )
      {
         MPI_Type_contiguous( sizeof(CheckPointHeader), MPI_BYTE, &type );
         MPI_Type_commit(&type);
      }
      return type;
   }
};


//  +--------------+
//  |  HindexType  |
//  +--------------+

/// Facilitates creation of hindexed types (used in ParticleIoTypeHelper)
struct HindexType : public ScopedMpiDataType 
{
   HindexType( size_t num_reserve, MPI_Datatype etype )
   : etype(etype)
   {
      blen.reserve(num_reserve);
      displ.reserve(num_reserve);
   }

   virtual int commit()
   {
      MPI_Type_create_hindexed(
            blen.size(),
            &blen[0],
            &displ[0],
            etype,
            ref()
      );
      return ScopedMpiDataType::commit();
   }

   std::vector<int> blen;
   std::vector<MPI_Aint> displ;
   MPI_Datatype etype;
};


//  +----------------------------+
//  |  Particle I/O Type Helper  |
//  +----------------------------+

/// Constructs MPI types for file and memory from carlo 

struct ParticleIoTypes
{
   boost::shared_ptr<HindexType> filetype;
   boost::shared_ptr<HindexType> memtype;

   // ------ 
   ParticleIoTypes( Iplmc& c, int max_particles_per_cell )
   {
      // count total number of home particles
      size_t num_particles = 0;
      for (size_t i=0; i<c.subdomain().nhomes(); ++i)
         num_particles += c.cells()[i].contained().size();

      // home cells 
      size_t num_cells = c.subdomain().nhomes();

      // fixed particle informational stuff
      const FixedParticleStuff& particle_stuff = FixedParticleStuff::get();

      // Begin creating memory and filetypes
      // File type is just dump of cells' contained particles (phi first, other fields next) 
      // Cells are dumped in global index order.
      // Memory type corresponds to a linked list of particles. addresses of 
      // each particles phi field and fixed fields members must be tracked separately

      filetype.reset( new HindexType( num_cells,  particle_stuff.etype ));
      memtype.reset( new HindexType( num_particles*2, mpi_real )); // 1 for fixed, 1 for phi arrays

      // Construct hindexed types
      BOOST_FOREACH( Subdomain::map_gid2lid_t::const_reference gidlid, c.subdomain().homes() )
      {
         // alias cell
         cell& mycell = c.cells()[gidlid.second];
         // Add particle address and block count for this cell
         MPI_Aint first_particle_address = gidlid.first * max_particles_per_cell * particle_stuff.size_bytes;
         filetype->displ.push_back( first_particle_address );
         filetype->blen.push_back( static_cast<int>( mycell.contained().size() ));

         // add particles one-by-one to the memory type
         BOOST_FOREACH( particle& p, mycell.contained() )
         {
            MPI_Aint address;

            // Add phi 
            MPI_Get_address( &p.phi[0], &address );
            memtype->displ.push_back(address);
            memtype->blen.push_back( particle_stuff.phi_size );

            // Add fixed struct fields
            MPI_Get_address( const_cast<real*>(p.first_fixed_member()), &address );
            memtype->displ.push_back(address);
            memtype->blen.push_back( particle_stuff.fixed_size );
         }
      }

      filetype->commit();
      memtype->commit();
   }

};

} // end unnamed namespace

//  +--------+
//  |  LOAD  |
//  +--------+
void Iplmc::load(MPI_File mpifile)
{
   using std::vector;
   const FixedParticleStuff& particle_stuff = FixedParticleStuff::get();

   //  +---------------+
   //  |  READ HEADER  |
   //  +---------------+
   mpi_file_reset_view(mpifile, CheckPointHeader::mpitype());

   CheckPointHeader hdr;
   CHECK(MPI_File_read_all( 
      mpifile,
      &hdr, 1, CheckPointHeader::mpitype(),
      MPI_STATUS_IGNORE
      ));

   // sanity checks
   CheckPointHeader hdr_expected( m_domain.npoints() );
   BOOST_ASSERT( hdr == hdr_expected );


   //  +-----------------------------+
   //  |  READ CELL PARTICLE COUNTS  |
   //  +-----------------------------+
   vector<int> particle_counts( m_subdomain.nhomes() );

   MPI_Offset offset_cell_particle_counts_begin = tell(mpifile);
   MPI_Aint intsize; MPI_File_get_type_extent(mpifile, MPI_INT, &intsize);
   MPI_Offset offset_cell_particle_counts_end   = offset_cell_particle_counts_begin
      + intsize* m_domain.npoints();

   // get gids into a vector
   vector<int> gids;
   gids.reserve( m_subdomain.nhomes() );
   BOOST_FOREACH( Subdomain::map_gid2lid_t::const_reference gidlid, m_subdomain.homes() )
      gids.push_back( gidlid.first );

   // create filetype using gids
   ScopedMpiDataType filetype;      
   MPI_Type_create_indexed_block(
      gids.size(),
      1,
      &gids[0],
      MPI_INT,
      filetype.ref()
      );
   filetype.commit();

   // Read 
   char amode[]="native"; // const char* to char* warning workaround
   CHECK(MPI_File_set_view(
      mpifile,
      offset_cell_particle_counts_begin,
      MPI_INT,
      filetype,
      amode,
      MPI_INFO_NULL
      ));

   CHECK(MPI_File_read_all( 
      mpifile,
      &particle_counts[0],
      static_cast<int>(particle_counts.size()),
      MPI_INT,
      MPI_STATUS_IGNORE
      ));


   //  +----------------------+
   //  |  Allocate Particles  |
   //  +----------------------+
   for (size_t i=0; i<m_subdomain.nhomes(); ++i)
      m_cells[i].contained().resize( particle_counts[i] );

   //  +-----------------+
   //  |  READ Particles |
   //  +-----------------+
   ParticleIoTypes particle_io_types(*this, hdr.max_particles_per_cell);

   CHECK(MPI_File_set_view(
      mpifile,
      offset_cell_particle_counts_end,
      particle_io_types.filetype->etype,
      *particle_io_types.filetype,
      amode,
      MPI_INFO_NULL));

   CHECK(MPI_File_read_all(
      mpifile,
      MPI_BOTTOM,
      1,
      *particle_io_types.memtype,
      MPI_STATUS_IGNORE
      ));

   // reset to end of particle data: uncomment if more data needs to be read
   //    mpi_file_reset_view(mpifile, particle_io_types.filetype->etype, particle_data_begin);
   //    MPI_File_seek(mpifile, num_particles_global, MPI_SEEK_SET);
}



//  +--------+
//  |  SAVE  |
//  +--------+
void Iplmc::save(MPI_File mpifile) const
{
   using std::vector; 
   const FixedParticleStuff& particle_stuff = FixedParticleStuff::get();

   int commrank, commsize;
   boost::tie(commrank, commsize) = mpi_get_comm_rank_size();
   int rootrank = 0;
   bool ioroot = commrank == rootrank;


   //  +----------------------+
   //  |  Get particle count  |
   //  +----------------------+
   int max_particles_per_cell=0; // for computing cell offsets (world)
   int num_particles_home = 0; // for ddiagnostic header info (root only)
   int num_particles; // to be MPI_Reduced
   vector<int> particle_counts( m_subdomain.nhomes() ); // for output (world)
   for ( size_t i=0; i<m_subdomain.nhomes(); ++i)
   {
      int nparticles = static_cast<int>( m_cells[i].contained().size() );
      max_particles_per_cell= (std::max)( max_particles_per_cell, nparticles );
      num_particles_home += nparticles;
      particle_counts[i] = nparticles;
   }


   MPI_Allreduce(MPI_IN_PLACE, &max_particles_per_cell, 1, MPI_INT, MPI_MAX, m_comm);
   MPI_Reduce(&num_particles_home, &num_particles, 1, MPI_INT, MPI_SUM, rootrank, m_comm);

   //  +----------------+
   //  |  WRITE HEADER  |
   //  +----------------+

   mpi_file_reset_view(mpifile, CheckPointHeader::mpitype());

   if ( ioroot ) 
   {
      CheckPointHeader hdr;
      hdr.num_cells     = m_domain.npoints();
      hdr.num_particles = num_particles;
      hdr.max_particles_per_cell = max_particles_per_cell;
      // write header
      CHECK(MPI_File_write(
            mpifile,
            &hdr,
            1,
            CheckPointHeader::mpitype(),
            MPI_STATUS_IGNORE
      ));
   }
   else 
      MPI_File_seek(mpifile, 1, MPI_SEEK_CUR);

   MPI_Offset offset_header_end = tell(mpifile);


   //  +---------------------+
   //  |  Compose datatypes  |
   //  +---------------------+

   // * Write actual cell particles counts in a global cell array of ints
   // * Write particles for each cell in global order, assuming maxmium particle count in each cell
   //   That way the offset for a given cell is fixed by its global id. 

   // Get home gid vector for indexed types
   vector<int> homegids;
   homegids.reserve( m_subdomain.nhomes() );
   BOOST_FOREACH( Subdomain::map_gid2lid_t::const_reference gidlid, m_subdomain.homes() )
   homegids.push_back( gidlid.first );

   // Cell particle counts file view
   ScopedMpiDataType pcounts_file_view;
   MPI_Type_create_indexed_block(
         homegids.size(),
         1,
         &homegids[0],
         MPI_INT,
         pcounts_file_view.ref()
   );
   pcounts_file_view.commit();

   //  +------------------------------------+
   //  |  WRITE CELL PARTICLE COUNTS ARRAY  |
   //  +------------------------------------+
   char amode[]="native"; // const char* to char* warning workaround
   MPI_File_set_view( 
         mpifile,
         offset_header_end,
         MPI_INT,
         pcounts_file_view,
         amode,
         MPI_INFO_NULL
   );
   CHECK(MPI_File_write_all(
         mpifile,
         &particle_counts[0],
         static_cast<int>(particle_counts.size()),
         MPI_INT,
         MPI_STATUS_IGNORE
   ));


   //  +-------------------+
   //  |  Write Particles  |
   //  +-------------------+

   // Create types
   ParticleIoTypes particle_io_types(const_cast<Iplmc&>(*this), max_particles_per_cell);
   MPI_Aint intsize; MPI_File_get_type_extent(mpifile, MPI_INT, &intsize);
   MPI_Offset offset_particles_begin = offset_header_end + intsize*m_domain.npoints();
   MPI_Offset offset_particles_end   = offset_particles_begin + FixedParticleStuff::get().size_bytes*m_domain.npoints();

   // Write
   CHECK(MPI_File_set_view(
      mpifile,
      offset_particles_begin,
      particle_io_types.filetype->etype,
      *particle_io_types.filetype,
      amode,
      MPI_INFO_NULL));

   CHECK(MPI_File_write_all(
      mpifile,
      MPI_BOTTOM,
      1,
      *particle_io_types.memtype,
      MPI_STATUS_IGNORE
      ));

   
   // reset to end of particle data: uncomment if more data needs to be written 
   //    mpi_file_reset_view(mpifile, particle_io_types.filetype->etype, particle_data_begin);
   //    MPI_File_seek(mpifile, num_particles_global, MPI_SEEK_SET);

   // write graph for Zoltan
   write_graph();
}

void Iplmc::write_graph() const
{
   // get particle size in bytes
   int psize = FixedParticleStuff::get().size_bytes;

   // open file
   std::ofstream os( "graph.txt" );
   
   // write in this format
   //  vtx_id  part_id     vtx_wgt  vtx_size   num_nbors  nbor_id  edge_wgt  nbor_id    edge_wgt
   //   0        0      	 5.0  	    100         2          1       2.0        5         2.0

   // output header line
   os << "## " 
      << "vtx_id "
      << "part_id "
      << "vtx_wgt "
      << "vtx_size "
      << "num_nbors "
      << "nbor_id "
      << "edge_wgt "
      << "..." << std::endl;

   BOOST_FOREACH( Subdomain::map_gid2lid_t::const_reference gidlid, m_subdomain.homes() )
   {
      // vertex ID is global ID
      os << gidlid.first << " ";

      // part ID is rank of partition
      os << m_subdomain.home_rank() << " ";

      // vertex weight is cpuload
      os << m_cells[gidlid.second].ctr.cputime << " ";

      // vertex size is number of particles in cell * psize
      os << m_cells[gidlid.second].contained().size() * psize << " ";

      // num neighbors is size of map
      int num_neigh = m_graph.nneighs( gidlid.second );
      os << num_neigh << " ";

      BOOST_FOREACH( Graph::pcomm_t::const_reference gidcount, m_graph.neighs( gidlid.second ) )
      {
         // neighbor ID is global ID of neighbor
         os << gidcount.first << " ";

         // edge weight is number of particles leaving * psize
         //os << gidcount.second * psize << " ";
         
         // lets assume edge weight is constant for now
         os << 1. << " ";
      }

      os << std::endl;

   }

}


} // namespace iplmcfd
