#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/cstdlib.hpp>
#include <boost/progress.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <mytools/mpi/environment.hpp>
#include <mytools/filesystem/change_directory.hpp>

#include <zoltan_cpp.h>

#include "defs.hpp"
#include "domain.hpp"
#include "simparam.hpp"
#include "subdomain.hpp"
#include "taylor_green.hpp"
#include "diagnostic/tracer.hpp"


// forward decl
boost::filesystem::path create_rank_dirs( int rank );

// +------+
// | main |
// +------+
int main( int argc, char* argv[] )
{

  using namespace iplmcfd;

  // print help
  if ( argc==2 )
    if ( std::string(argv[1]) == "--help" )
    {
      simparam::instance( argc, argv );
      return boost::exit_success;
    };

  // init MPI
  int rank, size;
  mytools::mpi::environment::init(argc,argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  // init Zoltan
  float version;
  if ( Zoltan_Initialize( argc, argv, &version ) )
    throw std::runtime_error( "Zoltan failed to Initialize." );
  if ( rank == 0 )
    std::cout << "Zoltan, version: " << version << std::endl;

  try {

    // parse simulation parameters
    const simparam& p = simparam::instance(argc, argv);
      
    // start with fresh tracer file
    tracer::remove_previous();

    // create and enter rank based subdirectories
    boost::filesystem::path root_dir = create_rank_dirs( rank );

    // construct domain
    Domain d( p );

    // construct subdomain
    boost::scoped_ptr< Subdomain > s;
    BoxPattern box(1);
    if ( ( p.decomp_type == "metis" ) || ( p.decomp_type == "wmetis" ) )
      s.reset( new MetisSubdomain( p, d, box, MPI_COMM_WORLD ) );
    else if ( p.decomp_type == "block" )
      s.reset( new BlockSubdomain( p, d, box, MPI_COMM_WORLD ) );
    else if ( p.decomp_type == "serial" )
      s.reset( new SerialSubdomain( p, d, box, MPI_COMM_WORLD ) );
    else
    {
      std::string error_str;
      error_str =  "Subdomain: Invalid decomposition type.\n";
      error_str += "decomp_type = " + p.decomp_type + "\n";
      throw std::invalid_argument( error_str );
    }

    TaylorGreen tg( p, d, *s );

    /*
    // DEBUG BLOCK
    if ( p.redecomp_type == "paragon" )
      tg.paragon_refinement();
    else
      tg.repartition();
    

    std::cout << s->my_rank() << " " << s->nhomes() << std::endl;

    MPI_Barrier( MPI_COMM_WORLD );
    if ( rank == 0 )
      std::cout << std::endl << "SUCCESS" << std::endl;
    return 0;
    // END DEBUG
    */

    // initialize
    tg.initialize();
    tg.output( 0 );
    
    // null stream 
    // useful for rank != 0 function calls
    std::ofstream null_stream;
    std::ostream* os = ( rank == 0 ) ? &std::cout : &null_stream;

    // init progress timer
    boost::progress_display prog( unsigned( p.iter_end ), *os );
    boost::progress_timer tmr;

    // open steps.out
    boost::filesystem::path steps = root_dir / "steps.out";
    std::ofstream sout( steps.c_str() );
    
    
    // main loop
    for ( int iter=1; iter<=p.iter_end; iter++, ++prog ) 
    {

      tg.step( iter );

      bool changes;
      if ( iter % p.repart_freq == 0 && iter != p.iter_end )
	if ( p.redecomp_type == "paragon" )  
	  changes = tg.paragon_refinement();
	else
	  changes = tg.repartition();
      
      if ( iter % p.iter_fdump == 0 )
      {
        if ( changes ) tg.reset_output();
        tg.output( iter );
      }

      if ( iter % p.iter_output_steps == 0 )
        tg.output_steps( sout, iter );

      // update time for next iteration
      tg.update_time();

    } // time loop
    
    if ( rank == 0 ) std::cout << "Success!" << std::endl;
    return boost::exit_success;

  }

  catch ( const std::exception& e )
  {
    std::cerr << "error: " << e.what() << std::endl;
  }
  
  return boost::exit_failure;

}

// +-------------------------+
// | CREATE RANK DIRECTORIES |
// +-------------------------+
boost::filesystem::path create_rank_dirs( int rank )
{
  iplmcfd::tracer::scope _( "create_rank_dirs" );

  // store root directory of the simulation, 
  // it will be changed below to rank enumerated dirs
  namespace bf = boost::filesystem;
  bf::path root_dir = bf::current_path();
  const bf::path& source_therm = root_dir / "therm.dat";
  const bf::path& source_tran  = root_dir / "tran.dat";
  const bf::path& target_therm = "therm.dat";
  const bf::path& target_tran  = "tran.dat";
  {
    bf::path rundir = bf::path("p") / 
      boost::lexical_cast< std::string >( rank );

    bf::create_directory("p");
    bf::create_directory(rundir);

    // Change directory for further computation
    mytools::filesystem::change_directory(rundir);

    // copy files
    if ( bf::exists( source_therm ) )
    {
      if ( bf::exists( target_therm ) ) bf::remove( target_therm );
      bf::copy_file( source_therm, target_therm );
    }

    if ( bf::exists( source_tran ) )
    {
      if ( bf::exists( target_tran ) ) bf::remove( target_tran );
      bf::copy_file( source_tran,  target_tran  );
    }
  }

  return root_dir;

}
