#include "simparam.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <boost/program_options.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/throw_exception.hpp>
#include <boost/filesystem/operations.hpp>

using boost::throw_exception;

// +---------------+
// | I/O OVERLOADS |
// +---------------+
// necessary for boost program options
namespace tvmet {

   // tvmet::Vector
   template< class T, size_t N >
   inline std::ostream& operator<<( std::ostream& os, const Vector< T, N >& vec )
   {
      for (size_t i=0; i<N-1; ++i)
         os << vec[i] << " ";
      os << vec[N-1];
      return os;
   }

   template< class T, size_t N >
   void validate( boost::any& v, const std::vector< std::string >& values, Vector< T, N >*, int )
   {
      Vector< T, N > vec;
      std::istringstream iss( values[0] );
      for (size_t i=0; i<N; ++i)
         iss >> vec[i] >> std::ws;
      v = vec;
   }

} // namespace tvmet

namespace std {

   // std::set
   inline ostream& operator<<( ostream& os, const set< int >& the_set )
   {
      BOOST_FOREACH( set<int>::const_reference iset, the_set )
         os << iset << " ";
      return os;
   }

   void validate( boost::any& v, const vector< string >& values, set< int >*, int )
   {
      set< int > the_set;
      istringstream iss( values[0] );
      while ( !iss.eof() )
      {
         int val;
         iss >> val >> ws;
         the_set.insert( val );
      }
      v = the_set;
   }

   // std::map
   typedef map< string, iplmcfd::real_t > name_value_map;
   inline ostream& operator<<( ostream& os, const name_value_map& the_map )
   {
      BOOST_FOREACH( name_value_map::const_reference imap, the_map )
         os << imap.first << " " << imap.second << " ";
      return os;
   }

   void validate( boost::any& v, const vector< string >& values, name_value_map*, int )
   {
      name_value_map the_map;
      istringstream iss( values[0] );
      while ( !iss.eof() )
      {
         string name;
         iss >> name >> ws;
         iss >> the_map[ name ] >> ws;
      }
      v = the_map;
   }

   // std::vector
   template< class T >
   inline ostream& operator<<( ostream& os, const vector< T >& vec )
   {
      for (size_t i=0; i<vec.size()-1; i++)
         os << vec[i] << " ";
      os << vec[ vec.size()-1 ];
      return os;
   }

   template< class T >
   inline void validate( boost::any& v, const vector< string >& values, vector< T >*, int )
   {
      vector< T > vec;
      istringstream iss( values[0] );
      while ( !iss.eof() )
      {
         T val;
         iss >> val >> ws;
         vec.push_back( val );
      }
      v = vec;
   }

} // namespace std

namespace iplmcfd {

   // +----------+
   // | SIMPARAM |
   // +----------+
   simparam::simparam( int argc, char* argv[] )
      : m_initial_path(boost::filesystem::current_path())
   {

      if ( argc<=0 ) 
         throw not_initialized();

      // max values
      real_t rmax  = 999999;
      int    imax  = 999999;
      size_t stmax = 999999;

      namespace po = boost::program_options;

      // +---------+
      // | general |
      // +---------+
      po::options_description general("General Options");
      general.add_options()
         ("help", "Print this help message." )
         ("write_default", "Writes a default simparam.in file." )
         ("input_file", po::value< std::string >(&filename)->default_value( "simparam.in" ),
         "Name of input file containing configuration options." );

      // +-----------------------+
      // | simulation parameters |
      // +-----------------------+
      po::options_description sim_param("Simulation Parameters");
      sim_param.add_options()
         ("iter_end", po::value< int >(&iter_end)->default_value( 0 ),
         "Maximum number of iterations to take."  )
         ("time_end", po::value< real_t >(&time_end)->default_value( rmax ),
         "Maximum nondimensional simulation time." )
         ("walltime", po::value< real_t >(&walltime)->default_value( rmax ),
         "Run for given walltime (in minutes)." )
         ("Re", po::value< real_t >(&Re)->default_value( 100 ),
         "Reynolds number for Taylor Green Vortex." )
         ("repart_freq", po::value< int >(&repart_freq)->default_value( 10 ),
         "Frequency to check the quality of the decomposition." )
         ("zoltan_debug_level", po::value< std::string >(&zoltan_debug_level)->default_value( "0" ),
         "Set Zoltan's debug level [0, 10]." );

      // +--------------------------+
      // | decomposition parameters |
      // +--------------------------+
      po::options_description decomp_param("Decomposition Parameters");
      decomp_param.add_options()
         ("save_cpuload", po::value< bool >(&save_cpuload)->default_value( true ),
         "[true,false]: Output cpu load for restart using wmetis decomposition type." )
         ("decomp_type", po::value< std::string >(&decomp_type)->default_value( "metis" ),
         "[metis]: Use Metis graph parition.\n"
         "[wmetis]: Use Metis weighted graph partition.\n"
         "\tRequires a cpuload file.\n"
         "[block]: Use block decomposition.  Specify ncores.")
         ("redecomp_type", po::value< std::string >(&redecomp_type)->default_value( "zoltan" ),
         "[parmetis]: Use ParMetis graph reparition.\n"
	 "[paragon]: Use Paragon graph repartition")
         ("degree_refine_parallelism", po::value< int >(&degree_refine_parallelism)->default_value( 8 ),
          "Number of MPI ranks used to do the refinement in parallel. It should be a value of [0, n/2]." )
         ("shuffle_refine_times", po::value< int >(&shuffle_refine_times)->default_value( 8 ),
         "Number of shuffle refinement times Paragon performs. It can be any value. The higher the value is, the better the resutling decomposition would be and the longer the refinement would take." )
	 ("imbalance_tolerance", po::value< real_t >(&imbalance_tolerance)->default_value( 1.02 ),
          "The amount of load imbalance we can have among partitions. It should be a value between 1 and 2. A value of 1.02 means we can have around 2% load imbalance among partitions."  )
	 ("degree_resource_contention", po::value< real_t >(&degree_resource_contention)->default_value( 0 ),
         "The degree of intra-node shared resource contention. It should be a value between 0 and 1. 0 means no contention and we should focus on minimizing the impact of network communication heterogeneity. 	1 means we should avoid intra-node shared resrouce contention at all cost. A value between 0 and 1 means we should consider both." )
	 ("ncores", po::value< int_v >(&ncores)->default_value( int_v(1,1,1) ),
         "Number of partitions in each direction to use for block decomposition." );

      // +---------------------+
      // | geometry parameters |
      // +---------------------+
      po::options_description geom_param("Geometry Parameters");
      geom_param.add_options()
         ("ppc", po::value< int_v >(&ppc)->default_value( int_v(5,10,15) ),
         "Particles per cell (min,mean,max): min and max are enforced through clustering and cloning." )
         ("use_clone_cluster", po::value< bool >(&use_clone_cluster)->default_value( true ),
         "[true,false]: Enforce ppc min and max through cloning and clustering." )
         ("xmin", po::value< real_v >(&xmin)->default_value( real_v(0,-1,-1) ),
         "Nondimensional minimum extents of cartesian domain." )
         ("xmax", po::value< real_v >(&xmax)->default_value( real_v(1,1,1) ),
         "Nondimensional maximum extents of cartesian domain." )
         ("ngrid", po::value< size_v >(&ngrid)->default_value( size_v(11,11,11) ),
         "Number of cells / FD points (nx, ny, nz)." )
         ("periods", po::value< bool_v >(&periods)->default_value( bool_v(false,false,false) ),
         "[true,false]: Periodicity of boundaries (x, y, z)." );

      // +------------------------------+
      // | constants & model parameters |
      // +------------------------------+
      po::options_description conmod_param("Constants and Model Parameters");
      conmod_param.add_options()
         ("Comega", po::value< real_t >(&Comega)->default_value( real_t(8) ),
         "LMSE mixing coefficent." )
         ("Sc", po::value< real_t >(&Sc)->default_value( real_t(.75), "0.75" ),
         "Schmidt number." )
         ("mu_power", po::value< real_t >(&mu_power)->default_value( real_t(.7), "0.7" ),
         "Sutherland's law exponent: mu = T^mu_power." )
         ("stable", po::value< real_t >(&stable)->default_value( real_t(.8), "0.8" ),
         "CFL safety factor: 0 < stable < 1." );

      // +--------------------------+
      // | I/O and trace parameters |
      // +--------------------------+
      // default trace_ranks contain root process only
      std::set<int> def_tr;
      def_tr.insert( 0 );
      //
      po::options_description iotrace_param("I/O and Trace Parameters");
      iotrace_param.add_options()
         ("iter_fdump", po::value< size_t >(&iter_fdump)->default_value( size_t(1000) ),
         "Iteration frequency to dump Eulerian data." )
         ("iter_output_steps", po::value< size_t >(&iter_output_steps)->default_value( size_t(1) ),
         "Iteration frequency to update steps.out." )
         ("trace_flush", po::value< bool >(&trace_flush)->default_value( true ),
         "[true,false]: flush output of tracer at every event." )
         ("trace_flush_iter", po::value< size_t >(&trace_flush_iter)->default_value( size_t(100) ),
         "Iteration frequency to flush tracer." )
         ("trace_rank", po::value< std::set<int> >(&trace_ranks)->default_value( def_tr ),
         "Space separated list of ranks to be flushed by tracer." )
         ("carlo_log_flush", po::value< bool >(&carlo_log_flush)->default_value( true ),
         "[true,false]: flush output of MC log at every event." )
         ("carlo_log_level", po::value< int >(&carlo_log_level)->default_value( 0 ),
         "MC activity log.\n"
         "[0]: disable logs.\n"
         "[1]: log summary at end of each timestep.\n"
         "[2]: log all activity.\n" );

      // +----------------------+
      // | chemistry parameters |
      // +----------------------+
      po::options_description chem_param("Chemistry Parameters");
      chem_param.add_options()
         ("chem", po::value< std::string >(&chem)->default_value( "mfc5"),
         "Chemistry mechanism: name of preinstalled mechanism or a chemkin mechanism filename.")
         ("chem_type", po::value< std::string>(&chem_type)->default_value( "inert"),
         "[full]: use finite rate kinetics, mix T.\n"
         "[fullmixh]: use finite rate kinetics, mix h.\n"
         "[fullmixrt]: use finite rate kinetics, mix RT.\n"
         "[inert]: just mixing, no reaction." )
         ("chem_solver", po::value< std::string >(&chem_solver)->default_value( "vode" ),
         "[isat, vode, lsoda, vodelsoda]: chemistry ODE integrator." )
         ("chem_reltol", po::value< real_t >(&chem_reltol)->default_value( real_t(1.e-6), "1.e-6" ),
         "Relative tolerence in ODE integration of kinetics." )
         ("chem_abstol", po::value< real_t >(&chem_abstol)->default_value( real_t(1.e-6), "1.e-6" ),
         "Absolute tolerance in ODE integration of kinetics." )
         ("chem_isatetol", po::value< real_t >(&chem_isatetol)->default_value( real_t(1.e-5), "1.e-5" ),
         "ISAT error tolerance." )
         ("chem_rxnstart", po::value< int >(&chem_rxnstart)->default_value( 0 ),
         "Iteration number to begin chemical reaction." )
         ("chem_rxnend", po::value< int >(&chem_rxnend)->default_value( imax ),
         "Iteration number to end chemical reaction." )
         ("chem_tran", po::value< std::string >(&chem_tran)->default_value( "tran.dat" ),
         "Chemkin transport coefficient database file." )
         ("chem_therm", po::value< std::string >(&chem_therm)->default_value( "therm.dat" ),
         "Chemkin thermal coefficient database file." );

      // options that can be used on the command line
      po::options_description cmdline_options("Command Line Options");
      cmdline_options.add( general );
      cmdline_options.add( sim_param );
      cmdline_options.add( decomp_param );
      cmdline_options.add( geom_param );
      cmdline_options.add( conmod_param );
      cmdline_options.add( iotrace_param );
      cmdline_options.add( chem_param );

      // options that can be set thru a config file
      po::options_description config_file_options("Configuration Options");
      config_file_options.add( sim_param );
      config_file_options.add( decomp_param );
      config_file_options.add( geom_param );
      config_file_options.add( conmod_param );
      config_file_options.add( iotrace_param );
      config_file_options.add( chem_param );

      // positional option for input filename
      po::positional_options_description p;
      p.add( "input_file", -1 );

      po::variables_map vm;

      // parse command line first it takes precedence over config file
      po::store( po::command_line_parser( argc, argv ).
         options(cmdline_options).positional(p).run(), vm);
      po::notify(vm);

      // parse config file
      std::ifstream ifs( filename.c_str() );
      if ( !ifs ) 
         throw std::ios_base::failure("Simparam: cannot open file: " + filename );
      po::store( po::parse_config_file( ifs, config_file_options ), vm );
      po::notify(vm);

      // PRINT HELP TO STDOUT AND RETURN
      if ( vm.count( "help" ) )
      {
         std::cout << cmdline_options << std::endl;
         return;
      }

      // +-----------------------+
      // | hidden and hard-coded |
      // +-----------------------+
      // eventually add this stuff to program_options
      // for now it is hard-coded
      debug_wait = false;
      restart = false;

      // +-----------+
      // | calculate |
      // +-----------+
      walltime *= 60; // convert min to sec
      Re_lambda = std::sqrt( 20 * Re / 3 );

      display_input_params();

   }

   // +----------------------+
   // | DISPLAY_INPUT_PARAMS |
   // +----------------------+
   void simparam::display_input_params() const
   {
     
     int rank;
     MPI_Comm_rank( MPI_COMM_WORLD, &rank );
     
     if ( rank == 0 )
     {
       // writes parameter values to info.out
       std::ofstream out( "info.out" );

       // unfortunately this is hard-coded
       out << "Simulation Parameters" << std::endl;
       out << "=====================" << std::endl;
       out << "iter_end " << iter_end << "\n"
         << "time_end " << time_end << "\n"
         << "walltime " << walltime << "\n"
         << "save_cpuload " << save_cpuload << "\n"
         << "decomp_type " << decomp_type << "\n"
         << "ncores " << ncores << "\n"
         << "ppc " << ppc << "\n"
         << "use_clone_cluster " << use_clone_cluster << "\n"
         << "xmin " << xmin << "\n"
         << "xmax " << xmax << "\n"
         << "ngrid " << ngrid << "\n"
         << "periods " << periods << "\n"
         << "Comega " << Comega << "\n"
         << "Sc " << Sc << "\n"
         << "mu_power " << mu_power << "\n"
         << "stable " << stable << "\n"
         << "iter_fdump " << iter_fdump << "\n"
         << "iter_output_steps " << iter_output_steps << "\n"
         << "trace_flush " << trace_flush << "\n"
         << "trace_flush_iter " << trace_flush_iter << "\n"
         << "trace_ranks " << trace_ranks << "\n"
         << "carlo_log_flush " << carlo_log_flush << "\n"
         << "carlo_log_level " << carlo_log_level << "\n"
         << "chem " << chem << "\n"
         << "chem_type " << chem_type << "\n"
         << "chem_solver " << chem_solver << "\n"
         << "chem_reltol " << chem_reltol << "\n"
         << "chem_abstol " << chem_abstol << "\n"
         << "chem_isatetol " << chem_isatetol << "\n"
         << "chem_rxnstart " << chem_rxnstart << "\n"
         << "chem_rxnend " << chem_rxnend << "\n"
         << "chem_tran " << chem_tran << "\n"
         << "chem_therm " << chem_therm << "\n";
     
     } // rank 0

   }

} //   namespace iplmcfd 
