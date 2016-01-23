#ifndef IPLMCFD_SIMPARAM_HPP_
#define IPLMCFD_SIMPARAM_HPP_

#include "defs.hpp"

#include <map>
#include <vector>
#include <string>
#include <set>
#include <boost/array.hpp>
#include <boost/filesystem/path.hpp>
#include <tvmet/Matrix.h>

namespace iplmcfd {

  /** \brief Process simulation parameters
  * \ingroup singleton
  * \todo For a restarting simulation implement parameter consistency checks
  * \todo Rename ngrid to ncell (this is how it is interpreted in the rest of the code!)
  */
  class simparam
  {
    // Using Singleton interface
    simparam(int argc, char* argv[]);
    ~simparam(){};

  public:
    static const simparam& instance(int argc=0, char* argv[] = 0) {
      static simparam singleton(argc, argv);
      return singleton;
    }

  public:
    class not_initialized : public std::exception 
    {
    public: 
      not_initialized () {} 
      virtual const char* what() const throw()
      { return "first call to simparam should include commandline arguments"; }
    };


    // DATA MEMBERS (all read-only)
  public:

    typedef std::map<std::string, real_t > name_value_map;

    std::string filename;

    // +-----------------------+
    // | simulation parameters |
    // +-----------------------+
    bool restart;
    int iter_end;
    real_t time_end;
    real_t walltime;
    real_t Re;
    real_t Re_lambda;
    int repart_freq;
    std::string zoltan_debug_level;

    // +--------------------------+
    // | decomposition parameters |
    // +--------------------------+
    bool save_cpuload;
    std::string decomp_type;
    int_v ncores;
      
    // +--------------------------+
    // | redecomposition parameters |
    // +--------------------------+
    std::string redecomp_type;      // can be either parmetis or zoltan
//    std::string refine_type;      // can be none, serial, or parallel
    int degree_refine_parallelism;  // can be a value between 0 and n/2
    int shuffle_refine_times;	    // can be any values. The higher the value is, the better the resulting decompositions will be and the longer the refinement took
    real_t imbalance_tolerance;	    // the amount of imblance we can have among partitions. 1.02 means that we allow 2% load imbalnace among partitions.
    real_t degree_resource_contention;	    // can be any value between 0 and 1: 0 means no intra-node resrouce contention, 1 means should avoid intra-node resrouce contention at all cost
    
    // +---------------------+
    // | geometry parameters |
    // +---------------------+
    int_v ppc;
    bool use_clone_cluster;
    real_v xmin;
    real_v xmax;
    size_v ngrid;
    bool_v periods;

    // +------------------------------+
    // | constants & model parameters |
    // +------------------------------+
    real_t Comega;
    real_t Sc;
    real_t mu_power;
    real_t stable;

    // +--------------------------+
    // | I/O and trace parameters |
    // +--------------------------+
    size_t iter_fdump;
    size_t iter_output_steps;
    bool  trace_flush;
    size_t trace_flush_iter;
    std::set<int> trace_ranks;
    bool carlo_log_flush;
    int  carlo_log_level;

    // +----------------------+
    // | chemistry parameters |
    // +----------------------+
    std::string chem;
    std::string chem_type;
    std::string chem_solver;
    real_t chem_reltol;
    real_t chem_abstol;
    real_t chem_isatetol;
    int chem_rxnstart;
    int chem_rxnend;
    std::string chem_tran;
    std::string chem_therm;

    // THIS STUFF REQUIRES SOME TLC
    // +-----------------------+
    // | hidden and hard-coded |
    // +-----------------------+
    // for now this is just a catch-all of things I don't want to deal with
    // work on this crap
    bool debug_wait;

    // save root path
    const boost::filesystem::path& initial_path() const { return m_initial_path; }

   private:

     void display_input_params() const;

     boost::filesystem::path m_initial_path; 

   };

} // namespace iplmcfd


#endif // IPLMCFD_SIMPARAM_HPP_
