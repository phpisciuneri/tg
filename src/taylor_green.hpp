#ifndef TAYLOR_GREEN_HPP_
#define TAYLOR_GREEN_HPP_

#include <fstream>
#include <boost/scoped_ptr.hpp>
#include <zoltan_cpp.h>
#include <paragon.h>

#include "defs.hpp"
#include "simparam.hpp"
#include "domain.hpp"
#include "subdomain.hpp"
#include "field.hpp"
#include "iplmc.hpp"
#include "migratezobjects.hpp"
#include "parallel_io.hpp"


namespace iplmcfd
{

  //!
  //! \file taylor_green.hpp
  //! \class TaylorGreen
  //! \brief
  //!
  class TaylorGreen
  {
  public:

    //!
    //! \brief
    //!
    TaylorGreen( const simparam& p, const Domain& d, Subdomain& s );

    //!
    //! \brief
    //!
    ~TaylorGreen();

    //!
    //! \brief
    //!
    void initialize();

    //!
    //! \brief
    //!
    void step( std::size_t iter );

    //!
    //! \brief
    //!
    void output( std::size_t iter );

    //!
    //!
    //!
    void output_steps( std::ostream& out, std::size_t iter );

    //!
    //! \brief
    //!
    void update_time() { m_mc.update_time(); }

    //!
    //! \brief
    //!
    bool repartition();

    bool paragon_refinement();
    //!
    //! \brief
    //!
    void reset_output();
      
  private:

    //!
    //! \brief
    //!
    void init_zoltan();
    void init_paragon();
    void build_local_graph(GraphStruct *localGraph);
  private:

    //!
    const Domain& m_domain;
    //!
    Subdomain& m_subdomain;
    //!
    const simparam& m_simparam;
    //!
    Iplmc m_mc;
    //!
    boost::scoped_ptr< Zoltan > m_zz;
    //!
    //boost::scoped_ptr< MigrateZObjects > m_zob;
    MigrateZObjects m_zob;
    //!
    boost::scoped_ptr< Field > m_fields;
    //!
    boost::scoped_ptr< Parallel_io > m_io;

    //! time step
    real_t m_dt;
    //! walltime spent checking partition / migrating to new decomposition 
    real_t m_repart_time;
    //! walltime for a single iteration, MPI_Wtime returns double
    double m_step_wtime;

    //!
    std::vector< float > m_cell_wgt;
    
    //!
    //! \brief rankCommCost[i * numProcs + j] denotes the relative network commu cost between rank/partition i and j
    //!
    float *m_rankCommCost;
    ParagonStruct_t m_paragon;
  };

}

#endif // TAYLOR_GREEN_HPP_
