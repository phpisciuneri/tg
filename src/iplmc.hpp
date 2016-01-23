#ifndef IPLMCFD_IPLMC_HPP_
#define IPLMCFD_IPLMC_HPP_

#include <vector>
#include <map>
#include <boost/scoped_ptr.hpp>

#include "domain.hpp"
#include "subdomain.hpp"
#include "simparam.hpp"
#include "cell.hpp"
#include "field.hpp"
#include "particle_p2p.hpp"
#include "taylor_green_vortex.hpp"

namespace iplmcfd {

  //!
  //! \file iplmc.hpp
  //! \class Iplmc
  //! \brief
  //! 
  class Iplmc
  {
    
    friend class MigrateZObjects;

  public: 
    
    //!
    typedef std::map< gid_t, Cell > cell_mp;

  public:
    
    //!
    //! \brief
    //!
    Iplmc( const simparam& p, const Domain& d, const Subdomain& s );

    //!
    //! \brief
    //!
    real_t update_dt();

    //!
    //! \brief
    //!
    void update_time() { m_time += m_dt; }

    //!
    //! \brief
    //!
    real_t time() const { return m_time; }

    //!
    //! \brief
    //!
    void initialize();

    //!
    //! \brief Integrate particles in physical space
    //!
    void position_step( real_t dt, std::size_t iter );

    //!
    //! \brief Integrate particles in composition space
    //!
    void scalar_step( real_t dt, std::size_t iter );

    //!
    //! \brief
    //!
    real_t ke() const { return m_ke; }

    //!
    //! \brief
    //!
    void register_fields( Field& fields );

    //!
    //! \brief
    //!
    void update_fields();

    //!
    //! \brief Get vector of cell weights
    //!
    void cell_wgts( std::vector< float >& wgts );

  private:

    //!
    //!
    //!
    void scalar_average();

    //!
    //!
    //!
    void integrate_particles( real_t dt, std::size_t iter );

    //!
    //!
    //!
    void scalar_integ_inert( real_t dt );

    //!
    //!
    //!
    void scalar_integ_full( real_t dt );

  private:

    //!
    //!
    //!
    real_t omega( const Cell& cell );

  private:

    //!
    boost::scoped_ptr< Particle_p2p > m_p2p;
    //!
    std::size_t m_min_ppc;
    //!
    std::size_t m_mean_ppc;
    //!
    std::size_t m_max_ppc;
    //!
    const simparam& m_simparam;
    //!
    const Domain& m_domain;
    //!
    const Subdomain& m_subdomain;
    //!
    boost::scoped_ptr< TaylorGreenVortex > m_tgv;
    //!
    boost::scoped_ptr< mix_t > m_chem;
    //!
    boost::scoped_ptr< cpr_t > m_cpr;
    //!
    cell_mp m_cells;
    //!
    cell_mp m_neigh_cells;
    //!
    mix_t::pressure m_pref;
    //!
    real_t m_nuref;
    //!
    real_t m_dt;
    //!
    real_t m_time;
    //! kinetic energy of entire system
    real_t m_ke;

    //! output buffers
    std::vector< real_t > m_colors;
    std::vector< real_t > m_u;
    std::vector< real_t > m_v;
    std::vector< real_t > m_w;
    std::vector< real_t > m_f;
    std::vector< real_t > m_T;
    std::vector< real_t > m_O2;
    std::vector< real_t > m_CH4;
    std::vector< real_t > m_OH;
    std::vector< real_t > m_CO2;
    std::vector< real_t > m_NP;
    std::vector< real_t > m_cputime;

  };

} // namespace iplmcfd


#endif // IPLMCFD_IPLMC_HPP_
