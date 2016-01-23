#include "../iplmc.hpp"

#include <cmath>
#include <fstream>
#include <boost/foreach.hpp>

#include "../diagnostic/cpublock.hpp"
#include "../diagnostic/tracer.hpp"

namespace iplmcfd {

  // +-------------+
  // | SCALAR_STEP |
  // +-------------+
  void Iplmc::scalar_step( real_t dt, std::size_t iter )
  {
    tracer::scope _( "Iplmc::scalar_step" );

    scalar_average();
    integrate_particles( dt, iter );

  }


  // +----------------------------+
  // | Iplmc::integrate_particles |
  // +----------------------------+
  void Iplmc::integrate_particles( real_t dt, std::size_t iter )
  {
    tracer::scope _( "Iplmc::integrate_particles" );

    if ( m_simparam.chem_type == "inert" || 
      iter < m_simparam.chem_rxnstart  || 
      iter >= m_simparam.chem_rxnend ) 
      scalar_integ_inert(dt);
    else if ( m_simparam.chem_type == "full" ) 
      scalar_integ_full(dt);
    else 
      throw( std::invalid_argument( "Unknown scalar integration strategy: " 
      + m_simparam.chem_type ) );

  }

  // +-----------------------+
  // | full_reaction_substep |
  // +-----------------------+
  inline void full_reaction_substep( 
    const mix_t& mix,
    cpr_t& cpr,
    cpr_t::phi& phi, 
    mix_t::pressure& press, 
    real_t dt )
  {
    static cpr_t::phi phi0( mix );
    phi0 = phi;
    if ( !cpr.react( phi, press, dt ) )
    {
      static int count = 0;
      static std::ofstream file;
      if ( !file.is_open() )
      {
        std::cerr << "Failed reaction encountered." << std::endl;
        file.open( "failed.reactions.bin" );
        file << phi0.size() << "\n";
      }
      if (++count < 1000) // don't overwrite
      {
        file.write( (const char*)&phi0[0], sizeof( real_t )*phi0.size() );
        file.write( (const char*)&dt, sizeof(dt) );
      }
    }
  }

  // +---------------------------+
  // | Iplmc::scalar_integ_inert | 
  // +---------------------------+
  void Iplmc::scalar_integ_inert( real_t dt )
  {
    tracer::scope _( "Iplmc::scalar_integ_inert" );

    BOOST_FOREACH( cell_mp::reference gidcell, m_cells )
    {
      Cell& cell = gidcell.second;
      real_t omega_mix = omega( cell );

      cpublock cpub( cell.ctr.cputime, cell.contained().size() );

      BOOST_FOREACH( Particle& p, cell.contained() )
      {
        // Mix particle (IEM)
        real_t decay = std::exp( -omega_mix * dt );
        for ( std::size_t j=0; j!=p.phi.size(); ++j )
          p.phi[j] = ( p.phi[j] - cell.scalars.phi[j] )*decay + cell.scalars.phi[j];
        p.f = ( p.f - cell.scalars.f )*decay + cell.scalars.f;

        p.srt = 0;
      }
    }

  }

  // +--------------------------+
  // | Iplmc::scalar_integ_full |
  // +--------------------------+
  void Iplmc::scalar_integ_full( real_t dt )
  {
    tracer::scope _( "Iplmc::scalar_integ_full" );

    BOOST_FOREACH( cell_mp::reference gidcell, m_cells )
    {
      Cell& cell = gidcell.second;
      real_t omega_mix = omega( cell );

      cpublock cpub( cell.ctr.cputime, cell.contained().size() );

      BOOST_FOREACH( Particle& p, cell.contained() )
      {
        // Mix particle (IEM)
        real_t decay = std::exp( -omega_mix * dt );
        for ( std::size_t j=0; j!=p.phi.size(); ++j )
          p.phi[j] = ( p.phi[j] - cell.scalars.phi[j] )*decay + cell.scalars.phi[j];
        p.f = ( p.f - cell.scalars.f )*decay + cell.scalars.f;
        
        // Get current RT
        real_t RTn = m_chem->R( p.phi.Y() ) * real_t( p.phi.T() );
        
        // react
        mix_t::pressure press( m_chem->atm2dynescc( 1 ) );
        full_reaction_substep( *m_chem, *m_cpr, p.phi, press, dt );
        cell.ctr.chemload += m_cpr->num_rhs();
	m_chem->realize( p.phi.Y() );
        
        // Set source RT
        real_t RTnp1 = m_chem->R( p.phi.Y() ) * real_t( p.phi.T() );
        p.srt = ( RTnp1 - RTn ) / dt;

      }
    }

  }

} // namespace iplmcfd
