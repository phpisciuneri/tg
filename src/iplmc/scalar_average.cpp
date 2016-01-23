#include "../iplmc.hpp"

#include <boost/foreach.hpp>

#include "../diagnostic/tracer.hpp"

namespace iplmcfd {

  // +-----------------------+
  // | Iplmc::scalar_average |
  // +-----------------------+
  void Iplmc::scalar_average()
  {
    tracer::scope _( "Iplmc::scalar_average" );

    // perform scalar averaging
    BOOST_FOREACH( cell_mp::reference gidcell, m_cells )
    {
      Cell& cell = gidcell.second;

      cell.scalars.reset();
      BOOST_FOREACH( const Particle& p, cell.contained() )
      {
        cell.scalars.mass += p.mass;
        cell.scalars.f    += p.mass * p.f;
        cell.scalars.h    += p.mass * m_chem->h( p.phi.T(), p.phi.Y() );
        cell.scalars.RT   += p.mass * m_chem->R( p.phi.Y() ) * real_t( p.phi.T() );
        for ( std::size_t j=0; j<p.phi.size(); j++ )
          cell.scalars.phi[j] += p.phi[j] * p.mass;
      }

    }

    // NORMALIZE GRIDS
    // keep an eye on zero mass grids
    std::stringstream zeromass;
    int nzeromass = 0;
    index_v sub;

    BOOST_FOREACH( cell_mp::reference gidcell, m_cells )
    {
      Cell& cell = gidcell.second;

      if ( cell.scalars.mass <= real_t( 0 ) )
      {
        m_domain.gid2ind( gidcell.first, sub );
        zeromass << m_subdomain.my_rank() << "/";
        for ( std::size_t i=0; i<sub.size(); ++i ) 
          zeromass << sub[i] << " ";
        zeromass << ", ";
        nzeromass++;
        continue;
      }

      cell.scalars.normalize();

    }

    //! \todo Maybe this ought to be fatal?
    if ( nzeromass > 0 )
    {
      std::cerr << "error: non-fatal: There were (" << nzeromass << ") "
        << "Zero-mass cells, see trace.out" << std::endl;
      tracer::instance().note(
        "!!! Zero mass grids Rank/Grid i j k: " + zeromass.str()
        );
      tracer::instance().flush();
    }

  }

} // namespace iplmcfd
