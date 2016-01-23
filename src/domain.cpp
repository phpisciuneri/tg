#include "domain.hpp"

#include <cmath>

#include "diagnostic/tracer.hpp"

namespace iplmcfd {

  // +----------------+
  // | Domain::Domain |
  // +----------------+
  Domain::Domain( const simparam& param ) 
  {
    tracer::scope _("Domain::Domain");
    
    m_periods = param.periods;
    m_shape   = param.ngrid;
    m_xmin    = param.xmin;
    m_xmax    = param.xmax;

    // calculate dx
    for (int i=0; i<NDIM; ++i)
    {
      if ( m_periods[i] )
        m_dx[i] = ( m_xmax[i] - m_xmin[i] ) / m_shape[i];
      else
        m_dx[i] = ( m_xmax[i] - m_xmin[i] ) / ( m_shape[i] - 1 );
    }

    // 2D hard-coded here
    //! \todo handle for the 3D case
    m_dV = m_dx[XDIM] * m_dx[YDIM];
    m_delta = 2 * std::sqrt( m_dV );

    m_h = 1 / m_dx;
    m_indsub.set_shape(m_shape);
  }

} // namespace iplmcfd
