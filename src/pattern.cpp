#include "pattern.hpp"

#include "diagnostic/tracer.hpp"

namespace iplmcfd {

  // +------------------------+
  // | BoxPattern::BoxPattern |
  // +------------------------+
  BoxPattern::BoxPattern( int depth )
  {
    tracer::scope _( "BoxPattern::BoxPattern" );

    for (int i=-depth; i<=depth; ++i)
      for (int j=-depth; j<=depth; ++j)
        for (int k=-depth; k<=depth; ++k)
        {
          int_v offset( i, j, k );
          if ( i==0 && j==0 && k==0 ) continue;
          m_elements.push_back( offset );
        }

  }

  // +--------------------------+
  // | StarPattern::StarPattern |
  // +--------------------------+
  StarPattern::StarPattern( int depth )
  {
    tracer::scope _( "StarPattern::StarPattern" );

    // we traverse +/- depth for each coordinate direction
    for (int idim=0; idim<NDIM; ++idim)
    {
      int_v offset( 0, 0, 0 );
      for (int idepth=1; idepth<=depth; ++idepth)
      {
        offset[idim] = idepth;
        m_elements.push_back( offset );
        offset[idim] = -idepth;
        m_elements.push_back( offset );
      }
    }

  }

} // namespace iplmcfd
