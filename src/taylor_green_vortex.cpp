#include "taylor_green_vortex.hpp"

#include <fstream>
#include <boost/filesystem/path.hpp>
#include <mpi.h>

#include "diagnostic/tracer.hpp"

namespace iplmcfd
{

  // +--------------------------------------+
  // | TaylorGreenVortex::TaylorGreenVortex |
  // +--------------------------------------+
  TaylorGreenVortex::TaylorGreenVortex( real_t Re, real_t a, real_t nu )
    : m_A( Re * a * nu ), m_a( a ), m_nu( nu )
  {
    tracer::scope _( "TaylorGreenVortex::TaylorGreenVortex" );

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    // append info about TGV to info.out
    if ( rank == 0 )
    {
      typedef boost::filesystem::path bfp;
      bfp info = bfp("..") / bfp("..") / "info.out";
      std::ofstream out( info.c_str(), std::ofstream::app );
      out << std::endl;
      out << "Taylor Green Vortex" << std::endl;
      out << "===================" << std::endl;
      out << "Re = " << Re << std::endl;
      out << "tau [s] = " << 1 / ( m_A * m_a ) << std::endl;
      out << "nu [cm^2 / s] = " << nu << std::endl;
    }

  };

  // +------------------------+
  // | TaylorGreenVortex::uvw |
  // +------------------------+
  void TaylorGreenVortex::uvw( real_t t, const real_v& xyz, real_v& vel ) const
  {
    // convenience aliases
    const real_t& A  = m_A;
    const real_t& a  = m_a;
    const real_t& nu = m_nu;

    const real_t& x = xyz[XDIM]; 
    const real_t& y = xyz[YDIM]; 
    const real_t& z = xyz[ZDIM]; 

    real_t& u = vel[XDIM];
    real_t& v = vel[YDIM];
    real_t& w = vel[ZDIM];

    using namespace std;

    real_t F = exp( -2 * a * a * nu * t );
    u =  A * sin( a * x ) * cos( a * y ) * F;
    v = -A * cos( a * x ) * sin( a * y ) * F;
    w = 0;

  }

} // namespace iplmcfd
