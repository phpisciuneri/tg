#ifndef TAYLORGREENVORTEX_HPP_
#define TAYLORGREENVORTEX_HPP_

#include "defs.hpp"

namespace iplmcfd
{

  //!
  //! \file taylor_green_vortex.hpp
  //! \class TaylorGreenVortex
  //! \brief Implements the classic Taylor Green Vortex
  //!
  //! The Taylor Green Vortex is based off of the original work of Taylor
  //!   and Green:
  //!   G. I. Taylor and A. E. Green, "Mechanism of the Production of Small
  //!     Eddies from Large Ones," Proc. Roy. Soc. A, vol. 158, no. 895, 
  //!     pp. 499-521, (1937).
  //!
  //! Every effort was made to keep the terms and equations similar to what
  //!   is presented in the reference.  Where it is helpful the equation 
  //!   number corresponding to the original work is noted.
  //!
  class TaylorGreenVortex
  {
  public:

    //!
    //! \brief Constructor, initializes reused terms
    //!
    //! \param[in] Re Reynolds number
    //! \param[in] a characteristic inverse length scale
    //! \param[in] nu constant viscosity
    //!
    //! From the input parameters, the magnitude of the velocity can be 
    //!   specified using the definition of the Reynolds number:
    //!   A = Re * a * nu
    //!
    TaylorGreenVortex( real_t Re, real_t a, real_t nu );

    //!
    //! \brief Given x, y and t calculate u and v 
    //! 
    //! \param[in] t time
    //! \param[in] xyz position vector
    //! \param[out] vel velocity vector
    //!
    void uvw( real_t t, const real_v& xyz, real_v& vel ) const;

  private:

    real_t m_A;  //!< magnitude of the velocity
    real_t m_a;  //!< characteristic inverse length scale
    real_t m_nu; //!< constant viscosity

  };

} // namespace iplmcfd

#endif // TAYLORGREENVORTEX_HPP_