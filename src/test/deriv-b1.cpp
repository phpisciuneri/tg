#include <iostream>
#include <boost/cstdlib.hpp>
#include <boost/exception/diagnostic_information.hpp>
#include <boost/tuple/tuple.hpp>
#include <mytools/mpi/environment.hpp>
#include <boost/foreach.hpp>

#include "../simparam.hpp"
#include "../domain.hpp"
#include "../subdomain.hpp"
#include "../grid.hpp"
#include "../localvector.hpp"
#include "../diff.hpp"
#include "../fd_family_b1.hpp"

int main( int argc, char* argv[] )
{

  using namespace iplmcfd;

  try 
  {

    mytools::mpi::environment::init(argc,argv);
    int commrank, commsize;
    boost::tuples::tie(commrank, commsize) = mpi_get_comm_rank_size();

    const simparam& p = simparam::instance(argc, argv);

    const real EPS = 1.e-12;
    const real THREE = 3;
    const real FOUR  = 4;
    const real FIVE  = 5;

    // initialize domain, subdomain and grid
    Domain d;
    Subdomain s( d, MPI_COMM_WORLD );
    Grid g( d, s );

    // allocate fcn
    Localvector< real > fcn;
    fcn.init( &(s.FDlinkage()), 1 );

    // initialize function
    BOOST_FOREACH(Subdomain::map_gid2lid_t::const_reference gidlid, s.homes() )
    {

      gid I = gidlid.first;
      lid i = gidlid.second;
      index_vector sub( d.ind2sub( I ) );
      real_vector xyz = d.xyz( sub );

      // f(x,y,z) = 3x + 4y + 5z
      // This is a linear function:
      //   first order derivatives should capture it exactly
      fcn[i] = THREE*xyz[XDIM] + FOUR*xyz[YDIM] + FIVE*xyz[ZDIM];

    }

    for (lid i=0; i< s.nhomes(); ++i)
    {
      typedef Diff< stencil_t > diff_t;
      diff_t d_dx( FD_family_b1::instance(), g.stencil[i], XDIM, d.dx() );
      real dfdx = d_dx( fcn.data() );

      // compare answer at each location with known result
      if ( std::abs( dfdx - THREE ) > EPS )
      {
        std::cout << "Error: Unexpected value for dfdx: "; 
        std::cout.precision( 16 );
        std::cout << dfdx << ".  Expected: 3" << std::endl;
        return 1;
      }

      diff_t d_dy( FD_family_b1::instance(), g.stencil[i], YDIM, d.dx() );
      real dfdy = d_dy( fcn.data() );
      if ( std::abs( dfdy - FOUR ) > EPS )
      {
        std::cout << "Error: Unexpected value for dfdy: " 
          << dfdy << ".  Expected: 4" << std::endl;
        return 1;
      }

      diff_t d_dz( FD_family_b1::instance(), g.stencil[i], ZDIM, d.dx() );
      real dfdz = d_dz( fcn.data() );
      if ( std::abs( dfdz - FIVE ) > EPS )
      {
        std::cout << "Error: Unexpected value for dfdz: " 
          << dfdz << ".  Expected: 5" << std::endl;
        return 1;
      }

    }


    std::cout << std::endl << "Success." << std::endl;

    return 0;

  }

  catch ( boost::exception& e )  {
    std::cerr << diagnostic_information(e) << std::endl;
  }
  catch ( const std::exception& e )  {
    std::cerr << "error: " << e.what() << std::endl;
  }
  // don't catch other errors
  return boost::exit_failure;

}
