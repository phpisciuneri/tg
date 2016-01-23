#ifndef IPLMCFD_SCALARAVERAGES_HPP_
#define IPLMCFD_SCALARAVERAGES_HPP_

#include <algorithm>
#include <boost/array.hpp>
#include <mpi.h>

#include "defs.hpp"

namespace iplmcfd {

  //!
  //! \file scalar_averages.hpp
  //! \struct ScalarAverages
  //! \brief
  //!
  struct ScalarAverages
  {
  
    //!
    cpr_t::phi phi;
    //! \todo what is up with the union here?
    union {
      struct {
        real_t mass;
        real_t f;
        real_t h;
        real_t RT;
      };
      boost::array< real_t, 4> scalars;
    };

    //!
    //! \brief
    //!
    ScalarAverages( const mix_t& mix ) : phi( mix ) {}

    //!
    //! \brief
    //!
    void reset()
    {
      std::fill( phi.begin(), phi.end(), real_t( 0 ) );
      std::fill( scalars.begin(), scalars.end(), real_t( 0 ) );
    }

    //!
    //! \brief
    //!
    void append( const ScalarAverages& s )
    {
      for( std::size_t i=0; i<phi.size(); i++)
        phi[i] += s.phi[i];
      for(size_t i=0; i<scalars.size(); i++)
        scalars[i] += s.scalars[i];
    }

    //!
    //! \brief
    //!
    void normalize()
    {
      for( std::size_t i = 1; i<scalars.size(); i++)
        scalars[i] /= mass;
      for(size_t i=0; i<phi.size(); i++)
        phi[i] /= mass;
    }

    //!
    //! \brief
    //!
    void append_to_message_real(
      std::vector< MPI_Aint >& displ
      , std::vector< int >& count
      )
    {
      MPI_Aint addr;
      // Append cell composition
      MPI_Address(&phi[0],&addr);
      displ.push_back(addr);
      count.push_back(phi.size());
      // Append other properties
      MPI_Address(&scalars,&addr);
      displ.push_back(addr);
      count.push_back(scalars.size());
    }

  };

} // namespace iplmcfd

#endif // IPLMCFD_SCALARAVERAGES_HPP_
