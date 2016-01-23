#ifndef IPLMCFD_DOMAIN_HPP_
#define IPLMCFD_DOMAIN_HPP_

#include <cassert>
#include <mytools/numeric/indsub.hpp>

#include "defs.hpp"
#include "simparam.hpp"

namespace iplmcfd {

  //! 
  //! \file domain.hpp
  //! \class Domain
  //! \brief Stores and provides information about the domain of the simulation
  //!
  //! \todo remove dependency on mytools/numeric/indsub.hpp
  //!
  class Domain
  {
  public:

    //! global index for "out of bounds" 
    static const gid_t OOB = static_cast<gid_t>(-1);

    //!
    //! \brief Constructs domain object
    //!
    //! \param[in] param Container of parsed simulation parameters
    //!
    Domain( const simparam& param );

    //!
    //! \brief Convert cartesian indices to global ID
    //!
    //! \param[in] sub cartesian i,j,k index (can be negative to accomodate 
    //!                periodic boundary conditions)
    //! \returns the corresponding global ID
    //!
    gid_t ind2gid( int_v ind ) const
    {
      for (int i=0; i<NDIM; ++i)
      {
        if ( m_periods[i] ) 
        { // wrap around
          ind[i] %= static_cast<int>( m_indsub.shape(i) );
          if ( ind[i] < 0 ) ind[i] += m_indsub.shape(i);				
        }
      }
      if (  m_indsub.subs_out_of_bounds(ind) ) return OOB;
      return m_indsub.ind(ind);	
    }

    //!
    //! \brief Converts global ID to cartesian i,j,k index
    //!
    //! \param[in] I Global ID
    //! \returns the corresponding i,j,k index
    //!
    index_v gid2ind( const gid_t& I ) const
    {
      index_v sub;
      gid2ind( I, sub );
      return sub;
    }
   
    //!
    //! \brief Converts global ID to cartesian i,j,k index (non-temp)
    //!
    //! \param[in] I Global ID
    //! \param[out] sub corresponding i,j,k index
    //!
   void gid2ind( const gid_t& I, index_v& ind ) const
   {
      assert( !m_indsub.ind_out_of_bounds(I) );
      m_indsub.sub( I, ind );
   }
    
   //!
   //! \brief Test if a global ID is part of the domain
   //!
   //! \param[in] I Global ID
   //! \returns True for out of bounds, false for in bounds
   //!
   bool out_of_bounds( const gid_t& I ) const { return I==OOB; }

   //!
   //! \brief Get the number of nodes in each coordinate direction
   //!
   //! \returns node extents
   //! \todo Better name than shape?
   //!
   const size_v& shape() const { return m_shape; }

   //!
   //! \brief Get the uniform grid spacing in each coordinate direction
   //!
   //! \returns grid spacing
   //!
   const real_v& dx() const { return m_dx; }

   //!
   //! \brief
   //!
   real_t delta() const { return m_delta; }

   //!
   //! \brief
   //!
   real_t dV() const { return m_dV; }

   //!
   //! \brief Get the inverse grid spacing in each coordinate direction
   //!
   //! \returns inverse grid spacing
   //!
   const real_v& h() const { return m_h; }

   //!
   //! \brief Get the maximum (dimensional) spatial extents in each coordinate 
   //!        direction
   //!
   //! \returns maximum (dimensional) spatial extents
   //!
   const real_v& xmax() const { return m_xmax; }

   //!
   //! \brief Get the minimum (dimensional) spatial extents in each coordinate 
   //!        direction
   //!
   //! \returns minimum (dimensional) spatial extents
   //!
   const real_v& xmin() const { return m_xmin; }

   //! 
   //! \brief Get the total number of nodes
   //!
   //! \returns product of node extents
   //!
   std::size_t npoints() const { return tvmet::product( m_shape ); }

   //!
   //! \brief Get the periodicity of the domain
   //!
   //! \returns periodicity of domain in each dim
   //!
   const bool_v& periodic() const { return m_periods; }
   
   //!
   //! \brief Get the cartesian position of a node
   //!
   //! \param[in] index cartesian (i,j,k) index
   //! \returns corresponding position (x,y,z)
   //!
   real_v xyz( const index_v& index ) const 
   { 
     return real_v( m_xmin + m_dx*index ); 
   }

   //!
   //! \brief Get the cartesian position of a node
   //!
   //! \param[in] I global index
   //! \returns corresponding position (x,y,z)
   //!
   real_v xyz( gid_t I ) const 
   { 
     index_v index;
     gid2ind( I, index );
     return real_v( m_xmin + m_dx*index ); 
   }

   //!
   //! \brief Get the cartesian position of a node (non-temp)
   //!
   //! \param[in] index cartesian (i,j,k) index
   //! \returns corresponding position (x,y,z)
   //!
   void xyz( const index_v& index, real_v& x ) const 
   { 
     x = m_xmin + m_dx*index;
   }

  private:

    //! ID to index converter class
    //! \todo Replace with local functionality
    mytools::numeric::indsub<NDIM> m_indsub;

    bool_v m_periods; //!< periodicity of domain in each dim
    size_v m_shape;   //!< number of nodes in each dim
    real_v m_dx;      //!< grid spacing in each dim
    real_v m_h;       //!< inverse grid spacing in each dim
    real_t m_dV;      //!< cell volume
    real_t m_delta;   //!<
    real_v m_xmax;    //!< maximum spatial extents (dimensional) in each dim
    real_v m_xmin;    //!< minimum spatial extents (dimensional) in each dim

  };

} // namespace iplmcfd

#endif // IPLMCFD_DOMAIN_HPP_
