#ifndef IPLMCFD_PATTERN_HPP_
#define IPLMCFD_PATTERN_HPP_

#include <iostream>

#include "defs.hpp"

namespace iplmcfd {


  //!
  //! \file pattern.hpp
  //! \class Pattern
  //! \brief Common communication patterns for a structured cartesian grid.
  //!
  class Pattern
  {
  public:

    //!
    //! \brief Get array of relative location of elements in a pattern
    //!
    virtual const std::vector< int_v >& elements() const = 0;

    //!
    //! \brief Print to std::out all elements in a pattern
    //!
    //! \param[out] out output stream
    //! \param[in] p Pattern object
    //! 
    //! __Note:__ Intended for info/debug only.
    //!
    friend std::ostream& operator<<( std::ostream& out, const Pattern& p )
    {
      for ( std::size_t i=0; i<p.elements().size(); ++i )
        out << p.elements()[i][XDIM] << " "
        << p.elements()[i][YDIM] << " "
        << p.elements()[i][ZDIM] << std::endl;
      return out;
    }

  };

  //!
  //! \class BoxPattern
  //! \brief Box communication pattern for a structured cartesian grid
  //!
  class BoxPattern : public Pattern
  {
  public:

    //!
    //! \brief Constructs the BoxPattern
    //!
    //! \param[in] depth level of box pattern replication
    //! 
    //! __Example:__
    //! \code
    //!              |    *  *  *  *  *
    //!   *  *  *    |    *  *  *  *  *
    //!   *  ^  *    |    *  *  ^  *  *
    //!   *  *  *    |    *  *  *  *  *
    //!              |    *  *  *  *  *
    //!              
    //!  depth = 1         depth = 2
    //! \endcode
    //!
    BoxPattern( int depth );

    //!
    //! \brief Get elements in the box pattern
    //!
    //! \returns array of relative location of elements in the box pattern
    //!
    const std::vector< int_v >& elements() const { return m_elements; }

  private:

    //! relative offsets of elements in the box pattern
    std::vector< int_v > m_elements; 

  };

  //!
  //! \class StarPattern
  //! \brief Star communication pattern for a structured cartesian grid
  //!
  class StarPattern: public Pattern
  {
  public:

    //!
    //! \brief Constructs the StarPattern
    //!
    //! \param[in] depth level of star pattern replication
    //! 
    //! __2D Examples:__
    //! \code
    //!              |          * 
    //!      *       |          * 
    //!   *  ^  *    |    *  *  ^  *  *
    //!      *       |          * 
    //!              |          * 
    //!              
    //!  depth = 1         depth = 2
    //! \endcode
    StarPattern( int depth );

    //!
    //! \brief Get elements in the star pattern
    //!
    //! \returns array of relative location of elements in the star pattern
    //!
    const std::vector< int_v >& elements() const { return m_elements; }

  private:

    //! relative offsets of elements in the star pattern
    std::vector< int_v > m_elements;

  };

} // namespace iplmcfd


#endif // IPLMCFD_PATTERN_HPP_
