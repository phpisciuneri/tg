#ifndef IPLMCFD_FIELD_HPP_
#define IPLMCFD_FIELD_HPP_

#include <vector>

#include "defs.hpp"

namespace iplmcfd {

  //!
  //! \file field.hpp
  //! \class Field
  //! \brief Registers and stores fields to be output
  //!
  class Field
  {
  public:

    //!
    //! \brief
    //!
    //! \todo handle null address
    //!
    void add( const std::vector< real_t >& f, const std::string& name )
    {
      if ( f.size() > 0 )
        m_data.push_back( &f[0] );
      else
        m_data.push_back( NULL );
      m_name.push_back( name );
    }

    //!
    //!
    //!
    real_t data( lid_t I, std::size_t field_index ) const 
    { 
      return *( m_data[field_index] + I );
    }

    //!
    //!
    //!
    const std::vector< std::string >& names() const { return m_name; }

    //!
    //!
    //!
    std::size_t size() const { return m_name.size(); }

  private:

    //!
    std::vector< const real_t* > m_data;
    //!
    std::vector< std::string >   m_name;

  };

} // namespace iplmcfd

#endif // IPLMCFD_FIELD_HPP_