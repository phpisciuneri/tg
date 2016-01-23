#include "version.hpp"
#include <boost/static_assert.hpp>
#include <boost/preprocessor/stringize.hpp>

/*
//! \file 
//! Keep revision information, and needs to be made everytime, always assumed to be out of date.

#ifndef LMCFD_REVISION
#warning No revision information available. 
#define LMCFD_REVISION 0
#endif

namespace iplmcfd {

   version::version()
   {
      const char r[] = BOOST_PP_STRINGIZE(LMCFD_REVISION);
      BOOST_STATIC_ASSERT( sizeof(r)<=LEN ); // in case LMCFD_REVISION macro is set to something senselessly long
      std::copy( r, r+sizeof(r), m_str );

      m_num = LMCFD_REVISION;
   }

} //  namespace iplmcfd
*/