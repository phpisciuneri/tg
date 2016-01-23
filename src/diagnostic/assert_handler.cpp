#include <boost/format.hpp>
#include <stdexcept>

//#define PAUSE_AT_ERROR

#ifdef PAUSE_AT_ERROR
#include <iostream>
#endif

#include "tracer.hpp"


namespace boost
{
   /** \brief Global assertion handler for LMCFD. 
    *  
    *  Instead of total halt due to assertion failure, this function will
    *  convert the assertion error to some useful text and throw a
    *  std::runtime_error for some higher level exception handler to manage
    *  the error. 
    *
    *  For this to work, BOOST_ENABLE_ASSERT_HANDLER should be globally defined, 
    *  that is defined with a compiler switch -DBOOST_ENABLE_ASSERT_HANDLER or 
    *  in defs.hpp. Otherwise, this function will never be called and the effect 
    *  of assertions will be the same as standard assert in \<cassert\>.
    *  
    *  For optimized code -DBOOST_DISABLE_ASSERTS to optimize away the assertion
    *  checks. 		 				 
    */
   void assertion_failed(char const * expr, char const * function, char const * file, long line)
   {
      using namespace std;
      std::string error_str = (boost::format(
         "IPLMCFD FATAL ERROR! Assertion failure.\n"
         "  %1%\n"
         "  in file `%2%:%3%'\n"
         "  in function %4%\n"         
         ) % expr % file % line % function         
         ).str(); 
      ::iplmcfd::tracer::instance().note("!!! ASSERTION FAILURE !!!");
#ifdef PAUSE_AT_ERROR
      std::cout << error_str;
      bool b = true;
      while(b);
#else
      throw runtime_error(error_str);
#endif
   }


   void assertion_failed_msg(char const * expr, char const* msg, char const * function, char const * file, long line)
   {
      using namespace std;
      std::string error_str = (boost::format(
         "IPLMCFD FATAL ERROR! Assertion failure, with message: %5%\n"
         "  %1%\n"
         "  in file `%2%:%3%'\n"
         "  in function %4%\n"         
         ) % expr % file % line % function % msg
         ).str(); 
      ::iplmcfd::tracer::instance().note( "!!! ASSERTION FAILURE !!!" );
#ifdef PAUSE_AT_ERROR
      std::cout << error_str;
      bool b = true;
      while(b);
#else
      throw runtime_error(error_str);
#endif
   }


}
