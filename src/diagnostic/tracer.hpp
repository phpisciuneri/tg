#ifndef IPLMCFD_TRACER_HPP_
#define IPLMCFD_TRACER_HPP_

#include <string>
#include <boost/scoped_ptr.hpp>
//
// not using the following, the name turns out ugly
//
//#include <boost/current_function.hpp>
//! Automated function scope naming utility -- relies on BOOST_CURRENT_FUNCTION being available
// #define TRACER_SCOPE_CURRENT_FUNCTION_MSG( ExtraMsg ) \
//    tracer::scope tracer_scope_var_ ## __LINE__ ( \
//    ::std::string(ExtraMsg).empty() ? BOOST_CURRENT_FUNCTION \
//    : BOOST_CURRENT_FUNCTION + ::std::string(": ") + ExtraMsg )
// 
// #define TRACER_SCOPE_CURRENT_FUNCTION TRACER_SCOPE_CURRENT_FUNCTION_MSG("")

namespace iplmcfd {



//! Singleton action timer/logger class.
/** \ingroup singleton

Three ways of usage: 

1. Scope
2. Range
3. Note

Scope traces from construction to the end of the encapsulating scope. 
Range marks arbitrary block of code within a scope. 
Call to endrange() or a start of a new range terminates active range. 
End of scope also terminates any active range.

Note creates a notice at the current location. It doesn't do interval timing. 

\code{.cpp}
void foo () {
  tracer::scope _("foo");
    
  // ... do work 

  tracer::instance().range("some block");
  // .. more work

  tracer::instance().range("some other block"); // terminates previous range
  // .. more work

  tracer::instance().endrange(); // terminates active range (no-op if no range exists)


}  // terminates foo scope. also terminates any active ranges
\endcode

*/
class tracer
{

   // Singleton
private:
   tracer();
    ~tracer();

public:
   static tracer& instance() {
      static tracer singleton;
      return singleton;
   }
   
   //! Remove existing trace out file
   static void remove_previous();

   //! Add a note within the block, no interval timing
   void note(const std::string& msg);
   //! Start a new range (ends previous in the same scope if any)
   void range(const std::string& msg);
   //! End the last range
   void endrange();

   //! Issue flush only if flush is set 
   void flush_soft();
   //! Force flush regardless of setting
   void flush();
   //! Set flush on
   void flush_on();
   //! Set flush off
   void flush_off();

public:
   
   // RAII helper for scope start()/end() logs
   class scope {
   public:
      explicit scope(const std::string& msg);
      ~scope();
   };


   
public:
   struct impl; // allow derivation
private:
   boost::scoped_ptr<impl> p;
};


} // namespace iplmcfd



#endif // IPLMCFD_TRACER_HPP_
