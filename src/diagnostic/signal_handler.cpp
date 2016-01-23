#ifdef LMCFD_ENABLE_SIGNAL_HANDLING

#include "signal_handler.hpp"
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "tracer.hpp"
#include <boost/lexical_cast.hpp>
#include <iostream>

namespace iplmcfd {


   volatile sig_atomic_t sigint_raised = 0;
   

   void 
   catch_sigint (int sig)
   {
      tracer::instance().note( "Caught signal " + boost::lexical_cast<std::string>(sig) );
      sigint_raised = 1;
      std::cout << "\ncaught signal\n" << std::endl;
   }
   
   bool 
   signal_handler::interrupted()
   {
      return sigint_raised == 1;
   }


   // register signal handler
   struct register_handler
   {
      register_handler()
      {
         signal( SIGTERM, catch_sigint );
         signal( SIGINT, catch_sigint );
      }      
   } instance;


} // namespace iplmcfd


#endif // LMCFD_ENABLE_SIGNAL_HANDLING
