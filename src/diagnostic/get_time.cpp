#include "get_time.hpp"

#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/date_time/local_time/local_time.hpp>


namespace iplmcfd {

   std::string get_current_time() 
   {
      // TODO Check this in new version of boost

      // fails in intel 11 with boost 1_3x
#if __INTEL_COMPILER == 1100
      return "unknown";
#else
      std::stringstream time;
      time << boost::posix_time::ptime(boost::posix_time::second_clock::local_time());
      return time.str();
#endif
   }

} // namespace iplmcfd
