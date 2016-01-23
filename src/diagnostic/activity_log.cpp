#include "activity_log.hpp"

#include <algorithm>
#include <boost/assign/list_of.hpp>
#include <boost/throw_exception.hpp>
#include <boost/filesystem/operations.hpp>

#include "../simparam.hpp"


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
namespace iplmcfd {
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



//==============================================================================
// CONSTRUCTOR
//==============================================================================
activity_log::activity_log()
{
   using namespace std;
   const simparam* prms = &simparam::instance();

   //
   // parameters
   //
   m_flush = prms->carlo_log_flush;
   set_log_level( prms->carlo_log_level );


   //
   // Open log files.
   // Truncate if new run, append if restarting
   //
   if ( m_loglevel == LOG_ALL )
   {
      m_logfile.open("mcdetail.log", prms->restart ? ios_base::app : ios_base::trunc );
      BOOST_ASSERT( m_logfile );
   }
   bool appending = false;
   if ( m_loglevel > LOG_NOT)
   {
      appending = boost::filesystem::exists("mcsummary.log") && prms->restart;
      m_sumfile.open("mcsummary.log", prms->restart ? ios_base::app : ios_base::trunc );
      BOOST_ASSERT( m_sumfile );
   }
   

   //
   // Compose formats
   //
   m_frogue.parse   ("[ROGUE    ] Cell %1%, Shift (%2%,%3%,%4%)\n" );
   m_fempty.parse   ("[EMPTY_C  ] Cell %1%\n");
   m_fbullied.parse ("[BULLIED_C] Cell %1%, Bullied %2%, Cloned %3%\n");
   m_fcrowded.parse ("[CROWDED_C] Cell %1%, Merged  %2%\n");
   m_funderpop.parse("[UNDERPP_C] Cell %1%, Cloned  %2%\n");

   // add header to summary file if not restarting
   if ( !appending && m_sumfile ) 
   {
      m_sumfile 
         << std::setw(10) << "iter"
         << std::setw(10) << "total"
         << std::setw(10) << "erased"
         << std::setw(10) << "created"
         << std::setw(10) << "incame"
         << std::setw(10) << "outgone"
         << std::setw(10) << "rogue"
         << std::setw(10) << "shifted"
         << std::setw(10) << "migrated"
         << std::setw(10) << "empty"
         << std::setw(10) << "bullied"
         << std::setw(10) << "underpop"
         << std::setw(10) << "crowded"
         << std::endl ;
   }
   // 
   // Reset iteration counter
   //
   m_iter = size_t(-1);
}



//==============================================================================
//==============================================================================
void activity_log::begin_iteration(size_t iter)
{
   if (m_iter == iter ) return;
   if (m_iter != size_t(-1) ) end_iteration();
   
   m_iter = iter;

   std::fill(m_counter.begin(), m_counter.end(), 0);
   if (m_logfile) 
      m_logfile << "*** " << m_iter << " ***\n";

}



void activity_log::end_iteration()
{
   if (m_logfile && m_flush ) m_logfile << std::flush;
   if (m_sumfile) 
   {
      m_sumfile << std::setw(10) << m_iter;
      for (int i=0; i!=N_ACTIVITY; ++i)
         m_sumfile << std::setw(10) << m_counter[i];
      m_sumfile << "\n";
      if ( m_flush ) m_sumfile << std::flush;
   }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
} // namespace iplmcfd
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
