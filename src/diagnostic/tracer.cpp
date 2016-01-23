#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <stack>

#include <boost/format.hpp>
#include <boost/progress.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp> 
#include <mpi.h>

#include "tracer.hpp"
#include "../simparam.hpp"

using namespace std;

namespace iplmcfd {


   //! No-op implementation (used in inactive ranks)
   struct tracer::impl
   {
      virtual void range(const string&)  {}
      virtual void endrange()  {}
      virtual void note(const std::string&)  {}
      virtual void flush_soft() {}
      virtual void flush() {}
      virtual void flush_on() {}
      virtual void flush_off() {}
      
      // stack ops
      virtual void push(const std::string&) {}
      virtual void pop() {}
   };


   // Interface to impl forwards
   void tracer::note(const std::string& msg)  { p->note(msg); }
   void tracer::range(const std::string& msg) { p->range(msg); }
   void tracer::endrange()   { p->endrange(); }
   void tracer::flush_soft() { p->flush_soft(); }
   void tracer::flush()      { p->flush(); }
   void tracer::flush_on()   { p->flush_on(); }
   void tracer::flush_off()  { p->flush_off(); }



/// Tracer implementation namespace (all details are local)
namespace {   

   // alternative implementations 
   struct tracer_impl; // will define, real definition (need a forward here)
   typedef tracer::impl tracer_impl_dummy; 
   
   //  +---------------+
   //  |  Definitions  | 
   //  +---------------+

typedef double time_t;
typedef double timediff_t;
const string START_TOKEN = "\\";
const string END_TOKEN   = "/";
const string FLOW_TOKEN  = "|";
const string TRACE_FILE_NAME = "trace.out";

//  +-------------+
//  |  Utilities  | 
//  +-------------+
time_t now()  { return MPI_Wtime(); }


//  +---------+
//  |  BLOCK  | 
//  +---------+

// timed block (stored in a stack)

struct block
{
   string msg;
   time_t begin_time;
   boost::scoped_ptr<block> active_range;
   tracer_impl* p;

   explicit block(const string& msg, tracer_impl* p);
   ~block();
};


//  +---------------+
//  |  TRACER_IMPL  | 
//  +---------------+

//! Real implementation of the tracer class
struct tracer_impl 
   : public tracer::impl
{

   //  +---------------+
   //  |  constructor  | 
   //  +---------------+
   tracer_impl( MPI_Comm comm )
   {   
      m_flush_on  = true;
      time_zero = now();
      blockcount = 0;

      // Open trace file 
      namespace bf = boost::filesystem;
      bf::path tracefile = simparam::instance().initial_path() / TRACE_FILE_NAME;
      log.open(tracefile.string().c_str(), std::ios_base::app );
      if (!log) 
         BOOST_THROW_EXCEPTION( std::runtime_error("cannot open " + TRACE_FILE_NAME + " for writing") );

      // Wait all to prevent writing before everyone opens it 
      // (This may not be necessary as everyone is opening it with ios::app)
      MPI_Barrier(comm);
      MPI_Comm_free(&comm); // not used beyond this point

      int rank,size;
      boost::tie(rank,size) = mpi_get_comm_rank_size(); // get ranks from WORLD, not the restricted tracer communicator `comm`. 

      // Compose time and rank formats (size depends on number of digits in comm_size)
      string ndigit    = boost::lexical_cast<string>(
         static_cast<int>( std::ceil(std::log10(float(size))) )
         );
      string formatstr = "%0" + ndigit + "d  ";
      rank_prefix = (boost::format(formatstr) % rank).str();
      //
      double wtick = MPI_Wtick();
      int time_digits = static_cast<int>( std::min(
         std::ceil(std::log10(1/MPI_Wtick())), double(6)
         ));

      time_format_long = boost::format(
         (boost::format( "%% %1%.%2%f") % (time_digits+7) % time_digits).str()
         );
      time_format_short = boost::format(
         (boost::format( "%% %1%.%2%f") % (time_digits+2) % time_digits).str()
         );
   }


   //  +---------+
   //  |  RANGE  | 
   //  +---------+
   void range( const string& msg )
   {
      if ( blockstack.empty() ) 
         blockstack.push( boost::make_shared<block>("NON_SCOPE: " + msg, this) );
      blockstack.top()->active_range.reset(0); // kill old
      blockstack.top()->active_range.reset(new block(msg,this)); // add this
   }
   void endrange()
   {
      if ( blockstack.empty() ) return;
      blockstack.top()->active_range.reset();
   }


   //  +---------+
   //  |  FLUSH  | 
   //  +---------+

   void flush()       {  log.flush(); }
   void flush_soft()  {  if ( m_flush_on ) flush(); }
   void flush_on()    {  m_flush_on = true; }
   void flush_off()   {  m_flush_on = false; }

   //  +------------+
   //  |  push/pop  | 
   //  +------------+

   void push(const std::string& msg) {   blockstack.push( boost::make_shared<block>(msg,this)); }
   void pop() { 
      if ( blockstack.empty() ) return;  
      blockstack.pop();  
   }
   
   //  +---------+
   //  |   ~()   |
   //  +---------+

   ~tracer_impl() {
      while( !blockstack.empty() ) blockstack.pop();
   }


   
   //  +---------+
   //  |  NOTE   | 
   //  +---------+
   void note(const std::string& msg)
   {
      log << prefix(now()) << flow() << "\t" << msg << "\n";
      flush_soft();
   }

   //  +---------------+
   //  |   START_MSG   |
   //  +---------------+
   void start_msg(const string& msg, time_t begin_time)
   {
      log << prefix(begin_time) << flow(START_TOKEN) << " "
         << msg 
         << "\n";
      flush_soft();
   }

   //  +-------------+
   //  |   END_MSG   |
   //  +-------------+
   void end_msg(const string& msg, time_t begin_time, time_t end_time)
   {
      log << prefix(end_time) << flow(END_TOKEN) << " "
         << msg << "\t" << timediff(begin_time,end_time)
         << "\n";
      flush_soft();
   }


   //  +----------+
   //  |   FLOW   |
   //  +----------+
   string flow (const std::string last_token = "") const 
   {
      string current_stack; 
      current_stack.reserve( blockcount*FLOW_TOKEN.size() );
      for(int i=0; i<blockcount; ++i)
         current_stack+=FLOW_TOKEN;
      return current_stack += last_token;
   }

   //  +------------+
   //  |   PREFIX   |
   //  +------------+
   string prefix(time_t t) 
   {
      return rank_prefix + (time_format_long % timediff(time_zero, t)).str() + "\t";
   }

   //  +--------------+
   //  |   TIMEDIFF   |
   //  +--------------+
   string timediff( time_t from, time_t to ) 
   {
      return (time_format_short % ( to - from )).str();
   }

      
   // data
   bool m_flush_on;
   std::ofstream log;
   string rank_prefix;
   time_t time_zero;
   boost::format time_format_long;
   boost::format time_format_short;
   std::stack<boost::shared_ptr<block> > blockstack;
   int blockcount; // all in stack plus all active ranges
}; 


//  +-------------------+
//  |  BLOCK Ctor/Dtor  | 
//  +-------------------+

block::block(const string& msg, tracer_impl* p)
   : msg(msg)
   , begin_time(now())
   , p(p)
{
   // do the output
   p->start_msg( msg, begin_time );
   p->blockcount++;
}  
block::~block() 
{
   // kill if any ranges blocks
   active_range.reset();
   p->blockcount--;
   // bye
   p->end_msg(msg, begin_time, now());
}


//  +----------------+
//  |   IS_ACTIVE()  |
//  +----------------+

//! Decides whether current rank is an active tracer
bool is_active()
{
   // Get rank/size 
   int rank, size;
   boost::tie(rank,size) = mpi_get_comm_rank_size();

   // Is my trace active or not
   const simparam& params = simparam::instance();
   // is my rank listed?
   bool my_rank_selected = params.trace_ranks.find(rank) != params.trace_ranks.end();
   // are the listed ranks to be turned off (if -1 is in the list), or on (if -1 is not in the list)?
   bool turn_on_selected = params.trace_ranks.find( -1 ) == params.trace_ranks.end();
   // on  and selected      -> active
   // off and not-selected  -> active
   // otherwise             -> not active
   // This is like invert xor (11 -> 1, 00 -> 1, o/w 0)
   return !(turn_on_selected ^ my_rank_selected);
}

}  // anon namespace



//  +------------------------+
//  |   TRACER Constructor   |
//  +------------------------+
tracer::tracer()
{
   bool active = is_active();
   int comm_color = active ? 1 : 2;
   MPI_Comm tracer_comm;
   // create a new communicator for active traces, because there will be a barrier 
   // in order sync access to the log file. 
   // key of 0 is passed to allow same order (but not the same rank ID) in the new
   // communicator. We will not use the new rank order anyway.
   MPI_Comm_split(MPI_COMM_WORLD, comm_color, 0, &tracer_comm);
   if( active )
      p.reset( new tracer_impl(tracer_comm) );
   else 
      p.reset( new tracer_impl_dummy );
}


//  +-----------------------+
//  |   TRACER Destructor   |
//  +-----------------------+

// this must be defined after impl is visible
tracer::~tracer() {}


//  +-------------------------------------+
//  |   REMOVE_PREVIOUS (static method)   |
//  +-------------------------------------+

void tracer::remove_previous()
{
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   if (rank==0)
   {
      namespace bf = boost::filesystem;
      bf::path tracefile = simparam::instance().initial_path() / TRACE_FILE_NAME;
      bf::remove(tracefile);
   }
   MPI_Barrier(MPI_COMM_WORLD);
}




//  +-----------------+
//  |  TRACER::SCOPE  | 
//  +-----------------+

tracer::scope::scope(const std::string& msg)
{
   instance().p->push(msg);
}

tracer::scope::~scope() try { // do not throw ever
   instance().p->pop();
} catch (...) {}

} // namespace iplmcfd

