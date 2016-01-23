#ifndef RANDOMSHUFFLER_HPP_
#define RANDOMSHUFFLER_HPP_


#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/utility.hpp>




class randomshuffler 
   : boost::noncopyable
{
public:
   //===========================================================================
   // TYPEDEFS -----------------------------------------------------------------
   typedef boost::uniform_int<size_t> uniform_dist;
   typedef boost::rand48 rng_engine;
   typedef void result_type;

   //===========================================================================
   // CONSTRUCTOR 
   //===========================================================================
   randomshuffler(size_t nCell)
      : rand( rng_engine(), uniform_dist(0,nCell-1) ) 
   {
      reset();
   }


   //===========================================================================
   // RESET    Reset counters (not rng)  
   //===========================================================================
   void reset() { total_moved = 0; }   



   //===========================================================================
   // OPERATOR()  Assign a random cell number
   //===========================================================================
   void operator()(particle& p)
   {
      size_t old = p.inCell;
      p.inCell = rand();
      if ( old!=p.inCell ) ++total_moved;
   }




   //===========================================================================
   // TOTAL_MOVES   Return number particles whose inCell changed.
   //===========================================================================
   size_t total_moves() const { return total_moved; }
   
   
   
   //===========================================================================
   // DATA MEMBERS -------------------------------------------------------------
private:
   size_t total_moved;
   boost::variate_generator< rng_engine, uniform_dist > rand;      
   //===========================================================================



};




#endif //RANDOMSHUFFLER_HPP_
