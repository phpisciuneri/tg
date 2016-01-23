#ifndef SORTEDVECTOR_HPP_
#define SORTEDVECTOR_HPP_

#include "particle.hpp"
#include "bin.hpp"
#include <list>
#include <vector>
#include <boost/progress.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

class sortedvector
{
public:
   typedef std::vector<particle> particle_array;
   typedef particle_array::iterator iterator;
   typedef std::vector<size_t> count_array;

   //===========================================================================
   // CONSTRUCTOR   Create vector particles.
   //    Assigns each created particle to cells in blocks
   //===========================================================================
   sortedvector( size_t nParticles, size_t nCells )
      : particles( nParticles ), cell_np(nCells)
   {
      size_t np_percell = nParticles / nCells;
   
      iterator p = particles.begin();
      iterator pend = particles.end();
      
      for (size_t ip = 0; p!=pend; ++p, ++ip)
         p->inCell = ip/np_percell;
      
      redistribute();
   }




   //===========================================================================
   // REDISTRIBUTE   Based on particle's inCell, sort them
   //===========================================================================
   void redistribute()
   {
      namespace l = boost::lambda;
      std::sort(particles.begin(), particles.end(), 
         l::bind(&particle::inCell, l::_1) 
         <
         l::bind(&particle::inCell, l::_2) 
         );
      update_counts();
   }


   //===========================================================================
   // UPDATE_COUNTS   Given sorted sequence of particles, update cell particle 
   //                 counts.
   //===========================================================================
   void update_counts()
   {
      namespace l = boost::lambda;
      iterator p = particles.begin();
      iterator pend = particles.end();
      //
      std::fill(cell_np.begin(), cell_np.end(), 0);      
      for (; p!=pend; ++p) ++cell_np[p->inCell];
   }




   //===========================================================================
   // FOR_EACH_P   Loop all particles and execute a callback on them
   //===========================================================================
   template < class Func_ >  void for_each_p( Func_ f )
   {
      iterator p,pend;
      p = particles.begin();
      pend = particles.end();
      for ( ; p!=pend; ++p)  f(*p);
   }



   //===========================================================================
   // DATA MEMBERS -------------------------------------------------------------
public:
   particle_array particles;
   count_array    cell_np;
   //===========================================================================


};





#endif // SORTEDVECTOR_HPP_
