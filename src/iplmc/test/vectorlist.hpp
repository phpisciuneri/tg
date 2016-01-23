#ifndef VECTORLIST_HPP_
#define VECTORLIST_HPP_

#include "particle.hpp"
#include "bin.hpp"
#include <list>
#include <vector>
#include <boost/progress.hpp>


class vectorlist
{
public:
   typedef bin< particle >       particle_bin;
   typedef std::vector< particle_bin >  bin_vector;

   //===========================================================================
   // CONSTRUCTOR   Create vector of lists of particles.
   //    Assigns each created particle its initial cell number. 
   //===========================================================================
   vectorlist( size_t nParticles, size_t nCells )
      : cells(nCells)
   {
      size_t npg = nParticles / nCells;
      size_t iCell = 0;
      bin_vector::iterator v    = cells.begin();
      bin_vector::iterator vend = cells.end();
      
      for (; v != vend; ++v, ++iCell)
      {
         v->contained.resize(npg);
         particle_bin::iterator p = v->contained.begin();
         particle_bin::iterator pend = v->contained.end();
         for (; p!=pend; ++p) p->inCell = iCell;
      }
      redistribute();
   }



   //===========================================================================
   // REDISTRIBUTE   Based on particle's inCell, put them in correct cell.
   //===========================================================================
   void redistribute()
   {
      size_t iCell = 0;
      bin_vector::iterator v = cells.begin();
      bin_vector::iterator vend = cells.end();
      boost::progress_display pg(cells.size());
      for ( ; v != vend; ++v, ++iCell, ++pg)
      {
         particle_bin::iterator p    = v->contained.begin();
         particle_bin::iterator pEnd = v->contained.end();
         for (; p!=pEnd; )
         {
            if (p->inCell == iCell) 
               ++p;
            else {
               assert(p->inCell<cells.size());
               transfer(p++, *v, cells[p->inCell] );
            }
         }
      }
      // complete transfer
      for (v = cells.begin(); v!=vend; ++v)  v->own();
   }
   



   //===========================================================================
   // FOR_EACH_P   Loop all particles and execute a callback on them
   //===========================================================================
   template < class Func_ >
   void for_each_p( Func_ f )
   {
      bin_vector::iterator v    = cells.begin();
      bin_vector::iterator vend = cells.end();
      for ( ; v != vend; ++v )
      {
         particle_bin::iterator p    = v->contained.begin();
         particle_bin::iterator pEnd = v->contained.end();
         for (; p!=pEnd; ++p)
            f ( *p );
      }
   }



   //===========================================================================
   // DATA MEMBERS -------------------------------------------------------------
public:
   bin_vector cells;
   //===========================================================================


};




#endif //VECTORLIST_HPP_
