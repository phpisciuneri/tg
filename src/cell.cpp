#include "cell.hpp"

#include <vector>
#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "simparam.hpp"
#include "diagnostic/activity_log.hpp"

namespace iplmcfd {

  // +-------------------------+
  // | Cell::adjust_population |
  // +-------------------------+
  void Cell::adjust_population(
    const size_t min_ppc,
    const size_t mean_ppc,
    const size_t max_ppc
    )
  {
    activity_log& log = activity_log::instance();

    // Get particle count
    size_t ppc = contained().size();
    enum {EMPTY, TOOMANY, TOOLITTLE } pop;
    if      ( ppc == 0     ) pop = EMPTY;
    else if ( ppc > max_ppc) pop = TOOMANY;
    else if ( ppc < min_ppc) pop = TOOLITTLE;
    else return; // healthy population; good work!

    // Empty cell is logged, nothing else can be done
    if ( pop == EMPTY ) 
    { 
      log.empty_cell(1); 
      return;
    }

    // Merge (cluster)
    if ( pop == TOOMANY )
    {
      // Merge the lighter ones.  
      namespace bl = boost::lambda;
      size_t excess =  contained().size() - mean_ppc;
      ctr.merged += excess;
      log.crowded_cell(excess);
      particle_lst temp;
      while (excess--)
      {
        // Find the minimum mass particle and move it into temporary list
        temp.splice( temp.end(), contained(), 
          std::min_element( contained().begin(), contained().end(), 
          bl::bind(&Particle::mass, bl::_1) < bl::bind(&Particle::mass, bl::_2) 
          ) );
        Particle& pmin = temp.back();
        // find the next minimum
        Particle& p = 
          *std::min_element( contained().begin(), contained().end(), 
          bl::bind(&Particle::mass, bl::_1) < bl::bind(&Particle::mass, bl::_2) 
          );
        // add these two
        real_t mass_merged = p.mass + pmin.mass;
#define merge(x) p.x = (p.x*p.mass + pmin.x*pmin.mass)/mass_merged
        for(int i=0; i!=p.phi.size(); ++i) merge(phi[i]);
        merge(X);
        merge(f);
#undef merge
        p.mass  = mass_merged;      
      }
      // get rid of merged particles
      temp.clear();
    }
    
    // CLONE
    else if ( pop == TOOLITTLE )
    {
      namespace bl = boost::lambda;
      //
      // get total mass
      real_t mass = 0;
      BOOST_FOREACH( const Particle& p, contained() )
        mass += p.mass;
      real_t mean_mass = mass / mean_ppc;

      // Divide heaviest ones
      int size_before = contained().size();
      while (contained().size() < mean_ppc)
      {
        // Find the maximum mass particle 
        Particle& p =         
          *std::max_element( contained().begin(), contained().end(), 
          bl::bind(&Particle::mass, bl::_1) < bl::bind(&Particle::mass, bl::_2) 
          );
        //
        // Divide its mass until less than mean_mass
        int div = static_cast<int>(std::ceil(p.mass / mean_mass));
        BOOST_ASSERT( div >= 2 );
        p.mass /= div;
        // Add div-1 of the same values and everything (+1 itself = div)
        contained().resize( contained().size()+div-1, p );

      }
      ctr.cloned += contained().size() - size_before;
      log.underpop_cell(contained().size() - size_before);
    }
  }

  // +------------------------------+
  // | Cell::append_to_message_real |
  // +------------------------------+
  void Cell::append_to_message_real(
    std::vector<MPI_Aint>& displ,
    std::vector<int>& count
    )
  {
    scalars.append_to_message_real( displ, count );
  }

  // +--------------------------------+
  // | Cell::append_to_message_counts |
  // +--------------------------------+
  void Cell::append_to_message_counts(
    std::vector<MPI_Aint>& chunk_addresses,
    std::vector<int>&      chunk_lengths
    )
  {
    MPI_Aint p;
    MPI_Address( &ctr, &p );
    chunk_addresses.push_back( p );
    chunk_lengths.push_back( Counters::NCOUNTS );
  }

} // namespace iplmcfd