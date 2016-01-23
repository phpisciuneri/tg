#include "../iplmc.hpp"

#include <cmath>
#include <boost/foreach.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_01.hpp>

#include "../diagnostic/cpublock.hpp"
#include "../diagnostic/activity_log.hpp"
#include "../diagnostic/tracer.hpp"

namespace iplmcfd {

  // +----------------------+
  // | Iplmc::position_step |
  // +----------------------+
  void Iplmc::position_step( real_t dt, std::size_t iter )
  {
    tracer::scope _( "Iplmc::position_step" );

    // Create random generators
    typedef boost::rand48 rngengine_t;
    typedef boost::normal_distribution< real_t > pdf_t;
    typedef boost::variate_generator< rngengine_t, pdf_t > normal_gen_t;
    typedef boost::uniform_01< rngengine_t, real_t > uniform_gen_t;
    static normal_gen_t gaussian(
      rngengine_t( rngengine_t::result_type( m_subdomain.my_rank() ) ), 
      pdf_t(0.,1.)
      );
    static uniform_gen_t uniform_rand(
      rngengine_t( rngengine_t::result_type( m_subdomain.my_rank() ) )
      );

    //
    real_t dtsqrt = std::sqrt( dt );
    real_v random_vector;
    real_v X0;
    int_v  shift;

    activity_log& log = activity_log::instance();
    log.begin_iteration( iter );

    BOOST_FOREACH( cell_mp::reference gidcell, m_cells )
    {
      Cell& cell = gidcell.second;

      log.attach_cell( gidcell.first );

      cell.ctr.reset();

      index_v cell_coordinates = m_domain.gid2ind( gidcell.first );
      cpublock cpub( cell.ctr.cputime, cell.contained().size() );

      // *** INTEGRATE POSITION ***
      // the loop is mutating the list, so do not use foreach
      Cell::particle_it p    = cell.contained().begin();
      Cell::particle_it pEnd = cell.contained().end();
      while ( p != pEnd )
      {
        log.attach_particle( *p );

        /* SDE coefficients are exact to the particle location and 
         * calculated as follows:
         *
         * Dx = u^+
         * Dy = v^+
         * Dz = 0
         * E  = sqrt( 2 * nu / Sc ) (nu = cnst)
         *
         */
        real_v dd( 0 );
        m_tgv->uvw( m_time, m_domain.xyz( gidcell.first ), dd );
        real_t ee = std::sqrt( 2 * m_nuref / m_simparam.Sc ); 

        // integrate position
        random_vector = gaussian(), gaussian(), 0;
        real_v drift, diffuse;
        drift   = dd * dt;
        diffuse = ee * dtsqrt * random_vector;
        // temporarily store current position (used if particle is rogue)
        X0 = p->X;
        p->X += m_domain.h()*( drift + diffuse );

        //
        // check new location
        shift = floor( p->X );

        //
        // if shift is all zero then we are in the same cell. relax and continue
        if ( all_elements(shift == 0) ) {
          ++p;
          continue;
        }

        // Before anything else see if this is a ROGUE particle. If a particle
        // moves further than dimensions of a single cell it is called ROGUE.
        // This may need to be dealt with properly if it crosses to another
        // subdomain. Also, too many rogues might be an indication of too long
        // timestep or excess diffusion or something. So, let's log it.
        bool is_rogue = any_elements( abs(shift) >= 2 );
        if ( is_rogue )
          log.rogue_particle( shift );
        else // the particle moved to one the adjacent cells
          ++(cell.ctr.wentout);

        //  Check if still inside the domain
        // now let's see where we landed 
        // (this handles periodicity inside domain ind2gid conversions)
        int_v new_sub( cell_coordinates + shift );
        gid_t new_gid = m_domain.ind2gid( new_sub );

        // everything is periodic, nothing should have left
        assert( !m_domain.out_of_bounds( new_gid ) );

        if ( m_domain.out_of_bounds( new_gid ) )
        {
          // no reflection of any sort. just erase the damn thing
          cell.contained().erase(p++);
          ++(cell.ctr.wentout);
          log.outgone_particle();
        }
        else
        {
          // in bounds, could be as follows
          // - lands on home
          // - lands on one the neighbor legs
          // - iff rogue, lands on beyond neighbor legs, which we can't handle but just log it
          // anything else is a programming error
          enum MoveCase { INHOME, INNEIGH, NOTLOCAL_ROGUE } movecase;

          cell_mp::iterator it;
          if ( ( it = m_cells.find( new_gid ) ) != m_cells.end() )
            movecase = INHOME;
          else if ( ( it = m_neigh_cells.find( new_gid ) ) != m_neigh_cells.end() )
            movecase = INNEIGH;
          else if (is_rogue)
            movecase = NOTLOCAL_ROGUE;
          else
            BOOST_ASSERT_MSG( false, "unknown MoveCase, something bad happened, come fix" );

          if ( movecase == INHOME )
          {
            Cell& new_cell = it->second;
            // offset particle position
            p->X -= shift;
            // put the particle in its new cell
            transfer( p++, cell, new_cell );
            log.shifted_particle();
          }
          else if ( movecase == INNEIGH )
          {
            Cell& new_cell = it->second;
            // offset particle position
            p->X -= shift;
            // put the particle in its new cell
            transfer( p++, cell, new_cell );
            log.migrated_particle();
          }
          else
          {
            // rogue
            // just put it back for now. If there are such particles too many
            // we'll deal with it
            p->X = X0;
            ++p;
          }

        }

      } // while

    } // foreach

    // message pass particles traveling to neighbor subdomain
    tracer::instance().range( "migrate");
    m_p2p->start_messaging();
    m_p2p->complete_messaging();

    // complete the move locally
    BOOST_FOREACH( cell_mp::reference gidcell, m_cells )
      gidcell.second.own_incoming();

    // ADJUST POPULATION
    if ( m_simparam.use_clone_cluster )
    {
      tracer::instance().range( "adjust_population" );
      BOOST_FOREACH( cell_mp::reference gidcell, m_cells )
      {
        Cell& cell = gidcell.second;
        log.attach_cell( gidcell.first );
        cell.adjust_population( m_min_ppc, m_mean_ppc, m_max_ppc );
      }

    }

  } // iplmc::position_step

} // namespace iplmcfd
