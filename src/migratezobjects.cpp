#include "migratezobjects.hpp"

#include <set>
#include <boost/foreach.hpp>

#include "cell.hpp"
#include "particle.hpp"
#include "diagnostic/tracer.hpp"

namespace iplmcfd
{

  // +----------------------------------+
  // | MigrateZObjects::MigrateZObjects |
  // +----------------------------------+
  MigrateZObjects::MigrateZObjects( Subdomain& s, Iplmc& mc ) 
    : ZObjects( s ), m_subdomain( s ), m_mc( mc ) 
  {
    tracer::scope _( "MigrateZObjects::MigrateZObjects" );

    // create sample particle to get the dynamic size of phi
    Particle p( *m_mc.m_chem );
    m_fix_bytes = p.fixed_num_bytes();
    m_phi_bytes = p.phi.size() * sizeof( real_t );
  
    // create distributed directory
    m_dd.reset( new Zoltan_DD( 
      MPI_COMM_WORLD, // communicator
      1, // num_gid_entries
      0, // num_lid_entries (ignore them)
      0, // length of user defined data
      0, // table length (0=default)
      0 // debug level
      ) );

    // update directory
    std::vector< ZOLTAN_ID_TYPE > gid( s.nhomes() );
      
    BOOST_FOREACH( Subdomain::gidlid_mp::const_reference gidlid, s.homes() )
      gid[ gidlid.second ] = ZOLTAN_ID_TYPE( gidlid.first );


    m_dd->Update( gid.data(), NULL, NULL, NULL, s.nhomes() );

  }

  // +---------------------------------+
  // | MigrateZObjects::obj_size_multi |
  // +---------------------------------+
  void MigrateZObjects::obj_size_multi( void* data,
    int num_gid_entries,
    int num_lid_entries,
    int num_ids,
    ZOLTAN_ID_PTR global_ids,
    ZOLTAN_ID_PTR local_ids,
    int* sizes,
    int* ierr )
  {
    tracer::scope _( "MigrateZObjects::obj_size_multi" );

    MigrateZObjects *ob = static_cast< MigrateZObjects* >( data );
    
    /* Here we describe the data, in terms of size that needs to migrate.
     * The migratable objects are portions of the local graph corresponding to
     * global IDs that are migrating, and the list of particles belonging to 
     * the corresponding cell.
     */

    Subdomain::adj_mp::const_iterator it_adj;
    Iplmc::cell_mp::const_iterator it_cell;
    for ( int i=0; i<num_ids; i++ )
    {
      
      it_adj = ob->m_subdomain.adj().find( global_ids[i] );
      assert( it_adj != ob->m_subdomain.adj().end() );

      it_cell = ob->m_mc.m_cells.find( global_ids[i] );
      assert( it_cell != ob->m_mc.m_cells.end() );

      *sizes++ = //assuming vertex local id equals i

        // 1. ADJACENCY LIST

        // buffer for the number of neighbors
        sizeof( std::size_t ) +                 
        // buffer for the neighbor list
        sizeof( int ) * it_adj->second.size() +
      
        // 2. PARTICLE LIST

        // buffer for number of particles in contained list
        sizeof( std::size_t ) +                 
        // buffer for static data of all particles
        ob->m_fix_bytes * it_cell->second.contained().size() +
        // buffer for dynamic data of all particles
        ob->m_phi_bytes * it_cell->second.contained().size();
    }

    *ierr = ZOLTAN_OK;

  }

  // +---------------------------------+
  // | MigrateZObjects::pack_obj_multi |
  // +---------------------------------+
  void MigrateZObjects::pack_obj_multi( void* data,
    int num_gid_entries,
    int num_lid_entries,
    int num_ids,
    ZOLTAN_ID_PTR global_ids,
    ZOLTAN_ID_PTR local_ids,
    int* dest,
    int* sizes,
    int* idx,
    char* buf,
    int* ierr )
  {
    tracer::scope _( "MigrateZObjects::pack_obj_multi" );

    MigrateZObjects *ob = static_cast< MigrateZObjects* >( data );

    Subdomain::adj_mp::iterator it_adj;
    Iplmc::cell_mp::iterator it_cell;

    // move data into buffer and remove from current owner
    for ( int i=0; i<num_ids; i++ )
    {
      // find objects
      it_adj = ob->m_subdomain._adj.find( global_ids[i] );
      assert( it_adj != ob->m_subdomain._adj.end() );

      it_cell = ob->m_mc.m_cells.find( global_ids[i] );
      assert( it_cell != ob->m_mc.m_cells.end() );

      // 1. ADJACENCY LIST

      // copy number of neighbors to buffer
      std::size_t* size_pbuf = reinterpret_cast< std::size_t* >( buf + idx[i] );
      *size_pbuf++ = it_adj->second.size();
            
      // copy neighbor list to buffer
      int* int_pbuf = reinterpret_cast< int* >( size_pbuf );
      for ( std::size_t j=0; j<it_adj->second.size(); j++ ) 
        *int_pbuf++ = it_adj->second[j];

      // 2. PARTICLE LIST

      // copy number of particles in contained list to buffer
      size_pbuf = reinterpret_cast< std::size_t* >( int_pbuf );
      *size_pbuf++ = it_cell->second.contained().size();
      
      // for all contained particles
      real_t* real_pbuf = reinterpret_cast< real_t* >( size_pbuf );
      BOOST_FOREACH( const Particle& p, it_cell->second.contained() )
      { 
        // copy static data of a particle to buffer
        *real_pbuf++ = p.mass;
        for ( std::size_t k=0; k<NDIM; k++ ) *real_pbuf++ = p.X[k];
        *real_pbuf++ = p.f;
        *real_pbuf++ = p.srt;

        // copy dynamic data of a particle to buffer
        for ( std::size_t k=0; k<p.phi.size(); k++ )
          *real_pbuf++ = p.phi[k];
      
      }

      // erase objects
      ob->m_subdomain._adj.erase( it_adj );
      ob->m_mc.m_cells.erase( it_cell );
    }

    *ierr = ZOLTAN_OK;

  }

  // +-----------------------------------+
  // | MigrateZObjects::unpack_obj_multi |
  // +-----------------------------------+
  void MigrateZObjects::unpack_obj_multi( void* data,
    int num_gid_entries,
    int num_ids,
    ZOLTAN_ID_PTR global_ids,
    int* sizes,
    int* idx,
    char* buf,
    int* ierr )
  {
    tracer::scope _( "MigrateZObjects::unpack_obj_multi" );

    MigrateZObjects *ob = static_cast< MigrateZObjects* >( data );

    for ( int i=0; i<num_ids; i++ )
    {
      // 1. ADJACENCY LIST

      // copy number of neighbors from buffer
      std::size_t* size_pbuf = reinterpret_cast< std::size_t* >( buf + idx[i] );
      std::vector< int > adj( *size_pbuf++, MAX_SIZE );

      // copy neighbor list from buffer
      int* int_pbuf = reinterpret_cast< int* >( size_pbuf );
      for ( std::size_t j=0; j<adj.size(); j++ ) 
        adj[j] = *int_pbuf++;

      // insert into graph
      ob->m_subdomain._adj.insert( 
        Subdomain::adj_mp::value_type( gid_t( global_ids[i] ), adj )
        );

      // 2. PARTICLE LIST
      
      // copy number of particles in contained list to buffer
      size_pbuf = reinterpret_cast< std::size_t* >( int_pbuf );
      std::size_t nparticles = *size_pbuf++;
      
      //std::cout << "nparticles: " << nparticles << std::endl;
      
      // create new cell
      Cell cell( *ob->m_mc.m_chem );
      
      // for all contained particles
      real_t* real_pbuf = reinterpret_cast< real_t* >( size_pbuf );
      for ( std::size_t j=0; j<nparticles; j++ )
      { 
        Particle p( *ob->m_mc.m_chem );
      
        // copy static data of a particle to buffer
        p.mass = *real_pbuf++;
        for ( std::size_t k=0; k<NDIM; k++ ) p.X[k] = *real_pbuf++;  
        p.f   = *real_pbuf++;
        p.srt = *real_pbuf++;
      
        // copy dynamic data of a particle to buffer
        for ( std::size_t k=0; k<p.phi.size(); k++ ) p.phi[k] = *real_pbuf++;
      
        // add particle to contained list of cell
        cell.contained().push_back( p );
      
      }
      
      // add cell to cell map
      ob->m_mc.m_cells.insert( 
        Iplmc::cell_mp::value_type( gid_t( global_ids[i] ), cell ) 
        );

    }

    *ierr = ZOLTAN_OK;

  }

  // +----------------------------------+
  // | MigrateZObjects::post_migrate_pp |
  // +----------------------------------+
  void MigrateZObjects::post_migrate_pp( void* data,
    int num_gid_entries,
    int num_lid_entries,
    int num_import,
    ZOLTAN_ID_PTR import_global_ids,
    ZOLTAN_ID_PTR import_local_ids,
    int* import_procs,
    int* import_to_part,
    int num_export,
    ZOLTAN_ID_PTR export_global_ids,
    ZOLTAN_ID_PTR export_local_ids,
    int* export_procs,
    int* export_to_part,
    int* ierr )
  {
    tracer::scope _( "MigrateZObjects::post_migrate_pp" );

    MigrateZObjects *ob = static_cast< MigrateZObjects* >( data );

    // remove exports from subdomain homes
    for ( int i=0; i<num_export; i++ )
      ob->m_subdomain._homes.erase( gid_t( export_global_ids[i] ) );

    // add imports to subdomain homes (don't know local id yet, use 0)
    for ( int i=0; i<num_import; i++ )
      ob->m_subdomain._homes[ gid_t( import_global_ids[i] ) ] = 0;

    // renumber local indices from 0 and fill update buffer
    lid_t I = 0;
    BOOST_FOREACH( Subdomain::gidlid_mp::reference gidlid, 
      ob->m_subdomain._homes )
      gidlid.second = I++;


    // update info in DD for new imports
    ob->m_dd->Update( import_global_ids, NULL, NULL, NULL, num_import );

  }

  // +-----------------------------------+
  // | MigrateZObjects::update_neighbors |
  // +-----------------------------------+
  void MigrateZObjects::update_neighbors()
  {
    tracer::scope _( "MigrateZObjects::update_neighbors" );

    /* At this point the adjacency graph and the distributed directory have
     * been updated.  Thus the neighbors, home_legs and neigh_legs can be 
     * reconstructed at this point. Let's rebuild them completely from scratch.
     * In the future I can explore some more efficient ways.
     */

    std::set< ZOLTAN_ID_TYPE > vals2find;
    BOOST_FOREACH( Subdomain::adj_mp::const_reference gid_vec,
      m_subdomain.adj() )
    {
      BOOST_FOREACH( int i, gid_vec.second )
        vals2find.insert( i );
    }

    // copy set to vector
    std::vector< ZOLTAN_ID_TYPE > adj_values( vals2find.size() );
    std::size_t j=0;
    BOOST_FOREACH( ZOLTAN_ID_TYPE i, vals2find )
      adj_values[j++] = i;

    std::vector< int > owners( vals2find.size() );

    int ierr = m_dd->Find( 
        adj_values.data(), NULL, NULL, NULL, vals2find.size(), owners.data() );
    if ( ierr != ZOLTAN_OK )
    {
      std::cout << "ierr: " << ierr << std::endl;
      throw std::runtime_error( "Zoltan DD Find operation failed." );
    }

    // rebuild neighbors
    m_subdomain._neighs.clear();
    m_mc.m_neigh_cells.clear();
    for ( std::size_t i=0; i<adj_values.size(); i++ )
    {
      if ( owners[i] != m_subdomain.my_rank() )
      {
        gid_t Ineigh = gid_t( adj_values[i] );
        // must be a neighbor
        m_subdomain._neighs.insert( 
          Subdomain::gidrank_mp::value_type( Ineigh, rank_t( owners[i] ) ) );

        // create the neighbor ghost cell
        Cell cell( *m_mc.m_chem );
        m_mc.m_neigh_cells.insert( 
          Iplmc::cell_mp::value_type( Ineigh, cell ) );
      }
    }

    // rebuild home and neighbor legs
    m_subdomain._padj.clear();
    m_subdomain._home_legs.clear();
    m_subdomain._neigh_legs.clear();
    
    BOOST_FOREACH( Subdomain::adj_mp::const_reference gid_vec,
      m_subdomain.adj() )
    {
      gid_t I = gid_vec.first;
      const std::vector< int >& adj = gid_vec.second;

      for ( std::size_t i=0; i<adj.size(); i++ )
      {
        gid_t Ineigh = adj[i];

        Subdomain::gidrank_mp::const_iterator it = 
          m_subdomain.neighs().find( Ineigh );
        
        if ( it != m_subdomain.neighs().end() )
        {
          // its an actual neighbor
          m_subdomain._home_legs[ it->second ].insert( I );
          m_subdomain._neigh_legs[ it->second ].insert( Ineigh );
          m_subdomain._padj[I].push_back( it->second );
        }
        else
          // neighbor rank is just my_rank
          m_subdomain._padj[I].push_back( m_subdomain.my_rank() );

      }
      
    }

  }


  // +----------------------------+
  // | MigrateZObjects::reset_p2p |
  // +----------------------------+
  void MigrateZObjects::reset_p2p()
  {
    tracer::scope _( "MigrateZObjects::reset_p2p" );

    // now that neighboring info has been rebuilt, particle_p2p needs reset
    m_mc.m_p2p.reset( 
      new Particle_p2p( m_subdomain, 
      m_mc.m_cells, 
      m_mc.m_neigh_cells, 
      *m_mc.m_chem ) 
      );


  }

} // namespace iplmcfd
