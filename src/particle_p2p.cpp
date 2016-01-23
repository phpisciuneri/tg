#include "particle_p2p.hpp"

#include <cassert>
#include <boost/format.hpp>
#include <boost/foreach.hpp>

#include "diagnostic/tracer.hpp"

namespace iplmcfd {

  // +----------------------------+
  // | Particle_p2p::Particle_p2p |
  // +----------------------------+
  Particle_p2p::Particle_p2p( 
    const Subdomain& s, 
    gidcell_mp& cells,
    gidcell_mp& neigh_cells,
    const mix_t& mix )
    : m_subdomain( s ), m_part_mpi_type( mix ), m_sample_particle( mix )
  {
    tracer::scope _( "Particle_p2p::Particle_p2p" );
    
    // store communicator for further use
    //MPI_Comm_dup( MPI_COMM_WORLD, &m_comm );
    m_comm = MPI_COMM_WORLD;

    // set tags (must be unique within the communicator)
    m_tags[PCOUNT_TAG] = 100;
    m_tags[PDATA_TAG]  = 101;

    // Create particle count send/receive requests
    create_pcount_requests();

    // store ranks of neighboring subdomain in std::vector
    BOOST_FOREACH( Subdomain::rank_gid_s_mp::const_reference rank_gidset,
      m_subdomain.home_legs() )
      m_neighbors_compv.push_back( rank_gidset.first );

    m_num_neigh = m_neighbors_compv.size();

    // Create incoming particle pointer lists 
    s_plist_p_cellv_compv.resize( m_num_neigh );
    r_plist_p_cellv_compv.resize( m_num_neigh );

    // Loop over all neighboring compute nodes and get bridged cells' incoming
    // lists' addresses. We will need this a lot later on, so for efficiency and
    // clarity let's take care of this once and for all.

    std::size_t i = 0;
    BOOST_FOREACH( Subdomain::rank_gid_s_mp::const_reference rank_gidset,
      m_subdomain.neigh_legs() )
    {
      BOOST_FOREACH( gid_t I, rank_gidset.second )
      {
        gidcell_mp::iterator it = neigh_cells.find( I );
        assert( it != neigh_cells.end() );
        s_plist_p_cellv_compv[i].push_back( &( it->second.incoming() ) );
      }
      i++;

    }

    // \todo CHECK ADDRESSES EXPLICITLY

    i = 0;
    BOOST_FOREACH( Subdomain::rank_gid_s_mp::const_reference rank_gidset,
      m_subdomain.home_legs() )
    {
      BOOST_FOREACH( gid_t I, rank_gidset.second )
      {
        gidcell_mp::iterator it = cells.find( I );
        assert( it != cells.end() );
        r_plist_p_cellv_compv[i].push_back( &( it->second.incoming() ) );
      }
      i++;
    }

    // Then store types of particle struct parts 
    m_part_mpi_type.get_section_types( m_prt_section_types );

  }

  // +-------------------------------+
  // | Particle_p2p::start_messaging |
  // +-------------------------------+
  void Particle_p2p::start_messaging()
  {
    //! short circuit if there is no neighbor nodes (single-carlo-node run)
    if (m_num_neigh==0) return;

    //! Create send requests and update particle counts
    reset_particle_send_requests();

    //! Synchronize particle counts

    MPI_Startall( m_num_neigh, s_pcount_req_compv[0].p() );
    MPI_Startall( m_num_neigh, r_pcount_req_compv[0].p() );
    //! Sit and wait for particle count sends to complete
    MPI_Waitall( m_num_neigh, s_pcount_req_compv[0].p(), MPI_STATUSES_IGNORE );
    //! Start particle sends
    MPI_Startall( m_num_neigh, s_particle_req_compv[0].p() );

    //! Wait count receives to complete
    MPI_Waitall( m_num_neigh, r_pcount_req_compv[0].p(), MPI_STATUSES_IGNORE );
    //! Add new particles based on refreshed and create recv requests
    reset_particle_recv_requests();

    //! Start particle receives 
    MPI_Startall( m_num_neigh, r_particle_req_compv[0].p() );

    //! Wait for particle sends to complete
    MPI_Waitall( m_num_neigh, s_particle_req_compv[0].p(), MPI_STATUSES_IGNORE);
  }

  // +----------------------------------+
  // | Particle_p2p::complete_messaging |
  // +----------------------------------+
  void Particle_p2p::complete_messaging() 
  {  
    //! short circuit if there is no neighbor nodes (single-carlo-node run)
    if ( m_num_neigh == 0 ) return;

    //! Wait for particle receives to complete
    MPI_Waitall( m_num_neigh, r_particle_req_compv[0].p(), MPI_STATUSES_IGNORE);

    // Purge incoming particles in neighbor-leg cells (they have been sent)
    for (plist_p_v_v::iterator // for all compute nodes
      iter1  = s_plist_p_cellv_compv.begin();
      iter1 != s_plist_p_cellv_compv.end(); ++iter1)
      for (plist_p_v::iterator // for all cells bridged with this compute node
        iter2  = iter1->begin();
        iter2 != iter1->end()  ;   ++iter2)
        // purge the list (*iter2 is pointer to that list)
        (*iter2)->clear();
  }


  // +--------------------------------------+
  // | Particle_p2p::create_pcount_requests |
  // +--------------------------------------+
  void Particle_p2p::create_pcount_requests()
  {

    // resize
    s_counts_cellv_compv.resize( m_subdomain.nneigh_subdomains() );
    s_pcount_req_compv.resize  ( m_subdomain.nneigh_subdomains() );
    r_counts_cellv_compv.resize( m_subdomain.nneigh_subdomains() );
    r_pcount_req_compv.resize  ( m_subdomain.nneigh_subdomains() );

    // send types
    std::size_t i = 0;
    BOOST_FOREACH( Subdomain::rank_gid_s_mp::const_reference rank_gidset, 
      m_subdomain.neigh_legs() )
    {

      // allocate space for counts
      int cell_count = static_cast<int>( rank_gidset.second.size() );
      s_counts_cellv_compv[i].resize( cell_count );

      // create persistent send
      MPI_Send_init(                 // send
        &s_counts_cellv_compv[i][0], // starting at this address
        cell_count,                  // this many
        MPI_INT,                     // of this datatype
        rank_gidset.first,           // to this node
        m_tags[PCOUNT_TAG],          // tag the message with this
        m_comm,                      // inside this communicator
        s_pcount_req_compv[i].p()    // assign request into this
        );

      i++;

    }

    // recv types
    i = 0;
    BOOST_FOREACH( 
      Subdomain::rank_gid_s_mp::const_reference rank_gidset,
      m_subdomain.home_legs()
      )
    {

      // allocate space for counts
      int cell_count = static_cast<int>( rank_gidset.second.size() );
      r_counts_cellv_compv[i].resize( cell_count );

      // create persistent receive
      MPI_Recv_init(                 // receive
        &r_counts_cellv_compv[i][0], // and write into this address
        cell_count,                  // this many elements
        MPI_INT,                     // of this datatype
        rank_gidset.first,           // from this compute node
        m_tags[PCOUNT_TAG],          // the message with this tag
        m_comm,                      // inside this communicator
        r_pcount_req_compv[i].p()    // assign receive request into this
        );

      i++;

    }

  }

  // +--------------------------------------------+
  // | Particle_p2p::reset_particle_send_requests |
  // +--------------------------------------------+
  /** \brief Clear previous particle send requests and issues new ones.
  * 
  * So, we will "send" the particles here ("send" in quotes, since we will take 
  * care of every detail of sending in here, but will not actually initiate the 
  * network message, i.e. we will create the send \e request). 
  * 
  * The "send" operation can be outlined as follows:
  * - Erase whatever request was allocated before and initialize new ones
  * - At this point we have the send counts per each cell. Using this 
  *
  * Implementation Notes: 
  * - Internal vectors are static to avoid destruction each time, with the hope
  *   that vector::resize() will optimize reallocation.
  * 
  * \todo Using MPI_Rsend_init leads to some runtime misbehavior. I do not know
  *       yet why. Investigate.
  */
  void Particle_p2p::reset_particle_send_requests()
  {
    s_particle_req_compv.clear(); // this automatically calls destructor for each request 
    s_particle_req_compv.resize( m_num_neigh );
    // datatype also needs to be stored. this is because in BigBen (CRAY MPICH) 
    // the persistent send can not start if the datatype is freed
    s_particle_dty_compv.clear(); 
    s_particle_dty_compv.resize( m_num_neigh );

    // Start looping for all compute nodes to send to
    for (size_t iNeigh = 0; iNeigh != m_num_neigh; ++iNeigh)
    {
      static address_v_a hindex_displ_a;   
      // get the particle structure part addresses
      count_gather_particle_addresses( iNeigh, hindex_displ_a );
      // Compose particle message, and select appropriate target
      // create_particle_message nullifies target if no particles are to be
      // sent. 
      int target = m_neighbors_compv[iNeigh];
      compose_particle_message( 
        hindex_displ_a,
        s_particle_dty_compv[iNeigh],
        target );

      // create persistent send
      MPI_Send_init(                      // send
        MPI_BOTTOM,                      // data at this address
        1,                               // this many
        s_particle_dty_compv[iNeigh],    // of this datatype
        target,                          // to this compute node
        //s_particle_tag_compv[iNeigh],  // message with this tag
        m_tags[PDATA_TAG],
        m_comm,                          // inside this communicator
        s_particle_req_compv[iNeigh].p() // store send request into this one
        );
    }
  }

  
  // +--------------------------------------------+
  // | Particle_p2p::reset_particle_recv_requests |
  // +--------------------------------------------+
  /** \brief Clear previous particle receive requests and issues new ones.
  *
  *   See implementation note of reset_particle_send_requests()
  */
  void Particle_p2p::reset_particle_recv_requests()
  {
    r_particle_req_compv.clear(); // this automatically calls destructor for each request 
    r_particle_req_compv.resize( m_num_neigh );
    r_particle_dty_compv.clear();
    r_particle_dty_compv.resize( m_num_neigh );

    // Start looping for all compute nodes to recv from
    for (size_t iNeigh = 0; iNeigh != m_num_neigh; ++iNeigh)
    {
      static address_v_a hindex_displ_a;
      // 
      // Add new particles and get their addresses
      //
      add_gather_particle_addresses( iNeigh, hindex_displ_a );
      //
      // Compose receive message using gathered addresses
      // 
      int target = m_neighbors_compv[iNeigh];
      compose_particle_message( 
        hindex_displ_a,
        r_particle_dty_compv[iNeigh],
        target); // nullifies target if necessary

      // create persistent recv
      MPI_Recv_init(                      // receive
        MPI_BOTTOM,                      // into this address
        1,                               // this many
        r_particle_dty_compv[iNeigh],    // of this datatype
        target,                          // from this compute node
        //r_particle_tag_compv[iNeigh],  // use this message tag
        m_tags[PDATA_TAG],
        m_comm,                          // inside this communicator
        r_particle_req_compv[iNeigh].p() // store receive request here
        );
    }
  }


  // +----------------------------------------+
  // | Particle_p2p::compose_particle_message |
  // +----------------------------------------+
  void Particle_p2p::compose_particle_message(
    address_v_a& addr,
    mytools::mpi::Datatype& message,
    int& target
    )
  {
    enum { NSECT = ParticleMpiType::NSECTION };
    // addr has NSECT set of addresses. Each set has identical number of address
    // values (corresponding to appropriate section of the particle struct see
    // particle_mpi_types class for more details on particle sections).
    //
    // Get number of particles
    int num_particles = addr[0].size();

    // Create message type:

    // Degenerate case of no particles  is dealt with as follows:
    // . Use a dummy type for the message 
    // . Use MPI_PROC_NULL for p2p communication 
    if (num_particles!=0) {
      // All block-length's are 1 (I think newer MPIs has this fixed;  
      //+ but there is no reason to risk portability)
      static std::vector<int> hindex_blen;
      hindex_blen.resize(num_particles,1);
      // one hindexed type per struct part will compose struct types
      mytools::mpi::Datatype struct_type[ NSECT ];
      for (int i=0; i!=NSECT; ++i) {
        MPI_Type_hindexed( num_particles, 
          &hindex_blen[0], &addr[i][0], m_prt_section_types[i],
          struct_type[i].p()
          );
      }
      int struct_blen[NSECT]; 
      std::fill_n(struct_blen , (size_t)NSECT, 1);
      MPI_Aint struct_displ[NSECT]; 
      std::fill_n(struct_displ, (size_t)NSECT, (MPI_Aint)MPI_BOTTOM);
      MPI_Type_struct( NSECT, struct_blen, struct_displ, struct_type[0].p(), 
        message.p() );
    } else { // zero particles:
      // Here, effectively, we won't be sending ANY message. However, to have a 
      // consistent interface with the non-degenerate cases we need to create a
      // dummy message type. There are things to consider while creating this
      // dummy type:
      // - We can not use elementary types (such as MPI_INT) since we need to 
      //   specify a message type that is a user defined type, i.e. a type that
      //   CAN be freed.
      // - The message should actually point to a valid memory location. Since
      //   sends will use MPI_BOTTOM as the buffer address, we need to pass a
      //   valid memory location. 
      // Given these issues, let's create a type that corresponds to the sample
      // particle's first struct section:
      target = MPI_PROC_NULL;
      int blen[] = {1};
      MPI_Aint addr[1]; MPI_Address( &m_sample_particle, &addr[0] );
      MPI_Type_hindexed(1, blen, addr, m_prt_section_types[0], message.p());
    }

    // commit eventually
    MPI_Type_commit( message.p() );
  }


  // +-----------------------------------------------+
  // | Particle_p2p::count_gather_particle_addresses |
  // +-----------------------------------------------+
  /** \brief Count the particles to be sent and gather their addresses.
  *
  * For a given neighbor node, look at the neighbor-leg cells and count their
  * incoming particles. Also, gather the addresses of those particles. 
  */
  void Particle_p2p::count_gather_particle_addresses(
    std::size_t iNeighbor,
    address_v_a& addrs
    )
  {
    // reset size
    for ( int i=0; i!=ParticleMpiType::NSECTION; addrs[i++].resize(0) );

    // loop for each neighbor-leg cells bridged to iNeighbor-th compute node
    plist_p_v::iterator incoming_particle_list = // loop their incoming particle lists
      s_plist_p_cellv_compv[iNeighbor].begin();
    pcount_v::iterator  incoming_particle_count = // and their particle counts 
      s_counts_cellv_compv[iNeighbor].begin();
    pcount_v::iterator  incoming_particle_count_end =
      s_counts_cellv_compv[iNeighbor].end();
    // 
    for (; incoming_particle_count!=incoming_particle_count_end; 
      ++incoming_particle_count, ++incoming_particle_list)
    {
      // Store number of particles in this grid
      *incoming_particle_count = (*incoming_particle_list)->size();
      // Store each incoming particle addresses
      gather_particle_addresses( 
        (*incoming_particle_list)->begin(),
        (*incoming_particle_list)->end(),
        addrs
        );
    }
  }


  // +---------------------------------------------+
  // | Particle_p2p::add_gather_particle_addresses |
  // +---------------------------------------------+
  /** \brief Create new particles and gather their addresses
  *
  *  For a given neighbor node, looks at how many particle are going to be 
  *  received into which home cells and then creates new particles inside that 
  *  cells (or rather allocates new particles). 
  *   
  *  Then for further stage of actually receiving the particle data, gathers the
  *  newly created particle addresses.
  */
  void Particle_p2p::add_gather_particle_addresses(
    std::size_t iNeighbor, address_v_a& addrs )
  {

    // reset size
    for ( int i=0; i!=ParticleMpiType::NSECTION; addrs[i++].resize(0) );

    // loop for each home-leg cell bridged to iNeighbor-th compute node
    plist_p_v::iterator incoming_particle_list = 
      r_plist_p_cellv_compv[iNeighbor].begin();
    pcount_v::iterator  incoming_particle_count  = 
      r_counts_cellv_compv[iNeighbor].begin();
    pcount_v::iterator  incoming_particle_count_end  = 
      r_counts_cellv_compv[iNeighbor].end();

    for (; incoming_particle_count != incoming_particle_count_end;
      ++incoming_particle_count, ++incoming_particle_list )
    {
      // Create new list of particles 
      Cell::particle_lst new_particles( *incoming_particle_count, 
        m_sample_particle );
      // gather the addresses of new particles
      gather_particle_addresses( 
        new_particles.begin(),
        new_particles.end(),
        addrs
        );
      // Splice new list into home-cell's incoming list
      (*incoming_particle_list)->splice( 
        (*incoming_particle_list)->begin(),
        new_particles
        );
    }

  }

  // +-----------------------------------------+
  // | Particle_p2p::gather_particle_addresses |
  // +-----------------------------------------+
  //! Loop through given iterator range and gather the addresses into last argument.
  void Particle_p2p::gather_particle_addresses(
    Cell::particle_it p,
    Cell::particle_it end,
    address_v_a& into
    )
  {
    // Loop each particle storing their addresses
    for (; p!=end; ++p)
    {
      ParticleMpiType::section_addresses pa;
      // get address of each particle struct part 
      m_part_mpi_type.get_section_addresses(*p, pa);
      // assign to proper locations 
      for ( int i=0; i!=ParticleMpiType::NSECTION; ++i )
        into[i].push_back(pa[i]);
    }
  }


} // namespace iplmcfd
