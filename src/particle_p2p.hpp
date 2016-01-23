#ifndef IPLMCFD_PARTICLE_P2P_HPP_
#define IPLMCFD_PARTICLE_P2P_HPP_

#include <mytools/mpi/opaque_objects.hpp>

#include "subdomain.hpp"
#include "cell.hpp"
#include "particle.hpp"
#include "particle_mpi_type.hpp"

namespace iplmcfd {


  /** \brief Responsible for particle inter-domain communication 
   * 
   * Every subdomain is surrounded by bunch of other subdomains and it is assumed
   * that there is one to one mapping between subdomains and compute nodes in the
   * particle communicator. So, when we refer to a subdomain we also mean a single
   * compute node. Let's call the current compute node the \e home subdomain (or
   * \e home compute node) and the surrounding subdomains the \e neighbor
   * subdomain (or \e neighbor compute node). The ranks of neighbor compute nodes 
   * and other graph connection details are provided by the port and domain 
   * classes. 
   *
   * For each neighbor subdomain we have a set of cells (or rather the particles
   * in that cell) to send and to receive. Thus, the communication entities
   * defined in this class will be arrays of size equal to the number of neighbor
   * subdomains. Let's create a naming scheme:
   *
   * \arg \<entity\>_gridv  
   * \arg \<entity\>_compv  Some set of entities for each neighboring 
   * \arg \<entity\>_prtcv  for some set of entities for each particle 
   * compute node.
   * 
   * For example, take the following entity and some related typedefs: 
   * 
   * \code
   *  typedef particle_station::particle_list* plist_ptr; typedef
   *  std::vector<plist_ptr> plist_ptr_vector; typedef
   *  std::vector<plist_ptr_vector> plist_prt_vector_vector; 
   * \endcode 
   * Then, 
   * \code
   *  plist_ptr_vector send_lists_gridv; 
   * \endcode 
   * will be set of particle list pointers for each grid. And \code plist_prt_vector_vector
   * send_lists_gridv_compv; \endcode will be collection of for each neighboring
   * compute node the set of particle lists for each grid.
   * 
   * The following naming scheme is also adopted for convenience typedefs:
   * \arg \<entity\>_p is pointer to entity
   * \arg \<entity\>_v is vector of entities
   * \arg \<entity\>_a is array (compile time size) of entities
   * \arg s_\<entity\> is related to send operations
   * \arg r_\<entity\> is related to recieve operations
   * 
   * 
   * \todo Something else to fix: Send nodes of port are receive nodes of
   *       particle_port. 
   * And vice versa. This requires clearer treatment. Currently ignoring.
   */

  //!
  //! \file particle_p2p.hpp
  //! \class Particle_p2p
  //! \brief
  //!
  class Particle_p2p : boost::noncopyable
  {

  public:

    //!
    typedef std::map< gid_t, Cell > gidcell_mp;

  public:

    //!
    //!
    //!
    Particle_p2p( 
      const Subdomain& s, 
      gidcell_mp& cells, 
      gidcell_mp& neigh_cells,
      const mix_t& mix
      );
    
    //!
    //!
    //!
    void start_messaging();
    
    //!
    //! \brief Finalize particle messaging
    //!
    void complete_messaging();
    

   
  private: 

    //!
    typedef Cell::particle_lst* plist_p;
    //!
    typedef std::vector< plist_p > plist_p_v;
    //!
    typedef std::vector< plist_p_v > plist_p_v_v;
    //!
    typedef std::vector< mytools::mpi::Request > request_v;
    //!
    typedef std::vector< MPI_Aint > address_v;
    //!
    typedef std::vector< MPI_Aint > address_v_a[ ParticleMpiType::NSECTION];
    //!
    typedef std::vector< address_v_a > address_v_a_v;
    //typedef std::vector<int> tag_v;
    //!
    typedef std::vector< int > pcount_v;
    //!
    typedef std::vector< rank_t > rank_v;
    //!
    typedef std::vector<pcount_v> pcount_v_v;
    //!
    typedef std::vector< mytools::mpi::Datatype > datatype_v;
    
    
  private:

    // note:
    // create_* are called only once during initialization
    // others are repeatedly called for every communication.
    //void create_message_tags(size_t,size_t);

    //!
    //! \brief Create persistent send/recvs for particle counts per cell 
    //!
    void create_pcount_requests();

    //void create_plist_p_vectors(const Subdomain::gid_vector&, grid::cell_array&, plist_p_v&);
    
    //!
    //! \brief Clear previous particle send requests and issues new ones.
    //!
    void reset_particle_send_requests();

    //!
    //! \brief Clear previous particle receive requests and issues new ones.
    //!
    void reset_particle_recv_requests();

    //!
    //! \brief Create message type given particle addresses
    //!
    void compose_particle_message( address_v_a&, mytools::mpi::Datatype&, int& );

    //!
    //! \brief Count the particles to be sent and gather their addresses.
    //!
    void count_gather_particle_addresses( std::size_t, address_v_a& );

    //!
    //! \brief Create new particles and gather their addresses
    //!
    void add_gather_particle_addresses( std::size_t, address_v_a& );

    //!
    //!
    //!
    void gather_particle_addresses( Cell::particle_it, Cell::particle_it, address_v_a& );
    
  private:
    // particle count related
    // s_* for send, r_* for receive
    
    //! 
    pcount_v_v s_counts_cellv_compv;
    //!
    pcount_v_v r_counts_cellv_compv;
    //!
    request_v  s_pcount_req_compv;
    //!
    request_v  r_pcount_req_compv;


    // particle message related

    ///@{ 
    //! addresses of particle lists for sending
    plist_p_v_v  s_plist_p_cellv_compv;
    //! addresses of particle lists for receiving
    plist_p_v_v  r_plist_p_cellv_compv;
    //!
    request_v    s_particle_req_compv;
    //!
    request_v    r_particle_req_compv;
    //!
    datatype_v   s_particle_dty_compv;
    //!
    datatype_v   r_particle_dty_compv;
    //!@}

    // other members


    //! list of ranks of neighboring subdomains
    rank_v m_neighbors_compv;
    //!
    MPI_Comm     m_comm; ///< particle only communicator
    //! total number of neighboring subdomains
    std::size_t m_num_neigh;
    //!
    ParticleMpiType m_part_mpi_type;
    //!
    ParticleMpiType::section_types m_prt_section_types;
    //!
    Particle m_sample_particle;
    //!
    const Subdomain& m_subdomain;

    //!
    enum { PCOUNT_TAG, PDATA_TAG, N_TAGS };
    //!
    int m_tags[N_TAGS];

  };

} // namespace iplmcfd

#endif // IPLMCFD_PARTICLE_P2P_HPP_
