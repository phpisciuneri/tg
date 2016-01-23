#ifndef IPLMCFD_SUBDOMAIN_HPP_
#define IPLMCFD_SUBDOMAIN_HPP_

#include <map>
#include <set>
#include <vector>
#include <mpi.h>

#include "defs.hpp"
#include "simparam.hpp"
#include "domain.hpp"
#include "pattern.hpp"

namespace iplmcfd
{

  //!
  //! \file subdomain.hpp
  //! \class Subdomain
  //! \brief Performs domain decomposition and stores information about the 
  //!        resulting local partition
  //!
  //! # Nomenclature:
  //!  - __homes__: The nodes assigned to this partition
  //!
  class Subdomain
  {

    friend class MigrateZObjects;

  public:

    //! Global ID, Local ID map type
    typedef std::map< gid_t, lid_t > gidlid_mp;
    //! Global ID, rank map type
    typedef std::map< gid_t, rank_t > gidrank_mp;
    //! Global ID set type
    typedef std::set< gid_t > gid_s;
    //! Neighboring subdomain rank, global ID set map
    typedef std::map< rank_t, gid_s > rank_gid_s_mp;
    
    //! \todo DOCUMENT
    //! \todo adj_mp should be composed of a set, not vector
    typedef std::map< gid_t, std::vector< int > > adj_mp;

  public:

    //!
    //! \brief Construct Subdomain
    //!
    //! \param[in] p Simulation parameters reference
    //! \param[in] d Domain reference 
    //! \param[in] pat Pattern reference 
    //! \param[in] comm MPI communicator
    //!
    Subdomain( 
      const simparam& p, 
      const Domain& d, 
      const Pattern& pat, 
      MPI_Comm comm )
      : _simparam( p ), _domain( d ), _comm( comm ), _pattern( pat )
    {
      MPI_Comm_rank( comm, &_rank );
    }
    
    //!
    //! \brief Get the number of nodes that call this subdomain home
    //!
    //! \returns the number of nodes
    //!
    std::size_t nhomes() const {  return _homes.size(); }

    //!
    //! \brief Get the global ID to local ID map for home nodes
    //!
    //! \returns the global ID to local ID map of home nodes
    //!
    const gidlid_mp& homes() const { return _homes; }

    //!
    //! \brief Get the number of neighboring nodes
    //!
    //! \returns the number of neighboring nodes
    //!
    std::size_t nneighs() const { return _neighs.size(); }

    //!
    //! \brief Get the global ID set of neighbor nodes
    //!
    //! \returns set of global IDs of neighboring nodes
    //!
    const gidrank_mp& neighs() const { return _neighs; }
    
    //!
    //! \brief Get the set of home node gids mapped to rank of each 
    //!        communication neighbor
    //!
    //! \returns
    //!
    const rank_gid_s_mp& home_legs() const { return _home_legs; }

    //!
    //! \brief Get the set of ghost (neighbor) node gids mapped to rank of 
    //!        owning process
    //!
    //! \returns
    //!
    const rank_gid_s_mp& neigh_legs() const { return _neigh_legs; }

    //!
    //! \brief Get the number of neighboring subdomains
    //!
    //! \returns
    //!
    std::size_t nneigh_subdomains() const { return _home_legs.size(); }

    //!
    //! \brief Get MPI rank that owns this subdomain instance
    //!
    //! \returns MPI rank
    //!
    rank_t my_rank() const { return _rank; }

    //!
    //! \brief Get the adjacency list
    //!
    //! \returns reference to adjacency list
    //!
    const adj_mp& adj() const { return _adj; }

    //!
    //! \brief Get processor for each neighbor in the adjacency list
    //!
    //! \returns reference to processor of each index in adjacency list
    //!
    const adj_mp& padj() const { return _padj; }

    //!
    //! \brief Output the adjacency graph (useful for debug)
    //!
    friend std::ostream& operator<<( std::ostream& out, const adj_mp& graph );

  protected:

    //!
    //! \brief Decomposes the domain
    //!
    //! \params[out] colors MPI rank for _every_ Global ID
    //!
    virtual void _decomp( std::vector< rank_t >& colors ) = 0;

    //!
    //! \brief Determines home and neighbor nodes.  Builds neighbor dependency
    //!        graph.
    //!
    //! \param[in] colors MPI rank for _every_ Global ID
    //!
    void _build_local_graph( const std::vector< rank_t >& colors );

  protected:

    //! Simulation parameters reference
    const simparam& _simparam;
    //! Domain reference
    const Domain& _domain;
    //! Reference to communication pattern
    const Pattern& _pattern;  
    //! MPI communicator
    MPI_Comm _comm;     
    //! MPI rank
    rank_t _rank;     
    //! global ID, local ID map of nodes
    gidlid_mp _homes; 
    //! global ID set of neighbor nodes
    gidrank_mp _neighs;
    //! set of home node gids mapped to rank of communication neighbor
    rank_gid_s_mp _home_legs;
    //! set of ghost (neighbor) node gids mapped to rank of owner
    rank_gid_s_mp _neigh_legs;


    //! \todo DOCUMENT
    adj_mp _adj;
    //!
    adj_mp _padj;

  };


  //!
  //! \class SerialSubdomain
  //! \brief Force creation of a serial subdomain
  //!
  //! This class mainly intended to aide debugging of repartitioning tools
  //!
  class SerialSubdomain : public Subdomain
  {
  public:

    friend class MigrateZObjects;

    //!
    //! \brief Construct BlockSubdomain
    //!
    //! \param[in] p Simulation parameters reference
    //! \param[in] d Domain reference 
    //! \param[in] pat Pattern reference 
    //! \param[in] comm MPI communicator
    //!
    SerialSubdomain( 
      const simparam& p, 
      const Domain& d, 
      const Pattern& pat, 
      MPI_Comm comm );

  protected:

    //!
    //! \brief Fill colors with 0.
    //!
    //! \param[out] colors MPI rank for _every_ Global ID
    //!
    void _decomp( std::vector< rank_t >& colors );

  };

  
  
  
  //!
  //! \class BlockSubdomain
  //! \brief Create a subdomain based on a block uniform decomposition
  //!
  class BlockSubdomain : public Subdomain
  {

  public:
    
    //!
    //! \brief Construct BlockSubdomain
    //!
    //! \param[in] p Simulation parameters reference
    //! \param[in] d Domain reference 
    //! \param[in] pat Pattern reference 
    //! \param[in] comm MPI communicator
    //!
    BlockSubdomain( 
      const simparam& p, 
      const Domain& d, 
      const Pattern& pat, 
      MPI_Comm comm );

  protected:

    //!
    //! \brief Calculate the block uniform decomposition
    //!
    //! \param[out] colors MPI rank for _every_ Global ID
    //!
    //! \todo clean up the block decomposition code
    //!
    void _decomp( std::vector< rank_t >& colors ); 

  };


  
  
  
  //!
  //! \class MetisSubdomain
  //! \brief Create a subdomain based on a (weighted) graph partition using
  //!        Metis.
  //!
  class MetisSubdomain : public Subdomain
  {

  public:

    //!
    //! \brief Constructs MetisSubdomain
    //!
    //! \param[in] p Simulation parameters reference
    //! \param[in] d Domain reference 
    //! \param[in] pat Pattern reference 
    //! \param[in] comm MPI communicator
    //!
    //! __Important Note:__ Construction must initialize private members of 
    //! the base class through the protected member functions.
    //!
    MetisSubdomain( 
      const simparam& param, 
      const Domain& d, 
      const Pattern& pat, 
      MPI_Comm comm );

  protected:

    //!
    //! \brief Calculate the (weighted) graph partitioned decomposition
    //!
    //! \param[out] colors MPI rank for _every_ Global ID
    //!
    void _decomp( std::vector< rank_t >& colors );      

  };

} // namespace iplmcfd

#endif // IPLMCFD_SUBDOMAIN_HPP_