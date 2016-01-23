#include "subdomain.hpp"

#include <algorithm>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <boost/foreach.hpp>
#include <mpi.h>
#include <metis.h>

#include "pattern.hpp"
#include "diagnostic/tracer.hpp"

namespace iplmcfd {

  std::ostream& operator<<( std::ostream& out, const Subdomain::adj_mp& graph )
  {

    BOOST_FOREACH( Subdomain::adj_mp::const_reference gid_neighlst, graph )
    {
      out << gid_neighlst.first << " ";
      for ( std::size_t i=0; i<gid_neighlst.second.size(); i++ )
        out << gid_neighlst.second[i] << " ";
      out << std::endl;
    }
    return out;

  }

  // +------------------------------+
  // | Subdomain::build_local_graph |
  // +------------------------------+
  void Subdomain::_build_local_graph( const std::vector< rank_t >& colors )
  {
    tracer::scope _( "Subdomain::build_local_graph" );

    // get homes
    for (gid_t I=0; I<colors.size(); ++I)
    {
      // if the vertex's color indicates it is not a home vertex, we are not
      // interested in it here.
      if ( colors[I] != _rank ) continue;

      _homes[I] = 0;
    }

    // enumerate local ids
    lid_t i = 0;
    BOOST_FOREACH( gidlid_mp::reference gidlid, _homes )
      gidlid.second = i++;

    // warn about empty partitions
    if ( _homes.size() == 0 )
      std::cerr << "Iplmcfd::Subdomain warning: " <<
      "Zero nodes in partition " << _rank << std::endl;

    // build neighbor dependency graph
    BOOST_FOREACH( gidlid_mp::const_reference gidlid, _homes )
    {
      gid_t I = gidlid.first;
      index_v index = _domain.gid2ind( I );

      // loop through stencil points
      std::vector< int_v > stencil = _pattern.elements();
      _adj[I].reserve( stencil.size() );
      _padj[I].reserve( stencil.size() );
      for ( std::size_t i=0; i<stencil.size(); i++ )
      {

        int_v neigh_index( index + stencil[i] );
        gid_t Ineigh = _domain.ind2gid( neigh_index );

        if ( _domain.out_of_bounds( Ineigh ) ) continue;

        _adj[I].push_back( Ineigh );

        // store the neighbor rank
        _padj[I].push_back( colors[Ineigh] );

        if ( colors[Ineigh] == _rank ) continue;

        // it is a neighbor
        _neighs.insert( gidrank_mp::value_type( Ineigh, colors[Ineigh] ) );
        _home_legs[ colors[Ineigh] ].insert( I );
        _neigh_legs[ colors[Ineigh] ].insert( Ineigh );
        
      }

    } // foreach

    assert( _home_legs.size() == _neigh_legs.size() );


    /*
    // output neighbor info
    if ( _rank == 0 )
    {

      // output colors
      for ( std::size_t i=0; i<colors.size(); i++ )
        std::cout << i << " " << colors[i] << std::endl;
      
      std::cout << " 0: HOME - NEIGH LEG PAIRS " << std::endl;
      BOOST_FOREACH( gidgid_mp::const_reference home_neigh, home_neigh_legs )
        std::cout << home_neigh.first << " " << home_neigh.second << std::endl;

      BOOST_FOREACH( rank_gid_s_mp::const_reference rank_set, _home_legs )
      {
        std::cout << " 0: HOME_LEGS: " << rank_set.first << "(rank): ";
        BOOST_FOREACH( gid_t homegid, rank_set.second )
          std::cout << homegid << " ";
        std::cout << std::endl;
      }

      BOOST_FOREACH( rank_gid_s_mp::const_reference rank_set, _neigh_legs )
      {
        std::cout << " 0: NEIGH_LEGS: " << rank_set.first << "(rank): ";
        BOOST_FOREACH( gid_t neighgid, rank_set.second )
          std::cout << neighgid << " ";
        std::cout << std::endl;
      }

    }
    MPI_Barrier( MPI_COMM_WORLD );
    if ( _rank == 1 )
    {

      std::cout << " 1: HOME - NEIGH LEG PAIRS " << std::endl;
      BOOST_FOREACH( gidgid_mp::const_reference home_neigh, home_neigh_legs )
        std::cout << home_neigh.first << " " << home_neigh.second << std::endl;

      BOOST_FOREACH( rank_gid_s_mp::const_reference rank_set, _home_legs )
      {
        std::cout << " 1: HOME_LEGS: " << rank_set.first << "(rank): ";
        BOOST_FOREACH( gid_t homegid, rank_set.second )
          std::cout << homegid << " ";
        std::cout << std::endl;
      }

      BOOST_FOREACH( rank_gid_s_mp::const_reference rank_set, _neigh_legs )
      {
        std::cout << " 1: NEIGH_LEGS: " << rank_set.first << "(rank): ";
        BOOST_FOREACH( gid_t neighgid, rank_set.second )
          std::cout << neighgid << " ";
        std::cout << std::endl;
      }

    }

    */


  }

  // +----------------------------------+
  // | SerialSubdomain::SerialSubdomain |
  // +----------------------------------+
  SerialSubdomain::SerialSubdomain( 
    const simparam& param, 
    const Domain& d,
    const Pattern& pat, 
    MPI_Comm comm ) 
    : Subdomain( param, d, pat, comm )
  {
    tracer::scope _( "SerialSubdomain::SerialSubdomain" );

    // allocate colors
    std::vector< rank_t > colors( d.npoints(), MAX_SIZE );

    // get decomposition
    _decomp( colors );

    // set the homes
    _build_local_graph( colors );

  }

  // +-------------------------+
  // | SerialSubdomain::decomp |
  // +-------------------------+
  void SerialSubdomain::_decomp( std::vector< rank_t >& colors )
  {
    tracer::scope _( "SerialSubdomain::decomp" );

    std::fill( colors.begin(), colors.end(), 0 );

  }

  // +--------------------------------+
  // | BlockSubdomain::BlockSubdomain |
  // +--------------------------------+
  BlockSubdomain::BlockSubdomain( 
    const simparam& param, 
    const Domain& d,
    const Pattern& pat, 
    MPI_Comm comm ) 
    : Subdomain( param, d, pat, comm )
  {
    tracer::scope _( "BlockSubdomain::BlockSubdomain" );

    // check args
    int npart;
    MPI_Comm_size( comm, &npart );
    const int_v& ncores = _simparam.ncores;
    if ( npart != tvmet::product( ncores ) )
    {
      std::ostringstream err;
      err << "MPI size is: " << npart << " Size requested by user is: ";
      err << tvmet::product( ncores ) << std::endl;
      throw std::invalid_argument( err.str() );
    }

    // allocate colors
    std::vector< rank_t > colors( d.npoints(), MAX_SIZE );

    // get decomposition
    _decomp( colors );

    // set the homes
    _build_local_graph( colors );

  }

  // +------------------------+
  // | BlockSubdomain::decomp |
  // +------------------------+
  void BlockSubdomain::_decomp( std::vector< rank_t >& colors )
  {
    tracer::scope _( "BlockSubdomain::decomp" );

    // determine number of FD points in each decomposition
    const int_v& ncores = simparam::instance().ncores;
    int local_remainder;

    // XDIM 
    std::vector<int> local_nx( ncores[XDIM] );
    for (int i=0; i<ncores[XDIM]; ++i)
      local_nx[i] = _domain.shape()[XDIM] / ncores[XDIM];
    // add any remainder to front
    local_remainder = _domain.shape()[XDIM] % ncores[XDIM];
    for (int i=0; i<local_remainder; ++i) local_nx[i]++;

    // YDIM 
    std::vector<int> local_ny( ncores[YDIM] );
    for (int i=0; i<ncores[YDIM]; ++i)
      local_ny[i] = _domain.shape()[YDIM] / ncores[YDIM];
    // add any remainder to back
    local_remainder = _domain.shape()[YDIM] % ncores[YDIM];
    for (int i=0; i<local_remainder; ++i) local_ny[ ncores[YDIM]-1-i ]++;

    // ZDIM 
    std::vector<int> local_nz( ncores[ZDIM] );
    for (int i=0; i<ncores[ZDIM]; ++i)
      local_nz[i] = _domain.shape()[ZDIM] / ncores[ZDIM];
    // add any remainder to front
    local_remainder = _domain.shape()[ZDIM] % ncores[ZDIM];
    for (int i=0; i<local_remainder; ++i) local_nz[i]++;

    // get start and end subscripts for each subdomain
    std::vector<int> istart( ncores[XDIM], 0 );
    std::vector<int> iend( ncores[XDIM], 0 );
    istart[0] = 0;
    iend[0]   = local_nx[0];
    for (int i=1; i<ncores[XDIM]; ++i)
    {
      istart[i] = iend[i-1];
      iend[i]   = istart[i] + local_nx[i];
    }

    std::vector<int> jstart( ncores[YDIM], 0 );
    std::vector<int> jend( ncores[YDIM], 0 );
    jstart[0] = 0;
    jend[0]   = local_ny[0];
    for (int j=1; j<ncores[YDIM]; ++j)
    {
      jstart[j] = jend[j-1];
      jend[j]   = jstart[j] + local_ny[j];
    }

    std::vector<int> kstart( ncores[ZDIM], 0 );
    std::vector<int> kend( ncores[ZDIM], 0 );
    kstart[0] = 0;
    kend[0]   = local_nz[0];
    for (int k=1; k<ncores[ZDIM]; ++k)
    {
      kstart[k] = kend[k-1];
      kend[k]   = kstart[k] + local_nz[k];
    }

    // populate colors array
    rank_t color = 0;
    for (int iproc=0; iproc<ncores[XDIM]; ++iproc) {
      for (int jproc=0; jproc<ncores[YDIM]; ++jproc) {
        for (int kproc=0; kproc<ncores[ZDIM]; ++kproc) {
          for (int i=istart[iproc]; i<iend[iproc]; ++i) {
            for (int j=jstart[jproc]; j<jend[jproc]; ++j) {
              for (int k=kstart[kproc]; k<kend[kproc]; ++k)
              {
                int_v sub(i,j,k);
                colors[ _domain.ind2gid(sub) ] = color;
              } // k
            } // j
          } // i
          ++color;
        } // kproc
      } // jproc
    } // iproc

  }

 
  // +--------------------------------+
  // | MetisSubdomain::MetisSubdomain |
  // +--------------------------------+
  MetisSubdomain::MetisSubdomain(
    const simparam& param, 
    const Domain& d,
    const Pattern& pat, 
    MPI_Comm comm )
    : Subdomain( param, d, pat, comm )
  {
    tracer::scope _("MetisSubdomain::MetisSubdomain");

    // allocate colors
    std::vector< rank_t > colors( d.npoints(), MAX_SIZE );

    // get decomposition 
    _decomp( colors );

    // set the homes
    _build_local_graph( colors );

  }

  // +-------------------------+
  // |  MetisSubdomain::decomp |
  // +-------------------------+
  void MetisSubdomain::_decomp( std::vector< rank_t >& colors )
  {
    tracer::scope _("MetisSubdomain::metis_decomp");
    tracer& log = tracer::instance();

    if ( my_rank() == 0 )
    {

      std::vector< metis::index_type > adj;
      std::vector< metis::index_type > xadj;

      adj.reserve( 2 * _domain.npoints() );
      xadj.reserve( _domain.npoints() );

      // for initial decomposition build graph with simple dependency
      BoxPattern box( 1 );
      const std::vector< int_v >& boxe = box.elements();
      for ( std::size_t I=0; I<_domain.npoints(); I++ )
      {

        index_v index = _domain.gid2ind( I );

        // Directly from Metis doc...
        /* The adjacency structure of the graph is stored as follows:
        *  the adjacency list of vertex i is stored in array adj starting 
        *  at index xadj[i] and ending at (but not including) index xadj[i+1]
        *  (i.e., adj[xadj[i]] through and including adj[xadj[i+1]-1]). That
        *  is, for each vertex i, its adjacency list is stored in consecutive
        *  locations in the array adj, and the array xadj is used to point to
        *  where it begins and where it ends.
        */

        // loop through stencil points
        xadj.push_back( metis::index_type( adj.size() ) );
        for ( std::size_t i=0; i<boxe.size(); i++ )
        {

          int_v neigh_index( index + boxe[i] );
          metis::index_type Ineigh = _domain.ind2gid( neigh_index );

          if ( _domain.out_of_bounds( Ineigh ) ) continue;

          adj.push_back( Ineigh );

        }

      }
      xadj.push_back( metis::index_type( adj.size() ) );

      int nVertex  = _domain.npoints();
      std::vector< int > wgts( _domain.npoints(), 0 );
      int wgtflag = metis::WEIGHT_NONE;
      int numflag = 0; // vertex indices start at 0
      int options[metis::N_OPTIONS]; options[0] = 0; // use defaults
      int edgecut;

      int npart;
      MPI_Comm_size( _comm, &npart );

      //! From metis documentation:
      //! "The number of parts that the graph will be partitioned into. 
      //! It should be greater than 1."
      if ( npart > 1 )
      {

        // See Metis documentation for suggested usage of 
        // PartGraphRecursive and PartGraphKway
        if ( npart <= 8 )
          METIS_PartGraphRecursive(
          &nVertex,     // this many vertices
          &xadj[0],     // adjacency array indices of vertices
          &adj[0],      // adjacency array
          &wgts[0],     // vertex weights
          NULL,         // edge weights
          &wgtflag,     // weight type selector
          &numflag,     // array origin index
          &npart,       // number of partitions
          options,      // options
          &edgecut,     // output: total edges cut with partition
          &colors[0] ); // output: partition colors
        else
          METIS_PartGraphKway(
          &nVertex,     // this many vertices
          &xadj[0],     // adjacency array indices of vertices
          &adj[0],      // adjacency array
          &wgts[0],     // vertex weights
          NULL,         // edge weights
          &wgtflag,     // weight type selector
          &numflag,     // array origin index
          &npart,       // number of partitions
          options,      // options
          &edgecut,     // output: total edges cut with partition
          &colors[0] ); // output: partition colors
      }
      else // serial sim - everything is owned by root
        std::fill( colors.begin(), colors.end(), 0 );

    } // my_rank() == 0

    log.range("Subdomain::metis_decomp: Bcast colors");
    MPI_Bcast(
      &colors[0],    // buffer
      colors.size(), // number of elements
      mytools::mpi::get_datatype( metis::index_type() ), // of this type
      0, // collective root
      _comm  // communicator
      );

  }

} // namespace iplmcfd
