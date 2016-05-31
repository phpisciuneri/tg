#ifndef ZOBJECTS_HPP_
#define ZOBJECTS_HPP_

#include <zoltan_cpp.h>
#include <vector>
#include <cassert>
#include <map>
#include <string>
#include <boost/foreach.hpp>
#include <paragon.h>

#include "defs.hpp"
#include "subdomain.hpp"


namespace iplmcfd {

  //!
  //!  \file zobjects.hpp
  //!  \class ZObjects
  //!  \brief Helper class for (re)partitioning using Zoltan.
  //!
  //! Usage:
  //!
  //!   -# call the constructor
  //!   -# user performs initial decomposition
  //!   -# call set_num_local_objs() to register number of local objects
  //!   -# call set_global_local_ids() to register global/local ID mapping of
  //!      local objects
  //!   -# call set_obj_wgts to register weights of each local object to be
  //!      considered in partitioning
  //!   -# user must set any Zoltan parameters: i.e.
  //!      \code
  //!        Zoltan *zz = new Zoltan( MPI_COMM_WORLD );
  //!        zz->Set_Param( "LB_METHOD", "BLOCK" );
  //!        zz->Set_Param( "NUM_GID_ENTRIES", "1" );
  //!        zz->Set_Param( "NUM_LID_ENTRIES", "1" );
  //!        zz->Set_Param( "OBJ_WEIGHT_DIM", "1" );
  //!        zz->Set_Param( "AUTO_MIGRATE", "TRUE" );
  //!      \endcode
  //!   -# user must register functions with Zoltan: i.e.
  //!      \code
  //!        ZObjects ob(N);
  //!        Zoltan *zz = new Zoltan( MPI_COMM_WORLD );
  //!        zz->Set_Num_Obj_Fn( ZObjects::get_num_local_objs, &ob );
  //!        zz->Set_Obj_List_Fn( ZObjects::get_obj_list, &ob );
  //!      \endcode
  //!   -# call Zoltan::LB_Partition()
  //!
  //!

  class ZObjects
  {

  public:

    //!
    //! \brief Constructor for ZObjects class.
    //!
    ZObjects( const Subdomain& s ) : m_subdomain( s ) 
    {
      m_obj_wgts.resize( s.nhomes(), 1. );
    }
    
    //!
    //! \brief Register the weights for local objects
    //!
    //! \param[in] wgts vector of wgts for each local ID
    //!
    void set_obj_wgts( const std::vector<float>& wgts )
    {
      assert( wgts.size() == m_subdomain.nhomes() );
      m_obj_wgts.resize( wgts.size() );
      for (std::size_t i=0; i<m_subdomain.nhomes(); i++)
        m_obj_wgts[i] = wgts[i];

    }
      
    //!
    //! \brief Function to register for ZOLTAN_NUM_OBJ_FN.
    //!
    //! \param[in] data point to user-defined data
    //! \param[out] ierr error code
    //! \return number of objects that are assigned to the processor
    //!
    static int get_num_local_objs( void* data, int* ierr )
    {

      ZObjects *ob = static_cast< ZObjects* >( data );
      *ierr = ZOLTAN_OK;
      int nVertex = ob->m_subdomain.nhomes();
      return nVertex;
    }

    //!
    //! \brief Function to register for ZOLTAN_OBJ_LIST_FN
    //!
    //! \param[in] data pointer to user-defined data
    //! \param[in] num_gid_entries number of array entries used to describe a
    //!            global ID
    //! \param[in] num_lid_entries number of array entries used to describe a
    //!            local ID
    //! \param[out] global_ids array of global IDs
    //! \param[out] local_ids array of local IDs
    //! \param[in] wgt_dim number of weights associated with an object
    //! \param[out] obj_wgts array of object weights
    //! \param[out] ierr error code
    //!
    static void get_obj_list( void* data, 
      int num_gid_entries, 
      int num_lid_entries, 
      ZOLTAN_ID_PTR global_ids,
      ZOLTAN_ID_PTR local_ids,
      int wgt_dim,
      float* obj_wgts,
      int* ierr )
    {

      ZObjects *ob = static_cast< ZObjects* >( data );

      BOOST_FOREACH( Subdomain::gidlid_mp::const_reference gidlid, 
        ob->m_subdomain.homes() )
      {
        global_ids[gidlid.second] = ZOLTAN_ID_TYPE( gidlid.first );
        local_ids[gidlid.second]  = ZOLTAN_ID_TYPE( gidlid.second );
        obj_wgts[gidlid.second]   = ob->m_obj_wgts[gidlid.second];
      }
      *ierr = ZOLTAN_OK;
    }

    //!
    //! \brief Returns a list of parts to which given objects are assigned
    //!
    //! \param[in] data pointer to user-defined data
    //! \param[in] num_gid_entries number of array entries used to describe a global ID
    //! \param[in] num_lid_entries number of array entries used to describe a local ID
    //! \param[in] num_obj number of object IDs in arrays global_ids and local_ids
    //! \param[in] global_ids Array of global IDs of objects whose part should be returned
    //! \param[in] local_ids Array of local IDs of objects whose part should be returned
    //! \param[out] parts Array of part numbers corresponding to the global and local IDs
    //! \param[out] ierr error code
    //!
    static void get_part_multi(
        void* data,
        int num_gid_entries,
        int num_lid_entries,
        int num_obj,
        ZOLTAN_ID_PTR global_ids,
        ZOLTAN_ID_PTR local_ids,
        int* parts,
        int* ierr )
    {
        ZObjects *ob = static_cast< ZObjects* >( data );
        
        BOOST_FOREACH( Subdomain::gidlid_mp::const_reference gidlid,
                      ob->m_subdomain.homes() )
        {
            parts[gidlid.second] = ob->m_subdomain.my_rank();
        }
        *ierr = ZOLTAN_OK;
    }
      

    //!
    //! \brief Provides the number of edges in the communication graph for each
    //!        object
    //!
    //! \param[in] data pointer to user-defined data
    //! \param[in] num_gid_entries number of array entries used to describe a 
    //!            single global ID
    //! \param[in] num_lid_entries number of array entries used to describe a 
    //!            single local ID
    //! \param[in] num_obj The number of object IDs in arrays global_ids and 
    //!            local_ids
    //! \param[in] global_ids Array of global IDs of objects whose number of 
    //!            edges should be returned.  
    //! \param[in] local_ids Array of local IDs of objects whose number of 
    //!            edges should be returned.
    //! \param[out] num_edges an array containing numbers of edges. For object 
    //!             i (specified by global_ids[i*num_gid_entries] and 
    //!             local_ids[i*num_lid_entries], i=0,1,...,num_obj-1), the 
    //!             number of edges should be stored in num_edges[i].
    //! \param[out] ierr error code
    //!
    static void get_num_edges_multi( 
      void* data,
      int num_gid_entries,
      int num_lid_entries,
      int num_obj,
      ZOLTAN_ID_PTR global_ids,
      ZOLTAN_ID_PTR local_ids,
      int* num_edges,
      int* ierr )
    {
      ZObjects *ob = static_cast< ZObjects* >( data );

      Subdomain::adj_mp::const_iterator 
        gid_adj = ob->m_subdomain.adj().begin();

      for ( int n=0; n<num_obj; n++ )
      {
        assert( gid_adj->first == global_ids[n] );
        num_edges[n] = int( gid_adj->second.size() );
        gid_adj++;
      }
      *ierr = ZOLTAN_OK;
    }

    //!
    //! \brief Provides lists of global IDs, processor IDs, and (optionally)
    //!        edge weights for objects sharing edges with objects specified 
    //!        in the global_ids input array
    //!
    //! \param[in] data pointer to user-defined data
    //! \param[in] num_gid_entries number of array entries used to describe a 
    //!            single global ID
    //! \param[in] num_lid_entries number of array entries used to describe a 
    //!            single local ID
    //! \param[in] num_obj The number of object IDs in arrays global_ids and 
    //!            local_ids
    //! \param[in] global_ids Array of global IDs of objects whose number of 
    //!            edges should be returned.  
    //! \param[in] local_ids Array of local IDs of objects whose number of 
    //!            edges should be returned.
    //! \param[in] num_edges an array containing numbers of edges. For object 
    //!             i (specified by global_ids[i*num_gid_entries] and 
    //!             local_ids[i*num_lid_entries], i=0,1,...,num_obj-1), the 
    //!             number of edges should be stored in num_edges[i].
    //! \param[out] nbor_global_id an array of global IDs of objects sharing 
    //!             edges with the objects specified in global_ids. For object 
    //!             i (specified by global_ids[i*num_gid_entries] and 
    //!             local_ids[i*num_lid_entries], i=0,1,...,num_obj-1), edges 
    //!             are stored in nbor_global_id[sum*num_gid_entries] to 
    //!             nbor_global_id[(sum+num_edges[i])*num_gid_entries-1], where
    //!             sum = the sum of num_edges[j] for j=0,1,...,i-1.
    //! \param[out] nbor_procs an array of processor IDs that identifies where 
    //!             the neighboring objects reside. For neighboring object i 
    //!             (stored in nbor_global_id[i*num_gid_entries]), the 
    //!             processor owning the neighbor is stored in nbor_procs[i].
    //! \param[in] wgt_dim The number of weights associated with an edge 
    //!            (typically 1), or 0 if edge weights are not requested
    //! \param[out] ewgts an array of edge weights, where 
    //!             ewgts[i*wgt_dim:(i+1)*wgt_dim-1] corresponds to the weights
    //!             for the ith edge. If wgt_dim=0, the return value of ewgts 
    //!             is undefined and may be NULL.
    //! \param[out] ierr error code
    //!
    static void get_edge_list_multi(
      void* data,
      int num_gid_entries,
      int num_lid_entries,
      int num_obj,
      ZOLTAN_ID_PTR global_ids,
      ZOLTAN_ID_PTR local_ids,
      int* num_edges,
      ZOLTAN_ID_PTR nbor_global_id,
      int* nbor_procs,
      int wgt_dim,
      float* ewgts,
      int *ierr )
    {
      ZObjects *ob = static_cast< ZObjects* >( data );

      Subdomain::adj_mp::const_iterator 
        gid_adj( ob->m_subdomain.adj().begin() ),
        gid_padj( ob->m_subdomain.padj().begin() );

      for ( int n=0; n<num_obj; n++ )
      {
        assert( gid_adj->first == global_ids[n] );
        assert( gid_adj->second.size() == num_edges[n] );
        for ( int e=0; e<num_edges[n]; e++ )
        {
          *nbor_global_id++ = ZOLTAN_ID_TYPE( gid_adj->second[e] );
          *nbor_procs++ = gid_padj->second[e];
        }
        gid_adj++;
        gid_padj++;
      }
      
      *ierr = ZOLTAN_OK;
    }
      
  private:

    std::vector<float> m_obj_wgts;    /*!< weights of objects */

    //! \todo DOCUMENT
    const Subdomain& m_subdomain;
  };

} // namespace iplmcfd

#endif // ZOBJECTS_HPP_
