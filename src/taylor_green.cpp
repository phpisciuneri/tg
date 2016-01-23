#include "taylor_green.hpp"

#include <mpi.h>
#include "diagnostic/tracer.hpp"

namespace iplmcfd
{

  // +--------------------------+
  // | TaylorGreen::TaylorGreen |
  // +--------------------------+
  TaylorGreen::TaylorGreen( const simparam& p, const Domain& d, Subdomain& s )
    : m_simparam( p ),
    m_domain( d ), 
    m_subdomain( s ), 
    m_mc( p, d, s ),
    m_zob( s, m_mc )
  {
    tracer::scope _( "TaylorGreen::TaylorGreen" );

    if ( p.redecomp_type == "paragon" )
      init_paragon();
    else
      init_zoltan();
  }

  // +---------------------------+
  // | TaylorGreen::~TaylorGreen |
  // +---------------------------+
  TaylorGreen::~TaylorGreen()
  {
    tracer::scope _( "TaylorGreen::~TaylorGreen" );

    if ( m_simparam.redecomp_type == "paragon" )
    {
      free( m_rankCommCost );
      ParagonDeinit( &m_paragon );
    }

  }

  // +-------------------------+
  // | TaylorGreen::initialize |
  // +-------------------------+
  void TaylorGreen::initialize()
  {
    tracer::scope _( "TaylorGreen::initialize" );

    m_mc.initialize();

    m_dt = m_mc.update_dt();

    /*
     * After initialization we can register fields and create our initial
     * parallel io class.  These need reset after each repartitioning step
     */
    reset_output();
  }

  // +-------------------+
  // | TaylorGreen::step |
  // +-------------------+
  void TaylorGreen::step( std::size_t iter )
  {
    tracer::scope _( "TaylorGreen::step" );

    m_repart_time = 0;
    static double last_time = MPI_Wtime();

    m_mc.position_step( m_dt, iter );

    m_mc.scalar_step( m_dt, iter );

    m_dt = m_mc.update_dt();
    /*
     * update_dt() ends with an all_reduce which is implicitly blocking.
     * Thus, this is a good location for iteration timing.
     */
    m_step_wtime = MPI_Wtime() - last_time;
    last_time    = MPI_Wtime();
  }

  // +---------------------+
  // | TaylorGreen::output |
  // +---------------------+
  void TaylorGreen::output( std::size_t iter )
  {
    tracer::scope _( "TaylorGreen::output" );

    m_mc.update_fields();
    m_io->output_paraview( *m_fields, iter );
  }

  // +---------------------------+
  // | TaylorGreen::reset_output |
  // +---------------------------+
  void TaylorGreen::reset_output()
  {
    tracer::scope _( "TaylorGreen::reset_output" );
    m_fields.reset( new Field );
    m_mc.register_fields( *m_fields );

    // create parallel io type
    m_io.reset( new Parallel_io( m_domain, m_subdomain, *m_fields ) );
  }

  // +---------------------------+
  // | TaylorGreen::output_steps |
  // +---------------------------+
  void TaylorGreen::output_steps( std::ostream& out, std::size_t iter )
  {

    static bool first_call = true;

    if ( m_subdomain.my_rank() != 0 ) 
      return;

    if ( first_call )
    {
      out << std::setw(6)  << "iter"  << " ";
      out << std::setw(12) << "time"  << " ";
      out << std::setw(12) << "dt"    << " ";
      out << std::setw(12) << "ke"    << " ";
      out << std::setw(12) << "wtime" << " ";
      out << std::setw(12) << "repart_time" << std::endl;

      first_call = false;
    }

    out << std::setw(6)  << iter << " ";
    out << std::setw(12) << m_mc.time() << " ";
    out << std::setw(12) << m_dt << " ";
    out << std::setw(12) << m_mc.ke() << " ";
    out << std::setw(12) << m_step_wtime << " ";
    out << std::setw(12) << m_repart_time << std::endl;
  }

  // +--------------------------+
  // | TaylorGreen::init_zoltan |
  // +--------------------------+
  void TaylorGreen::init_zoltan()
  {
    tracer::scope _( "TaylorGreen::init_zoltan" );

    m_zz.reset( new Zoltan( MPI_COMM_WORLD ) );

    m_zz->Set_Param( "DEBUG_LEVEL", m_simparam.zoltan_debug_level );
    m_zz->Set_Param( "LB_METHOD", "GRAPH" );
    m_zz->Set_Param( "NUM_GID_ENTRIES", "1" );
    m_zz->Set_Param( "NUM_LID_ENTRIES", "1" );
    m_zz->Set_Param( "AUTO_MIGRATE", "TRUE" );
    m_zz->Set_Param( "OBJ_WEIGHT_DIM", "1" );
    m_zz->Set_Param( "LB_APPROACH", "REPARTITION" );
    m_zz->Set_Param( "IMBALANCE_TOL", "1.01" );

    char buffer[64];
    sprintf(buffer, "%d", m_simparam.repart_freq);
    if( m_simparam.redecomp_type == "zoltan" )
    {
      m_zz->Set_Param( "PHG_OUTPUT_LEVEL", "0" );
      m_zz->Set_Param( "PHG_REFINEMENT_QUALITY", "2" );
      m_zz->Set_Param( "PHG_REPART_MULTIPLIER", buffer );
    }
    else if ( m_simparam.redecomp_type == "parmetis" )
    {
      m_zz->Set_Param("PARMETIS_OUTPUT_LEVEL", "0");
      m_zz->Set_Param("CHECK_GRAPH", "1" );
      m_zz->Set_Param("GRAPH_PACKAGE", "PARMETIS");
    }
    else
      throw std::runtime_error( "error redecomp type." );

    // graph repartitioning
    m_zz->Set_Num_Obj_Fn( MigrateZObjects::get_num_local_objs, 
      &m_zob );

    m_zz->Set_Obj_List_Fn( MigrateZObjects::get_obj_list, 
      &m_zob );

    m_zz->Set_Num_Edges_Multi_Fn( MigrateZObjects::get_num_edges_multi, 
      &m_zob );

    m_zz->Set_Edge_List_Multi_Fn( MigrateZObjects::get_edge_list_multi, 
      &m_zob);

    // TODO: is this being used anywhere??
    m_zz->Set_Part_Multi_Fn(MigrateZObjects::get_part_multi, 
      &m_zob );

    // manual/auto-migration
    m_zz->Set_Obj_Size_Multi_Fn( MigrateZObjects::obj_size_multi, 
      &m_zob );

    m_zz->Set_Pack_Obj_Multi_Fn( MigrateZObjects::pack_obj_multi, 
      &m_zob );

    m_zz->Set_Unpack_Obj_Multi_Fn( MigrateZObjects::unpack_obj_multi, 
      &m_zob );

    m_zz->Set_Post_Migrate_PP_Fn( MigrateZObjects::post_migrate_pp, 
      &m_zob );
  }

  // +---------------------------+
  // | TaylorGreen::init_paragon |
  // +---------------------------+
  void TaylorGreen::init_paragon()
  {
    tracer::scope _( "TaylorGreen::init_paragon()" );

    m_zz.reset( new Zoltan( MPI_COMM_WORLD ) );

    m_zz->Set_Param( "DEBUG_LEVEL", m_simparam.zoltan_debug_level );
    m_zz->Set_Param( "LB_METHOD", "GRAPH" );
    m_zz->Set_Param( "NUM_GID_ENTRIES", "1" );
    m_zz->Set_Param( "NUM_LID_ENTRIES", "1" );
    m_zz->Set_Param( "OBJ_WEIGHT_DIM", "1" );
    m_zz->Set_Param( "EDGE_WEIGHT_DIM", "1" );

    m_rankCommCost = hwloc_get_rank_comm_costsv2(MPI_COMM_WORLD);
    ParagonInit( &m_paragon, 
      MPI_COMM_WORLD, 
      m_simparam.repart_freq, 
      m_simparam.degree_refine_parallelism, 
      m_simparam.shuffle_refine_times, 
      m_simparam.imbalance_tolerance, 
      m_simparam.degree_resource_contention, 
      m_rankCommCost ); 

    // manual/auto-migration
    m_zz->Set_Obj_Size_Multi_Fn( MigrateZObjects::obj_size_multi, 
      &m_zob );

    m_zz->Set_Pack_Obj_Multi_Fn( MigrateZObjects::pack_obj_multi, 
      &m_zob );

    m_zz->Set_Unpack_Obj_Multi_Fn( MigrateZObjects::unpack_obj_multi, 
      &m_zob );

    m_zz->Set_Post_Migrate_PP_Fn( MigrateZObjects::post_migrate_pp, 
      &m_zob );

  } 

  // +--------------------------+
  // | TaylorGreen::repartition |
  // +--------------------------+
  bool TaylorGreen::repartition()
  {
    tracer::scope _( "TaylorGreen::repartition" );

    real_t st = MPI_Wtime();

    // update weights
    m_cell_wgt.resize( m_subdomain.nhomes() );
    m_mc.cell_wgts( m_cell_wgt );
    m_zob.set_obj_wgts( m_cell_wgt );

    // Repartitioning Phase
    int changes;
    int numGidEntries;
    int numLidEntries;
    int numImport;
    ZOLTAN_ID_PTR importGlobalIds;
    ZOLTAN_ID_PTR importLocalIds;
    int *importProcs;
    int *importToPart;
    int numExport;
    ZOLTAN_ID_PTR exportGlobalIds;
    ZOLTAN_ID_PTR exportLocalIds;
    int *exportProcs;
    int *exportToPart;

    int rc = m_zz->LB_Partition(
      changes,            /* 1 if partitioning was changed, 0 otherwise */
      numGidEntries,      /* Number of integers used for a global ID */
      numLidEntries,      /* Number of integers used for a local ID */
      numImport,          /* Number of vertices to be sent to me */
      importGlobalIds,    /* Global IDs of vertices to be sent to me */
      importLocalIds,     /* Local IDs of vertices to be sent to me */
      importProcs,        /* Process rank for source of each incoming vertex */
      importToPart,       /* New partition for each incoming vertex */
      numExport,          /* Number of vertices I must send to other processes*/
      exportGlobalIds,    /* Global IDs of the vertices I must send */
      exportLocalIds,     /* Local IDs of the vertices I must send */
      exportProcs,        /* Process to which I send each of the vertices */
      exportToPart );     /* Partition to which each vertex will belong */

    if ( rc != ZOLTAN_OK )
      throw std::runtime_error( "Zoltan LB_Partition failed." );

    rc = m_zz->LB_Free_Part( &exportGlobalIds, &exportLocalIds, &exportProcs,
      &exportToPart );
    rc = m_zz->LB_Free_Part( &importGlobalIds, &importLocalIds, &importProcs,
      &importToPart );

    if ( changes ) {
      m_zob.update_neighbors();

      // now that neighboring info has been rebuilt, particle_p2p needs reset
      m_zob.reset_p2p();
    }
    std::cout << m_subdomain.nhomes() << std::endl;

    m_repart_time = MPI_Wtime() - st;
    MPI_Allreduce( MPI_IN_PLACE, &m_repart_time, 1, mpi_real, MPI_MIN, MPI_COMM_WORLD );

    return bool( changes );
  }

  // +---------------------------------+
  // | TaylorGreen::paragon_refinement |
  // +---------------------------------+
  bool TaylorGreen::paragon_refinement()
  {
    tracer::scope _( "TaylorGreen::paragon_refinement" );

    if ( m_simparam.degree_refine_parallelism <= 0 ) 
      return 0; 

    real_t st = MPI_Wtime();

    // update weights
    m_cell_wgt.resize( m_subdomain.nhomes() );
    m_mc.cell_wgts( m_cell_wgt );
    m_zob.set_obj_wgts( m_cell_wgt );

    int changes;
    int numGidEntries;
    int numLidEntries;
    int *partVector = NULL;
    int numExport;
    int *exportProcs = NULL;
    int *exportToPart = NULL;
    ZOLTAN_ID_PTR exportGlobalIds;
    ZOLTAN_ID_PTR exportLocalIds;

    GraphStruct localGraph;
    build_local_graph( &localGraph );

    // Architecture-Aware Refinement Phase
    ParagonRefinement(
      &m_paragon, 
      localGraph, 
      &partVector, 
      &changes,
      &numExport,
      &exportLocalIds,
      &exportGlobalIds,
      &exportToPart );
    free(partVector);

    std::cout << changes << " " << numExport << std::endl;

    if ( changes )
    {

      /* We must invert the lists here.  The import lists are used during the 
       * migration phase to update the distributed directory.
       */
      int numImport;
      ZOLTAN_ID_PTR importGlobalIds;
      ZOLTAN_ID_PTR importLocalIds;
      int *importProcs;
      int *importToPart;
      
      m_zz->Invert_Lists( numExport, exportGlobalIds, exportLocalIds, 
			  exportToPart, exportToPart, numImport, importGlobalIds,
			  importLocalIds, importProcs, importToPart );

      std::cout << m_subdomain.my_rank() << " " << numImport 
		<< " " << numExport << std::endl;

      // Manual Migration Phase
      m_zz->Migrate( numImport,
		     importGlobalIds,
		     importLocalIds,
		     importProcs,
		     importToPart,
		     numExport,
		     exportGlobalIds,
		     exportLocalIds,
		     exportToPart,      //rank id = part number
		     exportToPart );

      m_zob.update_neighbors();
      // now that neighboring info has been rebuilt, particle_p2p needs reset
      m_zob.reset_p2p();    

    }

    free( exportGlobalIds );
    free( exportLocalIds );
    free( exportToPart );
    graphDeinit( &localGraph );

    m_repart_time = MPI_Wtime() - st;
    MPI_Allreduce( MPI_IN_PLACE, &m_repart_time, 1, mpi_real, MPI_MIN, MPI_COMM_WORLD );

    return bool( changes );
  }

  // +--------------------------------+
  // | TaylorGreen::build_local_graph |
  // +--------------------------------+
  void TaylorGreen::build_local_graph( GraphStruct *localGraph )
  {
    tracer::scope _( "TaylorGreen::build_local_graph" );

    int i;
    int error;
    int nVertex = m_zob.get_num_local_objs( &m_zob, &error );
    graphInitVertexLists( localGraph, nVertex );
    
    m_zob.get_obj_list( &m_zob, 1, 1, localGraph->vertexGIDs, 
      localGraph->vertexLIDs, 1, localGraph->vertexWgts, &error );
    
    m_zob.get_part_multi( &m_zob, 1, 1, nVertex, localGraph->vertexGIDs, 
      localGraph->vertexLIDs, localGraph->partVector, &error );
    
    m_zob.obj_size_multi( &m_zob, 1, 1, nVertex, localGraph->vertexGIDs, 
      localGraph->vertexLIDs, localGraph->vertexSize, &error );
    
    m_zob.get_num_edges_multi( &m_zob, 1, 1, nVertex, localGraph->vertexGIDs, 
      localGraph->vertexLIDs, localGraph->nborIndex, &error );
    
    int numNbors = 0;
    for( i=0; i<nVertex; i++)
      numNbors += localGraph->nborIndex[i];
    
    graphInitEdgeLists( localGraph, numNbors );

    m_zob.get_edge_list_multi( &m_zob, 1, 1, nVertex, localGraph->vertexGIDs, 
      localGraph->vertexLIDs, localGraph->nborIndex, localGraph->nborGIDs, 
      localGraph->nborParts, 1, localGraph->edgeWgts, &error );
    
    for(i=1; i<=nVertex; i++)
      localGraph->nborIndex[i] += localGraph->nborIndex[i-1];

    localGraph->nborIndex[0] = 0;
  }

}
