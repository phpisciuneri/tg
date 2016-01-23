#ifndef IPLMCFD_MIGRATEZOBJECTS_HPP_
#define IPLMCFD_MIGRATEZOBJECTS_HPP_

#include <vector>
#include <boost/foreach.hpp>
#include <boost/scoped_ptr.hpp>
#include <cassert>

#include "defs.hpp"
#include "zobjects.hpp"
#include "iplmc.hpp"

namespace iplmcfd
{

//!
//! \file migratezobjects.hpp
//! \class MigrateZObjects
//! \brief Helper class for data migration due to repartitioning with Zoltan
//!
class MigrateZObjects : public ZObjects
{

public:

  //!
  //! \brief Constructor for ZMigrateObjects class
  //!
  //! \param[in] s Subdomain reference
  //!
  //! Additionally creates and populates the distributed directory (DD).  The
  //!   DD is used to quickly update neighbor ownership after migration.
  //!
  MigrateZObjects( Subdomain& s, Iplmc& mc );
	
  //!
  //! \brief Sets the size (in bytes) of each migrate-able object
  //!
  //! \param[in] data pointer to user-defined data
  //! \param[in] num_gid_entries number of array entries used to describe a
  //!            global ID
  //! \param[in] num_lid_entries number of array entries used to describe a
  //!            local ID
  //! \param[in] num_ids number of exporting objects whose sizes are to be 
  //!            returned
  //! \param[in] global_ids array of exported global IDs
  //! \param[in] local_ids array of exported local IDs
  //! \param[out] sizes array of sizes (in bytes) for each object in ID list
  //! \param[out] ierr error code
  //!
  static void obj_size_multi( void* data,
    int num_gid_entries,
    int num_lid_entries,
    int num_ids,
    ZOLTAN_ID_PTR global_ids,
    ZOLTAN_ID_PTR local_ids,
    int* sizes,
    int* ierr );

  //!
  //! \brief Packs and removes data to be exported for migration
  //!
  //! \param[in] data pointer to user-defined data
  //! \param[in] num_gid_entries number of array entries used to describe a
  //!            global ID
  //! \param[in] num_lid_entries number of array entries used to describe a
  //!            local ID
  //! \param[in] num_ids number of objects with data to pack
  //! \param[in] global_ids array of global IDs with data to pack
  //! \param[in] local_ids array of local IDs with data to pack
  //! \param[in] dest array of destination part numbers
  //! \param[in] sizes per-object sizes (bytes) of each object in the
  //!            communication buffer
  //! \param[in] idx array of indices into the buf array giving the starting
  //!            location of that objects data
  //! \param[out] buf communication buffer into which objects are packed
  //! \param[out] ierr error code
  //!
  static void pack_obj_multi( void* data,
    int num_gid_entries,
    int num_lid_entries,
    int num_ids,
    ZOLTAN_ID_PTR global_ids,
    ZOLTAN_ID_PTR local_ids,
    int* dest,
    int* sizes,
    int* idx,
    char* buf,
    int* ierr );

  //!
  //! \brief Unpacks and inserts imported data.
  //!
  //! \param[in,out] data pointer to user-defined data
  //! \param[in] num_gid_entries number of array entries used to describe a
  //!            global ID
  //! \param[in] num_ids number of objects with data to import
  //! \param[in] global_ids array of global IDs with data to import
  //! \param[in] sizes per-object sizes (bytes) of each object in the
  //!            communication buffer
  //! \param[in] idx array of indices into the buf array giving the starting
  //!            location of that objects data
  //! \param[in] buf communication buffer into which objects are packed
  //! \param[out] ierr error code
  //!
  static void unpack_obj_multi( void* data,
    int num_gid_entries,
    int num_ids,
    ZOLTAN_ID_PTR global_ids,
    int* sizes,
    int* idx,
    char* buf,
    int* ierr );

  //!
  //! \brief Performs any post-processing at the end of data Migration
  //!
  //! \param[in,out] data pointer to user-defined data
  //! \param[in] num_gid_entries number of array entries used to describe a
  //!            global ID
  //! \param[in] num_ids number of objects with data to import
  //! \param[in] global_ids array of global IDs with data to import
  //! \param[in] sizes per-object sizes (bytes) of each object in the
  //!            communication buffer
  //! \param[in] idx array of indices into the buf array giving the starting
  //!            location of that objects data
  //! \param[in] buf communication buffer into which objects are packed
  //! \param[out] ierr error code
  //!
  static void post_migrate_pp( void* data,
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
    int* ierr );

  void update_neighbors();

  void reset_p2p();

//  int get_m_fix_bytes() { return m_fix_bytes; }
//  int get_m_phi_bytes() { return m_fix_bytes; }

private:

  //!
  Subdomain& m_subdomain;
  //!
  Iplmc& m_mc;
  //!
  boost::scoped_ptr< Zoltan_DD > m_dd;
  //!
  int m_fix_bytes;
  //!
  int m_phi_bytes;

};

} // namespace iplmcfd

#endif // IPLMCFD_MIGRATEZOBJECTS_HPP_
