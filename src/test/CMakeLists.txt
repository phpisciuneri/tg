#
# Unit tests for IPLMCFD modules/classes/overall
#

#  +--------------+
#  |  Test Lists  |
#  +--------------+

set( tests  mpiio ) #subdomain ipfd_intersolver )
set( testsrc  test-mpiio.cpp )
#test-subdomain.cpp )
#test-ipfd_intersolver.cpp )

include_directories(${PROJECT_SOURCE_DIR}/src)

#  +------------+
#  |  Packages  |
#  +------------+

find_package( Boost 1.34 REQUIRED  )
include_directories(${Boost_INCLUDE_DIRS})
add_definitions( -DBOOST_TEST_NO_MAIN )
remove_definitions( -DBOOST_ENABLE_ASSERT_HANDLER )

# MPI header is included throught the tested units
find_package( MPI REQUIRED )
include_directories( ${MPI_INCLUDE_DIRS} )

#  +----------+
#  |  driver  |
#  +----------+

add_executable( serial_test_driver serial_test_driver.cpp ${testsrc})   
# even though the serial driver don't use MPI, some libraries 
# require linkage when the <mpi.h> header is included 
target_link_libraries( serial_test_driver ${MPI_LIBRARIES} ipfdlib)


#  +------------+
#  |  Add tests |
#  +------------+
foreach( x ${tests} )
   add_test(${x} serial_test_driver --run_test=${x})
endforeach( x )

