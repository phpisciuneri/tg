# +-----------------+
# | PROJECT IPLMCFD |
# +-----------------+
project( IPLMCFD CXX C )
cmake_minimum_required( VERSION 2.8.8 )

# +------------------------------+
# | CHECK NATIVE FORTRAN SUPPORT |
# +------------------------------+
include( CheckLanguage )
check_language( Fortran )
if( CMAKE_Fortran_COMPILER )
  enable_language( Fortran )
endif()

# +---------+
# | OPTIONS |
# +---------+
option( IPLMCFD_HANDLE_ASSERTIONS "Catch and cleanup on assertion errors instead of execution abort" ON )
option( IPLMCFD_UNIT_TESTS        "Enable unit tests?" OFF )
option( IPLMCFD_DEV_DOX           "Build developer Doxygen documentation?" OFF )

# +----------+
# | PACKAGES |
# +----------+

# where to look for custom cmake modules
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/site-macros/Modules)

find_package( Boost 1.39 REQUIRED COMPONENTS filesystem regex system date_time program_options )
find_package( MPI REQUIRED )
find_package( Zoltan REQUIRED )
find_package( hwloc REQUIRED )
find_package( SIGNALHANDLING )
find_package( Subversion )

set( INCS ${INCS}
  ${MPI_C_INCLUDE_PATH}
  ${Boost_INCLUDE_DIRS}
  ${Zoltan_INCLUDE_DIRS}
  ${hwloc_INCLUDE_DIRS} )

set( LIBS ${LIBS} 
  ${Boost_LIBRARIES}
  ${MPI_C_LIBRARIES}
  ${Zoltan_LIBRARIES}
  ${hwloc_LIBRARIES} )

if (WIN32)
   # otherwise boost auto link fails eventho ${Boost_LIBRARIES} is using full path names
   link_directories( "${Boost_LIBRARY_DIRS}" )
   set( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -Dunix" )
   set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Dunix" )
endif (WIN32)

# +-------------------+
# | EXTERNAL PROJECTS |
# +-------------------+
add_subdirectory( ext/chemkinpp )
add_subdirectory( ext/metis )
add_subdirectory( ext/paragon )

# +---------------------+
# | SOURCES AND TARGETS |
# +---------------------+

# Add all source files
#file( GLOB src src/*.cpp src/ipfd/*.cpp src/iplmc/*.cpp )
file( GLOB src
	src/cell.cpp
	src/domain.cpp
	#src/iplmc/ensemble.cpp
	src/iplmc/iplmc.cpp
	src/iplmc/position_step.cpp
	src/iplmc/scalar_average.cpp
	src/iplmc/scalar_step.cpp
	src/migratezobjects.cpp
	src/parallel_io.cpp
	src/particle_p2p.cpp
	src/pattern.cpp
	src/simparam.cpp
	src/subdomain.cpp
	src/taylor_green.cpp
	src/taylor_green_vortex.cpp
	)

# Add all header files
#file( GLOB hdr src/*.hpp src/ipfd/*.hpp src/iplmc/*.hpp )
file( GLOB hdr 
	src/cpuloadheader.hpp
	src/cell.hpp
	src/defs.hpp
	src/domain.hpp
	src/field.hpp
	src/iplmc.hpp
	src/migratezobjects.hpp
	src/parallel_io.hpp
	src/particle.hpp
	src/particle_mpi_type.hpp
	src/particle_p2p.hpp
	src/pattern.hpp
	src/scalar_averages.hpp
	src/simparam.hpp
	src/subdomain.hpp
	src/taylor_green.hpp
	src/taylor_green_vortex.hpp
	src/zobjects.hpp
	)
	
add_library( diagnosticlib
	src/diagnostic/activity_log.cpp
	src/diagnostic/activity_log.hpp
	src/diagnostic/assert_handler.cpp
	src/diagnostic/binary_cstdio.hpp
	src/diagnostic/binary_fstream.hpp
	src/diagnostic/binary_stream.hpp
	src/diagnostic/binary_stream.cpp
	src/diagnostic/cpublock.cpp
	src/diagnostic/cpublock.hpp
	src/diagnostic/get_time.hpp
	src/diagnostic/get_time.cpp
	src/diagnostic/signal_handler.hpp
	src/diagnostic/signal_handler.cpp
	src/diagnostic/tracer.hpp
	src/diagnostic/tracer.cpp
	src/diagnostic/version.hpp
	src/diagnostic/version.cpp
	)
  
# remove stuff
#foreach(x main program_monitor)
#  list(REMOVE_ITEM src ${PROJECT_SOURCE_DIR}/src/${x}.cpp)
#endforeach(x)

# create iplmcfd library
add_library( iplmcfdlib ${src} ${hdr} )

# add niceties for MS Visual Studio
if( WIN32 )
  add_definitions( -D_CRT_SECURE_NO_WARNINGS )
  source_group( src REGULAR_EXPRESSION "\\.cpp$" )
  source_group( hpp REGULAR_EXPRESSION "\\.hpp$" )
  source_group( iplmc-src "src/iplmc/.*\\.cpp" )  
  source_group( ipfd-src "src/ipfd/.*\\.cpp" )  
endif( WIN32)

# add executable
add_executable( tg src/main.cpp )

# +------------------+
# | VERSION AND HOST |
# +------------------+
if(Subversion_FOUND)
  Subversion_WC_INFO(${PROJECT_SOURCE_DIR} "IPLMCFD")
  message(STATUS "Current revision is ${IPLMCFD_WC_REVISION}")
ELSE( Subversion_FOUND )
  set( _WC_REVISION "unknown" )
ENDIF(Subversion_FOUND)

# site-macros/Modules/hostname.cmake
include( hostname )

# +-------------------------------+
# | LINK AND COMPILATION SETTINGS |
# +-------------------------------+
include_directories( 
  ${INCS}
  ext/tvmet/include
  ext/chemkinpp
  ext/chemkinpp/ext/mytools/include
  ext/metis/include
  ext/paragon/src
  )

target_link_libraries( tg 
  iplmcfdlib
  diagnosticlib
  chemkinpp
  metis
  paragon
  ${LIBS} )
  
  # test executables
#  enable_testing()
#  add_executable( deriv-f1 src/test/deriv-f1.cpp )
#  target_link_libraries( deriv-f1 
#	iplmcfdlib
#	metis
#	${Boost_LIBRARIES}
#	${MPI_C_LIBRARIES} )
#  add_test( DERIV-F1 deriv-f1 ${PROJECT_SOURCE_DIR}/src/test/simparam.in )
#  
#  add_executable( deriv-f4 src/test/deriv-f4.cpp )
#  target_link_libraries( deriv-f4 
#	iplmcfdlib
#	metis
#	${Boost_LIBRARIES}
#	${MPI_C_LIBRARIES} )
#  add_test( DERIV-F4 deriv-f4 ${PROJECT_SOURCE_DIR}/src/test/simparam.in )
#  
#  add_executable( deriv-b1 src/test/deriv-b1.cpp )
#  target_link_libraries( deriv-b1 
#	iplmcfdlib
#	metis
#	${Boost_LIBRARIES}
#	${MPI_C_LIBRARIES} )
#  add_test( DERIV-B1 deriv-b1 ${PROJECT_SOURCE_DIR}/src/test/simparam.in )
#  
#  add_executable( deriv-b4 src/test/deriv-b4.cpp )
#  target_link_libraries( deriv-b4 
#	iplmcfdlib
#	metis
#	${Boost_LIBRARIES}
#	${MPI_C_LIBRARIES} )
#  add_test( DERIV-B4 deriv-b4 ${PROJECT_SOURCE_DIR}/src/test/simparam.in )
#  
#  add_executable( deriv-c2 src/test/deriv-c2.cpp )
#  target_link_libraries( deriv-c2 
#	iplmcfdlib
#	metis
#	${Boost_LIBRARIES}
#	${MPI_C_LIBRARIES} )
#  add_test( DERIV-C2 deriv-c2 ${PROJECT_SOURCE_DIR}/src/test/simparam.in )

# Some MPICH implementations require this. 
# Should not hurt setting it on other implementations.
set_property( 
  SOURCE ${src} 
  APPEND PROPERTY 
  COMPILE_DEFINITIONS MPICH_IGNORE_CXX_SEEK )


set_property(
  SOURCE ${src}
  APPEND PROPERTY 
  COMPILE_FLAGS ${MPI_COMPILE_FLAGS} )

if( CMAKE_SYSTEM_NAME MATCHES "Catamount" )
  add_definitions( -DLMCFD_INEFFICIENT_FSTREAM )
endif( CMAKE_SYSTEM_NAME MATCHES "Catamount" )

set_property(
  SOURCE src/version.cpp
  APPEND
  PROPERTY COMPILE_DEFINITIONS LMCFD_REVISION=${IPLMCFD_WC_REVISION} )

set_property( 
  SOURCE src/info.cpp
  APPEND
  PROPERTY COMPILE_DEFINITIONS LMCFD_HOSTNAME=${HOSTNAME} )

if ( SIGNALHANDLING_FOUND )
  add_definitions(-DLMCFD_ENABLE_SIGNAL_HANDLING )
endif ( SIGNALHANDLING_FOUND )

if ( IPLMCFD_HANDLE_ASSERTIONS )
  add_definitions( -DBOOST_ENABLE_ASSERT_HANDLER )
endif ( IPLMCFD_HANDLE_ASSERTIONS )

# Prevent bogus warning on VS, harmless on other systems
add_definitions( -D_SCL_SECURE_NO_WARNINGS )

# We need the following to link Boost Program Options on 
# windows when our environment is MSVC and MinGW Fortran
if( WIN32 AND ( NOT CMAKE_Fortran_COMPILER ) )
  add_definitions( -DBOOST_ALL_DYN_LINK )
endif()

# +------+
# | MISC |
# +------+

# site-macros/Modules/utils.cmake
include( utils )
set_build_type( Release )

# sample input
configure_file( test/taylor_green/simparam.in . COPY_ONLY )
configure_file( test/taylor_green/therm.dat   . COPY_ONLY )
configure_file( test/taylor_green/tran.dat    . COPY_ONLY )

# Testing
if ( IPLMCFD_UNIT_TESTS )
  enable_testing()
  add_subdirectory( src/test )
endif ( IPLMCFD_UNIT_TESTS ) 

# Documentation
if ( IPLMCFD_DEV_DOX )
  add_subdirectory( doc )
endif ( IPLMCFD_DEV_DOX ) 

# Other utilities
option( IPLMCFD_POSTPROCESSING_TOOLS "Build post-processing utilities?" OFF )
if ( IPLMCFD_POSTPROCESSING_TOOLS )
  add_subdirectory( post/ppp )
endif()

# See ${PROJECT_SOURCE_DIR}/ext/tvmet/include/tvmet/local/README.txt
# for an explanation of what is going on here.
if( CMAKE_CXX_COMPILER_ID MATCHES "Intel" )

  # use intel configuration and default CommaInitializer.h
  configure_file( ${PROJECT_SOURCE_DIR}/ext/tvmet/include/tvmet/local/config_intel.h 
    ${PROJECT_SOURCE_DIR}/ext/tvmet/include/tvmet/config.h COPY_ONLY )
  configure_file( ${PROJECT_SOURCE_DIR}/ext/tvmet/include/tvmet/local/CommaInitializer.h 
    ${PROJECT_SOURCE_DIR}/ext/tvmet/include/tvmet/CommaInitializer.h  COPY_ONLY )
	
elseif( CMAKE_CXX_COMPILER_ID MATCHES "MSVC" )

  # use windows vs configuration and custom CommaInitializer.h
  configure_file( ${PROJECT_SOURCE_DIR}/ext/tvmet/include/tvmet/local/config_VS.h 
    ${PROJECT_SOURCE_DIR}/ext/tvmet/include/tvmet/config.h COPY_ONLY )
  configure_file( ${PROJECT_SOURCE_DIR}/ext/tvmet/include/tvmet/local/CommaInitializer_VS.h 
    ${PROJECT_SOURCE_DIR}/ext/tvmet/include/tvmet/CommaInitializer.h  COPY_ONLY )
  
else()

  # else use GNU and default CommaInitializer.h
  configure_file( ${PROJECT_SOURCE_DIR}/ext/tvmet/include/tvmet/local/config_gnu.h 
    ${PROJECT_SOURCE_DIR}/ext/tvmet/include/tvmet/config.h COPY_ONLY )
  configure_file( ${PROJECT_SOURCE_DIR}/ext/tvmet/include/tvmet/local/CommaInitializer.h 
    ${PROJECT_SOURCE_DIR}/ext/tvmet/include/tvmet/CommaInitializer.h  COPY_ONLY )
	
endif( CMAKE_CXX_COMPILER_ID MATCHES "Intel" )

# +------------+
# | POST-BUILD |
# +------------+

# For msvc, with no integrated Fortran support (use gfortran)
# all libraries are shared.  Therefore, all .dlls
# must be in the path, or in the same directory of the
# executables.  Below we copy the generated .dlls to the
# directory of the executables.
IF( WIN32 AND ( NOT CMAKE_Fortran_COMPILER ) )

  # copy mytools dll
  ADD_CUSTOM_COMMAND(
    TARGET iplmcfd
	POST_BUILD
	COMMAND ${CMAKE_COMMAND} ARGS -E copy 
      ${PROJECT_BINARY_DIR}/ext/chemkinpp/ext/mytools/src/fortran/libmytools.dll
      ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR} )
	  
  # copy chemkin dll
  ADD_CUSTOM_COMMAND(
    TARGET iplmcfd
	POST_BUILD
	COMMAND ${CMAKE_COMMAND} ARGS -E copy 
      ${PROJECT_BINARY_DIR}/ext/chemkinpp/ext/chemkin/libchemkin.dll
      ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR} )
	  
  # copy odepack dll
  ADD_CUSTOM_COMMAND(
    TARGET iplmcfd
	POST_BUILD
	COMMAND ${CMAKE_COMMAND} ARGS -E copy 
      ${PROJECT_BINARY_DIR}/ext/chemkinpp/ext/odepack/libodepack.dll
      ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR} )
	  
endif( WIN32 AND ( NOT CMAKE_Fortran_COMPILER ) )
