#
# Developer docmentation using Doxygen
#

# Find Doxygen
find_package(Doxygen REQUIRED)
#  this module will define
#     DOXYGEN_EXECUTABLE     = The path to the doxygen command.
#     DOXYGEN_FOUND          = Was Doxygen found or not?
#     DOXYGEN_DOT_EXECUTABLE = The path to the dot program used by doxygen.
#     DOXYGEN_DOT_FOUND      = Was Dot found or not?
#     DOXYGEN_DOT_PATH       = The path to dot not including the executable
#   Some of these variables are used in Doxyfile, and replaced for build

# Set configuration variables
#  SRCROOT will be stripped from filenames
set( SRCROOT ${PROJECT_SOURCE_DIR}/src )
#set( ENV{PATH} "c:/Program Files/gs/gs8.63/bin";ENV{PATH} )

# Configure
configure_file( Doxyfile . @ONLY ) 
#configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY ) 

# Create documentation target
add_custom_target(doc ${DOXYGEN_EXECUTABLE})
