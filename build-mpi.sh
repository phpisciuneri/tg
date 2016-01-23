# source this script to build
# $ . build-frank.sh

# set environment
module purge
spack load boost
spack load hwloc%intel

# uncomment the next two lines to use mvapich2
#spack load mvapich2
#spack load trilinos^mvapich2

# uncomment the next two lines to use mpich
#spack load mpich
#module load trilinos@12.0.1%intel@15.0.3-mj35adm

# uncomment the next two lines to use openmpi
spack load openmpi@1.8.6
module load trilinos@12.0.1%intel@15.0.3-736pgl4

spack load cmake@3.3.0

# create build directory
mkdir -p build

cd build

# configure and generate makefile
CC=icc CXX=icpc FC=ifort cmake \
    -DCMAKE_BUILD_TYPE=release \
    -DCHEMKINPP_IMPLICIT_BLAS=off \
    ..



