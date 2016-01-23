# source this script to build
# $ . build-frank.sh

# set environment
module purge
module load env/intel-2011.1-openmpi-1.4
module load parmetis/4.0.2-intel12-openmpi
module load zoltan/3.8.1-intel12
module load boost/1.48.0-intel12-openmpi
module load cmake/2.8.11.2

# create build directory
mkdir -p build
rm -rf build
mkdir -p build
cp submit.qsub build/
cd build

# configure and generate makefile
cmake -C ../site-macros/caches/intel.cmake \
    -DCMAKE_BUILD_TYPE=release \
    -DIPLMCFD_IMPLICIT_MPI=off \
    -DCHEMKINPP_IMPLICIT_BLAS=off \
    ..



