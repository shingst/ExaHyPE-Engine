#!/bin/bash

cd $(dirname "$0")

# note that this Application uses the MakefileFORTRAN which does not
# yet use the COMPILER=GNU.

export CC=gcc
export SHAREDMEM="TBB" # None
export TBB_INC=/usr/include/tbb
MPI_LDFLAGS="$(mpicc -showme:link)"
export TBB_SHLIB="-L/usr/lib -ltbb $MPI_LDFLAGS"

# Debugging mode
export MODE="DEBUG"

set -e

cd ../../

# Vasco: do not run toolkit on SRHD right now, as this is the new Fortran prototype.

#java -jar ExaHyPE.jar  --not-interactive Applications/srhd3dfortran.exahype

cd -

#make clean

# broken build system, do this by hand:
gfortran -c Parameters.f90 

make -j $(nproc)

