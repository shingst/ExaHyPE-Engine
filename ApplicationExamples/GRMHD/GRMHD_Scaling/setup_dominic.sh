# helper script for Dominic; will be removed later
mkdir -p output
rm -r *.o cipofiles.mk cfiles.mk ffiles.mk kernels  2> /dev/null
../../../Toolkit/toolkit.sh -d -s --format=json GRMHD_Limiting_dominic.exahype2 --mesh_info
ln -sf ../../../ExternalLibraries/TOVSolver/src/ tovsolver_lib
