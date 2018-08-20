module switch python/3.5_intel
module switch java/1.8
module switch intel/18.0
module switch tbb/2018
module switch gcc/5

module load git

module load ddt

export COMPILER=Intel
export MODE=Release
export EXAHYPE_CC=mpicc
export EXAHYPE_FC=ifort
export DISTRIBUTEDMEM=MPI
export COMPILER_CFLAGS=" -DnoMultipleThreadsMayTriggerMPICalls -DMPIProgressionReliesOnMPITest -DUsePeanosAggregationBoundaryExchanger -DUsePeanosAggregationBoundaryExchangerForMetaData "
export PROJECT_CFLAGS=" -Ilib/ -g -DTBB_USE_DEBUG -O0 "
export PROJECT_LFLAGS=" -Llib/ -ltecio "
export USE_IPO=on
export SHAREDMEM=TBB
