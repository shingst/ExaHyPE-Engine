source /etc/profile.d/modules.sh

module unload gcc
module switch python/3.5_intel
module load gcc/4.9
module switch java/1.8
module switch intel/17.0 
module switch tbb/2017
module switch mpi.intel/2017

export EXAHYPE_CC="mpicc -DnoParallelExchangePackedRecordsAtBoundary -DnoParallelExchangePackedRecordsBetweenMasterAndWorker -DnoParallelExchangePackedRecordsInHeaps -DnoParallelExchangePackedRecordsThroughoutJoinsAndForks"

export MODE=Release
export COMPILER=Intel
export DISTRIBUTEDMEM=MPI
export GPROF=off
export USE_IPO=off
