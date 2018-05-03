export I_MPI_FABRICS="tmi"
export COMPILER_CFLAGS="-g -DnoParallelExchangePackedRecordsAtBoundary -DnoParallelExchangePackedRecordsBetweenMasterAndWorker -DnoParallelExchangePackedRecordsInHeaps -DnoParallelExchangePackedRecordsThroughoutJoinsAndForks"

#export EXAHYPE_CC="scorep --noonline-access --nocompiler --mpp=mpi --thread=none --user mpiicpc -DUseScoreP"

module load java/1.8
module load python/3.5_intel
module unload intel
module load intel/17.0
module unload mpi.intel
module load mpi.intel/2017
module unload gcc
module load gcc/5
module unload tbb
module load tbb/2017

