export COMPILER=Intel
export MODE=Release
export SHAREDMEM=TBB
export DISTRIBUTEDMEM=MPI
export OMP_NUM_THREADS=28
export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS
export I_MPI_PIN=on
export I_MPI_PIN_DOMAIN=omp:compact     #   i.e. =OMP_NUM_THREADS
