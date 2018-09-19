module load intel/analyze/2018.2

export TBB_SHLIB="-L/ddn/apps/Cluster-Apps/intel/xe_2018.2/tbb/lib/intel64/gcc4.7 -ltbb_debug"
export COMPILER_CFLAGS=${COMPILER_CFLAGS}:"-debug all -DTBB_USE_THREADING_TOOLS"
