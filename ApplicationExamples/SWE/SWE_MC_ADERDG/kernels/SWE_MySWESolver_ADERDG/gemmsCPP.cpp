
#include "kernels/SWE_MySWESolver_ADERDG/gemmsCPP.h"
#include "kernels/SWE_MySWESolver_ADERDG/Kernels.h" //for the libxsmm flop counter

// amrRoutines gemms
#include "kernels/SWE_MySWESolver_ADERDG/asm_amrRoutines.c"


// fusedSpaceTimePredictorVolumeIntegral gemms
#include "kernels/SWE_MySWESolver_ADERDG/asm_fstpvi.c"

// Limiter gemms
#include "kernels/SWE_MySWESolver_ADERDG/asm_limiter.c"
