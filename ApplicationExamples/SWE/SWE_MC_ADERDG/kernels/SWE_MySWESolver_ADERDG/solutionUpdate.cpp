
// update the elements 

#include "kernels/SWE_MySWESolver_ADERDG/Kernels.h"
#include "kernels/SWE_MySWESolver_ADERDG/Quadrature.h"

void SWE::MySWESolver_ADERDG_kernels::aderdg::solutionUpdate( 
  double* restrict luh,
  const double* restrict const luhOld,
  const double* restrict const lduh, 
  const double dt
) {
#ifdef __INTEL_COMPILER
  __assume_aligned(iweights3, ALIGNMENT);
  __assume_aligned(luh,       ALIGNMENT); //luh    should be aligned, see Solver.h
  __assume_aligned(luhOld,    ALIGNMENT); //luhOld should be aligned, see Solver.h
  __assume_aligned(lduh,      ALIGNMENT); //lduh   should be aligned, see Solver.h
#endif

  for (int xyz = 0; xyz < 4; xyz++) {
    const double coeff = dt*iweights3[xyz];
    #pragma omp simd aligned(luh,luhOld,lduh:ALIGNMENT)
    for (int n = 0; n < 4; n++) { //update only the variables, lduh contains no parameters
      luh[xyz*4+n] = luhOld[xyz*4+n] + coeff*lduh[xyz*4+n]; //simd+fma
    }
  }
 
}