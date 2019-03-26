
#include "kernels/SWE_MySWESolver_ADERDG/Kernels.h"
#include "kernels/SWE_MySWESolver_ADERDG/DGMatrices.h"
#include "kernels/SWE_MySWESolver_ADERDG/Quadrature.h"


void SWE::MySWESolver_ADERDG_kernels::aderdg::surfaceIntegral( 
  double* restrict lduh, 
  const double* restrict const lFhbnd, 
  const double inverseDx //Assume dx[0] == dx[1] == dx[2]
) {

#ifdef __INTEL_COMPILER
  __assume_aligned(FRCoeff,  ALIGNMENT);
  __assume_aligned(FLCoeff,  ALIGNMENT);
  __assume_aligned(weights2, ALIGNMENT);
  __assume_aligned(lFhbnd,   ALIGNMENT);
  __assume_aligned(lduh,     ALIGNMENT);
#endif

  // x faces
  for (int yz = 0; yz < 2; yz++) {
    const double weight = weights2[yz] * inverseDx; //Assume dx[0] == dx[1] == dx[2]
    for (int x = 0; x < 2; x++) {
      #pragma omp simd aligned(lduh,lFhbnd:ALIGNMENT)
      for (int n = 0; n < 4; n++) {
        lduh[n+4*(x+2*yz)] -= weight *
            (lFhbnd[n+4*yz+8] * FRCoeff[x] -              lFhbnd[n+4*yz+0] * FLCoeff[x]);
      }
    }
  }

  // y faces
  for (int xz = 0; xz < 2; xz++) {
    const double weight = weights2[xz] * inverseDx; //Assume dx[0] == dx[1] == dx[2]
    const int xzLuhIndex = xz*4;
    for (int y = 0; y < 2; y++) {
      #pragma omp simd aligned(lduh,lFhbnd:ALIGNMENT)
      for (int n = 0; n < 4; n++) {
        lduh[n+xzLuhIndex+y*8] -= weight *
            (lFhbnd[n+4*xz+24] * FRCoeff[y] -
             lFhbnd[n+4*xz+16] * FLCoeff[y]);
        }
      }
  }
}
