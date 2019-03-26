
#include "kernels/SWE_MySWESolver_ADERDG/Kernels.h"
#include "kernels/SWE_MySWESolver_ADERDG/DGMatrices.h"
#include "kernels/SWE_MySWESolver_ADERDG/Quadrature.h"


void SWE::MySWESolver_ADERDG_kernels::aderdg::faceIntegral(
                              double *lduh, 
                              const double *const lFhbnd,
                              const int direction, 
                              const int orientation,
                              const double inverseDxDirection
) {

  const double* const FCoeff = (orientation == 0 ? FLCoeff : FRCoeff);
  
#ifdef __INTEL_COMPILER
  __assume_aligned(FRCoeff,  ALIGNMENT);
  __assume_aligned(FLCoeff,  ALIGNMENT);
  __assume_aligned(FCoeff,   ALIGNMENT);
  __assume_aligned(weights2, ALIGNMENT);
  __assume_aligned(lFhbnd,   ALIGNMENT);
  __assume_aligned(lduh,     ALIGNMENT);
#endif
  
  const double scaling = (2.0 * orientation - 1.0) * inverseDxDirection;
  
  switch (direction) {
    case 0:   // x faces, left-right flux
      for (int zy = 0; zy < 2; zy++) { // zy
        const double scaledWeight = scaling * weights2[zy];
        for (int x = 0; x < 2; x++) { // x
          #pragma omp simd aligned(lduh,FCoeff,lFhbnd:ALIGNMENT)
          for (int n = 0; n < 4; n++) {
            lduh[(zy*2+x)*4+n] -= scaledWeight * FCoeff[x] * lFhbnd[zy*4+n];
          }
        }
      }
      break;
    case 1: // y faces
      for (int z = 0; z < 1; z++) { // z
      for (int x = 0; x < 2; x++) { // x
        const double scaledWeight = scaling * weights2[z*2+x];
        for (int y = 0; y < 2; y++) { // y
          #pragma omp simd aligned(lduh,FCoeff,lFhbnd:ALIGNMENT)
          for (int n = 0; n < 4; n++) {
            lduh[((z*2+y)*2+x)*4+n] -= scaledWeight * FCoeff[y] * lFhbnd[(z*2+x)*4+n];
          }
        }
      }
      }
      break;
    case 2: // z faces
    
      for (int yx = 0; yx < 4; yx++) { // yx
        const double scaledWeight = scaling * weights2[yx];
        for (int z = 0; z < 2; z++) { // x
          #pragma omp simd aligned(lduh,FCoeff,lFhbnd:ALIGNMENT)
          for (int n = 0; n < 4; n++) {
            lduh[(z*4+yx)*4+n] -= scaledWeight * FCoeff[z] * lFhbnd[yx*4+n];
          }
        }
      }
      break;
  }
}