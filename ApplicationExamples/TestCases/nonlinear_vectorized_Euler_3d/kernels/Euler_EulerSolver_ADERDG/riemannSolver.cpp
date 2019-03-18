
#include <algorithm>
#include <cstring>
#include <cmath>

#include "kernels/Euler_EulerSolver_ADERDG/Kernels.h"
#include "kernels/Euler_EulerSolver_ADERDG/DGMatrices.h"
#include "kernels/Euler_EulerSolver_ADERDG/Quadrature.h"

#include "EulerSolver_ADERDG.h"

void Euler::EulerSolver_ADERDG_kernels::aderdg::riemannSolver( 
  Euler::EulerSolver_ADERDG& solver,
  double* restrict FL, double* restrict FR,
  const double* restrict const QL, const double* restrict const QR,
  const double dt,
  const int direction
) {

#ifdef __INTEL_COMPILER
  __assume_aligned(weights2, ALIGNMENT);
  __assume_aligned(FL, ALIGNMENT);
  __assume_aligned(FR, ALIGNMENT);
  __assume_aligned(QL, ALIGNMENT);
  __assume_aligned(QR, ALIGNMENT);
#endif

  double* tmp_bnd = ((double *) _mm_malloc(sizeof(double)*392, ALIGNMENT));
  double QavL[8] __attribute__((aligned(ALIGNMENT))) = {0.0};
  double QavR[8] __attribute__((aligned(ALIGNMENT))) = {0.0};
  
  for (int xy = 0; xy < 49; xy++) { //xy = 1 or 2 spatial dim (1 in 2D, 2 in 3D)
    #pragma omp simd aligned(QavL,QavR,QL,QR:ALIGNMENT)
    for (int n = 0; n < 8; n++) {
      QavL[n] += weights2[xy] * QL[xy*8+n];
      QavR[n] += weights2[xy] * QR[xy*8+n];
    }
  }
  
  double lambdaL[8] __attribute__((aligned(ALIGNMENT)));
#ifdef USE_IPO
  #pragma forceinline recursive
#endif
  solver.Euler::EulerSolver_ADERDG::eigenvalues(&QavL[0], direction, &lambdaL[0]);
  double lambdaR[8] __attribute__((aligned(ALIGNMENT)));
#ifdef USE_IPO
  #pragma forceinline recursive
#endif
  solver.Euler::EulerSolver_ADERDG::eigenvalues(&QavR[0], direction, &lambdaR[0]);
  
  double smax = 0.;
  for (int ivar = 0; ivar < 5; ivar++) {
    smax = std::max(smax, std::max(fabs(lambdaL[ivar]), fabs(lambdaR[ivar])));
  }
  
  for (int xy = 0; xy < 49; xy++){
    #pragma omp simd aligned(tmp_bnd,QL,QR:ALIGNMENT)
    for (int n = 0; n < 5; n++) { //skip parameters
      tmp_bnd[xy*8+n] = smax * (QL[xy*8+n]-QR[xy*8+n]);
    }
  }
  
  #pragma omp simd aligned(FL,FR,tmp_bnd:ALIGNMENT)
  for (int xyn = 0; xyn < 392; xyn++) {
    FL[xyn] = 0.5 * (FL[xyn] + FR[xyn] + tmp_bnd[xyn]);
  }
  std::copy_n(FL, 392, FR);
  _mm_free(tmp_bnd);
  

}