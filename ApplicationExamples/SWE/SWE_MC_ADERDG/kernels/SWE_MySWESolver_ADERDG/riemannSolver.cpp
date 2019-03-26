
#include <algorithm>
#include <cstring>
#include <cmath>

#include "kernels/SWE_MySWESolver_ADERDG/Kernels.h"
#include "kernels/SWE_MySWESolver_ADERDG/DGMatrices.h"
#include "kernels/SWE_MySWESolver_ADERDG/Quadrature.h"

#include "MySWESolver_ADERDG.h"

void SWE::MySWESolver_ADERDG_kernels::aderdg::riemannSolver( 
  SWE::MySWESolver_ADERDG& solver,
  double* restrict FL, double* restrict FR,
  const double* restrict const QL, const double* restrict const QR,
  const double t,
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

  double tmp_bnd[8] __attribute__((aligned(ALIGNMENT)));
  double QavL[4] __attribute__((aligned(ALIGNMENT))) = {0.0};
  double QavR[4] __attribute__((aligned(ALIGNMENT))) = {0.0};
  
  for (int xy = 0; xy < 2; xy++) { //xy = 1 or 2 spatial dim (1 in 2D, 2 in 3D)
    #pragma omp simd aligned(QavL,QavR,QL,QR:ALIGNMENT)
    for (int n = 0; n < 4; n++) {
      QavL[n] += weights2[xy] * QL[xy*4+n];
      QavR[n] += weights2[xy] * QR[xy*4+n];
    }
  }
  
  double lambdaL[4] __attribute__((aligned(ALIGNMENT)));
#ifdef USE_IPO
  #pragma forceinline recursive
#endif
  solver.SWE::MySWESolver_ADERDG::eigenvalues(&QavL[0], direction, &lambdaL[0]);
  double lambdaR[4] __attribute__((aligned(ALIGNMENT)));
#ifdef USE_IPO
  #pragma forceinline recursive
#endif
  solver.SWE::MySWESolver_ADERDG::eigenvalues(&QavR[0], direction, &lambdaR[0]);
  
  double smax = 0.;
  for (int ivar = 0; ivar < 4; ivar++) {
    smax = std::max(smax, std::max(fabs(lambdaL[ivar]), fabs(lambdaR[ivar])));
  }
  
  for (int xy = 0; xy < 2; xy++){
    #pragma omp simd aligned(tmp_bnd,QL,QR:ALIGNMENT)
    for (int n = 0; n < 4; n++) { //skip parameters
      tmp_bnd[xy*4+n] = smax * (QL[xy*4+n]-QR[xy*4+n]);
    }
  }
  
  #pragma omp simd aligned(FL,FR,tmp_bnd:ALIGNMENT)
  for (int xyn = 0; xyn < 8; xyn++) {
    FL[xyn] = 0.5 * (FL[xyn] + FR[xyn] + tmp_bnd[xyn]);
  }
  std::copy_n(FL, 8, FR);
    
  //add non-conservative product part

  double Qavg[4] __attribute__((aligned(ALIGNMENT))) = {0.0};
  double gradQ[8] __attribute__((aligned(ALIGNMENT))) = {0.0};
  double ncp[4] __attribute__((aligned(ALIGNMENT))) = {0.0};
  for (int xy = 0; xy < 2; xy++) {
    #pragma omp simd aligned(Qavg,QL,QR,gradQ:ALIGNMENT)
    for (int n = 0; n < 4; n++) {
       Qavg[n] = 0.5 * (QL[xy*4+n] + QR[xy*4+n]);
       gradQ[direction*4+n] = QR[xy*4+n] - QL[xy*4+n];
    }
#ifdef USE_IPO
    #pragma forceinline recursive
#endif
    solver.SWE::MySWESolver_ADERDG::nonConservativeProduct(Qavg, gradQ, ncp);
    #pragma omp simd aligned(FL,FR,ncp:ALIGNMENT)
    for (int n = 0; n < 4; n++) {
      FR[xy*4+n] -= 0.5*ncp[n];
      FL[xy*4+n] += 0.5*ncp[n];
    }
  }

}