
#include <limits>

#include "kernels/SWE_MySWESolver_ADERDG/Kernels.h"

#include "MySWESolver_ADERDG.h"

double SWE::MySWESolver_ADERDG_kernels::aderdg::stableTimeStepSize(
    SWE::MySWESolver_ADERDG& solver,
    const double* restrict const luh,
    const double inverseDx //Assume dx[0] == dx[1] == dx[2]
) {
  constexpr double cflFactor       = SWE::MySWESolver_ADERDG::CFL;
  constexpr double PNPM[10]        = {1.0,   0.33,  0.17, 0.1,  0.069,
                                      0.045, 0.038, 0.03, 0.02, 0.015};
  
  double lambda[4] __attribute__((aligned(ALIGNMENT))) = {0.0};
  double dt = std::numeric_limits<double>::max();
  
  for (int xyz = 0; xyz < 4; xyz++) { //xyz
    double denominator = 0.0;
    for (int d = 0; d < 2; d++) {
#ifdef USE_IPO
      #pragma forceinline recursive
#endif
      solver.SWE::MySWESolver_ADERDG::eigenvalues(luh+(4*xyz), d, lambda);
      
      double maxEigenvalue = 0.0;
      for (int n = 0; n < 4; n++) {
        maxEigenvalue = std::max(fabs(lambda[n]), maxEigenvalue);
      }
      denominator += maxEigenvalue * inverseDx; //Assume dx[0] == dx[1] == dx[2]
    }

    dt = std::min(dt, cflFactor * PNPM[1] / denominator);  // order = nDof-1

  }

  return dt;
}