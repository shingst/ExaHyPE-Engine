
#include "kernels/SWE_MySWESolver_ADERDG/Kernels.h"
#include "kernels/SWE_MySWESolver_ADERDG/Quadrature.h"

#include "MySWESolver_ADERDG.h"

void SWE::MySWESolver_ADERDG_kernels::aderdg::solutionAdjustment(
  SWE::MySWESolver_ADERDG& solver,
  double* luh,
  const double* const center,
  const double dx, //Assume dx[0] == dx[1] == dx[2]
  const double t,
  const double dt
) {
  double x0[2];
  int pos = 0;

    for (int y = 0; y < 2; y++) {
      x0[1] = center[1] + dx * (nodes[y] - 0.5);
      for (int x = 0; x < 2; x++) {
        x0[0] = center[0] + dx * (nodes[x] - 0.5);

        solver.SWE::MySWESolver_ADERDG::adjustPointSolution(x0, t, dt, luh+(pos*4));
        pos++;
      }
    }

}