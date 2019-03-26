
#include <algorithm>

#include "kernels/SWE_MySWESolver_ADERDG/Kernels.h"
#include "kernels/SWE_MySWESolver_ADERDG/Quadrature.h"
#include "kernels/SWE_MySWESolver_ADERDG/DGMatrices.h"

#include "MySWESolver_ADERDG.h"

void SWE::MySWESolver_ADERDG_kernels::aderdg::boundaryConditions(
                        SWE::MySWESolver_ADERDG& solver,
                        double* fluxOut,
                        double* stateOut,
                        const double* const fluxIn,
                        const double* const stateIn,
                        const double* const cellCentre,
                        const double* const cellSize,
                        const double t,const double dt,
                        const int faceIndex,
                        const int normalNonZero
) {                         
  // Compute if face is "left" (0,2,4) or "right" face (1,2,3).
  const int f = faceIndex-2*normalNonZero;

  // Determine the free directions from the non-zero normal entry.
  const int d1 = (3 ^ normalNonZero) / 3; //0->1, 1->0, 2->0 , ^ is bitwise XOR
  
  double x[2];
  x[normalNonZero] = cellCentre[normalNonZero] + (-0.5 + f)*cellSize[normalNonZero];
  
  for (int jj = 0; jj < 1; jj++) {  // loop over dof, loop removed by compiler if 2D
    for (int ii = 0; ii < 2; ii++) {  // loop over dof
      x[d1] = cellCentre[d1] + cellSize[d1] * (nodes[ii] - 0.5); 
#ifdef USE_IPO
    #pragma forceinline recursive
#endif
    // TODO(JMG): Pass gradient here in case of viscous flux
    solver.SWE::MySWESolver_ADERDG::boundaryValues(x,t,dt,faceIndex,normalNonZero,
                          &fluxIn[(jj*1+ii)*4], &stateIn[(jj*1+ii)*4],
					  nullptr,
                          &fluxOut[(jj*1+ii)*4],&stateOut[(jj*1+ii)*4]);   
    }
  }
  
}