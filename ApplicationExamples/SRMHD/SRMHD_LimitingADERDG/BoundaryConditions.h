#ifndef __BOUNDARY_CONDITIONS_SRMHD__
#define __BOUNDARY_CONDITIONS_SRMHD__

#include <algorithm> // string lower
#include "AbstractMHDSolver_ADERDG.h" // Todo: Include a make file variable
// like -DSRMHD_LimitingADERDG so we now about the application name and can
// include "AbstractMHDSolver.h" when this is inside the SRMHD (no limiting)
// application

#include <memory>
#include <cstring>
#include "kernels/GaussLegendreBasis.h"
#include "kernels/KernelUtils.h" // matrix indexing

#include "InitialData.h"
#include "PDE.h"


/**
 * Note that the function signatures here differ from the signature of
 * 
 * void MHD::MHDSolver****::boundaryValues(
 *     const double* const x,const double t,
 *     const double dt, const int faceIndex, const int normalNonZero,
 *     const double * const fluxIn, const double* const stateIn,
 *     double *fluxOut, double* stateOut)
 * 
 * This is because the signatures here are for the FORTRAN functions which
 * pass everything by value and also follow the normalVector DOUBLE(3)
 * convention instead of a int normalNonZero.
 * 
 **/

extern "C" { // Fortran functions:

void boundaryoutflow_(
	const double* const x, const double* const t, const double* const dt,
	const int* const faceIndex, double* nv,
	const double* const fluxIn, const double* const stateIn,
	double* fluxOut, double* stateOut);

void boundaryjet_(
	const double* const x, const double* const t, const double* const dt,
	const int* const faceIndex, double* nv,
	const double* const fluxIn, const double* const stateIn,
	double* fluxOut, double* stateOut);

} /* extern "C" */

typedef void (*BoundaryConditionHandler)(
	const double* const x, const double* const t, const double* const dt,
	const int* const faceIndex, double* nv,
	const double* const fluxIn, const double* const stateIn,
	double* fluxOut, double* stateOut);
extern BoundaryConditionHandler bcfunc;

// The following C++ functions are inlined for convenience.

// Exact Boundary conditions in a generic way, in C:
// actually not used.
inline void BoundaryExact(
	const double* const x, const double* const t, const double* const dt,
	const int* const faceIndex, double* nv,
	const double* const fluxIn, const double* const stateIn,
	double* fluxOut, double* stateOut) {
  // Time integrate the idfunc. This needs idfunct to get the time passed...
  // determine the fluxes for each.
}

// instead, for the time being, we do it only for AlfenWave:

// We could easily do this in Fortran, but this would require that
// kernels::legendre::weights and kernels::legendre::nodes
// can be accessed from Fortran.
inline void BoundaryAlfenWave(
	const double* const x, const double* const t, const double* const dt,
	const int* const faceIndex, double* nv,
	const double* const fluxIn, const double* const stateIn,
	double* fluxOut, double* stateOut) {

  const int nVar = MHD::AbstractMHDSolver_ADERDG::NumberOfVariables;
  const int order = MHD::AbstractMHDSolver_ADERDG::Order;
  const int basisSize = order + 1;
  const int nDim = DIMENSIONS;

  double Qgp[nVar];
  std::memset(stateOut, 0, nVar * sizeof(double));
  std::memset(fluxOut, 0, nVar * sizeof(double));

  double F[3 * nVar]; // Fortran needs continous storage!
                      // Use always 3 dimensions here since the kernels works with those internally; see nDim in PDE.f90;

  kernels::idx2 F_idx(nDim, nVar);

  // Outflow BC conditions.
  // TOOD(Sven): Why don't we use reflection BC?
  if (*faceIndex==2 || *faceIndex==3) {
     for(int m=0; m < nVar; m++) {
        stateOut[m] = stateIn[m];
        fluxOut[m]  = fluxIn[m];
     }
     return;
  }
  
  // the lazy way, determine normalNonZero again from nv.
  int normalNonZero = -1;
  for(int i=0; i<3; i++) if(nv[i]) { normalNonZero=i; break; }
  if(normalNonZero<0) exit(-123); // failure.

  // Integrate solution in gauss points (Qgp) in time
  for(int i=0; i < basisSize; i++)  { // i == time
     const double weight = kernels::legendre::weights[order][i];
     const double xi = kernels::legendre::nodes[order][i];
     double ti = t[0] + xi * dt[0];

     alfenwave_(x, Qgp, &ti);
     pdeflux_(F, Qgp);
     for(int m=0; m < nVar; m++) {
  //if(m==checkm) printf("fluxOut[%d] += %.20e\n", m, weight * F[normalNonZero][m]);
        stateOut[m] += weight * Qgp[m];
        fluxOut[m] += weight * F[F_idx(normalNonZero, m)];
     }
  }
}

#endif /* __BOUNDARY_CONDITIONS_SRMHD__ */
