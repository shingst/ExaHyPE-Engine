#include "MHDSolver.h"

#include <memory>

#include "InitialDataAdapter.h"
#include "PDE.h"
#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/KernelUtils.h" // matrix indexing
#include <cstring> // memset


// Fortran functions:
extern "C" {
void adjustedsolutionvalues_(const double* const x,const double* const w,const double* const t,const double* const dt,double* const Q);
}/* extern "C" */


void MHDSolver::MHDSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // This function is called by the constructor.
  // You can access spec file parameters as well as command line arguments (argv as std::vector).
  // @todo Please implement/augment if required.
}

//************************************************* 
//for FORTRAN kernels the fluxes and eigenvalues 
//have to be implemented in the file ./PDE.f90 and ./typesDef.f90 
//
//You have to create these files yourself
//and follow the sample applications in the wiki
//************************************************* 






void MHDSolver::MHDSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  // Dimensions             = 3
  // Number of variables    = 9 (#unknowns + #parameters)


  const int basisSize = order + 1;

  double Qgp[nVar];
  std::memset(stateOut, 0.0, nVar * sizeof(double));
  std::memset(fluxOut, 0.0, nVar * sizeof(double));

  double F[nDim * nVar]; // Fortran needs continous storage!
  kernels::idx2 F_idx(nDim, nVar);

  // Integrate solution in gauss points (Qgp) in time
  for(int i=0; i < basisSize; i++)  { // i == time
     const double weight = kernels::gaussLegendreWeights[order][i];
     const double xi = kernels::gaussLegendreNodes[order][i];
     double ti = t + xi * dt;

     alfenwave_(x, Qgp, &ti);
     pdeflux_(F, Qgp);
     for(int m=0; m < nVar; m++) {
	//if(m==checkm) printf("fluxOut[%d] += %.20e\n", m, weight * F[normalNonZero][m]);
	stateOut[m] += weight * Qgp[m];
        fluxOut[m] += weight * F[F_idx(normalNonZero, m)];
     }
  }
  const int statem=6; // the interesting component
}






bool MHDSolver::MHDSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t, double dt) {
  // @todo Please implement
  return (t < 1e-10);
}



void MHDSolver::MHDSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* const Q) {
  // Dimensions             = 3
  // Number of variables    = 9 (#unknowns + #parameters)
  // @todo Please implement

  adjustedsolutionvalues_(x, &w, &t, &dt, Q);
}



exahype::solvers::Solver::RefinementControl MHDSolver::MHDSolver::refinementCriterion(const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}









