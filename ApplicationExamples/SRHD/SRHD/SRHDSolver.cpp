#include "SRHDSolver.h"
#include "GeneratedConstants.h"

#include <memory>

using std::endl;
using std::cout;

extern "C" {
void hastoadjustsolution_(double* const time, bool* refine);
void adjustedsolutionvalues_(const double* const x,const double* const w,const double* const t,const double* const dt,double* const Q);
void pdeflux_(double* const F, const double* const Q);
void pdeeigenvalues_(double* const lambda, const double* const Q, const double* const nv);
}

void SRHD::SRHDSolver::init(){
  // implement if wanted
}

void SRHD::SRHDSolver::flux(const double* const Q, double** const F) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  pdeflux_(F[0], Q);
}



void SRHD::SRHDSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* const lambda) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  // normal vector: Allocate for 3 dimensions for convenience
  double nv[3] = {0.};
  nv[normalNonZeroIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}



bool SRHD::SRHDSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t, double dt) {
  return (t < 1e-10);

  // This would be the alternative invocation via Fortran.
  // However this crashes for some reason and is also unneccessary.
  bool refine;
  hastoadjustsolution_(&t, &refine);
  return refine;
}



void SRHD::SRHDSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* const Q) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  adjustedsolutionvalues_(x, &w, &t, &dt, Q);
}

void SRHD::SRHDSolver::algebraicSource(const double* const Q, double* const S){
  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;
}


exahype::solvers::Solver::RefinementControl SRHD::SRHDSolver::refinementCriterion(const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  return exahype::solvers::Solver::RefinementControl::Keep;
}

void SRHD::SRHDSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)


  // @todo Please implement
  // fluxOut
  fluxOut[0] = fluxIn[0];
  fluxOut[1] = fluxIn[1];
  fluxOut[2] = fluxIn[2];
  fluxOut[3] = fluxIn[3];
  fluxOut[4] = fluxIn[4];
  // stateOut
  // @todo Please implement
  stateOut[0] = stateIn[0];
  stateOut[1] = stateIn[1];
  stateOut[2] = stateIn[2];
  stateOut[3] = stateIn[3];
  stateOut[4] = stateIn[4];
}


