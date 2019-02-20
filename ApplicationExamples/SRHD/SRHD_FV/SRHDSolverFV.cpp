#include "SRHDSolverFV.h"

#include <memory>

using std::endl;
using std::cout;

extern "C" {
void hastoadjustsolution_(double* const time, bool* refine);
void adjustedsolutionvalues_(const double* const x,const double* const w,const double* const t,const double* const dt,double* const Q);
void pdeflux_(double* const F, const double* const Q);
void pdeeigenvalues_(double* const lambda, const double* const Q, const double* const nv);
}


void SRHD::SRHDSolverFV::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* const Q) {
  if ( tarch::la::equals(t,0.0) ) {
  // @todo Please implement and set initial conditions
    adjustedsolutionvalues_(x, &w, &t, &dt, Q);
  }
  // @todo Feel free to add further conditions
}


bool SRHD::SRHDSolverFV::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, double dt) {
  // @todo Please implement
  if ( tarch::la::equals(t,0.0) ) {
    // Tell kernel that you want to set initial conditions 
    return true;
  }
  else {
    // @todo Please implement
    return false; 
  }
}


exahype::solvers::Solver::RefinementControl SRHD::SRHDSolverFV::refinementCriterion(const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void SRHD::SRHDSolverFV::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* const lambda) {
  // @todo Please implement
  double nv[3] = {0.};
  nv[normalNonZeroIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}


void SRHD::SRHDSolverFV::flux(const double* const Q, double** const F) {
  // @todo Please implement
  pdeflux_(F[0], Q);
}


void SRHD::SRHDSolverFV::algebraicSource(const double* const Q, double* const S) {
  S[0] = 0.;
  S[1] = 0.;
  S[2] = 0.;
  S[3] = 0.;
  S[4] = 0.;
}

void SRHD::SRHDSolverFV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int normalNonZero,
    const double* const stateIn,
    double* const stateOut) {

  // The problem with these definitions is that in a simulation
  // with a global nonzero velocity (as in MovingGauss2D), errnous
  // values move into the simulation domain very quickly. So these
  // boundary conditions are not good at all. Instead, we should
  // have per default "vacuum" boundary conditions, so that vacuum
  // values enter the grid as soon as matter moves away.

  //  // stateOut
  //  // @todo Please implement
  stateOut[0] = stateIn[0];
  stateOut[1] = stateIn[1];
  stateOut[2] = stateIn[2];
  stateOut[3] = stateIn[3];
  stateOut[4] = stateIn[4];
}
