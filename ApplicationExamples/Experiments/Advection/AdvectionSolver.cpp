#include "AdvectionSolver.h"

#include "AdvectionSolver_Variables.h"
#include <cmath>

tarch::logging::Log Trivial::AdvectionSolver::_log( "Trivial::AdvectionSolver" );

double f(const double* const x, const double t) {
  return 1.0 + 0.4 * std::sin(x[0] - 0.2*t);
}

void Trivial::AdvectionSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

void Trivial::AdvectionSolver::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
  if (tarch::la::equals(t,0.0)) {
    Q[0] = f(x, t);
  }
}

void Trivial::AdvectionSolver::eigenvalues(const double* const Q,const int d,double* const lambda) {
  lambda[0] = 1.0;
}


void Trivial::AdvectionSolver::flux(const double* const Q,double** const F) {
  F[0][0] = 1.0;
  F[1][0] = 0.0;
  F[2][0] = 0.0;
}


void Trivial::AdvectionSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double* const fluxOut,double* const stateOut) {
  stateOut[0] = f(x, t);
  fluxOut[0] = 1.0;
}


exahype::solvers::Solver::RefinementControl Trivial::AdvectionSolver::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  return exahype::solvers::Solver::RefinementControl::Keep;
}
