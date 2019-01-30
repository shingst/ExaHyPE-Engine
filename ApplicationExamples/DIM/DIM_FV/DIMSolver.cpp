#include "DIMSolver.h"

#include "DIMSolver_Variables.h"

#include "PDE.h"

#include "InitialData.h"


tarch::logging::Log DIM::DIMSolver::_log( "DIM::DIMSolver" );


void DIM::DIMSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

bool DIM::DIMSolver::useAdjustSolution(const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, const double t, const double dt) const {
  // @todo Please implement/augment if required
  return tarch::la::equals(t,0.0);
}

void DIM::DIMSolver::adjustSolution(const double* const x,const double w,const double t,const double dt, double* Q) {
  // Dimensions             = 2
  // Number of variables    = 14 + #parameters
  
  // Fortran
  initialdata_(x, &t, Q);
}

exahype::solvers::Solver::RefinementControl DIM::DIMSolver::refinementCriterion(const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void __attribute__((optimize("O0"))) DIM::DIMSolver::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
  // Dimensions             = 2
  // Number of variables    = 14 + #parameters
  
  double nv[3] = {0.};
  nv[dIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
    
}

/* void DIM::DIMSolver::flux(const double* const Q, double** F) {
  // Dimensions             = 2
  // Number of variables    = 14 + #parameters
  
  // @todo Please implement/augment if required
  F[0][ 0] = 0.0;
  F[0][ 1] = 0.0;
  F[0][ 2] = 0.0;
  F[0][ 3] = 0.0;
  F[0][ 4] = 0.0;
  F[0][ 5] = 0.0;
  F[0][ 6] = 0.0;
  F[0][ 7] = 0.0;
  F[0][ 8] = 0.0;
  F[0][ 9] = 0.0;
  F[0][10] = 0.0;
  F[0][11] = 0.0;
  F[0][12] = 0.0;
  F[0][13] = 0.0;

  F[1][ 0] = 0.0;
  F[1][ 1] = 0.0;
  F[1][ 2] = 0.0;
  F[1][ 3] = 0.0;
  F[1][ 4] = 0.0;
  F[1][ 5] = 0.0;
  F[1][ 6] = 0.0;
  F[1][ 7] = 0.0;
  F[1][ 8] = 0.0;
  F[1][ 9] = 0.0;
  F[1][10] = 0.0;
  F[1][11] = 0.0;
  F[1][12] = 0.0;
  F[1][13] = 0.0;
}

*/ 

void DIM::DIMSolver::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* stateOutside) {
  // Dimensions             = 2
  // Number of variables    = 14 + #parameters

   initialdata_(x, &t, stateOutside);
  
}

void  __attribute__((optimize("O0"))) DIM::DIMSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
    pdencp_(BgradQ, Q, gradQ);
}

void  __attribute__((optimize("O0"))) DIM::DIMSolver::algebraicSource(const double* const Q, double* S) {
  pdesource_(S, Q);

}


void DIM::DIMSolver::coefficientMatrix(const double* const Q,const int d,double* Bn) {
  double nv[3] = {0.};
  nv[d] = 1;
  pdematrixb_(Bn, Q, nv);
}

