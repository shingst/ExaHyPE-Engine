#include "DIMSolver.h"

#include "DIMSolver_Variables.h"

#include "PDE.h"

#include "InitialData.h"

#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreBasis.h"

tarch::logging::Log DIM::DIMSolver::_log( "DIM::DIMSolver" );


void DIM::DIMSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

exahype::solvers::ADERDGSolver::AdjustSolutionValue DIM::DIMSolver::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
  // @todo Please implement/augment if required
  return tarch::la::equals(t,0.0) ? exahype::solvers::ADERDGSolver::AdjustSolutionValue::PointWisely : exahype::solvers::ADERDGSolver::AdjustSolutionValue::No;
}

void  __attribute__((optimize("O0"))) DIM::DIMSolver::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
  // Dimensions             = 2
  // Number of variables    = 14 + #parameters
  
  // @todo Please implement/augment if required
  // State variables:
 /*
  Q[ 0] = 0.0;
  Q[ 1] = 0.0;
  Q[ 2] = 0.0;
  Q[ 3] = 0.0;
  Q[ 4] = 0.0;
  Q[ 5] = 0.0;
  Q[ 6] = 0.0;
  Q[ 7] = 0.0;
  Q[ 8] = 0.0;
  Q[ 9] = 0.0;
  Q[10] = 0.0;
  Q[11] = 0.0;
  Q[12] = 0.0;
  Q[13] = 0.0;
*/
  // Fortran
  if (tarch::la::equals(t,0.0)) {
    initialdata_(x, &t, Q);
  }
}

void  __attribute__((optimize("O0"))) DIM::DIMSolver::eigenvalues(const double* const Q,const int d,double* const lambda) {
  // Dimensions             = 2
  // Number of variables    = 14 + #parameters
  
  // @todo Please implement/augment if required
/*
  lambda[ 0] = 1.0;
  lambda[ 1] = 1.0;
  lambda[ 2] = 1.0;
  lambda[ 3] = 1.0;
  lambda[ 4] = 1.0;
  lambda[ 5] = 1.0;
  lambda[ 6] = 1.0;
  lambda[ 7] = 1.0;
  lambda[ 8] = 1.0;
  lambda[ 9] = 1.0;
  lambda[10] = 1.0;
  lambda[11] = 1.0;
  lambda[12] = 1.0;
  lambda[13] = 1.0;
*/
  double nv[3] = {0.};
  nv[d] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}


void  __attribute__((optimize("O0"))) DIM::DIMSolver::flux(const double* const Q,double** const F) {
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


void  __attribute__((optimize("O0"))) DIM::DIMSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  // Dimensions             = 2
  // Number of variables    = 14 + #parameters

  const int nVar = DIM::AbstractDIMSolver::NumberOfVariables;
  const int order = DIM::AbstractDIMSolver::Order;
  const int basisSize = order + 1;
  const int nDim = DIMENSIONS;

  double Qgp[nVar], F[nDim][nVar];

  std::memset(stateOut, 0, nVar * sizeof(double));
  std::memset(fluxOut, 0, nVar * sizeof(double));
  
  for(int i=0; i < basisSize; i++)  { // i == time
     const double weight = kernels::legendre::weights[order][i];
     const double xi = kernels::legendre::nodes[order][i];
     double ti = t + xi * dt;

     initialdata_(x, &ti, Qgp);
     pdeflux_(F[0], F[1], (nDim==3)?F[2]:nullptr, Qgp);
     for(int m=0; m < nVar; m++) {
        stateOut[m] += weight * Qgp[m];
        fluxOut[m] += weight * F[d][m];
     }
  }
//initialdata_(x, &t, stateOut);
  // @todo Please implement/augment if required

/*
  stateOut[ 0] = 0.0;
  stateOut[ 1] = 0.0;
  stateOut[ 2] = 0.0;
  stateOut[ 3] = 0.0;
  stateOut[ 4] = 0.0;
  stateOut[ 5] = 0.0;
  stateOut[ 6] = 0.0;
  stateOut[ 7] = 0.0;
  stateOut[ 8] = 0.0;
  stateOut[ 9] = 0.0;
  stateOut[10] = 0.0;
  stateOut[11] = 0.0;
  stateOut[12] = 0.0;
  stateOut[13] = 0.0;
*/


// initialdata_(x, &t, stateOut);
/*
  fluxOut[ 0] = 0.0;
  fluxOut[ 1] = 0.0;
  fluxOut[ 2] = 0.0;
  fluxOut[ 3] = 0.0;
  fluxOut[ 4] = 0.0;
  fluxOut[ 5] = 0.0;
  fluxOut[ 6] = 0.0;
  fluxOut[ 7] = 0.0;
  fluxOut[ 8] = 0.0;
  fluxOut[ 9] = 0.0;
  fluxOut[10] = 0.0;
  fluxOut[11] = 0.0;
  fluxOut[12] = 0.0;
  fluxOut[13] = 0.0;
*/
// initialdata_(x, &t, fluxOut);
}


exahype::solvers::Solver::RefinementControl DIM::DIMSolver::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  if( (center[0]*center[0]+center[1]*center[1]) < 0.25*0.25 )
  {
    if(dx[0]>getMaximumMeshSize()/9)
    {
      return exahype::solvers::Solver::RefinementControl::Refine;
    }
  }
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void  __attribute__((optimize("O0"))) DIM::DIMSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {

  pdencp_(BgradQ, Q, gradQ);


/*
  BgradQ[ 0] = 0.0;
  BgradQ[ 1] = 0.0;
  BgradQ[ 2] = 0.0;
  BgradQ[ 3] = 0.0;
  BgradQ[ 4] = 0.0;
  BgradQ[ 5] = 0.0;
  BgradQ[ 6] = 0.0;
  BgradQ[ 7] = 0.0;
  BgradQ[ 8] = 0.0;
  BgradQ[ 9] = 0.0;
  BgradQ[10] = 0.0;
  BgradQ[11] = 0.0;
  BgradQ[12] = 0.0;
  BgradQ[13] = 0.0;
*/  
}


void DIM::DIMSolver::coefficientMatrix(const double* const Q,const int d,double* const Bn) {

  double nv[3] = {0.};
  nv[d] = 1;
  pdematrixb_(Bn, Q, nv);

/*
  Bn[ 0] = 0.0;
  Bn[ 1] = 0.0;
  Bn[ 2] = 0.0;
  Bn[ 3] = 0.0;
  Bn[ 4] = 0.0;
  Bn[ 5] = 0.0;
  Bn[ 6] = 0.0;
  Bn[ 7] = 0.0;
  Bn[ 8] = 0.0;
  Bn[ 9] = 0.0;
  Bn[10] = 0.0;
  Bn[11] = 0.0;
  Bn[12] = 0.0;
  Bn[13] = 0.0;
*/
}

