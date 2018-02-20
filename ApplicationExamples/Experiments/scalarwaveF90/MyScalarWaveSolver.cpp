#include "MyScalarWaveSolver.h"

void ScalarWave::MyScalarWaveSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

bool ScalarWave::MyScalarWaveSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) {
  // @todo Please implement/augment if required
   // @todo Please implement/augment if required
   if (tarch::la::equals(t, 0.0, 1e-15)) {  // @todo precision
    return true;
  }
 
}

void ScalarWave::MyScalarWaveSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Dimensions             = 3
  // Number of variables    = 1 (#unknowns + #parameters)
  
  // @todo Please implement/augment if required
  // State variables:
  Q[0] =  1.0 +
        std::exp(-sqrt((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5) +
                   (x[2] - 0.5) * (x[2] - 0.5)) /((0.25)*(0.25)));

  //0.5*x[0] + 1.0*x[1] + 2.0*x[2];
}

exahype::solvers::Solver::RefinementControl ScalarWave::MyScalarWaveSolver::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


// void ScalarWave::MyScalarWaveSolver::eigenvalues(const double* const Q,const int normalNonZeroIndex,double* lambda) {
//   // Dimensions             = 3
//   // Number of variables    = 1 (#unknowns + #parameters)
  
//   // @todo Please implement/augment if required
//   lambda[0] = 0.0;
// }

// void ScalarWave::MyScalarWaveSolver::flux(const double* const Q,double** F) {
//   // Dimensions             = 3
//   // Number of variables    = 1 (#unknowns + #parameters)
  
//   // @todo Please implement/augment if required
//   F[0][0] = 0.0;

//   F[1][0] = 0.0;

//   F[2][0] = 0.0;
// }


// void ScalarWave::MyScalarWaveSolver::algebraicSource(const double* const Q,double* S) {
//   // Dimensions             = 3
//   // Number of variables    = 1 (#unknowns + #parameters)
  
//   // @todo Please implement/augment if required
//   S[0] = 0.0;
// }


void ScalarWave::MyScalarWaveSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
  // Dimensions             = 3
  // Number of variables    = 1 (#unknowns + #parameters)

  // @todo Please implement/augment if required
  stateOut[0] = stateIn[0];

  fluxOut[0] = fluxIn[0];
}


// void ScalarWave::MyScalarWaveSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
//   // Dimensions             = 3
//   // Number of variables    = 1 (#unknowns + #parameters)

//   // @todo Please implement/augment if required
//   BgradQ[0] = 0.0;
// }

    
// void ScalarWave::MyScalarWaveSolver::coefficientMatrix(const double* const Q,const int normalNonZero,double* Bn) {
//   // Dimensions             = 3
//   // Number of variables    = 1 (#unknowns + #parameters)
  
//   // @todo Please implement/augment if required
//   Bn[0] = 0.0;
// }
