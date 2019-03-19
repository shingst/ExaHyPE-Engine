#include "MyAcousticWaveSolver.h"

void AcousticWave::MyAcousticWaveSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

bool AcousticWave::MyAcousticWaveSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) {
  // @todo Please implement/augment if required
   if (tarch::la::equals(t, 0.0, 1e-15)) {  // @todo precision
    return true;
  }
  return false;
 
}

void AcousticWave::MyAcousticWaveSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* const Q) {
  // Dimensions             = 3
  // Number of variables    = 4 (#unknowns + #parameters)
  
  // @todo Please implement/augment if required
  // State variables:
   Q[0] =  1.0 +
        std::exp(-((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5) +
		   (x[2] - 0.5) * (x[2] - 0.5)) /((0.25)*(0.25)));
  Q[1] = 0.0;
  Q[2] = 0.0;
  Q[3] = 0.0;
}

exahype::solvers::Solver::RefinementControl AcousticWave::MyAcousticWaveSolver::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


// void ScalarWave::MyScalarWaveSolver::eigenvalues(const double* const Q,const int normalNonZeroIndex,double* const lambda) {
//   // Dimensions             = 3
//   // Number of variables    = 4 (#unknowns + #parameters)
  
//   // @todo Please implement/augment if required
//   lambda[0] = 0.0;
//   lambda[1] = 0.0;
//   lambda[2] = 0.0;
//   lambda[3] = 0.0;
// }

// void ScalarWave::MyScalarWaveSolver::flux(const double* const Q,double** const F) {
//   // Dimensions             = 3
//   // Number of variables    = 4 (#unknowns + #parameters)
  
//   // @todo Please implement/augment if required
//   F[0][0] = 0.0;
//   F[0][1] = 0.0;
//   F[0][2] = 0.0;
//   F[0][3] = 0.0;

//   F[1][0] = 0.0;
//   F[1][1] = 0.0;
//   F[1][2] = 0.0;
//   F[1][3] = 0.0;

//   F[2][0] = 0.0;
//   F[2][1] = 0.0;
//   F[2][2] = 0.0;
//   F[2][3] = 0.0;
// }


// void ScalarWave::MyScalarWaveSolver::algebraicSource(const double* const Q,double* const S) {
//   // Dimensions             = 3
//   // Number of variables    = 4 (#unknowns + #parameters)
  
//   // @todo Please implement/augment if required
//   S[0] = 0.0;
//   S[1] = 0.0;
//   S[2] = 0.0;
//   S[3] = 0.0;
// }


void AcousticWave::MyAcousticWaveSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  // Dimensions             = 3
  // Number of variables    = 4 (#unknowns + #parameters)

  // @todo Please implement/augment if required
  stateOut[0] = stateIn[0];
  stateOut[1] = stateIn[1];
  stateOut[2] = stateIn[2];
  stateOut[3] = stateIn[3];

  fluxOut[0] = fluxIn[0];
  fluxOut[1] = fluxIn[1];
  fluxOut[2] = fluxIn[2];
  fluxOut[3] = fluxIn[3];
}


// void ScalarWave::MyScalarWaveSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
//   // Dimensions             = 3
//   // Number of variables    = 4 (#unknowns + #parameters)

//   // @todo Please implement/augment if required
//   BgradQ[0] = 0.0;
//   BgradQ[1] = 0.0;
//   BgradQ[2] = 0.0;
//   BgradQ[3] = 0.0;
// }

    
// void ScalarWave::MyScalarWaveSolver::coefficientMatrix(const double* const Q,const int normalNonZero,double* const Bn) {
//   // Dimensions             = 3
//   // Number of variables    = 4 (#unknowns + #parameters)
  
//   // @todo Please implement/augment if required
//   Bn[0] = 0.0;
//   Bn[1] = 0.0;
//   Bn[2] = 0.0;
//   Bn[3] = 0.0;
//   Bn[4] = 0.0;
//   Bn[5] = 0.0;
//   Bn[6] = 0.0;
//   Bn[7] = 0.0;
//   Bn[8] = 0.0;
//   Bn[9] = 0.0;
//   Bn[10] = 0.0;
//   Bn[11] = 0.0;
//   Bn[12] = 0.0;
//   Bn[13] = 0.0;
//   Bn[14] = 0.0;
//   Bn[15] = 0.0;
// }
