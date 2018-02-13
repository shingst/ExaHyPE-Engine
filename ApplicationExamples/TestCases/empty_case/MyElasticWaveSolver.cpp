#include "MyElasticWaveSolver.h"

#include "MyElasticWaveSolver_Variables.h"


tarch::logging::Log Elastic::MyElasticWaveSolver::_log( "Elastic::MyElasticWaveSolver" );


void Elastic::MyElasticWaveSolver::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

void Elastic::MyElasticWaveSolver::adjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 9 + 3
  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
    Q[0] = 0.0;
    Q[1] = 0.0;
    Q[2] = 0.0;
    Q[3] = 0.0;
    Q[4] = 0.0;
    Q[5] = 0.0;
    Q[6] = 0.0;
    Q[7] = 0.0;
    Q[8] = 0.0;
    Q[9] = 1.0; // test parameter
    Q[10] = 0.0;
    Q[11] = 0.0;
  }
}

void Elastic::MyElasticWaveSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 9 + 3

  // @todo Please implement/augment if required
  stateOut[0] = 0.0;
  stateOut[1] = 0.0;
  stateOut[2] = 0.0;
  stateOut[3] = 0.0;
  stateOut[4] = 0.0;
  stateOut[5] = 0.0;
  stateOut[6] = 0.0;
  stateOut[7] = 0.0;
  stateOut[8] = 0.0;
  stateOut[9] = 1.0;

  fluxOut[0] = 0.0;
  fluxOut[1] = 0.0;
  fluxOut[2] = 0.0;
  fluxOut[3] = 0.0;
  fluxOut[4] = 0.0;
  fluxOut[5] = 0.0;
  fluxOut[6] = 0.0;
  fluxOut[7] = 0.0;
  fluxOut[8] = 0.0;
}

exahype::solvers::Solver::RefinementControl Elastic::MyElasticWaveSolver::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void Elastic::MyElasticWaveSolver::eigenvalues(const double* const Q,const int d,double* lambda) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 9 + 3
  
  // @todo Please implement/augment if required
  lambda[0] = 1.0;
  lambda[1] = 1.0;
  lambda[2] = 1.0;
  lambda[3] = 1.0;
  lambda[4] = 1.0;
  lambda[5] = 1.0;
  lambda[6] = 1.0;
  lambda[7] = 1.0;
  lambda[8] = 1.0;
}


void Elastic::MyElasticWaveSolver::flux(const double* const Q,double** F) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 9 + 3
  
  // @todo Please implement/augment if required

  if((Q[9]-1.0)*(Q[9]- 1.0)> 1.0e-5){
    std::cout<< "Test parameter not read correctly" << std::endl;
    std::terminate();
  }
  F[0][0] = 0.0;
  F[0][1] = 0.0;
  F[0][2] = 0.0;
  F[0][3] = 0.0;
  F[0][4] = 0.0;
  F[0][5] = 0.0;
  F[0][6] = 0.0;
  F[0][7] = 0.0;
  F[0][8] = 0.0;
  
  F[1][0] = 0.0;
  F[1][1] = 0.0;
  F[1][2] = 0.0;
  F[1][3] = 0.0;
  F[1][4] = 0.0;
  F[1][5] = 0.0;
  F[1][6] = 0.0;
  F[1][7] = 0.0;
  F[1][8] = 0.0;
  
  F[2][0] = 0.0;
  F[2][1] = 0.0;
  F[2][2] = 0.0;
  F[2][3] = 0.0;
  F[2][4] = 0.0;
  F[2][5] = 0.0;
  F[2][6] = 0.0;
  F[2][7] = 0.0;
  F[2][8] = 0.0;
  
}



void  Elastic::MyElasticWaveSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  // @todo Please implement/augment if required
  BgradQ[0] = 0.0;
  BgradQ[1] = 0.0;
  BgradQ[2] = 0.0;
  BgradQ[3] = 0.0;
  BgradQ[4] = 0.0;
  BgradQ[5] = 0.0;
  BgradQ[6] = 0.0;
  BgradQ[7] = 0.0;
  BgradQ[8] = 0.0;
}


