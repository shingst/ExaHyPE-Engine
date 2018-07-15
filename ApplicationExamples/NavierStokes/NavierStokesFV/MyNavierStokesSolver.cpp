#include "MyNavierStokesSolver.h"

#include "MyNavierStokesSolver_Variables.h"


tarch::logging::Log NavierStokesFV::MyNavierStokesSolver::_log( "NavierStokesFV::MyNavierStokesSolver" );

void NavierStokesFV::MyNavierStokesSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

void NavierStokesFV::MyNavierStokesSolver::adjustSolution(const double* const x,const double t,const double dt, double* Q) {
  // Dimensions             = 2
  // Number of variables    = 5 + #parameters
  
  // @todo Please implement/augment if required
  Q[0] = 0.0;
  Q[1] = 0.0;
  Q[2] = 0.0;
  Q[3] = 0.0;
  Q[4] = 0.0;
}

void NavierStokesFV::MyNavierStokesSolver::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
  // Dimensions             = 2
  // Number of variables    = 5 + #parameters
  
  // @todo Please implement/augment if required
  lambda[0] = 1.0;
  lambda[1] = 1.0;
  lambda[2] = 1.0;
  lambda[3] = 1.0;
  lambda[4] = 1.0;
}

void NavierStokesFV::MyNavierStokesSolver::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* stateOutside) {
  // Dimensions             = 2
  // Number of variables    = 5 + #parameters

  // @todo Please implement/augment if required
  stateOutside[0] = stateInside[0];
  stateOutside[1] = stateInside[1];
  stateOutside[2] = stateInside[2];
  stateOutside[3] = stateInside[3];
  stateOutside[4] = stateInside[4];
}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit




void NavierStokesFV::MyNavierStokesSolver::parabolicFlux(const double* const Q,const double* const gradQ, double** F) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 5 + 0
  
  // @todo Please implement/augment if required
  F[0][0] = 0.0;
  F[0][1] = 0.0;
  F[0][2] = 0.0;
  F[0][3] = 0.0;
  F[0][4] = 0.0;
  
  F[1][0] = 0.0;
  F[1][1] = 0.0;
  F[1][2] = 0.0;
  F[1][3] = 0.0;
  F[1][4] = 0.0;
  
}




