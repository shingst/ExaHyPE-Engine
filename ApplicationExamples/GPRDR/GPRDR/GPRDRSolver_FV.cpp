#include "GPRDRSolver_FV.h"
#include "kernels/finitevolumes/musclhancock/c/musclhancock.h"
#include "GPRDRSolver_FV_Variables.h"


tarch::logging::Log GPRDR::GPRDRSolver_FV::_log( "GPRDR::GPRDRSolver_FV" );

void GPRDR::GPRDRSolver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // Tip: You find documentation for this method in header file "GPRDR::GPRDRSolver_FV.h".
  
  // @todo Please implement/augment if required
}

void GPRDR::GPRDRSolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {
  // if (tarch::la::equals(t,0.0)) {
  //   int md = exahype::solvers::Solver::getMaximumAdaptiveMeshDepth();
  //   double cms = exahype::solvers::Solver::getCoarsestMeshSize();
  //   const int order = 0;
  //   initialdata_(x, &t, Q,&md,&cms,&order);
  // } 

  Q[0] = 0.0;
  Q[1] = 0.0;
  Q[2] = 0.0;
  Q[3] = 0.0;
  Q[4] = 0.0;
  Q[5] = 0.0;
  Q[6] = 0.0;
  Q[7] = 0.0;
  Q[8] = 0.0;
  Q[9] = 0.0;
  Q[10] = 0.0;
  Q[11] = 0.0;
  Q[12] = 0.0;
  Q[13] = 0.0;
  Q[14] = 0.0;
  Q[15] = 0.0;
  Q[16] = 0.0;
  Q[17] = 0.0;
  Q[18] = 0.0;
  Q[19] = 0.0;
  Q[20] = 0.0;
  Q[21] = 0.0;
  Q[22] = 0.0;
  Q[23] = 0.0;
}

void GPRDR::GPRDRSolver_FV::eigenvalues(const double* const Q, const int dIndex, double* const lambda) {
  // Tip: You find documentation for this method in header file "GPRDR::GPRDRSolver_FV.h".
  // Tip: See header file "GPRDR::AbstractGPRDRSolver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
  
  // @todo Please implement/augment if required
  // double nv[3] = {0.};
  // nv[d] = 1;
  // pdeeigenvalues_(lambda, Q, nv);
  lambda[0] = 1.0;
  lambda[1] = 1.0;
  lambda[2] = 1.0;
  lambda[3] = 1.0;
  lambda[4] = 1.0;
  lambda[5] = 1.0;
  lambda[6] = 1.0;
  lambda[7] = 1.0;
  lambda[8] = 1.0;
  lambda[9] = 1.0;
  lambda[10] = 1.0;
  lambda[11] = 1.0;
  lambda[12] = 1.0;
  lambda[13] = 1.0;
  lambda[14] = 1.0;
  lambda[15] = 1.0;
  lambda[16] = 1.0;
  lambda[17] = 1.0;
  lambda[18] = 1.0;
  lambda[19] = 1.0;
  lambda[20] = 1.0;
  lambda[21] = 1.0;
  lambda[22] = 1.0;
  lambda[23] = 1.0;

}

void GPRDR::GPRDRSolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int direction,
    const double* const stateInside,
    double* const stateOutside) {
  // Tip: You find documentation for this method in header file "GPRDR::GPRDRSolver_FV.h".
  // Tip: See header file "GPRDR::AbstractGPRDRSolver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.

  // @todo Please implement/augment if required
  // const int nVar = GPRDR::GPRDRSolver_FV::NumberOfVariables;	
  // double Qgp[nVar];
  // int md=0;
  // double cms=0;
  // const int order=0;
  
  // double ti = t + 0.5 * dt;
  // // Compute the outer state according to the initial condition
  // initialdata_(x, &ti, Qgp,&md,&cms,&order);
  // // Assign the proper outer state
  // for(int m=0; m < nVar; m++) {
  //   stateOutside[m] = Qgp[m];
  // }

  stateOutside[0] = stateInside[0];
  stateOutside[1] = stateInside[1];
  stateOutside[2] = stateInside[2];
  stateOutside[3] = stateInside[3];
  stateOutside[4] = stateInside[4];
  stateOutside[5] = stateInside[5];
  stateOutside[6] = stateInside[6];
  stateOutside[7] = stateInside[7];
  stateOutside[8] = stateInside[8];
  stateOutside[9] = stateInside[9];
  stateOutside[10] = stateInside[10];
  stateOutside[11] = stateInside[11];
  stateOutside[12] = stateInside[12];
  stateOutside[13] = stateInside[13];
  stateOutside[14] = stateInside[14];
  stateOutside[15] = stateInside[15];
  stateOutside[16] = stateInside[16];
  stateOutside[17] = stateInside[17];
  stateOutside[18] = stateInside[18];
  stateOutside[19] = stateInside[19];
  stateOutside[20] = stateInside[20];
  stateOutside[21] = stateInside[21];
  stateOutside[22] = stateInside[22];
  stateOutside[23] = stateInside[23];
}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit


void GPRDR::GPRDRSolver_FV::flux(const double* const Q,double** const F) {
  // Tip: You find documentation for this method in header file "GPRDR::GPRDRSolver_FV.h".
  // Tip: See header file "GPRDR::AbstractGPRDRSolver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
  
  // @todo Please implement/augment if required
  // const int nVar = GPRDR::GPRDRSolver_FV::NumberOfVariables;
  // if(DIMENSIONS == 2){
  //   double F_3[nVar];
  //   pdeflux_(F[0], F[1],F_3, Q);
  // }else{
  //   pdeflux_(F[0], F[1],F[2], Q);
  // }
  F[0][0] = 0.0;
  F[0][1] = 0.0;
  F[0][2] = 0.0;
  F[0][3] = 0.0;
  F[0][4] = 0.0;
  F[0][5] = 0.0;
  F[0][6] = 0.0;
  F[0][7] = 0.0;
  F[0][8] = 0.0;
  F[0][9] = 0.0;
  F[0][10] = 0.0;
  F[0][11] = 0.0;
  F[0][12] = 0.0;
  F[0][13] = 0.0;
  F[0][14] = 0.0;
  F[0][15] = 0.0;
  F[0][16] = 0.0;
  F[0][17] = 0.0;
  F[0][18] = 0.0;
  F[0][19] = 0.0;
  F[0][20] = 0.0;
  F[0][21] = 0.0;
  F[0][22] = 0.0;
  F[0][23] = 0.0;
  
  F[1][0] = 0.0;
  F[1][1] = 0.0;
  F[1][2] = 0.0;
  F[1][3] = 0.0;
  F[1][4] = 0.0;
  F[1][5] = 0.0;
  F[1][6] = 0.0;
  F[1][7] = 0.0;
  F[1][8] = 0.0;
  F[1][9] = 0.0;
  F[1][10] = 0.0;
  F[1][11] = 0.0;
  F[1][12] = 0.0;
  F[1][13] = 0.0;
  F[1][14] = 0.0;
  F[1][15] = 0.0;
  F[1][16] = 0.0;
  F[1][17] = 0.0;
  F[1][18] = 0.0;
  F[1][19] = 0.0;
  F[1][20] = 0.0;
  F[1][21] = 0.0;
  F[1][22] = 0.0;
  F[1][23] = 0.0;
  
}




//You can either implement this method or modify fusedSource
void GPRDR::GPRDRSolver_FV::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
  // Tip: You find documentation for this method in header file "GPRDR::GPRDRSolver_FV.h".
  // Tip: See header file "GPRDR::AbstractGPRDRSolver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
  // @todo Please implement/augment if required
  //pdesource_(S, Q);
  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;
  S[5] = 0.0;
  S[6] = 0.0;
  S[7] = 0.0;
  S[8] = 0.0;
  S[9] = 0.0;
  S[10] = 0.0;
  S[11] = 0.0;
  S[12] = 0.0;
  S[13] = 0.0;
  S[14] = 0.0;
  S[15] = 0.0;
  S[16] = 0.0;
  S[17] = 0.0;
  S[18] = 0.0;
  S[19] = 0.0;
  S[20] = 0.0;
  S[21] = 0.0;
  S[22] = 0.0;
  S[23] = 0.0;
}

void  GPRDR::GPRDRSolver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  // Tip: You find documentation for this method in header file "GPRDR::GPRDRSolver_FV.h".
  // Tip: See header file "GPRDR::AbstractGPRDRSolver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
  
  // @todo Please implement/augment if required
  //pdencp_(BgradQ, Q, gradQ);
  BgradQ[0] = 0.0;
  BgradQ[1] = 0.0;
  BgradQ[2] = 0.0;
  BgradQ[3] = 0.0;
  BgradQ[4] = 0.0;
  BgradQ[5] = 0.0;
  BgradQ[6] = 0.0;
  BgradQ[7] = 0.0;
  BgradQ[8] = 0.0;
  BgradQ[9] = 0.0;
  BgradQ[10] = 0.0;
  BgradQ[11] = 0.0;
  BgradQ[12] = 0.0;
  BgradQ[13] = 0.0;
  BgradQ[14] = 0.0;
  BgradQ[15] = 0.0;
  BgradQ[16] = 0.0;
  BgradQ[17] = 0.0;
  BgradQ[18] = 0.0;
  BgradQ[19] = 0.0;
  BgradQ[20] = 0.0;
  BgradQ[21] = 0.0;
  BgradQ[22] = 0.0;
  BgradQ[23] = 0.0;
}


void GPRDR::GPRDRSolver_FV::solutionUpdate(double* luh,const tarch::la::Vector<DIMENSIONS,double>& cellCenter,const tarch::la::Vector<DIMENSIONS,double>& cellSize,const double t, const double dt,double& maxAdmissibleDt) {
  maxAdmissibleDt = kernels::finitevolumes::musclhancock::c::solutionUpdate<
    true, true, true, false, false,
    kernels::finitevolumes::commons::c::minmod,
    GPRDRSolver_FV
    >(*static_cast<GPRDRSolver_FV*>(this),luh,cellCenter,cellSize,t,dt);
}

