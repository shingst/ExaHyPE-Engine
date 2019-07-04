#include "GPRDRSolver_FV.h"
#include "kernels/finitevolumes/musclhancock/c/musclhancock.h"
#include "GPRDRSolver_FV_Variables.h"
#include "InitialData.h"
#include "PDE.h"
#include "ODE.h"


tarch::logging::Log GPRDR::GPRDRSolver_FV::_log( "GPRDR::GPRDRSolver_FV" );

void GPRDR::GPRDRSolver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
}

void GPRDR::GPRDRSolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {
    if (tarch::la::equals(t,0.0)) {
    int md = exahype::solvers::Solver::getMaximumAdaptiveMeshDepth();
    double cms = exahype::solvers::Solver::getCoarsestMeshSize();
    const int order = GPRDR::GPRDRSolver_FV::PatchSize;
    std::fill_n(Q,24,0.0);
    
    //    initialdata_(x, &ti, Qgp,&md,&cms,&order);
    double x_3[3];
    x_3[2]=0;
    std::copy_n(&x[0],DIMENSIONS,&x_3[0]);
    
    initialdata_(x_3, &t, Q);
    for(int i = 0; i< 24 ; i++){
      assert(std::isfinite(Q[i]));
    }
  }
}

void GPRDR::GPRDRSolver_FV::eigenvalues(const double* const Q, const int dIndex, double* const lambda) {
  // Tip: You find documentation for this method in header file "GPRDR::GPRDRSolver_FV.h".
  // Tip: See header file "GPRDR::AbstractGPRDRSolver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
  
  // @todo Please implement/augment if required
  double nv[3] = {0.};
  nv[dIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}

void GPRDR::GPRDRSolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int direction,
    const double* const stateInside,
    double* const stateOutside) {
  const int nVar = GPRDR::GPRDRSolver_FV::NumberOfVariables;
  std::copy_n(stateInside,nVar,stateOutside);
}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit


void GPRDR::GPRDRSolver_FV::flux(const double* const Q,double** const F) {
  const int nVar = GPRDR::GPRDRSolver_FV::NumberOfVariables;
  if(DIMENSIONS == 2){
    double F_3[nVar];
    pdeflux_(F[0], F[1],F_3, Q);
  }else{
    pdeflux_(F[0], F[1],F[2], Q);
  }
}




//You can either implement this method or modify fusedSource
void GPRDR::GPRDRSolver_FV::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
  pdesource_(S, Q);
}

void  GPRDR::GPRDRSolver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  pdencp_(BgradQ, Q, gradQ);
}


void GPRDR::GPRDRSolver_FV::solutionUpdate(double* luh,const tarch::la::Vector<DIMENSIONS,double>& cellCenter,const tarch::la::Vector<DIMENSIONS,double>& cellSize,const double t, const double dt,double& maxAdmissibleDt) {
  GPRDR::AbstractGPRDRSolver_FV::solutionUpdate(luh, cellCenter, cellSize, t, dt, maxAdmissibleDt);
  constexpr int patchSize          = GPRDR::GPRDRSolver_FV::PatchSize;
  constexpr int ghostLayerWidth    = GPRDR::GPRDRSolver_FV::GhostLayerWidth;
  constexpr int patchBegin         = ghostLayerWidth; // patchBegin cell is inside domain
  constexpr int patchEnd           = patchBegin+patchSize; // patchEnd cell is outside domain
#if DIMENSIONS==3
  kernels::idx4 idx(patchEnd,patchEnd,patchEnd,NumberOfVariables);
  for (int k = patchBegin; k < patchEnd; k++) {
#else
  kernels::idx3 idx(patchEnd,patchEnd,NumberOfVariables);
#endif    
    for (int j = patchBegin; j < patchEnd; j++) {
      for (int i = patchBegin; i < patchEnd; i++) {
#if DIMENSIONS==3
	double* luh_cell = luh + idx(k,j,i,0);
#else
	double* luh_cell = luh + idx(j,i,0);
#endif

	double luh_cell_new[NumberOfVariables];
	updatesolutionode_(&luh_cell_new[0],luh_cell,&dt);
	std::copy_n(luh_cell_new,NumberOfVariables,luh_cell);
      }
    }
#if DIMENSIONS==3
  }
#endif  
  
}

