#include "GPRSolver_FV.h"

#include "GPRSolver_FV_Variables.h"
#include "PDE.h"
#include "InitialData.h"
#include "C2P-GRGPR.h"

#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing

tarch::logging::Log GRGPR::GPRSolver_FV::_log( "GRGPR::GPRSolver_FV" );

void GRGPR::GPRSolver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

void GRGPR::GPRSolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* Q) {
  // Dimensions             = 2
  // Number of variables    = 30 + #parameters
  
  // @todo Please implement/augment if required
  initialdata_(x, &t, Q);
}

void GRGPR::GPRSolver_FV::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
  // Dimensions             = 2
  // Number of variables    = 30 + #parameters
  
  // @todo Please implement/augment if required
  double nv[3] = {0.};
  nv[dIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}

void GRGPR::GPRSolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* stateOutside) {
		
	const int nVar = GRGPR::AbstractGPRSolver_FV::NumberOfVariables;	
	double Qgp[nVar];
  // Dimensions             = 2
  // Number of variables    = 30 + #parameters

  // @todo Please implement/augment if required
  double ti = t + 0.5 * dt;
  initialdata_(x, &ti, Qgp);
  for(int m=0; m < nVar; m++) {
        stateOutside[m] = Qgp[m];
  }
}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit


void GRGPR::GPRSolver_FV::flux(const double* const Q,double** F) {
	const int nVar = GRGPR::AbstractGPRSolver_FV::NumberOfVariables;
  // Dimensions                        = 3
  // Number of variables + parameters  = 17 + 0
  
  // @todo Please implement/augment if required
  // pdeflux_(F[0], F[1], (DIMENSIONS==3)?F[2]:nullptr, Q);
    if(DIMENSIONS == 2){
		double F_3[nVar];
		pdeflux_(F[0], F[1],F_3, Q);
	}else{
		pdeflux_(F[0], F[1],F[2], Q);
	}
}


//You can either implement this method or modify fusedSource
void GRGPR::GPRSolver_FV::algebraicSource(const double* const Q,double* S) {
  // @todo Please implement/augment if required
  pdesource_(S, Q);
}

void  GRGPR::GPRSolver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  // @todo Please implement/augment if required
  pdencp_(BgradQ, Q, gradQ);
}

