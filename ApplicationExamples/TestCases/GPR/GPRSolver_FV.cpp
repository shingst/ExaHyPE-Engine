#include "GPRSolver_FV.h"

#include "GPRSolver_FV_Variables.h"
// User defined calls
#include "PDE.h"
#include "InitialData.h"
#include "C2P-GPR.h"

#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing

tarch::logging::Log GPR::GPRSolver_FV::_log( "GPR::GPRSolver_FV" );


void GPR::GPRSolver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

void GPR::GPRSolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {
  // Dimensions             = 3
  // Number of variables    = 17 + #parameters
  
  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
  initialdata_(x, &t, Q);
  }
}

void GPR::GPRSolver_FV::eigenvalues(const double* const Q, const int dIndex, double* const lambda) {
  // Dimensions             = 3
  // Number of variables    = 17 + #parameters
  
  // @todo Please implement/augment if required
  double nv[3] = {0.};
  nv[dIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}

void GPR::GPRSolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* const stateOutside) {
	const int nVar = GPR::AbstractGPRSolver_FV::NumberOfVariables;	
	double Qgp[nVar];
  // Dimensions             = 3
  // Number of variables    = 17 + #parameters

  // @todo Please implement/augment if required
  /*
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
  */
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


void GPR::GPRSolver_FV::flux(const double* const Q,double** const F) {
	const int nVar = GPR::AbstractGPRSolver_FV::NumberOfVariables;
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
void GPR::GPRSolver_FV::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
  // @todo Please implement/augment if required
  pdesource_(S, Q);
}

void  GPR::GPRSolver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  // @todo Please implement/augment if required
  pdencp_(BgradQ, Q, gradQ);
}

