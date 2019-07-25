#include "DIMSolver_FV.h"

#include "DIMSolver_FV_Variables.h"

#include "DIMSolver_FV_Variables.h"
#include "PDE.h"
#include "C2P-DIM.h"
#include "InitialData.h"

#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreBasis.h"

tarch::logging::Log DIM::DIMSolver_FV::_log( "DIM::DIMSolver_FV" );


void DIM::DIMSolver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
    std::cout << " ==================================================================================" << std::endl;
	//readcgfile_(&_domainOffset[0],&_domainSize[0]);
}

void DIM::DIMSolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {
  // Dimensions             = 3
  // Number of variables    = 14 + #parameters
  const int nVar = DIM::AbstractDIMSolver_FV::NumberOfVariables;

  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
  initialdata_(x, &t, Q);
  //Q[nVar]=Q[9];
  //Q[nVar+1]=Q[10];
  //Q[nVar+2]=Q[11];
  }

	Q[13]=1.0;
  //Q[nVar+3]=1.0;
}
/*
void DIM::DIMSolver_FV::algebraicSource(const double* const Q,double* const S) {
	const int nVar = DIM::AbstractDIMSolver_FV::NumberOfVariables;
  // @todo Please implement/augment if required
  for(int m=0; m < nVar; m++) {
	S[m]=0.0;  
  }
}
*/
void DIM::DIMSolver_FV::eigenvalues(const double* const Q, const int dIndex, double* const lambda) {
  // Dimensions             = 3
  // Number of variables    = 14 + #parameters
  
  // @todo Please implement/augment if required
  double nv[3] = {0.};
  nv[dIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}

void DIM::DIMSolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* const stateOutside) {
	const int nVar = DIM::AbstractDIMSolver_FV::NumberOfVariables;	
	const int nPar = DIM::AbstractDIMSolver_FV::NumberOfParameters;
	double Qgp[nVar];
  // Dimensions             = 3
  // Number of variables    = 14 + #parameters

  // @todo Please implement/augment if required
  double ti = t + 0.5 * dt;
  initialdata_(x, &ti, Qgp);
  for(int m=0; m < nVar; m++) {
        stateOutside[m] = Qgp[m];
  }
  for(int m=nVar; m < nVar+nPar; m++) {
        stateOutside[m] = stateInside[m];
  }
}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit


void DIM::DIMSolver_FV::flux(const double* const Q,double** const F) {
	const int nVar = DIM::AbstractDIMSolver_FV::NumberOfVariables;
  // Dimensions                        = 3
  // Number of variables + parameters  = 14 + 0
  
  // @todo Please implement/augment if required
    if(DIMENSIONS == 2){
		double F_3[nVar];
		pdeflux_(F[0], F[1],F_3, Q);
	}else{
		pdeflux_(F[0], F[1],F[2], Q);
	}
}



void  DIM::DIMSolver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  // @todo Please implement/augment if required
  pdencp_(BgradQ, Q, gradQ);
  //pdencp_lk_(BgradQ, Q, gradQ);
}

double DIM::DIMSolver_FV::riemannSolver(double* fL, double *fR, const double* qL, const double* qR, const double* gradQL, const double* gradQR, const double* cellSize, int normalNonZero) {
  const int numberOfVariables  = DIM::AbstractDIMSolver_FV::NumberOfVariables;
  const int numberOfParameters = DIM::AbstractDIMSolver_FV::NumberOfParameters;
  const int numberOfData       = numberOfVariables+numberOfParameters;
  const int order              = 0;
  const int basisSize          = order+1;
  // Compute the average variables and parameters from the left and the right
  double QavL[numberOfData] = {0.0}; // ~(numberOfVariables+numberOfParameters)
  double QavR[numberOfData] = {0.0}; // ~(numberOfVariables+numberOfParameters)
  
  // std::cout << "opened ---------------------"<< std::endl;
  
    kernels::idx2 idx_QLR(basisSize, numberOfData);
    for (int j = 0; j < basisSize; j++) {
      const double weight = kernels::legendre::weights[order][j];

      for (int k = 0; k < numberOfData; k++) {
        QavL[k] += weight * qL[idx_QLR(j, k)];
        QavR[k] += weight * qR[idx_QLR(j, k)];
      }
    }
	
// Call the Fortran routine
// std::cout << "normalNonZero=" << normalNonZero << std::endl;
//std::cout << "Means done ---------------------"<< std::endl;	
//std::cout << "numberOfVariables=" << numberOfVariables << std::endl;
//std::cout << "numberOfParameters=" << numberOfParameters << std::endl;
//std::cout << "numberOfData=" << numberOfData << std::endl;
//std::cout << "QavR=" << QavR[0] << std::endl;
//hllemriemannsolver_(&basisSize, &normalNonZero, fL,fR,qL, qR,QavL, QavR);
//std::cout << "Fl_before=" << fL[0] << "||" << std::endl;
hllemriemannsolver_(&basisSize, &normalNonZero, fL,fR,qL, qR,QavL, QavR);
//testriemannsolver_(&basisSize, &normalNonZero, fL,fR,qL, qR,QavL, QavR);
return 2;
//std::cout << "closed ---------------------"<< std::endl;	
}

void DIM::DIMSolver_FV::pointSource(const double* const Q,const double* const x,const double t,const double dt, double* const forceVector,int n){
	
}

//void DIM::DIMSolver_FV::multiplyMaterialParameterMatrix(const double* const Q, double* const rhs){}

