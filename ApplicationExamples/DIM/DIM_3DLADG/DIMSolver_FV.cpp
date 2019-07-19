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
	readcgfile_();
}

void DIM::DIMSolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {
  // Dimensions             = 3
  // Number of variables    = 14 + #parameters
  
  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
  initialdata_(x, &t, Q);
  }
}

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
  // Dimensions             = 3
  // Number of variables    = 14 + #parameters

  // @todo Please implement/augment if required
  initialdata_(x, &t, stateOutside);
}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit


void DIM::DIMSolver_FV::flux(const double* const Q,double** const F) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 14 + 0
  
  // @todo Please implement/augment if required
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
  
  F[2][0] = 0.0;
  F[2][1] = 0.0;
  F[2][2] = 0.0;
  F[2][3] = 0.0;
  F[2][4] = 0.0;
  F[2][5] = 0.0;
  F[2][6] = 0.0;
  F[2][7] = 0.0;
  F[2][8] = 0.0;
  F[2][9] = 0.0;
  F[2][10] = 0.0;
  F[2][11] = 0.0;
  F[2][12] = 0.0;
  F[2][13] = 0.0;
  
}



void  DIM::DIMSolver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  // @todo Please implement/augment if required
  pdencp_(BgradQ, Q, gradQ);
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
