#include "DIMSolver_FV.h"

#include "DIMSolver_FV_Variables.h"
#include "PDE.h"

#include "InitialData.h"

#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreQuadrature.h"


tarch::logging::Log DIM::DIMSolver_FV::_log( "DIM::DIMSolver_FV" );


void DIM::DIMSolver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

void DIM::DIMSolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {
  // Dimensions             = 2
  // Number of variables    = 14 + #parameters
  
  // @todo Please implement/augment if required
  // Fortran
  if (tarch::la::equals(t,0.0)) {
    initialdata_(x, &t, Q);
  }
}

void DIM::DIMSolver_FV::eigenvalues(const double* const Q, const int dIndex, double* const lambda) {
  // Dimensions             = 2
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
  // Dimensions             = 2
  // Number of variables    = 14 + #parameters

  // @todo Please implement/augment if required
initialdata_(x, &t, stateOutside);
}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit

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
      const double weight = kernels::gaussLegendreWeights[order][j];

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
hllemriemannsolver_(&basisSize, &normalNonZero, fL,fR,qL, qR,QavL, QavR);
//std::cout << "closed ---------------------"<< std::endl;

	return 1.0; // if you don't want to return a proper maximum eigenvalue
}

void  DIM::DIMSolver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  // @todo Please implement/augment if required
    pdencp_(BgradQ, Q, gradQ);
}

