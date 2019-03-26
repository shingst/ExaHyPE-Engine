#include "DIMSolver_FV.h"

#include "DIMSolver_FV_Variables.h"
#include "PDE.h"
#include "C2P-DIM.h"
#include "InitialData.h"

#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreQuadrature.h"

tarch::logging::Log GPRDIM::DIMSolver_FV::_log( "GPRDIM::DIMSolver_FV" );

void GPRDIM::DIMSolver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
  // Place here some initialization functions (like read DTM file)
	//std::cout << " ==================================================================================" << std::endl;
	//std::cout << " ==================================================================================" << std::endl;
	//std::cout << " ==================================================================================" << std::endl;
	//std::cout << _maximumMeshSize << std::endl;
	//std::cout << _coarsestMeshLevel << std::endl;
	//std::cout << _coarsestMeshSize << std::endl;
	//std::cout << _maximumAdaptiveMeshDepth << std::endl;
	//std::cout << _maxLevel << std::endl;
	std::cout << " ==================================================================================" << std::endl;
	const int md = exahype::solvers::Solver::getMaximumAdaptiveMeshDepth();
	std::cout << md << std::endl;
}

void GPRDIM::DIMSolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {
  // Dimensions             = 2
  // Number of variables    = 24 + #parameters
  
  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
		int md = exahype::solvers::Solver::getMaximumAdaptiveMeshDepth();
		double cms = exahype::solvers::Solver::getCoarsestMeshSize();
		const int order = 0;
		initialdata_(x, &t, Q,&md,&cms,&order);
  } 
  dynamicrupture_(x, &t, Q);
  // Place here the code for the dynamic rupture
}

void GPRDIM::DIMSolver_FV::eigenvalues(const double* const Q, const int dIndex, double* const lambda) {
  // Dimensions             = 2
  // Number of variables    = 24 + #parameters
  
  // @todo Please implement/augment if required
  double nv[3] = {0.};
  nv[dIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}

void GPRDIM::DIMSolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* const stateOutside) {
	// Local variables
	const int nVar = GPRDIM::AbstractDIMSolver_FV::NumberOfVariables;	
	double Qgp[nVar];
	int md=0;
	double cms=0;
	const int order=0;
	
	double ti = t + 0.5 * dt;
	// Compute the outer state according to the initial condition
	initialdata_(x, &ti, Qgp,&md,&cms,&order);
	// Assign the proper outer state
	for(int m=0; m < nVar; m++) {
        stateOutside[m] = Qgp[m];
	}
}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit


void GPRDIM::DIMSolver_FV::flux(const double* const Q,double** const F) {
	const int nVar = GPRDIM::AbstractDIMSolver_FV::NumberOfVariables;
	if(DIMENSIONS == 2){
		double F_3[nVar];
		pdeflux_(F[0], F[1],F_3, Q);
	}else{
		pdeflux_(F[0], F[1],F[2], Q);
	}
  
}


//You can either implement this method or modify fusedSource
void GPRDIM::DIMSolver_FV::algebraicSource(const double* const Q,double* const S) {
	//S[0] = 0.0;
	//S[1] = 0.0;
	//S[2] = 0.0;
	//S[3] = 0.0;
	//S[4] = 0.0;
	//S[5] = 0.0;
	//S[6] = 0.0;
	//S[7] = 0.0;
	//S[8] = 0.0;
	//S[9] = 0.0;
	//S[10] = 0.0;
	//S[11] = 0.0;
	//S[12] = 0.0;
	//S[13] = 0.0;
	//S[14] = 0.0;
	//S[15] = 0.0;
	//S[16] = 0.0;
	//S[17] = 0.0;
	//S[18] = 0.0;
	//S[19] = 0.0;
	//S[20] = 0.0;
	//S[21] = 0.0;
	//S[22] = 0.0;
	//S[23] = 0.0;
	pdesource_(S, Q);
}

void  GPRDIM::DIMSolver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
	pdencp_(BgradQ, Q, gradQ);
}

double GPRDIM::DIMSolver_FV::riemannSolver(double* fL, double *fR, const double* qL, const double* qR, const double* gradQL, const double* gradQR, const double* cellSize, int normalNonZero) {
  const int numberOfVariables  = GPRDIM::AbstractDIMSolver_FV::NumberOfVariables;
  const int numberOfParameters = GPRDIM::AbstractDIMSolver_FV::NumberOfParameters;
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
hllemriemannsolver_(&basisSize, &direction, fL,fR,qL, qR,QavL, QavR);
//testriemannsolver_(&basisSize, &normalNonZero, fL,fR,qL, qR,QavL, QavR);
return 1;
}



