#include "GPRSolver_ADERDG.h"

#include "GPRSolver_ADERDG_Variables.h"

// User defined calls
#include "PDE.h"
#include "InitialData.h"
#include "C2P-GPR.h"

#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreQuadrature.h"

tarch::logging::Log GPR::GPRSolver_ADERDG::_log( "GPR::GPRSolver_ADERDG" );


void GPR::GPRSolver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
    const int order = GPR::AbstractGPRSolver_ADERDG::Order;
    //inittecplot_(&order,&order);
}

void GPR::GPRSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 17 + 0
  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
    initialdata_(x, &t, Q);
  }
}

void GPR::GPRSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
 const double * const fluxIn,const double* const stateIn, const double* const gradStateIn,
					   double *fluxOut,double* stateOut) {
  const int nVar = GPR::AbstractGPRSolver_ADERDG::NumberOfVariables;
  const int order = GPR::AbstractGPRSolver_ADERDG::Order;
  const int basisSize = order + 1;
  const int nDim = DIMENSIONS;

  double Qgp[nVar],*F[nDim], Fs[nDim][nVar];

  std::memset(stateOut, 0, nVar * sizeof(double));
  std::memset(fluxOut, 0, nVar * sizeof(double));
  
  for(int dd=0; dd<nDim; dd++) F[dd] = Fs[dd];
  
  for(int i=0; i < basisSize; i++)  { // i == time
     const double weight = kernels::gaussLegendreWeights[order][i];
     const double xi = kernels::gaussLegendreNodes[order][i];
     double ti = t + xi * dt;

     initialdata_(x, &ti, Qgp);
    //pdeflux_(F[0], F[1], F[2], Qgp);
	flux(Qgp, F);
     for(int m=0; m < nVar; m++) {
        stateOut[m] += weight * Qgp[m];
        fluxOut[m] += weight * Fs[normalNonZero][m];
     }
  }
  /*
	for(int m=0; m < nVar; m++) {
	stateOut[m] = stateIn[m];
	fluxOut[m] = fluxIn[m];
	}
	*/
}

exahype::solvers::Solver::RefinementControl GPR::GPRSolver_ADERDG::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void GPR::GPRSolver_ADERDG::eigenvalues(const double* const Q,const int d,double* const lambda) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 17 + 0
  
  // @todo Please implement/augment if required
  double nv[3] = {0.};
  nv[d] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}


void GPR::GPRSolver_ADERDG::flux(const double* const Q,double** const F) {
	const int nVar = GPR::AbstractGPRSolver_ADERDG::NumberOfVariables;
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
void GPR::GPRSolver_ADERDG::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
  // @todo Please implement/augment if required
  pdesource_(S, Q);
}

void  GPR::GPRSolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  // @todo Please implement/augment if required
  pdencp_(BgradQ, Q, gradQ);
}

void GPR::GPRSolver_ADERDG::mapDiscreteMaximumPrincipleObservables(double* const observables, const double* const Q) const {
  assertion(NumberOfDMPObservables==1);
  ReadOnlyVariables vars(Q);

  observables[0]=Q[0]; //extract alpha
}

bool GPR::GPRSolver_ADERDG::isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,
      const bool wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t) const {
  int limvalue;
  //pdelimitervalue_(&limvalue,&center[0]);
  int NumberOfDMPObservables=1;
  pdelimitervalue_(&limvalue,&center[0],&NumberOfDMPObservables, observablesMin, observablesMax);
  if(limvalue>0){
	  return false;
  }else{
	  return true;
  };
}
