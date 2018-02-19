#include "DIMSolver_ADERDG.h"

#include "DIMSolver_ADERDG_Variables.h"

// User defined calls
#include "PDE.h"
#include "InitialData.h"
#include "C2P-DIM.h"
// Used for the rieman-solver part
#include "peano/utils/Loop.h"
#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreQuadrature.h"

tarch::logging::Log DIM::DIMSolver_ADERDG::_log( "DIM::DIMSolver_ADERDG" );


void DIM::DIMSolver_ADERDG::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
  //  std::cout << " ==================================================================================" << std::endl;
  //	readcgfile_(&_domainOffset[0],&_domainSize[0]);
}

void DIM::DIMSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 14 + 0
  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
    initialdata_(x, &t, Q);
  }
}

void DIM::DIMSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
const int nVar = DIM::AbstractDIMSolver_ADERDG::NumberOfVariables;
  const int order = DIM::AbstractDIMSolver_ADERDG::Order;
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
}

exahype::solvers::Solver::RefinementControl DIM::DIMSolver_ADERDG::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void DIM::DIMSolver_ADERDG::eigenvalues(const double* const Q,const int d,double* lambda) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 14 + 0
  
  // @todo Please implement/augment if required
  double nv[3] = {0.};
  nv[d] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}


void DIM::DIMSolver_ADERDG::flux(const double* const Q,double** F) {
	const int nVar = DIM::AbstractDIMSolver_ADERDG::NumberOfVariables;
  // Dimensions                        = 3
  // Number of variables + parameters  = 14 + 0
  
  // @todo Please implement/augment if required
    if(DIMENSIONS == 2){
		double F_3[nVar];
		pdeflux_(F[0], F[1],F_3, Q);
	}else{
		pdeflux_(F[0], F[1],F[3], Q);
	}
}

void DIM::DIMSolver_ADERDG::mapDiscreteMaximumPrincipleObservables(
    double* observables,const int numberOfObservables,
    const double* const Q) const {
  assertion(numberOfObservables==1);
  ReadOnlyVariables vars(Q);

  observables[0]=Q[12]; //extract alpha
}

bool DIM::DIMSolver_ADERDG::isPhysicallyAdmissible(
  const double* const solution,
  const double* const observablesMin,const double* const observablesMax,const int numberOfObservables,
  const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx,
  const double t, const double dt) const {
  int limvalue;
  // Variant 1 (cheapest, currently works only in 2D)
  //  double outerRadius = 1.25*0.25;
  //  double innerRadius = 0.75*0.25;
  //  double radiusSquared = (center[0])*(center[0])+(center[1])*(center[1])+(center[2])*(center[2]);
  
  //  if (
  //    radiusSquared<outerRadius*outerRadius &&
  //    radiusSquared>=innerRadius*innerRadius
  //  ) {
  //    return false;
  //  }  
//	return true;
// if ((center[0]-0.5)*(center[0]-0.5)+(center[1]-0.5)*(center[1]-0.5)<0.25*dx[0]*dx[0]) return false;
 //  Worked for the sphere
 //if ((center[0])*(center[0])+(center[1])*(center[1])+(center[2])*(center[2])<0.5*0.5) return false;
 // if (observablesMin[0] <= 0.0) return false;
  //if (observablesMax[0] >= 1.0) return false;
 // return true;
  // Slow bug has to works
  //pdelimitervalue_(&limvalue,xx);
  pdelimitervalue_(&limvalue,&center[0],&numberOfObservables, observablesMin, observablesMax);
  if(limvalue>0){
	  return false;
  }else{
	  return true;
  };
}

void  DIM::DIMSolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  // @todo Please implement/augment if required
   pdencp_(BgradQ, Q, gradQ);
/*    BgradQ[0] = 0.0;
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
  */
}


