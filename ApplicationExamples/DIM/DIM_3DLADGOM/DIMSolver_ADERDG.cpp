#include "DIMSolver_ADERDG.h"

#include "DIMSolver_ADERDG_Variables.h"

// User defined calls
#include "PDE.h"
#include "InitialData.h"
#include "C2P-DIM.h"
#include "TECPLOTinterface.h"
// Used for the rieman-solver part
#include "peano/utils/Loop.h"
#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreBasis.h"

tarch::logging::Log DIM::DIMSolver_ADERDG::_log( "DIM::DIMSolver_ADERDG" );


void DIM::DIMSolver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
  //  std::cout << " ==================================================================================" << std::endl;
  //	readcgfile_(&_domainOffset[0],&_domainSize[0]);
  const int order = DIM::AbstractDIMSolver_ADERDG::Order;
  inittecplot_(&order,&order);
}

void DIM::DIMSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 14 + 0
  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
    initialdata_(x, &t, Q);
	//Q[13]=1.0;
  }
  Q[13]=1.0;
}

void DIM::DIMSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut)
const int nVar = DIM::AbstractDIMSolver_ADERDG::NumberOfVariables;
  const int order = DIM::AbstractDIMSolver_ADERDG::Order;
  const int basisSize = order + 1;
  const int nDim = DIMENSIONS;

  double Qgp[nVar],*F[nDim], Fs[nDim][nVar];

  std::memset(stateOut, 0, nVar * sizeof(double));
  std::memset(fluxOut, 0, nVar * sizeof(double));
  
  for(int dd=0; dd<nDim; dd++) F[dd] = Fs[dd];
  
  for(int i=0; i < basisSize; i++)  { // i == time
     const double weight = kernels::legendre::weights[order][i];
     const double xi = kernels::legendre::nodes[order][i];
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
/*
void DIM::DIMSolver_ADERDG::algebraicSource(const double* const Q,double* const S) {
	const int nVar = DIM::AbstractDIMSolver_ADERDG::NumberOfVariables;
  // @todo Please implement/augment if required
  for(int m=0; m < nVar; m++) {
	S[m]=0.0;  
  }
}
*/
exahype::solvers::Solver::RefinementControl DIM::DIMSolver_ADERDG::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  if ( level > getCoarsestMeshLevel() )
    return exahype::solvers::Solver::RefinementControl::Erase;
  else return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void DIM::DIMSolver_ADERDG::eigenvalues(const double* const Q,const int d,double* const lambda) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 14 + 0
  
  // @todo Please implement/augment if required
  double nv[3] = {0.};
  nv[d] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}


void DIM::DIMSolver_ADERDG::flux(const double* const Q,double** const F) {
	const int nVar = DIM::AbstractDIMSolver_ADERDG::NumberOfVariables;
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

void DIM::DIMSolver_ADERDG::mapDiscreteMaximumPrincipleObservables(double* const observables, const double* const Q) const {
  assertion(NumberOfDMPObservables==2);
  ReadOnlyVariables vars(Q);

  observables[0]=Q[12]; //extract alpha
  //observables[1]=Q[14];
}

bool DIM::DIMSolver_ADERDG::isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,
      const bool wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t) const {
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
  //if (tarch::la::equals(t,0.0)) {
  
  // Another variant: Just check solution at GLob nodes
  //  double obsMin[(NumberOfVariables+NumberOfParameters)];
  //  double obsMax[(NumberOfVariables+NumberOfParameters)];
  //  kernels::limiter::generic::c::computeMinimumAndMaximumValueAtGaussLobattoNodes<Order+1,NumberOfVariables+NumberOfParameters>(
  //      solution,obsMin,obsMax);
  //  pdelimitervalue_(&limvalue,&center[0],&NumberOfDMPObservables, obsMin, obsMax);

  //pdelimitervalue_(&limvalue,&center[0],&NumberOfDMPObservables, observablesMin, observablesMax);
  if (tarch::la::equals(t,0.0)) {
	  //pdelimitervalue_(&limvalue,&center[0],&NumberOfDMPObservables, observablesMin, observablesMax);
	  pdegeometriclimitervalue_(&limvalue,&center[0]);
	  if(limvalue>0){
		  return false;
	  }else{
		  return true;
	  };
  }else{
		return !wasTroubledInPreviousTimeStep;
  };
  /*}else{
	if(!tarch::la::equals(observablesMax[1],observablesMin[1])){
		std::cout << "different max an min limiter values " << observablesMax[1] << " " <<observablesMin[1] << std::endl;
	}
	if (tarch::la::equals(observablesMax[1],0.0)) {  
		  return true;
	}else{ 
		return false;
	}
  }*/
}

void  DIM::DIMSolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
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


