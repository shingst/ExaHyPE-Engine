#include "GRMHDSolver_ADERDG.h"

#include "GRMHDSolver_ADERDG_Variables.h"
#include "InitialData/InitialData.h"
#include "Fortran/PDE.h"

#include <cstring> // memset
#include <cmath> //min,max
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreBasis.h"

#include "exahype/disableOptimization.h" // we experience compiler bugs sometimes.

#include "BoundaryConditions/BoundaryConditions_ADERDG.h"


const double excision_radius = 1.0;
constexpr int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;
constexpr int order = GRMHD::AbstractGRMHDSolver_ADERDG::Order;
constexpr int basisSize = order + 1;
constexpr int nDim = DIMENSIONS;

tarch::logging::Log GRMHD::GRMHDSolver_ADERDG::_log("GRMHDSolver_ADERDG");

typedef BoundaryConditions<GRMHD::GRMHDSolver_ADERDG> ADERDG_BC;
ADERDG_BC* abc;

// enable nan tracker
#include <fenv.h>

void GRMHD::GRMHDSolver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // NAN checker
//  feenableexcept(FE_INVALID | FE_OVERFLOW);  // Enable all floating point exceptions but FE_INEXACT
	
  // Todo: Move this to specfile once we have working constants.
   std::string id_default = "Fortran";
   std::string bc_default = "left:exact,right:exact,top:exact,bottom:exact,front:exact,back:exact";
  //std::string id_default = "TOVSolver";
  //std::string bc_default = "left:exact,right:exact,top:exact,bottom:exact,front:exact,back:exact";

  // alternatives:
  //std::string id_default = "RNSID";
  //std::string bc_default = "left:refl,right:exact,bottom:refl,top:exact,front:refl,back:exact";

  // try to obtain requested initial data and boundary conditions from the
  // environment variables, as the specfile parameter system is still broken.
  std::string tid = getenv("EXAHYPE_ID") ? getenv("EXAHYPE_ID") : id_default;
  std::string tbc = getenv("EXAHYPE_BC") ? getenv("EXAHYPE_BC") : bc_default;

  if(!prepare_id(tid)) {
	  logError("prepare_id", "Could not setup Initial Data '" << tid << "', probably misspelled.");
	  std::abort();
  }

  abc = new ADERDG_BC(this);
  //if(!abc->setFromSpecFile<exahype::parser::ParserView>(constants)) {
  if(!abc->setFromSpecFile<StringMapView>(StringMapView(tbc))) {
	logError("boundaryValues", "Some Boundary faces are missing in Specfile. Need: left,right,top,bottom,front,back. Got:" << tbc);
	std::abort();
  }
}

void __attribute__((optimize("O0"))) initialData(const double* const x,const double t,const double dt,double* const Q) {
  id->Interpolate(x, t, Q);
  //printf("Interpoalted at x=[%f,%f,%f], t=%f, Q0=%f\n", x[0],x[1],x[2], t, Q[0]);
  for(int i=0; i<nVar; i++) {
    if(!std::isfinite(Q[i])) {
      printf("NAN in i=%d at t=%f, x=[%f,%f,%f], Q[%d]=%f\n", i, t, x[0],x[1],x[2], i, Q[i]);
    }
  }
}


void __attribute__((optimize("O0"))) GRMHD::GRMHDSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
  bool insideExcisionBall = false;
  bool hastoadjust = tarch::la::equals(t,0.0) || insideExcisionBall;

  using namespace tarch::la;
  if(equals(t,0.0)) {
    initialData(x, t, dt, Q); 
//    if( (x[1] > -0.1) && (x[1]  < 0.1)) {
//      if( (x[2] > -0.1) && (x[2]  < 0.1)) {
//        printf("adjusting ADERDG:  Q[0]=%.5e , x,y,z = %f,%f,%f, t= %f\n",Q[0],x[0],x[1],x[2],t);
//      }
//    }
  }
}

void __attribute__((optimize("O0"))) GRMHD::GRMHDSolver_ADERDG::eigenvalues(const double* const Q,const int d,double* const lambda) {
  double nv[3] = {0.};
  nv[d] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}


void GRMHD::GRMHDSolver_ADERDG::flux(const double* const Q,double** const F) {
  pdeflux_(F[0], F[1], (DIMENSIONS==3)?F[2]:nullptr, Q);
}

// Source is zero
/*
void __attribute__((optimize("O0"))) GRMHD::GRMHDSolver_ADERDG::algebraicSource(const double* const Q,double* const S) {
  pdesource_(S, Q);
}
*/


void GRMHD::GRMHDSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
	 // for debugging, to make sure BC are set correctly
	double snan = std::numeric_limits<double>::signaling_NaN();
	double weird_number = -1.234567;
	std::memset(stateOut, weird_number, nVar * sizeof(double));
	std::memset(fluxOut,  weird_number, nVar * sizeof(double));
	
	//abc->apply(ADERDG_BOUNDARY_CALL);
	//abc->exact(ADERDG_BOUNDARY_CALL);
	
	/////// EXACT
	  // employ time-integrated exact BC for AlfenWave.
  
  double Qgp[nVar], Fs[nDim][nVar], *F[nDim];
  for(int dd=0; dd<nDim; dd++) F[dd] = Fs[dd];
  // zeroise stateOut, fluxOut
  for(int m=0; m<nVar; m++) {
    stateOut[m] = 0;
    fluxOut[m] = 0;
  }
  for(int i=0; i < basisSize; i++)  { // i == time
    const double weight = kernels::legendre::weights[order][i];
    const double xi = kernels::legendre::nodes[order][i];
    double ti = t + xi * dt;

    initialData(x, ti, dt, Qgp);
    flux(Qgp, F);
    
    for(int m=0; m < nVar; m++) {
      stateOut[m] += weight * Qgp[m];
      fluxOut[m] += weight * Fs[normalNonZero][m];
    }
  }
  ///// EXACT
	
	
  for(int i=0; i<NumberOfVariables; i++) {
	if(!std::isfinite(stateOut[i])) {
		printf("BoundaryValues stateOut NAN in i=%d at t=%f, x=[%f,%f,%f], stateOut[%d]=%f\n", i, t, x[0],x[1],x[2], i, stateOut[i]);
	}
	if(!std::isfinite(fluxOut[i])) {
		printf("BoundaryValues NAN in i=%d at t=%f, x=[%f,%f,%f], fluxOut[%d]=%f\n", i, t, x[0],x[1],x[2], i, fluxOut[i]);
	}
  }
}

bool isInRefinementZone(const tarch::la::Vector<DIMENSIONS,double>& center){
	double radius = 8.12514;
	// lower left, upper right radius of cell
	double cen = tarch::la::norm2(center);
	double dr = 0.5;
  bool shouldRefine = (cen > (radius -dr) ) && ( cen  <= (radius+dr) ); 
  return shouldRefine;
}


exahype::solvers::Solver::RefinementControl GRMHD::GRMHDSolver_ADERDG::refinementCriterion(const double* const luh,
        const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
    // @todo Please implement/augment if required

    if(isInRefinementZone(center))
        return exahype::solvers::Solver::RefinementControl::Refine;
    return exahype::solvers::Solver::RefinementControl::Keep;
}

// only evaluated in Limiting context
/*
void GRMHD::GRMHDSolver_ADERDG::mapDiscreteMaximumPrincipleObservables(double* const observables, const double* const Q) const {
  assertion(NumberOfDMPObservables==2);

  observables[0] = Q[0]; // rho
  observables[1] = Q[4]; // dens
}
*/


/*
void GRMHD::GRMHDSolver_ADERDG::mapDiscreteMaximumPrincipleObservables(
  double* const observables, const int NumberOfVariables,
  const double* const Q) const {
  for (int i = 0; i < NumberOfVariables; ++i) {
    observables[i] = Q[i];
  }
}
*/



bool GRMHD::GRMHDSolver_ADERDG::isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,
      const bool wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t) const {

	double radius = 8.12514;
	double cen = tarch::la::norm2(center);
	double dr = 1.5;

  bool shouldLimit = (cen > (radius -dr) ) && ( cen  <= (radius+dr) ); 
  // return TRUE if the cell does not need limited
	return !shouldLimit;

}

int getGeometricLoadBalancingWeight(
        const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS,double>& cellSize) {
  double cen = tarch::la::norm2(cellCentre);
  double dr = std::max(0.1, tarch::la::norm2(cellSize));
  if(cen < 8.12514+dr)
    return 27;
  return 1;
}





void __attribute__((optimize("O0"))) GRMHD::GRMHDSolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  pdencp_(BgradQ, Q, gradQ);
  
//  for(int i=0; i<NumberOfVariables; i++) {
//	if(!std::isfinite(BgradQ[i])) {
//		printf("NCP NAN in BgradQ[%d]=>%f\n", i, BgradQ[i]);
//		for(int j=0; j<NumberOfVariables; j++) {
//			printf("Q[%d]=%f\n", j, Q[j]);
//			printf("BgradQ[%d]=%f\n", j, BgradQ[j]);
//		}
//	}
//  }
}

/*
void GRMHD::GRMHDSolver_ADERDG::coefficientMatrix(const double* const Q,const int d,double* const Bn) {
  // new scheme has no coefficient matrix
  static tarch::logging::Log _log("GRMHDSolver");
  logError("coefficientMatrix()", "ADERDG Coefficient Matrix invoked");
  exit(-2);

  double nv[3] = {0.};
  nv[d] = 1;
  pdematrixb_(Bn, Q, nv);
}
*/
