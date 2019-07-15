#include "MHDSolver_ADERDG.h"
//#include "fortran.h" _ltob

#include "InitialData.h"
#include "BoundaryConditions.h"
#include "PDE.h"

#include <memory>
#include <cstring>
#include <stdio.h>
#include "kernels/GaussLegendreBasis.h"
#include "kernels/KernelUtils.h" // matrix indexing

// storage
InitialDataHandler idfunc;
BoundaryConditionHandler bcfunc;

void MHD::MHDSolver_ADERDG::init(std::vector<std::string>& cmdargs) {
	/**
	 * InitialData with consistent BC can be set by environment:
	 *    export EXAHYPE_INITIALDATA="Blast"
	 * or by command line:
	 *    ./ExaHyPE-MHD ../MHD_Limiting.exahype Blast
	 * You can also set the default_id to trigger another ID.
	 *
	 * Obey that you need the matching grid for your simulation.
	 **/

	static tarch::logging::Log _log("MHDSolver_ADERDG::init");
	// Don't rely on exahype specfile parameters, use Environment
	// variables or command line variables instead.
	const char* env_id = std::getenv("EXAHYPE_INITIALDATA");
	const char* cmd_id = cmdargs.size()>1 ? cmdargs[1].c_str() : nullptr;
	std::string default_id = "Jet";

	std::string id, bc;	
	std::map<std::string, InitialDataHandler> ids;
	std::map<std::string, BoundaryConditionHandler> bcs;
	
	ids["Jet"] = initialjet_;
	bcs["Jet"] = boundaryjet_;

	ids["AlfenWave"] = initialalfenwave_;
	bcs["AlfenWave"] = BoundaryAlfenWave;
	
	ids["Rotor"] = initialrotor_;
	bcs["Rotor"] = boundaryoutflow_;
	
	ids["Blast"] = initialblast_;
	bcs["Blast"] = boundaryoutflow_;
	
	ids["OrsagTang"] = initialorsagtang_;
	bcs["OrsagTang"] = boundaryoutflow_;
	
	ids["ShockTube"] = initialshocktube_;
	bcs["ShockTube"] = boundaryoutflow_;
	
	if(env_id) id = env_id;
	else if(cmd_id) id = cmd_id;
	else id = default_id;
	
	if(id==default_id)
		logError("ID", "Loading default Initial Data: " << id);
	if(!ids.count(id) && !bcs.count(id)) {
		logError("ID", "Cannot understand requested ID: " << id);
		exit(-1);
	} else {
		logInfo("ID", "Will load requested Initial Data: " << id);
	}
	
	// TODO: Lowercase key lookup
	idfunc = ids[id];
	bcfunc = bcs[id];
}

void MHD::MHDSolver_ADERDG::flux(const double* const Q, double** const F) {
  // Caveats: Fortran accepts a uniform array of size (nVar*nDim), however C passes an array of pointers.
  // This Fortran interface works only if F is a continous array and F[1]==F[nDim+1] etc!
  pdeflux_(F[0], Q);
}

void MHD::MHDSolver_ADERDG::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* const lambda) {
  double nv[3] = {0.};
  nv[normalNonZeroIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}

void MHD::MHDSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
  if (tarch::la::equals(t, 0.0)) {
    idfunc(x, Q);
  }
}

void MHD::MHDSolver_ADERDG::algebraicSource(const double* const Q,double* const S) {
  pdesource_(S, Q);
}

exahype::solvers::Solver::RefinementControl MHD::MHDSolver_ADERDG::refinementCriterion(const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}

void MHD::MHDSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  // pass to Fortran/generic function pointer
  double nv[3] = {0.};
  nv[normalNonZero] = 1;
  bcfunc(x, &t, &dt, &faceIndex, nv, fluxIn, stateIn, fluxOut, stateOut);
}


// only evaluated in Limiting context
void MHD::MHDSolver_ADERDG::mapDiscreteMaximumPrincipleObservables(double* const observables, const double* const Q) const {
  assertion(NumberOfDMPObservables==2);
  observables[0] = Q[0];
  observables[1] = Q[4];
}


// only evaluated in Limiting context
bool MHD::MHDSolver_ADERDG::isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,
      const bool wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t) const {
  if (observablesMin[0] < 0.0) return false;
  if (observablesMin[1] < 0.0) return false;

  return true;
}



