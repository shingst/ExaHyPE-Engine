// new BC code written 2017-12-13 for GRMHD by SvenK.
// Superseding BoundaryConditions_{ADERDG,FV}.h.
// Implementing the new parameter code.
#ifndef __GRMHD_BOUNDARY_CONDITIONS__
#define __GRMHD_BOUNDARY_CONDITIONS__

#include "GRMHDSolver_ADERDG.h"
#include "GRMHDSolver_FV.h"
#include "Parameters/mexa.h"

namespace GRMHD{
  class BoundaryConditions;
  class GlobalBoundaryConditions; // Registration
}

// The difference between the ADERDG and FV boundary signature is that in ADERDG, also
// the fluxIn and fluxOut pointers are given. Therefore, the ADERDG signature is the
// broader one and we use it all the time.
#define BOUNDARY_SIGNATURE \
  const double* const x,const double t,const double dt,const int faceIndex,const int d,\
  const double * const fluxIn,const double* const stateIn,\
  double *fluxOut, double* stateOut

// Calls: Needed for calling a function with a specific signature. This is what you
// should use when using the GRMHD_BoundaryConditions class.
#define ADERDG_BOUNDARY_CALL  x, t, dt, faceIndex, d, fluxIn,  stateIn, fluxOut, stateOut
#define FV_BOUNDARY_CALL      x, t, dt, faceIndex, d, nullptr, stateIn, nullptr, stateOut

// And this is what we pass internally:
#define BOUNDARY_CALL ADERDG_BOUNDARY_CALL

/**
 * A unified BoundaryConditions class which in principle does everything for the DG solver
 * but also can work for the FV solver.
 **/
class GRMHD::BoundaryConditions {
	tarch::logging::Log _log;
public:	
	typedef GRMHD::GRMHDSolver_ADERDG  Solver_ADERDG;
	typedef GRMHD::BoundaryConditions me;
	typedef void (me::*boundarymethod)(BOUNDARY_SIGNATURE);
	
	/// The solver class reference is only needed in the ADERDG case to set
	/// the fluxes.
	Solver_ADERDG* aderdg_solver;
	
	// this could be improved by using boundarymethod side[6]
	// and having apply defined just as  side[faceIndex](...);.
	boundarymethod left, right, front, back, bottom, top;
	
	/// Constructor for the FV solver
	BoundaryConditions();
	
	/// Constructor for the ADERDG solver.
	BoundaryConditions(Solver_ADERDG* _aderdg_solver);
	
	bool allFacesDefined();
	
	/// Apply by looking up the saved boundarymethods
	/// This is the main method how to use an instance.
	void apply(BOUNDARY_SIGNATURE);
	
	/***************************************************************************
	 * Boundary Methods.
	 * 
	 * To add your own, do the following:
	 * 
	 *  1) Implement it with a method void yourName(ADERDG_BOUNDARY_SIGNATURE)
	 *  2) Add yourName to the parseFromString thing
	 * 
	 ***************************************************************************/
	
	bool setFromParameters(const mexa::mexafile& constants, bool raiseOnFailure=true);

	/**
	 * assign a boundary method by some string value
	 * @returns pinter if value could be correctly parsed, in case of error nullptr.
	 **/
	static boundarymethod parseFromString(const std::string& value);
	
	bool fluxesRequested(BOUNDARY_SIGNATURE) {
		return (fluxOut != nullptr); // could also check for fluxIn != nullptr.
	}
	/// Just an alias.
	bool workingForADERDG(BOUNDARY_SIGNATURE) {
		return fluxesRequested(ADERDG_BOUNDARY_CALL);
	}
	
	/**
	 * Compute consistent fluxes.
	 **/
	void deriveAderdgFlux(BOUNDARY_SIGNATURE);

	/**
	 * The vacuum is problem dependent, thus this is really MHD vacuum in
	 * Minkowski spacetime for GRMHD.
	 **/
	void vacuum(BOUNDARY_SIGNATURE);
	
	/// Outflow / Copy BC
	void outflow(BOUNDARY_SIGNATURE);
	
	/// Reflection/Hard wall BC
	void reflective(BOUNDARY_SIGNATURE);

	/// Exact BC
	void exact(BOUNDARY_SIGNATURE);
};


/**
 * A singleton BoundaryConditions registration. This is useful when you only
 * have one solver and don't want to store BoundaryConditions instances on
 * your own.
 **/
class GRMHD::GlobalBoundaryConditions {
	tarch::logging::Log _log;
public:
	
	GRMHD::BoundaryConditions bc;
	
	GlobalBoundaryConditions() : _log("GlobalBoundaryConditions"), bc() {}
	
	/// Gives a singleton instance
	static GlobalBoundaryConditions& getInstance();
	
	/// initialize. Chainable.
	GlobalBoundaryConditions& initializeDG(GRMHD::GRMHDSolver_ADERDG* aderdg_solver) {
		bc.aderdg_solver = aderdg_solver; return *this; }
	GlobalBoundaryConditions& initializeFV(GRMHD::GRMHDSolver_FV* fv_solver) { return *this; }
	
	/// Read the parameters for the bc.
	void readParameters(const mexa::mexafile& parameters) { bc.setFromParameters(parameters); }
	
	// shorthands
	void apply(BOUNDARY_SIGNATURE) { bc.apply(BOUNDARY_CALL); }
};

#endif /* __GRMHD_BOUNDARY_CONDITIONS__ */
