#ifndef __EXAHYPE_Z4_EXACT_F_ID__
#define __EXAHYPE_Z4_EXACT_F_ID__

#include "AbstractCCZ4Solver_ADERDG.h"
#include <stdio.h>
#include "Parameters/mexa.h"

/**
 * Once we move back again from Fortran Initial Data handling to C++, this structure
 * plays an even more important role than e.g. it's GRMHD counterpart.
 **/
struct InitialDataCode {
	virtual void readParameters(const mexa::mexafile& parameters) {} // by default do nothing.
	virtual void prepare() {} // Run a ID code. Parameters should be set before.

	virtual void adjustPointSolution(const double* const x,const double t,double* Q) {
		printf("ERROR: Not implemented pointSolution for ID\n");
		exit(-2);
	}
	virtual void adjustPatchSolution(
		const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
		const tarch::la::Vector<DIMENSIONS, double>& dx,
		const double t, const double dt, double* luh) {
		printf("ERROR: Not implemented pointSolution for ID\n");
		exit(-2);
	}
	
	/// This is the actual initial data registration
	static InitialDataCode* getInstanceByName(std::string name);
};

// a shorthand for pointSolutions.
void InitialData(const double* x, double t, double* Q);

// global initial data
class GlobalInitialData {
	tarch::logging::Log _log;
public:
	/// The pointer is by default null, meaning no ID can be used. The setIdByName method
	/// sets this id pointer, finally.
	InitialDataCode* id;
	std::string name; ///< set by setIdByName
	
	GlobalInitialData(InitialDataCode* _id, std::string _name) : _log("GlobalInitialData"), id(_id), name(_name) {}
	
	/// Returns true in case of success
	bool setIdByName(const std::string& name);
	
	/// Tries to read ID from parameters. In case of problems, raises exceptions.
	void setByParameters(const mexa::mexafile& parameters);
	
	/// Gives a singleton instance
	static GlobalInitialData& getInstance();
	
	// as a shorthand
	static InitialDataCode& getInitialDataCode() { return *(getInstance().id); }
};

/*
/// This is hacky
#include "kernels/aderdg/generic/c/3d/solutionAdjustment.cpph"
class PointwiseInitialData : public InitialDataCode {
public:
	constexpr int NumberOfVariables  = AbstractCCZ4Solver_ADERDG::NumberOfVariables;
	constexpr int NumberOfParameters = AbstractCCZ4Solver_ADERDG::NumberOfParameters;
	constexpr int Order              = AbstractCCZ4Solver_ADERDG::Order;
  
	virtual void adjustPatchSolution(
	  const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
	  const tarch::la::Vector<DIMENSIONS, double>& dx,
	  const double t, const double dt, double* luh) {
		kernels::aderdg::generic::c::solutionAdjustment<PointwiseInitialData>(
			*this, luh, cellCentre, dx, t, dt, luh);
	}
}

// should be:
template <void PDEInitialData(const double* const x, double* Q, const double* t)>
struct FortranPointwiseID : public PointwiseInitialData {
	virtual void adjustPointSolution(const double* const x, double* Q) { double t=0; PDEInitialData(x,Q,&t); }
};

typedef FortranPointwiseID<initialfield_> Z4_InitialFortranID;
*/


/**
 * The interface to the Fortran Initial Data Functions
 **/
extern "C" {

// This is a Fortran function which allows to use any METRIC_ fortran function
// to set basically the ADM quantities and computes then the CCZ4 quantities
// from that. That is, see Fortran/Init.f90 for the switch for different ID.
void initialfield_(double* u0, const double* const xGP, const double* tGP);

// should be:
struct FortranPointwiseID : public InitialDataCode {
	FortranPointwiseID(int ic_type);
	virtual void adjustPointSolution(const double* const x, const double t, double* Q) {
		initialfield_(Q, x, &t); }
};

}/* extern "C" */



#ifdef TWOPUNCTURES_AVAILABLE

// C interface to TwoPunctures:
namespace TP { class TwoPunctures; } // Forward declaration
#include "tarch/la/Vector.h"
#include "peano/utils/Globals.h" // DIMENSIONS
class ImportedTwoPunctures : public InitialDataCode {
public:
	static constexpr int numberOfVariables = CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables;
	static constexpr int order = CCZ4::AbstractCCZ4Solver_ADERDG::Order;
	static constexpr int basisSize = order + 1;
	static constexpr int basisSize3 = basisSize*basisSize*basisSize;
	TP::TwoPunctures* tp;
	
	ImportedTwoPunctures();

	// Find a consistent solution of Einsteins Equations on a spectral domain
	virtual void prepare();
	
        // this sets the InitialData on points using finite differencing to determine
        // the derivatives. This is also used at the boundaries, right...
	virtual void adjustPointSolution(const double* const x,const double t,double* Q);
        
	// Caveats when using this function: It does only set gij, kij, alpha.
	void Interpolate(const double* const x, double* Q);
	void dInterpolate(const double* const xc, double *Qx, double *Qy, double *Qz);

	// This function computes the derived quantities AA, BB, DD
	void ComputeAuxillaries(const double* const x, double* Q, double *Qx, double *Qy, double *Qz);
	
	// This is the actual InitialData function for a whole patch, computing
	// the derivatives Qx,Qy,Qz and using Interpolate(x,Q,Qx,Qy,Qz) on all DOF.
	virtual void adjustPatchSolution(
		const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
		const tarch::la::Vector<DIMENSIONS, double>& dx,
		const double t, const double dt, double* luh);
};

// Fortran interface to TwoPunctures, using a singleton
extern "C" {
void twopunctures_interpolate_(const double* x, double* Q);
}
#endif


#endif /* __EXAHYPE_Z4_EXACT_F_ID__ */
