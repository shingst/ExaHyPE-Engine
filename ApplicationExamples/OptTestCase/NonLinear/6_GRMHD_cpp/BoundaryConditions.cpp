
#include "BoundaryConditions.h"

#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreBasis.h"

#include "InitialData/InitialData.h"
#include "PDE/PDE.h"
#include "Parameters/mexa.h"

/// REGISTRATION

GRMHD::GlobalBoundaryConditions& GRMHD::GlobalBoundaryConditions::getInstance() {
	static GRMHD::GlobalBoundaryConditions* me = nullptr;
	if(!me) me = new GRMHD::GlobalBoundaryConditions();
	return *me;
}

/// END OF REGISTRATION

// face indices: 0 x=xmin 1 x=xmax, 2 y=ymin 3 y=ymax 4 z=zmin 5 z=zmax
// corresponding 0-left, 1-right, 2-front, 3-back, 4-bottom, 5-top
static constexpr int EXAHYPE_FACE_LEFT = 0;     // Dimension 1
static constexpr int EXAHYPE_FACE_RIGHT = 1;    // Dimension 1
static constexpr int EXAHYPE_FACE_FRONT = 2;    // Dimension 2
static constexpr int EXAHYPE_FACE_BACK = 3;     // Dimension 2
static constexpr int EXAHYPE_FACE_BOTTOM = 4;   // Dimension 3
static constexpr int EXAHYPE_FACE_TOP = 5;      // Dimension 3

// solver-unspecific number of Variables
#include "AbstractGRMHDSolver_FV.h"
typedef GRMHD::AbstractGRMHDSolver_FV  Solver_FV;
static constexpr int nVar = Solver_FV::NumberOfVariables;
static constexpr int nDim = DIMENSIONS;

// ADERDG specific, for time integration:
#include "AbstractGRMHDSolver_ADERDG.h"
typedef GRMHD::AbstractGRMHDSolver_ADERDG  Solver_ADERDG;
static constexpr int basisSize = Solver_ADERDG::Order + 1;
static constexpr int order = Solver_ADERDG::Order;

GRMHD::BoundaryConditions::BoundaryConditions() :
	_log("BoundaryCondition"),
	aderdg_solver(nullptr),
	left(nullptr), right(nullptr), front(nullptr), back(nullptr), bottom(nullptr), top(nullptr) {}
	

GRMHD::BoundaryConditions::BoundaryConditions(Solver_ADERDG* _aderdg_solver) : 
	GRMHD::BoundaryConditions() {
	aderdg_solver = _aderdg_solver;
	}

bool GRMHD::BoundaryConditions::allFacesDefined() {
	return (
		left != nullptr && right != nullptr &&
		front != nullptr && back != nullptr &&
		(DIMENSIONS==3?(bottom != nullptr && top != nullptr):true)
	);
}

/// Apply by looking up the saved boundarymethods
void GRMHD::BoundaryConditions::apply(BOUNDARY_SIGNATURE) {
	switch(faceIndex) {
		// Dim 1
		case EXAHYPE_FACE_LEFT:   (this->*left  )(BOUNDARY_CALL); break;
		case EXAHYPE_FACE_RIGHT:  (this->*right )(BOUNDARY_CALL); break;
		
		// Dim 2
		case EXAHYPE_FACE_FRONT:  (this->*front )(BOUNDARY_CALL); break;
		case EXAHYPE_FACE_BACK:   (this->*back  )(BOUNDARY_CALL); break;
		
		#if DIMENSIONS==3
		case EXAHYPE_FACE_TOP:    (this->*top   )(BOUNDARY_CALL); break;
		case EXAHYPE_FACE_BOTTOM: (this->*bottom)(BOUNDARY_CALL); break;
		#endif
		default: throw std::runtime_error("Inconsistent face index");
	}
}

#define SET_BC(name) \
	if(constants.contains(#name)) { \
		std::string val = constants.get(#name).get_string();\
		name = parseFromString(val);\
		if(name==nullptr) {\
			logError("setFromSpecFile", "Boundary condition at " #name ": Invalid method: '" << val << "'");\
		} else {\
			logInfo("setFromSpecFile", #name " boundary: " << val);\
		}\
	} else {\
		logError("setFromSpecFile", "Boundary condition " #name " is not defined.");\
	}

bool GRMHD::BoundaryConditions::setFromParameters(const mexa::mexafile& constants, bool raiseOnFailure) {
	static tarch::logging::Log _log("BoundaryConditions");
	SET_BC(left);
	SET_BC(right);
	SET_BC(front);
	SET_BC(back);
	if(DIMENSIONS==3) {
		SET_BC(top);
		SET_BC(bottom);
	}
	if(!allFacesDefined() && raiseOnFailure) {
		logError("setFromparameters()", "Not all boundary faces have been defined. Cannot continue. I got the constants definition " << constants.toString());
		if(left==nullptr) logError("setFromparameters()", "left boundary is not defined");
		if(right==nullptr) logError("setFromparameters()", "right boundary is not defined");
		if(front==nullptr) logError("setFromparameters()", "front boundary is not defined");
		if(back==nullptr) logError("setFromparameters()", "back boundary is not defined");
		if(DIMENSIONS==3) {
			if(top==nullptr) logError("setFromparameters()", "top boundary is not defined");
			if(bottom==nullptr) logError("setFromparameters()", "bottom boundary is not defined");
		}
		std::abort();
	}
	return allFacesDefined();
}

GRMHD::BoundaryConditions::boundarymethod
  GRMHD::BoundaryConditions::parseFromString(const std::string& value) {
	GRMHD::BoundaryConditions::boundarymethod target=nullptr;
	if(value == "zero")
		target = &me::vacuum;
	if(value == "reflective" || value == "refl" || value=="wall")
		target = &me::reflective;
	if(value == "copy" || value == "outflow")
		target = &me::outflow;
	if(value == "exact")
		target = &me::exact;
	return target;
	//return !(target==nullptr);
	//throw std::exception("Invalid boundary method"); // Todo: return a usable error, or std::optional<BoundaryMethod> or so.
}

/**
 * Compute consistent fluxes.
 **/
void GRMHD::BoundaryConditions::deriveAderdgFlux(BOUNDARY_SIGNATURE) {
	if(fluxesRequested(BOUNDARY_CALL)) {
		// compute all fluxes but extract the one needed in direction "d".
		// Construct dummy storage Fs:
		double Fs[nDim][nVar], *F[nDim];
		for(int e=0; e<nDim; e++) F[e] = Fs[e];
		// Only point the direction "d" to the real outgoing storage.
		F[d] = fluxOut;
		// Call the fluxes, afterwards throwing away the fluxes in other
		// direction.
		aderdg_solver->flux(stateOut, F); 
	}
}

void GRMHD::BoundaryConditions::vacuum(BOUNDARY_SIGNATURE) {
	VacuumInitialData foo(stateOut);
	deriveAderdgFlux(BOUNDARY_CALL);
}

/// Outflow / Copy BC
void GRMHD::BoundaryConditions::outflow(BOUNDARY_SIGNATURE) {
	NVARS(i) stateOut[i] = stateIn[i];
	deriveAderdgFlux(BOUNDARY_CALL);
}

void GRMHD::BoundaryConditions::reflective(BOUNDARY_SIGNATURE) {
	// First, copy state
	NVARS(i) stateOut[i] = stateIn[i];
	
	SVEC::GRMHD::Shadow Q(stateOut);
	
	// vectors: flip sign
	Q.Si.lo(d) *= -1;
	Q.Bmag.up(d) *= -1;
	
	// Note: ADM BC should never be used anywhere as
	// ADM variables are parameters in this setup anyway!
	
	Q.beta.up(d) *= -1;
	// Symmetric Tensors: flip sign of off-diagonal
	DFOR(i) if(i != d) Q.gam.lo(i,d) *= -1;
	
	deriveAderdgFlux(BOUNDARY_CALL);
}

/// Exact BC for ADERDG
void GRMHD::BoundaryConditions::exact(BOUNDARY_SIGNATURE) {
	if(workingForADERDG(BOUNDARY_CALL)) {
		double Qgp[nVar], Fs[nDim][nVar], *F[nDim];
		for(int di=0; di<nDim; di++)
			F[di] = Fs[di];
		// zeroise stateOut, fluxOut
		for(int m=0; m<nVar; m++) {
			stateOut[m] = 0;
			if(fluxesRequested(BOUNDARY_CALL))
				fluxOut[m] = 0;
		}
		// Time integrate exact BC
		for(int i=0; i < basisSize; i++)  { // i == time
			const double weight = kernels::legendre::weights[order][i];
			const double xi = kernels::legendre::nodes[order][i];
			double ti = t + xi * dt;

			// do *not* use adjustPointSolution, it only works if ti~0.
			//solver->adjustPointSolution(x, ti, dt, Qgp);
			InitialData(x, ti, Qgp);
			aderdg_solver->flux(Qgp, F);
			
			for(int m=0; m < nVar; m++) {
				stateOut[m] += weight * Qgp[m];
				fluxOut[m] += weight * Fs[d][m];
			}
		}
		
		// no need to call deriveAderdgFlux here because fluxes have been
		// set exactly by time integration.
	} else {
		// FV: No time integration
		InitialData(x, t, stateOut);
	}
}
