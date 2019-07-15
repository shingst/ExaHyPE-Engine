/* Outsourcing from GRMHDSolver_*.cpp */
#ifndef GENERIC_BOUNDARY_CONDITIONS_WHICH_SHOULD_GO_SOMEWHERE_CENTRAL
#define GENERIC_BOUNDARY_CONDITIONS_WHICH_SHOULD_GO_SOMEWHERE_CENTRAL

#include "GRMHDSolver_ADERDG_Variables.h" // just for GRMHD-specific variable names

/** Helper class to parse key value strings */
class StringMapView {
   std::string base;
   public:
    StringMapView(std::string _base) : base(_base) {}
    std::string getValueAsString(const std::string& key) const {
	std::size_t startIndex = base.find(key);
	startIndex = base.find(":", startIndex);
	std::size_t endIndex = base.find_first_of("}, \n\r", startIndex + 1);
	return base.substr(startIndex + 1, endIndex - startIndex - 1);
    }
    bool isValueValidString(const std::string& key) const {
	std::size_t startIndex = base.find(key);
	return (startIndex != std::string::npos);
    }
  };

#define ADERDG_BOUNDARY_SIGNATURE \
  const double* const x,const double t,const double dt,const int faceIndex,const int d,\
  const double * const fluxIn,const double* const stateIn,\
  double *fluxOut, double* stateOut
#define ADERDG_BOUNDARY_CALL \
  x, t, dt, faceIndex, d, fluxIn, stateIn, fluxOut, stateOut

// The difference between the ADERDG and FV boundary signature is that in ADERDG, also
// the fluxIn and fluxOut pointers are given.
// Thus the ADERDG BC signature is a real superset of the FV BC signature and thus we
// stick to the ADERDG signature as the more generic one.

#define BOUNDARY_SIGNATURE   ADERDG_BOUNDARY_SIGNATURE
#define BOUNDARY_CALL        ADERDG_BOUNDARY_CALL


// Note that passing a number of variables like this is an ugly way to circumvent a proper
// object oriented structure which is passed instead.

template <typename SolverType>
class BoundaryConditions {
public:
	// face indices: 0 x=xmin 1 x=xmax, 2 y=ymin 3 y=ymax 4 z=zmin 5 z=zmax
	// corresponding 0-left, 1-right, 2-front, 3-back, 4-bottom, 5-top
	static constexpr int EXAHYPE_FACE_LEFT = 0;
	static constexpr int EXAHYPE_FACE_RIGHT = 1;
	static constexpr int EXAHYPE_FACE_FRONT = 2;
	static constexpr int EXAHYPE_FACE_BACK = 3;
	static constexpr int EXAHYPE_FACE_BOTTOM = 4;
	static constexpr int EXAHYPE_FACE_TOP = 5;
	
	SolverType* solver;
	
	// a callback to SolverType. Typically this is fulfilled by ADERDG::adjustPointSolution
	// or FV::adjustSolution (for inconsistency-reasons wrongly labeled).
	typedef void (SolverType::*adjustSolutionType)(double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,
		const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt);
	
	// constexpr SolverType::VariableShortcuts vars; // this does not work
	//const GRMHD::AbstractGRMHDSolver_ADERDG::VariableShortcuts vars;
	
	typedef BoundaryConditions<SolverType> me;
	typedef void (me::*boundarymethod)(BOUNDARY_SIGNATURE);
	
	// this could be improved by using boundarymethod side[6]
	// and having apply defined just as  side[faceIndex](...);.
	boundarymethod left, right, front, back, bottom, top;
	
	BoundaryConditions(SolverType* _solver) : solver(_solver),
		left(nullptr), right(nullptr), front(nullptr), back(nullptr), bottom(nullptr), top(nullptr) {}
	
	bool allFacesDefined() {
		return (left != nullptr && right != nullptr && front != nullptr && back != nullptr && bottom != nullptr && top != nullptr);
	}
	
	/// Apply by looking up the saved boundarymethods
	void apply(BOUNDARY_SIGNATURE) {
		switch(faceIndex) {
			case EXAHYPE_FACE_LEFT:   (this->*left  )(BOUNDARY_CALL); break;
			case EXAHYPE_FACE_FRONT:  (this->*front )(BOUNDARY_CALL); break;
			case EXAHYPE_FACE_BOTTOM: (this->*bottom)(BOUNDARY_CALL); break;
			case EXAHYPE_FACE_RIGHT:  (this->*right )(BOUNDARY_CALL); break;
			case EXAHYPE_FACE_BACK:   (this->*back  )(BOUNDARY_CALL); break;
			case EXAHYPE_FACE_TOP:    (this->*top   )(BOUNDARY_CALL); break;
			default: throw std::runtime_error("Inconsistent face index");
		}
	}

	
	// just to highlight that apply() is the main method:
	// void operator()(ADERDG_BOUNDARY_SIGNATURE) { apply(ADERDG_BOUNDARY_CALL); }
	
	/***************************************************************************
	 * Boundary Methods.
	 * 
	 * To add your own, do the following:
	 * 
	 *  1) Implement it with a method void yourName(ADERDG_BOUNDARY_SIGNATURE)
	 *  2) Add yourName to the parseFromString thing
	 * 
	 ***************************************************************************/
	
	#define BC_FROM_SPECFILE_SIDE(name) \
		if(constants.isValueValidString(#name)) { \
			std::string val = constants.getValueAsString(#name);\
			name = parseFromString(val);\
			if(name==nullptr) {\
				logError("setFromSpecFile", "Boundary condition at " #name ": Invalid method: '" << val << "'");\
			} else {\
				logInfo("setFromSpecFile", #name " boundary: " << val);\
			}\
		} else {\
			logError("setFromSpecFile", "Boundary condition at " #name " is not valid according to specfile parser.");\
		}

	/**
	 * MapView shall have methods
	 *    bool isValueValidString(std::string)
	 *    std:.string getvalueAsString(std::string)
	 * such as the Parser::ParserView class has.
	 * 
	 * @returns true in case of success, false otherwise
	 **/
	template<typename MapView>
	bool setFromSpecFile(MapView constants) {
		tarch::logging::Log _log("GRMHD::BoundaryCondition");
		BC_FROM_SPECFILE_SIDE(left);
		BC_FROM_SPECFILE_SIDE(right);
		BC_FROM_SPECFILE_SIDE(front);
		BC_FROM_SPECFILE_SIDE(back);
		BC_FROM_SPECFILE_SIDE(bottom);
		BC_FROM_SPECFILE_SIDE(top);
		return allFacesDefined();
	}

	/**
	 * assign a boundary method by some string value
	 * @returns pinter if value could be correctly parsed, in case of error nullptr.
	 **/
	static boundarymethod parseFromString(const std::string& value) {
		boundarymethod target=nullptr;
		if(value == "zero")
			target = &me::vacuum;
		if(value == "reflective" || value == "refl")
			target = &me::reflective;
		if(value == "copy" || value == "outflow")
			target = &me::outflow;
		if(value == "exact")
			target = &me::exact;
		return target;
		//return !(target==nullptr);
		//throw std::exception("Invalid boundary method"); // Todo: return a usable error, or std::optional<BoundaryMethod> or so.
	}

	static constexpr int nVar = SolverType::NumberOfVariables;
	static constexpr int nDim = DIMENSIONS;
	// ADERDG specific:
	static constexpr int basisSize = SolverType::Order + 1;
	static constexpr int order = SolverType::Order;
	
	bool fluxesRequested(BOUNDARY_SIGNATURE) {
		return (fluxOut != nullptr); // could also check for fluxIn != nullptr.
	}
	
	/**
	 * Compute consistent fluxes.
	 **/
	void deriveAderdgFlux(BOUNDARY_SIGNATURE) {
		if(fluxesRequested(BOUNDARY_CALL)) {
			// compute all fluxes but extract the one needed in direction "d".
			// Construct dummy storage Fs:
			double Fs[nDim][nVar], *F[nDim];
			for(int e=0; e<nDim; e++) F[e] = Fs[e];
			// Only point the direction "d" to the real outgoing storage.
			F[d] = fluxOut;
			// Call the fluxes, afterwards throwing away the fluxes in other
			// direction.
			solver->flux(stateOut, F);
		}
	}

	/**
	 * The vacuum is problem dependent, thus this is really MHD vacuum in
	 * Minkowski spacetime for GRMHD.
	 **/
	void vacuum(BOUNDARY_SIGNATURE) {
		using namespace GRMHD::GRMHDSolver_ADERDG_Variables::shortcuts;
		stateOut[rho] = 1e-10; // Some vacuum value
		stateOut[vel+0] = 0;
		stateOut[vel+1] = 0;
		stateOut[vel+2] = 0;
		stateOut[E] = 1.0; // or some reasonably vacuum value
		stateOut[psi] = 0;

		// Minkowski spacetime:
		stateOut[lapse] = 1;
		stateOut[shift+0] = 0;
		stateOut[shift+1] = 0;
		stateOut[shift+2] = 0;
		stateOut[gij+0] = 1; // gxx
		stateOut[gij+1] = 0; // gxy
		stateOut[gij+2] = 0; // gxz
		stateOut[gij+3] = 0; // gyx
		stateOut[gij+4] = 1; // gyy
		stateOut[gij+5] = 0; // gyz
		stateOut[gij+6] = 1; // gzz
		
		deriveAderdgFlux(BOUNDARY_CALL);
	}
	
	/// Outflow / Copy BC
	void outflow(BOUNDARY_SIGNATURE) {
		for(int m=0; m<nVar; m++) {
			stateOut[m] = stateIn[m];
		}

		deriveAderdgFlux(BOUNDARY_CALL);
	}
	
	/// Reflection/Hard wall BC, Currently WRONG.
	void reflective(BOUNDARY_SIGNATURE) {
		// First, copy state
		for(int m=0; m<nVar; m++) {
			stateOut[m] = stateIn[m];
		}
		
		using namespace GRMHD::GRMHDSolver_ADERDG_Variables::shortcuts;
		
		// second, flip the sign of the vectors in normal direction
		stateOut[vel+d] *= -1;
		stateOut[shift+d] *= -1;
		
		// Third, flip the sign of the tensors in normal direction
		// TODO: Look up how this is done correctly!!!!
		
		deriveAderdgFlux(BOUNDARY_CALL);
	}
		
	/// Exact BC for ADERDG
	void exact(BOUNDARY_SIGNATURE) {
		double Qgp[nVar], Fs[nDim][nVar], *F[nDim];
		for(int d=0; d<nDim; d++)
			F[d] = Fs[d];
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

			solver->adjustPointSolution(x, ti, dt, Qgp);
			//id->Interpolate(x, ti, Qgp);
			solver->flux(Qgp, F);
			
			for(int m=0; m < nVar; m++) {
				stateOut[m] += weight * Qgp[m];
				fluxOut[m] += weight * Fs[d][m];
			}
		}
		
		// no need to call deriveAderdgFlux here because fluxes have been
		// set exactly by time integration.
	}
};

#endif /* GENERIC_BOUNDARY_CONDITIONS_WHICH_SHOULD_GO_SOMEWHERE_CENTRAL */
