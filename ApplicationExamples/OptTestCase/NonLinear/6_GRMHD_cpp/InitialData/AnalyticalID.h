#ifndef ANALYTICAL_ID_GRMHD_OCT17
#define ANALYTICAL_ID_GRMHD_OCT17

#include "InitialData.h"

// Adapter to functions or classes which create ID on the fly ("analytical"
template<typename IDCreator>
struct AnalyticalID : public InitialDataCode {
	void Interpolate(const double* x, double t, double* Q) {
		IDCreator(x, t, Q);
	}
};

// Initial State: Masking the State Vector
struct InitialState : public SVEC::GRMHD::Shadow, public SVEC::Hydro::Primitives::Stored {
	InitialState(double* Q) : SVEC::GRMHD::Shadow(Q), SVEC::Hydro::Primitives::Stored() {
		// by default, set NaNs for uninitialized values to catch mistakes when setting ID
		Dens = tau = phi = rho = press = alpha = NAN;
		DFOR(i) { Bmag.up(i) = Si.lo(i) = vel.up(i) = beta.up(i) = NAN; }
		DFOR(i) DFOR(j) gam.lo(i,j) = NAN;
	}
};

/// (Analytical) vacuum initial data
struct VacuumInitialData : public InitialState {
	static constexpr double rho0 = 1.;
	static constexpr double p0   = 1.;
	
	void flatSpace() {
		// ADM base: Flat space
		alpha = 1.0;
		DFOR(i) beta.up(i) = 0;
		SYMFOR(i,j) gam.lo(i,j) = SVEC::sym::delta(i,j);
	}
	
	VacuumInitialData(double* Q) : InitialState(Q) {
		// zero first all variables. This is useful for a 2D workaround as 3D variables
		// are stored but not initialized (only 2D storage is actually used)
		// NVARS(i) Q[i] = 0;
		
		// For test, do not set anything
		//return;
		
		// c2p invariants:
		DFOR(i) Bmag.up(i) = 0;
		phi = 0; // damping term
		
		// primitives:
		rho = rho0;
		press = p0;
		DFOR(i) vel.up(i) = 0;
		
		// ADM
		flatSpace();
		
		// Don't forget to call the Cons2Prim afterwards.
	}
};

/// AlfenWave analytical initial data
struct AlfenWave : public VacuumInitialData {
	AlfenWave(const double* const x, const double t, double* Q);
};

/// Conserved AlfenWave Initial Data
struct AlfenWaveCons : public AlfenWave {
	double V[100];
	AlfenWaveCons(const double* const x, const double t, double* Q) : AlfenWave(x,t,V) {
			SVEC::GRMHD::Prim2Cons(Q, V).copyFullStateVector();
		}
};

struct GRMHD_Shocktube : public InitialDataCode {
	tarch::logging::Log _log;
	typedef SVEC::GRMHD::Primitives::Stored state;
	struct stateFromParameters : public state {
		void readParameters(const mexa::mexafile& parameters);
		std::string toString() const;
		
	};

	stateFromParameters prims_right, prims_left;
	double sep_x;
	
	GRMHD_Shocktube(); // start with reasonable defaults
	void readParameters(const mexa::mexafile& parameters) override;
	void Interpolate(const double* x, double t, double* Q);
};

#endif /* ANALYTICAL_ID_GRMHD_OCT17 */
