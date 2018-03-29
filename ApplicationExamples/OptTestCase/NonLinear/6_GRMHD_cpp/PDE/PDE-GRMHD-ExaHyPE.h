#ifndef __GRMHD_FOR_EXAHYPE_ADAPTER__
#define __GRMHD_FOR_EXAHYPE_ADAPTER__

/**
 * Coupling of the SVEC GRMHD PDE to ExaHyPE.
 * You can actually use it like
 *
 *  using namespace SVEC::GRMHD::ExaHyPEAdapter;
 *
 * in your ExaHyPE solvers and then just pass all PDE calls.
 **/

#include "PDE/PDE.h"

#if !defined(TEST_NEW_PDE_AUTONOMOUSLY)
#include "peano/utils/Dimensions.h" // Defines DIMENSIONS
#include "AbstractGRMHDSolver_ADERDG.h" // Defines: nVar
#endif

#include "peano/utils/Dimensions.h" // Defines DIMENSIONS


namespace SVEC {
namespace GRMHD {
namespace ExaHyPEAdapter {

	#ifdef TEST_NEW_PDE_AUTONOMOUSLY
		// for an autonomous test distinct from ExaHyPE
		constexpr int nVar = 19;
	#else
		constexpr int nVar = /*global*/::GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;  // ensure this is 19 or so
	#endif

	#if defined(NVARS)
	#error NVARS was already defined. Maybe you try to load several different PDEs at a time?
	#endif
	
	#if !defined(DIMENSIONS)
	#error DIMENSIONS needs to be defined. Maybe include tarch/la/Dimensions.h
	#endif

	// a very handy macro.
	#define NVARS(i) for(int i=0; i < SVEC::GRMHD::ExaHyPEAdapter::nVar; i++)


	/**
	 * A generic storage class for the high level shadows of the fluxes and gradients.
	 * Typically,
	 *  - Base is an actual shadow class which uses this storage
	 *  - T is the high level shadow type (i.e. GRMHD::Shadow)
	 *  - I is a double pointer, either "const double*const" or "double*const".
	 * 
	 * As this is for the interface between ExaHyPE and Tensish, this class does all the
	 * dimensional fits.
	 **/
	template<class Base, class T, typename I>
	struct Storage : Base {
		T data[TDIM];
		
		/**
		* In a linearized matrix (as today gradQ), give the position of the j. vector inside
		* the matrix.
		* Attention: Here we implement a HACK to cover the situation that tensish works in
		*   3 dimensions but ExaHyPE only works in 2 dimensions (which is typical for GRMHD).
		*   We then give back a dummy storage which is zeroized at each invocation
		**/
		static I idxlin(I p, int j) {
			if((TDIM == 3 && DIMENSIONS == 2) && j==2) {
				static double dummy[nVar];
				std::fill_n(dummy,nVar,0);
				return dummy;
			}
			return p+j*nVar;
		}

		
		#if TDIM == 3
			/// individual vectors
			Storage(I Ax, I Ay, I Az) : Base(data), data{Ax,Ay,Az} {}
			/// similar form as the fluxes (matrix)
			Storage(I*AA) : Storage(AA[0],AA[1],AA[2]) {}
			/// linear array
			Storage(I AA) : Storage(idxlin(AA,0),idxlin(AA,1),idxlin(AA,2)) {}
		#elif TDIM == 2
			/// individual vectors
			Storage(I Ax, I Ay) : Base(data), data{Ax,Ay} {}
			/// Similar form as the fluxes
			Storage(I*AA) : Storage(AA[0],AA[1]) {}
			/// Linear array
			Storage(I AA) : Storage(idxlin(AA,0),idxlin(AA,1)) {}
			/// For convenience a similar interface as what TDIM==3 provides:
			Storage(I Ax, I Ay, I Az) : Storage(Ax,Ay) {}
		#endif
	};


	/********* CONSERVATIVE PART ************/
	
	typedef Storage<GRMHD::PDE::Fluxes, GRMHD::PDE::Flux, double*const> BaseFluxes;
	struct Fluxes : public BaseFluxes {
		using BaseFluxes::BaseFluxes;

		void zeroMaterialFluxes() {
			for(int d=0; d<DIMENSIONS; d++) up(d).zero_adm();
		}
	};
	
	inline Fluxes flux(const double* const Q, double** F) {
		Fluxes f(F);
		GRMHD::PDE(Q).flux(f);
		return f;
	}

	// we also can work with splitted fluxes:
	inline Fluxes flux(const double* const Q, double *Fx, double *Fy, double *Fz) {
		Fluxes f(Fx,Fy,Fz);
		GRMHD::PDE(Q).flux(f);
		return f;
	}
	
	/********* NONCONSERVATIVE PART ************/
	
	typedef Storage<GRMHD::PDE::Gradients, const GRMHD::PDE::Gradient, const double* const> Gradients;

	inline GRMHD::Shadow nonConservativeProduct(const double* const Q, const double* const gradQ, double* BgradQ) {
		GRMHD::Shadow n(BgradQ);
		Gradients g(gradQ);
		GRMHD::PDE(Q).nonConservativeProduct(g, n);
		return n;
	}
	
	// for a future API
	inline GRMHD::Shadow nonConservativeProduct(const double* const Q, const double* const Qx, const double* const Qy, const double* const Qz, double* BgradQ) {
		GRMHD::Shadow n(BgradQ);
		Gradients g(Qx,Qy,Qz);
		GRMHD::PDE(Q).nonConservativeProduct(g, n);
		return n;
	}
	
	inline GRMHD::Shadow algebraicSource(const double* const Q, double* S) {
		// GRMHD::Densitied state(Q); // GRMHD::Shadow not sufficient, need densitied + ADM for gam.det
		// subject to further improvements: More efficient source without C2P.
		GRMHD::Shadow s(S);
		GRMHD::PDE(Q).algebraicSource(s);
		return s;
	}
	
	inline GRMHD::Shadow fusedSource(const double* const Q, const double* const gradQ, double* S) {
		GRMHD::Shadow s(S);
		Gradients g(gradQ);
		GRMHD::PDE(Q).fusedSource(g, s);
		return s;
	}

	// for a future API
	inline GRMHD::Shadow fusedSource(const double* const Q, const double* const Qx, const double* const Qy, const double* const Qz, double* S) {
		GRMHD::Shadow s(S);
		Gradients g(Qx,Qy,Qz);
		GRMHD::PDE(Q).fusedSource(g, s);
		return s;
	}
	
	inline GRMHD::Shadow eigenvalues(const double* const Q, const int dIndex, double* lambda) {
		GRMHD::Shadow L(lambda);
		GRMHD::PDE(Q).eigenvalues(L,dIndex);
		return L;
	}
	
	// Also provide the Prim2Cons and Cons2Prim for convenience.
	typedef SVEC::GRMHD::Prim2Cons Prim2Cons;
	typedef SVEC::GRMHD::Cons2Prim Cons2Prim;
} // ns ExaHyPEAdapter
} // ns GRMHD
} // ns SVEC

#endif /* __GRMHD_FOR_EXAHYPE_ADAPTER__ */
