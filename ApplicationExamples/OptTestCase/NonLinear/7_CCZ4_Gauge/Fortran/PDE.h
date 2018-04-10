#ifndef __EXAHYPE_CCZ4_PDE_FORTRAN__
#define __EXAHYPE_CCZ4_PDE_FORTRAN__

// forward declaration instead of #include "Parameters/mexa.h"
namespace mexa { class mexafile; }

namespace CCZ4Fortran {

	/**
	* The CCZ4 PDE system runtime parameters.
	* This structure has to been kept in sync with the Fortran user type in typesDef.f90
	**/
	struct tEquations {
		double k1;        // CCZ4 damping parameter k1 
		double k2;        // CCZ4 damping parameter k2  
		double k3;        // CCZ4 damping parameter k3 
		double eta;       // CCZ4 damping parameter for the PDE of b^i in the gamma driver 
		double itau;      // inverse relaxation time for the relaxed approach in order to enforce the unit determinant of the conformal metric and the trace-free matrix A_ij
		double f;         // set f=0.75 or f=1.0 for the gamma driver. Typical BSSNOK value: f=0.75. Set f=0 to avoid shift evolution, i.e. to have d/dt beta^i=0.  
		double g;         // not used at the moment: reserved for the generalized harmonic shift   
		double xi;        // set to zero to switch off the gamma driver and to avoid shift evolution (i.e. to have d/dt b^i = 0) 
		double e;         // cleaning speed e>=1 for the Hamiltonian constraint. Typical value for CCZ4 is e=1. However, e>1 gives better constraints and better hyperbolicity//  
		double c;         // set c=0 to remove the algebraic source terms of the type -2*Theta 
		double mu;        // mu=1 adds second order ordering constraints. Has shown to be important for the gamma driver, but too large values can also become a problem... 
		double ds;        // set this value always to ds=1, unless you know what you are doing. It allows to increase the cleaning speed for the momentum constraints, but in CCZ4 this does not seem to do what it should do...  
		double sk;        // setting sk=0 removes the contribution of the shift also in the auxiliary variables. If you want to evolve the shift, set sk=1. 
		double bs;        // set bs=1 if you want to activate the shift convection for beta, b and B (standard CCZ4 formulation). set it to bs=0 to switch off shift convection for those quantities 
		int LapseType;    // LapseType = 0 is harmonic lapse, LapseType = 1 is 1+log slicing.
		
		void setByParameters(const mexa::mexafile& parameters);
		std::string toString() const;
	};

	struct GlobalPDEParameters {
		tEquations* const parameters;
		GlobalPDEParameters();

		/// Tries to parameters from mexa. In case of problems, raises exceptions.
		void setByParameters(const mexa::mexafile& parameters);
	
		/// Gives a singleton instance
		static GlobalPDEParameters& getInstance();
	};

	// it is pointless to have this being a class as the methods only mask Fortran functions.
	namespace PDE {
		void nonConservativeProduct(const double* const Q, const double* const gradQ, double* BgradQ);
		void algebraicSource(const double* const Q, double* S);
		void fusedSource(const double* const Q, const double* const gradQ, double* S);
		void eigenvalues(const double* const Q, const int dIndex, double* lambda);
	}

	void Cons2Prim(double* const V, const double* const Q);
	void Prim2Cons(double* const Q, const double* const V);

	void AdjustPointSolution(const double* const x,const double t,const double dt,double* Q);
	void DeriveADMConstraints(double* constraints, const double* const Q, const double* const gradQ);
	void EnforceCCZ4Constraints(double *Q);

} // ns CCZ4Fortran

#endif /* __EXAHYPE_CCZ4_PDE_FORTRAN__ */
