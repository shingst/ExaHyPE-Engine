/**
 * A C++ variant of the GRMHD equations for ExaHyPE.
 * 
 * This formulation uses the Tensish library for the Tensor abstraction.
 * We mimic the ExaHyPE notation of state vector abstractions with structures.
 * We copy all neccessary parts for practical purposes (scalar, vector, tensor
 * notation) as well as distinguation between conserved, primitive, material
 * parameters.
 *
 * Written at 2017-08-06 by SvenK.
 **/
#ifndef GRMHD_PDE_CPP
#define GRMHD_PDE_CPP

/**
 * TDIM defines the dimensionality of the PDE (as well as the underlying tensish
 * library). It can be independent of the ExaHyPE DIMENSIONS. For the GRMHD,
 * with certain initial data, only 3 dimensions are senseful.
 **/
#define TDIM 3

/* This header file needs the definition of a constexpr int GRMHD::nVar
 * as well as a preprocessor variable TDIM.
 * You should provide this in the including file, for instance:
 *
 *   namespace SVEC { constexpr int nVar = 10; }
 *   #define TDIM 3
 *   #include "PDE.h"
 *
 * TDIM has to be a C preprocessor variable for C macro ifdefs.
 * As ExaHyPE, there are some places where this code only works for 2 and 3
 * dimensions.
 */

#include "tensish.cpph"
#include <cstring> // memcpy
// see also at the end for further includes!


// TODO: Collect documentation at some central place.

/**
* Explanation of the different state vectors at the example
* of Hydro::Conserved::StateVector's ancestors and typedefs:
* 
* # ConstShadow:
* 
* Access the Conserved Variables Q fully read-only without access to anything
* beyond what's inside Q (no further storage allocated except pointers to Q).
* That is, you cannot access S^i or B_i as we don't compute it. Therefore, we
* say this structure is just "shadowing" Q. You might want to use it like
* 
*    const ConstShadow Q(Q_);
*    double foo = Q.Bmag.up(1) * Q.Dens;
* 
* # Shadow:
* The variable version of the ConstShadow: Allows to change all
* entries of Q by their aliases, as well as directly changing the vectors, i.e.
*
*    ConservedVariableShadow Q(Q_);
*    Q.Bmag.up(2) = Q.Dens + 15;
* 
* But not:
* 
*    Q.Bmag.lo(2) = 15; 
* 
* That is, we dont' allow setting (or even getting) B_i or S^j. This can be
* helpful to avoid errors in the formulation.
*
* # ConstShadowExtendable:
* Access to the full set of conserved variables in a read only fashion which
* however allows you to retrieve B_i and S^j, i.e. we have storage for these
* variables. These two vectors are the only non-constant members of the
* data structure as you need to fill them after construction.
* 
* # ShadowExtendable
* The variable version of ConservedConstFull: You can change all conserved
* variables as well access the B_i and S^j. Note that
* 
*   Q.Bmag.lo(2) = 15;
* 
* will not change the shadowed double* Q_.
* 
* # Stored
* A version of shadow with actual storage, i.e. not shadowing. Can be used without
* a state vector reference.
* 
**/

namespace SVEC {
	using namespace tensish; 
	using sym::delta; // Kronecker symbol
	
	// this could be the actual interface to a state vector, quite similar
	// to what vec::shadow<N> does, just with another name.
	// Actually this object could also take over ownership (std::move semantics)
	// or exploit caching better.
	/*
	template<int N, typename scalar> struct Shadow {
		scalar* const Q;
		constexpr Shadow(scalar* const Q_) : Q(Q_) {}
		double operator=(double val) { std::fill_n(data,N,val); return val; }
	};
	
	StateVector<Shadow>(QShadow)
	*/
	
	// geht:
	//template<int N>
	//using WhateverThisWorksAtleast = vec::shadow<N>;
	
	/**
	 * Hydrodynamics (fluid) state vectors in their conservative and primitive
	 * variants. The state vectors always store size = 5 elements and one could
	 * do classical, special relativistic and general relativistic hydrodynamics
	 * with these elements.
	 **/
	namespace Hydro {
		// Hydro vars are always with offset = 0.
		constexpr int size = 5;
		
		/**
		* Storage for the conserved Hydro variables, i.e. the tuple
		*   Q = (Dens, S_i, tau).
		* The structure allows addressing the vector Q with names.
		* The template allows to specify types for the scalars and 3-vectors inside Q.
		* See the typedefs below for usable instances of the "Conserved" structure.
		**/
		namespace Conserved {
			/**
			* These structs store the index positions of fields, i.e. Q.Dens=0.
			* They are similiar to the VariableShortcuts which are generated
			* by the ExaHyPE java toolkit. However, by defining our own we are
			* independent and can make use of the functional notation.
			**/
			namespace Indices {
				// Conserved hydro variables
				constexpr int Dens = 0;
				constexpr int Si_lo = 1;
				constexpr int tau = 4;
			}
			
			/// The Conserved hydro variables state vector
			template<typename state_vector, typename vector_lo, typename scalar>
			struct StateVector {
				state_vector Q;
				scalar &Dens, &tau;
				vector_lo Si; /* sic */
				
				constexpr StateVector(state_vector Q_) :
					Q(Q_),
					Dens (Q[Indices::Dens]), // Conserved scalars
					tau  (Q[Indices::tau]),
					Si   (Q+Indices::Si_lo) // Conserved vectors
					{}
				
				/// Multiply a factor onto the conserved variables
				void multiply_conserved(double factor) {
					Dens *= factor;
					tau *= factor;
					DFOR(i) Si.lo(i) *= factor;
				}
			};
			
			typedef StateVector<const double* const, const ConstLo<vec::const_shadow_D>, const double> ConstShadow;
			typedef StateVector<double* const, Lo<vec::shadow_D>, double> Shadow;				
			typedef StateVector<const double* const, UpLo<vec::stored_D, vec::const_shadow_D>::ConstLo, const double> ConstShadowExtendable;
			typedef StateVector<double* const, UpLo<vec::stored_D, vec::shadow_D>::InitLo, double> ShadowExtendable;
			struct Stored : public Shadow { double QStorage[size]; Stored() : Shadow(QStorage) {} };
		} // ns ConservedHydro
		
		/**
		* Prmitive hydro Variables, i.e. the tuple
		*   V = (rho, vel^i, press).
		* They only store the hydro variables, not Maxwell or ADM variables.
		**/
		namespace Primitives {
			namespace Indices {
				// The positions of the primitive hydro variables where there exists
				// a 1:1 mapping to the conserved hydro variables.
				constexpr int rho = 0;
				constexpr int vec_up = 1; /* sic */
				constexpr int press = 4;
			}
			
			/// The Primitive hydro variables state vector
			template<typename state_vector, /*typename vector_lo, */typename vector_up, typename scalar>
			struct StateVector {
				state_vector V;
				
				// Primitive Scalars
				scalar &rho, &press;
				
				// Primitive vectors
				vector_up vel;
				
				constexpr StateVector(state_vector V_) :
					V(V_),
					// Primitive Scalars
					rho  (V[Indices::rho]),
					press(V[Indices::press]),
					// Primitive vectors
					vel  (V+Indices::vec_up)
					{}
			};
		
			typedef StateVector<double*, Up<vec::shadow_D>, double> Shadow;
			typedef StateVector<const double* const, UpLo<vec::const_shadow_D,vec::stored_D>::ConstUp, const double> ConstShadowExtendable;
			typedef StateVector<double*, UpLo<vec::shadow_D, vec::stored_D>::InitUp, double> ShadowExtendable;
			struct Stored : public Shadow { double VStorage[size]; Stored() : Shadow(VStorage) {} };
		} // ns Primitives
	} // ns Hydro

	/**
	 * Magnetic field for a MHD formulation with constraint damping term Phi.
	 **/
	namespace Magneto {
		constexpr int size = 4; // length of Magneto/Maxwell-vector + constraint damping
		
		namespace RelativeIndices {
			constexpr int Bmag_up = 0;
			constexpr int Phi = 4;
		}
		namespace AbsoluteIndices {
			// Conserved (C2P invariant) Magnetic/Maxwell variables
			constexpr int offset = 5;
			constexpr int Bmag_up = 5;    /* length: 3 */
			constexpr int Phi = 8;
		}
		
		/// The Magneto state vector (magnetic part of MHD)
		template<typename state_vector, typename vector_up, typename scalar>
		struct StateVector {
			state_vector Q_Magneto;
			scalar &phi;
			vector_up Bmag; /* sic */
			
			constexpr StateVector(state_vector Q) :
				Q_Magneto(Q),
				phi  (Q[AbsoluteIndices::Phi]),   // Conserved scalars
				Bmag (Q+AbsoluteIndices::Bmag_up) // Conserved vectors
				{}
			
			/// Copy the magnetic state variables to another state vector
			void copy_magneto(double* target) {
				std::memcpy(target + AbsoluteIndices::offset, Q_Magneto + AbsoluteIndices::offset, size*sizeof(double));
			}
			
			/// Multiply a factor onto the magneto variables
			void multiply_magneto(double factor) {
				phi *= factor;
				DFOR(i) Bmag.up(i) *= factor;
			}
		};
		
		typedef StateVector<const double* const, const ConstUp<vec::const_shadow_D>, const double> ConstShadow;
		typedef StateVector<double* const, Up<vec::shadow_D>, double> Shadow;
		typedef StateVector<const double* const, UpLo<vec::const_shadow_D, vec::stored_D>::ConstUp, const double> ConstShadowExtendable;
		typedef StateVector<double* const, UpLo<vec::shadow_D, vec::stored_D>::InitUp, double> ShadowExtendable;
		struct Stored : public Shadow { double QMStored[size]; Stored() : Shadow(QMStored) {} };
		
		/// A zero-mimicro magnetic state vector for effectively turning off magnetic
		/// contributions at compile time.
		struct Zero {
			scalar::zero phi;
			UpLo<vec::zero, vec::zero> Bmag;
			constexpr Zero(const double* const Q) {}
			void copy_magneto(double* target) { std::abort(); }
			void multiply_magneto(double* target) { std::abort(); }
		};
	} // ns Magneto

	/**
	 * The ADM variables (alpha,beta^i,gamma_ij).
	 * TODO: We lack an inclusion for the extrinsic curvature tensor.
	 **/
	namespace ADMBase {
		constexpr int size = 10;
		
		namespace RelativeIndices {
			constexpr int lapse = 0;
			constexpr int shift_lo = 1;  /* length: 3 */
			constexpr int gam_lo = 4;    /* length: 6 */
		}
		namespace AbsoluteIndices {
			// The positions of the material parameters in the GRMHD state vector. There is an offset
			// due to the size of the hydro and maxwell variables before.
			constexpr int offset = 9;
			constexpr int lapse = 9;
			constexpr int shift_lo = 10;  /* length: 3 */
			constexpr int gam_lo = 13;    /* length: 6 */
			//constexpr int detg = 19;    // not yet
		}
		
		/**
		 * In pure GRMHD, the ADM variables are "material" parameters, i.e. they are constant after
		 * initial data setting. In GRMHD, the ADM set is just (alpha,beta,gamma). With a dynamical
		 * spacetime, it is (alpha,beta,gamma,kextr).
		 **/
		template<typename state_vector, typename scalar, typename vector_up, typename metric_lo, typename extrinsic_lo>
		struct StateVector {
			state_vector Q_ADM;
			scalar &alpha; ///< Lapse
			// &detg; ///< Determinant of g_ij as material parameter.
			vector_up beta; // Shift vector: (Conserved) Material parameter vector
			metric_lo gam;  // 3-Metric: (Conserved) Material parameter tensor
			extrinsic_lo Kext; ///< Extrinsic curvature tensor
			
			constexpr StateVector(state_vector Q) :
				Q_ADM(Q),
				alpha(Q[AbsoluteIndices::lapse]),
				//detg (Q[AbsoluteIndices::detg]),
				beta (Q+AbsoluteIndices::shift_lo),
				gam  (Q+AbsoluteIndices::gam_lo)
				{}
				
			void copy_adm(double* target) {
				std::memcpy(target + AbsoluteIndices::offset, Q_ADM + AbsoluteIndices::offset, size*sizeof(double));
			}
			
			// set all ADM components to zero. This is unphysical (not identity) for real metric but
			// practical if gradients of the metric are stored.
			void zero_adm() {
				alpha = beta.up = gam.lo = 0;
			}
		};
		
		// A version of the ADMVariables where beta_lo and the upper metric
		// are not recovered.
		typedef StateVector<
			const double* const,
			const double,                       // lapse
			const ConstUp<vec::const_shadow_D>, // shift
			const ConstLo<sym::const_shadow_D>, // metric
			sym::illegal                        // extrinsic curvature
			> ConstShadow;
			
		// A writable shadowed ADMBase. Writeable is only useful for setting the initial data.
		// You can only set the lower components of the metric.
		typedef StateVector<
			double* const,
			double,
			Up<vec::shadow_D>,
			Lo<sym::shadow_D>,
			sym::illegal                        // extrinsic curvature
			> Shadow;
		
		// This is what the PDE needs: A working metric (upper and lower) but only upper shift.
		// This is read only (especially inside metric3).
		typedef StateVector<
			const double* const,
			const double,
			// ConstUp_Lo<vec::const_shadow, vec::stored_D>, // if you ever need beta_lo
			ConstUp<vec::const_shadow_D>, // if you never need beta_lo
			metric3, // could store only parts of the metric here, too.
			sym::illegal                        // extrinsic curvature
			> Full;
		
	} // ns ADM
	
	/// a generic way putting the state vectors of two theories next to each other
	namespace TwoTheories {
		template<class constructor, class A, class B>
		struct StateVector : public A, public B {
			StateVector(constructor X) : A(X), B(X) {}
		};
	} // ns TwoTheories
	
	/**
	 * Stores linearly or shadows the Hydro+Magneto variables in a common MHD.
	 * The linear storage can be advanteous at some points.
	 **/
	namespace MHD {
		constexpr int size = Hydro::size + Magneto::size;
		
		namespace Conserved {
			typedef TwoTheories::StateVector<double* const, Hydro::Conserved::Shadow, Magneto::Shadow> Shadow;
			typedef TwoTheories::StateVector<const double* const, Hydro::Conserved::ConstShadow, Magneto::ConstShadow> ConstShadow;
			typedef TwoTheories::StateVector<const double* const, Hydro::Conserved::ConstShadowExtendable, Magneto::ConstShadowExtendable> ConstShadowExtendable;
			struct Stored : public Shadow { double MHDStorage[size]; Stored() : Shadow(MHDStorage) {} };
		} // ns Conserved
		
		namespace Primitives {
			typedef TwoTheories::StateVector<double* const, Hydro::Primitives::Shadow, Magneto::Shadow> Shadow;
			typedef TwoTheories::StateVector<const double* const, Hydro::Primitives::ConstShadowExtendable, Magneto::ConstShadowExtendable> ConstShadowExtendable;
			struct Stored : public Shadow { double MHDStorage[size]; Stored() : Shadow(MHDStorage) {} };
		} // ns Primitives
		
		
		// TODO: Do not use these abbreviations
		
		typedef Conserved::Shadow Shadow;
		typedef Conserved::ConstShadow ConstShadow;
		typedef Conserved::ConstShadowExtendable ConstShadowExtendable;
		typedef Conserved::Stored Stored;
	} // ns MHD
	
	// Generic PDE:
	/*
	template<class System>
	struct GenericPDE : public System {
		// the System defines the types of the constituents:
		typedef typename System::State State;
		typedef typename System::Fluxes Fluxes;
		typedef typename System::Gradients Gradients;
		
		GenericPDE(const double* const Q) : System(Q) {}
	
		/// Conserved fluxes
		void flux(Fluxes& flux);
	
		/// Not to be used today.
		void eigenvalues(State& lambda, const int d);
	
		/// The BgradQ
		void nonConservativeProduct(const Gradients& grad, State& ncp);
	
		/// Sets the algebraic source
		void algebraicSource(State& Source_data);
		
		/// Adds the algebraic source (assumes Source=0 or similiar)
		void addAlgebraicSource(State& Source_data);
		
		/// Computes the fusedSource = algebraicSource - NCP. This is the classical RightHandSide.
		void fusedSource(const Gradients& grad, State& source);
	};
	*/
	

	// Naming suggestion: "GRMHD" for this namespace and "SVEC" for the parent.
	namespace GRMHD { // need a better name
		constexpr int size = MHD::size + ADMBase::size;

		namespace Conserved {
			typedef TwoTheories::StateVector<double* const, MHD::Conserved::Shadow, ADMBase::Shadow> Shadow;
			typedef TwoTheories::StateVector<const double* const, MHD::Conserved::ConstShadow, ADMBase::ConstShadow> ConstShadow;
			struct Stored : public Shadow { double GRMHDStorage[size]; Stored() : Shadow(GRMHDStorage) {} };
			// The advantage of a single total system storage is a linear storage
			// which can be copied, mapped and addressed easier.
		} // ns Conserved

		namespace Primitives {
			typedef TwoTheories::StateVector<double* const, MHD::Primitives::Shadow, ADMBase::Shadow> Shadow;
			struct Stored : public Shadow { double GRMHDStorage[size]; Stored() : Shadow(GRMHDStorage) {} };
		} // ns Primitives

		
		// TODO: Do not use these abbreviations
		
		typedef Conserved::Shadow Shadow;
		typedef Conserved::ConstShadow ConstShadow;
		typedef Conserved::Stored Stored;
		
	
		/**
		 * The link between tensor densities and tensors: QDensity is a real conserved
		 * vector and we remove the sqrt(gam.det).
		 * This structure uses the Conserved quantities.
		 **/
		template<class MHDtype>
		struct DensitiedState : public MHDtype, public ADMBase::Full {
			//double MHDStorage[MHD::size]; // what we used in the past
			//double MHDStorage[GRMHD::size]; // what we use now for simplicty for ADM
			
			constexpr static int debuggingGRMHDsize = 23;
			double MHDStorage[debuggingGRMHDsize]; // what we use even beyond to store the debugging stuff
			DensitiedState(const double* const QDensity) : MHDtype(MHDStorage), ADMBase::Full(QDensity) {
				// For the quick and dirty, assume all MHD parts to be next to each other.
				//for(int i=0; i<MHD::size; i++)
				//	MHDStorage[i] = QDensity[i] / gam.sqdet;
				
				// we copy everything for simple ADM access in the primitives
				std::copy_n(QDensity, debuggingGRMHDsize, MHDStorage); // DEBU
				for(int i=0; i<MHD::size; i++)
					MHDStorage[i] /= gam.sqdet;
			}
		};
				
		typedef DensitiedState<MHD::ConstShadow> Densitied;
		typedef DensitiedState<MHD::ConstShadowExtendable> DensitiedExtendable;
		
					
		/// a temporary space for global system parameters
		///
		/// Later, when we want to have runtime parameters, these parameters should probably
		/// be replaced by a pointer to a central parameter class.
		struct Parameters {
			// Ideal EOS:
			// 4/3 used in ADERDG3D-SRMHD-Z4 by MD, 01.Nov2016
			// 2.0 used for TOV stars
			static constexpr double gamma = 2.0; //4./3.;
	
			// Divergence cleaning:
			// 1.0 used in ADERDG3D-SRMHD-Z4 by MD, 01.Nov2016
			static constexpr double DivCleaning_a = 1.0;
			
			static constexpr double atmo_rho = 1e-11;
			static constexpr double atmo_press = 1e-7;
		};
		
		
		/**
		* The PDE structure wraps ("shadows") the read-only conserved vector Q in order to provide
		* the possibility to compute the flux and the source terms. It inherits the conserved
		* variables, primitive variables and ADM variables for convenient access when implementing
		* the source and the flux.
		*
		* It used to provide storage for the 3-Energy momentum tensor in the past, but now this is
		* stored directly in the methods because they need different forms (S_ij, S^ij or S^i_j)
		* which we address differently in the current formalism.
		**/
		// OLD:
		/*struct Cons2Prim : public Hydro::Conserved::ConstShadowExtendable, public Hydro::Primitives::ShadowExtendable, public Magneto::ConstShadowExtendable, public ADMBase::Full {
			Cons2Prim(double* const V, const double*const Q_) :
				Hydro::Conserved::ConstShadowExtendable(Q_),
				Hydro::Primitives::ShadowExtendable(V),
				Magneto::ConstShadowExtendable(Q_),
				ADMBase::Full(Q_)
		*/
		// NEW:
		struct Cons2Prim : public GRMHD::DensitiedExtendable, public Hydro::Primitives::ShadowExtendable, public Parameters {
			Cons2Prim(double* const V, const double* const Q_, bool _crash_on_failure=true) :
				GRMHD::DensitiedExtendable(Q_),
				Hydro::Primitives::ShadowExtendable(V),
				crash_on_failure(_crash_on_failure)
				{ prepare(); perform(); followup(); }

			// Quantities for computing the energy momentum tensor
			// as well as the C2P.
			double	WLorentz,	///< Lorentz factor (Gamma)
				WW,		///< Squared Lorentz factor (Gamma^2)
				RhoEnthWW,	///< rho*h*Gamma^2 as given by rtsafe
				SconScon,	///< S^2 = S_i*S^i: Squared conserved momentum
				BmagBmag,	///< B^2 = B_i*B^i: Squared magnetic field
				BmagVel,	///< B_i*v^i: Magn field times three velocity
				BmagScon,	///< B_i*S^i: Magn field times squared cons momentum
				VelVel,		///< V^2 = v_i*v^i: Squared three velocity
				ptot;		///< ptot = p_hydro + p_mag: Total pressure
				
			bool	crash_on_failure, ///< whether to abort the program in case of c2p failure
				failed; ///< Stores whether the C2P failed or not
			
			/// Reports about failed c2p. Pass your calling context as string.
			void report_c2p(const char* const context);
			
			/// Prepares the Conserved B_i, S^i as well as the quantities neccessary to compute
			/// The energy momentum tensor.
			void prepare();
			
			/// Computes the primitive variables with knowledge of all conserved/adm/local class variables.
			void perform();
			
			/// A modern interface to the RTSAFE C2P root finding procedure
			bool rtsafe(double& x, double& y);
			
			/// Compute the magnetic pressure and thelike. 
			void followup();
			
			/**
			* Copy the full state vector V, not only the hydro variables. You want this if you
			* want a "traditional" Cons2Prim in the Trento sense. Just call it as
			*    Cons2Prim(Q,V).copyFullStateVector()
			**/
			void copyFullStateVector();
			
			struct Stored;
		};

		/// A version which does not write to a shadowed storage but a local one
		struct Cons2Prim::Stored : public Cons2Prim {
			double V[Hydro::size];
			Stored(const double* const Q_) : Cons2Prim(V, Q_) {}
		};
		
		/**
		 * The raw PDE makes use of sqrt(det(g_ij)) stripping of the input state and
		 * adding on the output state. It also gets the full set of Primitives, the full
		 * ADM variables and all variables extended (i.e. vel_lo and vel^up).
		 * Therefore it is ready to do something.
		 **/
		struct RawPDE : public Cons2Prim::Stored {
			// the System defines the types of the constituents:
			typedef GRMHD::Shadow State;
			typedef const GRMHD::ConstShadow Gradient;
			typedef GenericLo<generic::shadow<Gradient, TDIM>, Gradient*> Gradients;
			typedef GRMHD::Shadow Flux;
			typedef GenericUp<generic::shadow<Flux, TDIM>, Flux*> Fluxes;
			
			RawPDE(const double* const Q) : Cons2Prim::Stored(Q) {}
		
			/// Conserved fluxes
			void flux(Fluxes& flux);
		
			/// Not to be used today.
			void eigenvalues(State& lambda, const int d);
		
			/// The BgradQ
			void nonConservativeProduct(const Gradients& grad, State& ncp);
		
			/// Sets the algebraic source
			void algebraicSource(State& Source_data);
			
			/// Adds the algebraic source (assumes Source=0 or similiar)
			void addAlgebraicSource(State& Source_data);
			
			/// Computes the fusedSource = algebraicSource - NCP. This is the classical RightHandSide.
			void fusedSource(const Gradients& grad, State& source);
		};
		
		/**
		 * Add the sqrt(det(g_ij)) factor to all MHD variables on every PDE output.
		 **/
		template<class P>
		struct DensitiedPDE : public P {
			typedef typename P::State State;
			typedef typename P::Fluxes Fluxes;
			typedef typename P::Gradients Gradients;
			
			DensitiedPDE(const double* const Q) : P(Q) {}

			void weight(State& state) {
				// For the quick and dirty, assume all MHD parts to be next to each other
				for(int i=0; i<MHD::size; i++) state.Q[i] *= this->gam.sqdet;
			}
			
			void flux(Fluxes& flux) {
				P::flux(flux); DFOR(i) weight(flux.up(i)); }
			void nonConservativeProduct(const Gradients& grad, State& ncp) {
				P::nonConservativeProduct(grad,ncp); weight(ncp); }
			void addAlgebraicSource(State& source) {
				P::addAlgebraicSource(source); weight(source); }
		};
		
		typedef DensitiedPDE<RawPDE> PDE;
		
		/**
		* The analytic GRHydro Primitive to Conserved conversion.
		* 
		* We shadow the prim input V onto
		*  1) Hydro::Primitives::ConstShadowExtendable , extendable for computing vel^2
		*  2) Magneto::ConstShadowExtendable , extendable for computing B^2
		*  3) ADMBase::Full because we need the full gij.
		* 
		* And we shadow the cons output Q onto
		*  1) Hydro::Conserved: Shadowing output Q
		*
		* Use `copyFullStateVector` if you want to copy the Magneto and ADM from V to Q.
		**/
		struct Prim2ConsRaw : public Hydro::Conserved::Shadow, public Magneto::ConstShadowExtendable, public Hydro::Primitives::ConstShadowExtendable, public ADMBase::Full {
			double W, enth, BmagBmag, BmagVel, VelVel;
			
			/**
			* Attention: We map {Hydro::Primitives,ADMBase,Magneto} <-> V
			*            and only the Hydro::Conserved <-> Q.
			* This is good if you set everything to the primitives.
			* However, if you have a mixed V-Q setting (i.e. B in Q, etc)
			* then you need a different mapping.
			**/
			Prim2ConsRaw(double*const Q_, const double* const V) :
				Hydro::Conserved::Shadow(Q_),
				Magneto::ConstShadowExtendable(V),
				Hydro::Primitives::ConstShadowExtendable(V),
				ADMBase::Full(V)
				{ prepare(); perform(); }
			
			void prepare();
			void perform();
			void copyFullStateVector(); ///< See Cons2Prim copyFullStateVector
		};
		
		struct Prim2ConsDensified : public Prim2ConsRaw {
			Prim2ConsDensified(double*const Q_, const double* const V) : Prim2ConsRaw(Q_,V) {
				// For the quick and dirty, assume all MHD parts to be next to each other
				for(int i=0; i<Hydro::size; i++)
					Q[i] = Q[i] * gam.sqdet;
			}
			
			void copyFullStateVector(); ///< Copy state vector incl Bmag*sqrt(detgam).
		};
		
		// this is the actual Prim2Cons you should use
		typedef Prim2ConsDensified Prim2Cons;
	} // ns GRMHD
	
	/**
	 * The GRHD (General-relativistic Hydrodynamics) implementation here is a cheap thing
	 * which fell of the GRMHD formulation. It just makes sure GRMHD sees Bmag=0 and does
	 * not store anything for Bmag.
	 **/
	namespace GRHD {
		struct empty {
		};
	} // ns GRHD (not GRMHD)
} // namespace SVEC

#endif /* GRMHD_PDE_CPP */
