#include "GRMHDSolver_ADERDG.h"
#include "GRMHDSolver_ADERDG_Variables.h"

#include <algorithm> // fill_n
#include <cstring> // memset

#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreBasis.h"
#include "peano/utils/Dimensions.h" // Defines DIMENSIONS

//namespace GRMHD { constexpr int nVar = GRMHDSolver_FV::NumberOfVariables; } // ensure this is 19 or so
#include "InitialData/InitialData.h"
#include "PDE/PDE-GRMHD-ExaHyPE.h"
#include "DebuggingHelpers.h"
#include "BoundaryConditions.h"

//namespace SVEC::GRMHD::ExaHyPEAdapter = ExaGRMHD;
#define ExaGRMHD SVEC::GRMHD::ExaHyPEAdapter
using namespace tensish;



constexpr int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;
constexpr int order = GRMHD::AbstractGRMHDSolver_ADERDG::Order;
constexpr int basisSize = order + 1;
constexpr int nDim = DIMENSIONS;

tarch::logging::Log GRMHD::GRMHDSolver_ADERDG::_log( "GRMHD::GRMHDSolver_ADERDG" );

void GRMHD::GRMHDSolver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
	// feenableexcept(FE_INVALID | FE_OVERFLOW);  // Enable all floating point exceptions but FE_INEXACT
	
	mexa::mexafile mf = mexa::fromSpecfile(constants.getAllAsOrderedMap(), constants.toString());
	
	std::cout << "Mexa configuration: \n" << mf.toString();
	std::cout << "ID configuration: \n" << mf("initialdata").toString();
	std::cout << "ID NAME: '" << mf.get("initialdata/name").get_string() << "'\n";
	std::cout << "ID subquery NAME: '" << mf("initialdata").get("name").get_string() << "'\n";
	
	GlobalInitialData::getInstance().setByParameters(mf).prepare();
	GlobalBoundaryConditions::getInstance().initializeDG(this).readParameters(mf("boundaries"));
}


/**
 * This is a modified riemannSolver which turns off the numerical diffusion
 * for material parameters. Since we need the derivatives of these parameters
 * which are not computed by the engine, we don't treat them as parameters in
 * ExaHyPE and instead do this patching here.
 **/
#include "kernels/aderdg/generic/Kernels.h"
void GRMHD::GRMHDSolver_ADERDG::riemannSolver(double* const FL,double* const FR,const double* const QL,const double* const QR,const double dt,const int normalNonZeroIndex,bool isBoundaryFace, int faceIndex) {
	kernels::aderdg::generic::c::riemannSolverNonlinear<true,GRMHDSolver_ADERDG>(*static_cast<GRMHDSolver_ADERDG*>(this),FL,FR,QL,QR,dt,normalNonZeroIndex);
	
	#if DIMENSIONS == 3
	kernels::idx3 idx_FLR(basisSize, basisSize, NumberOfVariables);
	for (int i = 0; i < basisSize; i++) {
	for (int j = 0; j < basisSize; j++) {
		resetNumericalDiffusionOnADM(FL+idx_FLR(i, j, 0));
		resetNumericalDiffusionOnADM(FR+idx_FLR(i, j, 0));
	}}
	#elif DIMENSIONS == 2
	kernels::idx2 idx_FLR(basisSize, NumberOfVariables);
	for (int i = 0; i < basisSize; i++) {
		resetNumericalDiffusionOnADM(FL+idx_FLR(i, 0));
		resetNumericalDiffusionOnADM(FR+idx_FLR(i, 0));
	}
	#endif
}

void GRMHD::GRMHDSolver_ADERDG::mapDiscreteMaximumPrincipleObservables(double* const observables, const double* const Q) const {
	// ensure NumberOfDMPObservables == 2
	observables[0] = Q[0]; // conserved density
	observables[1] = Q[4]; // conserved tau
	
	// Q: Shall we make a C2P here and map the primitive RMD and pressure?
	// observables[2] = prim.rho;
	// observables[3] = prim.press;
	// Can we use the DMP to check for superluminant velocities if we store a
	// observables[3] = prim.VelVel;
}

bool pointwiseisPhysicallyAdmissible(const double* const Q) {
	// pre-C2P plausibility check on the conserved quantities
	// This is *probably* already caught by the DMP.
	if(Q[0]<0) return false;

	// try to do a C2P and check the primitives for correctness
	SVEC::GRMHD::Cons2Prim::Stored prim(Q);
	return (!prim.failed) &&
	       (prim.rho > 0) && (prim.press > 0) && (prim.VelVel < 1);
}

bool GRMHD::GRMHDSolver_ADERDG::isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,
      const bool wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t) const {

	// Q: Do we have to check the observables{Min,Max} whether they are useful?

	if(observablesMin[0] < 0) return false; // conserved density
	if(observablesMin[1] < 0) return false; // conserved tau

	// delegate to pointwise evaluation
	constexpr int basisSize3D = (DIMENSIONS==3) ? basisSize : 1;
	kernels::idx4 idx(basisSize3D,basisSize,basisSize,NumberOfVariables);
	for(int iz = 0; iz < basisSize3D; iz++) {
	for(int iy = 0; iy < basisSize;   iy++) {
	for(int ix = 0; ix < basisSize;   ix++) {
		if(!pointwiseIsPhysicallyAdmissible(solution+idx(iz,iy,ix,0)))
			return false;
	}}}

	return true;
}


void GRMHD::GRMHDSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
  if (tarch::la::equals(t,0.0)) {
    InitialData(x,t,Q);
  } else {
	constexpr bool overwriteADMalways = true;
	if(overwriteADMalways) {
		overwriteADM(x,t,Q);
	}
  }
}

void GRMHD::GRMHDSolver_ADERDG::eigenvalues(const double* const Q,const int d,double* const lambda) {
	// Provide NVARS eigenvalues
	//PDE::eigenvalues(Q, d, lambda);
	
	// Just set all eigenvalues to one.
	NVARS(i) lambda[i] = 1.0;
}


void GRMHD::GRMHDSolver_ADERDG::flux(const double* const Q,double** const F) {
	
	// as TDIM is 3 and DIMENSIONS is 2, come up with this dirty wrapper:
	#if DIMENSIONS == 2
	double *FT[3], Fz[nVar];
	FT[0] = F[0];
	FT[1] = F[1];
	FT[2] = Fz;
	ExaGRMHD::Fluxes f(FT);
	#else
	ExaGRMHD::Fluxes f(F);
	#endif
	
	// Debugging: set all fluxes to zero
	for(int d=0;d<DIMENSIONS;d++) NVARS(i) F[d][i] = 0;
	
	SVEC::GRMHD::PDE pde(Q);
	pde.flux(f);
	f.zeroMaterialFluxes();
	//ExaGRMHD::flux(Q, FT).zeroMaterialFluxes();
	
	//DFOR(d) { zero2Din3D(F[d]); zeroHelpers(F[d]); }
	for(int d=0;d<DIMENSIONS;d++) zeroHelpers(F[d]);
	
	for(int d=0;d<DIMENSIONS;d++) NVARS(i) { if(!std::isfinite(F[d][i])) { printf("F[%d][%d] = %e\n", d, i, F[d][i]); std::abort(); } }
}


void GRMHD::GRMHDSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
	// for debugging, to make sure BC are set correctly
	/*
	double snan = std::numeric_limits<double>::signaling_NaN();
	double weird_number = -1.234567;
	std::memset(stateOut, weird_number, nVar * sizeof(double));
	std::memset(fluxOut,  weird_number, nVar * sizeof(double));
	*/
	
	// Let it be managed by user parameters
	GlobalBoundaryConditions::getInstance().apply(ADERDG_BOUNDARY_CALL);
	
	/*
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

		InitialData(x, ti, Qgp);
		flux(Qgp, F);

		for(int m=0; m < nVar; m++) {
			stateOut[m] += weight * Qgp[m];
			fluxOut[m] += weight * Fs[d][m];
		}
	}
	*/
	/*
	// Neutron star Reflective + Outflow BC
	setNeutronStarBoundaryConditions(faceIndex, d, stateIn, stateOut);
	
	double Fs[nDim][nVar], *F[nDim];
	for(int dd=0; dd<nDim; dd++) F[dd] = Fs[dd];
	F[d] = fluxOut;
	flux(stateOut, F);
	*/

	// ONLY FOR DEBUGGING:
	
	//NVARS(i) printf("stateOut[%d]=%e\n", i, stateOut[i]);
	//NVARS(i) printf("fluxOut[%d]=%e\n", i, fluxOut[i]);
	//std::abort();

}

//#include "kernels/aderdg/generic/c/computeGradients.cpph"
exahype::solvers::Solver::RefinementControl GRMHD::GRMHDSolver_ADERDG::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
	/*
	// TODO: Continue here.
	// Compute the gradient of the Density
	using namespace kernels::aderdg::generic::c;
	typedef GRMHD::AbstractGRMHDSolver_ADERDG self;
	double gradDens[basisSizeD(self::Order+1)*DIMENSIONS];
	cpmstexpr int positiondens = 0;
	computeGradQi<self>(gradDens, luh, positiondens, dx);
	
	// Apply the Loehner scheme criterion
	// cf. kernels/aderdg/generic/c/loehnerScheme.cpph for prelimiinary work
	*/

	// @todo Please implement/augment if required
	return exahype::solvers::Solver::RefinementControl::Keep;
}

void GRMHD::GRMHDSolver_ADERDG::algebraicSource(const double* const Q,double* const S) {
	// No algebraic source
	for(int i=0;i<nVar;i++) S[i] = 0;
	
	return;	
	ExaGRMHD::algebraicSource(Q, S).zero_adm();
	zeroHelpers(S);
	//zero2Din3D(S_);
}

void GRMHD::GRMHDSolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ, double* const BgradQ) {
	// debugging: set BgradQ initially to zero:
	NVARS(i) BgradQ[i] = 0;
	
	// check at this place
	/*
	constexpr double x=0.023143948067657905, y=0.68981061473432448;
	if(tarch::la::equals(Q[19],x,1e-8) && tarch::la::equals(Q[20],y,1e-8)) {
		printf("Debug guy\n");
	}
	}
	*/
	
	// as we only have DIMENSIONS == 2 but TDIM == 3, construct a zero gradient
	#if DIMENSIONS == 2
	double gradZ[nVar] = {0.};
	// deconstruct the gradient matrix because we use more than the 19 variables of GRMHD.
	ExaGRMHD::nonConservativeProduct(Q,gradQ+0,gradQ+nVar,gradZ,BgradQ).zero_adm();
	#else
	ExaGRMHD::nonConservativeProduct(Q,gradQ,BgradQ).zero_adm();
	#endif
	
	zeroHelpers(BgradQ);
	//zero2Din3D(BgradQ);
	
	// the NCP must be zero in flat space
	/*
	constexpr double eps = 1e-10;
	NVARS(i) { if(BgradQ[i]>eps) { printf("BgradQ[%d] = %e\n", i, BgradQ[i]); std::abort(); } }
	*/
}

// The optimized fusedSource
/*
void GRMHD::GRMHDSolver_ADERDG::fusedSource(const double* const Q, const double* const gradQ, double* const S_) {
	PDE::Source S(S_);
	PDE pde(Q);
	GRMHD::AbstractGRMHDSolver_ADERDG::ReadOnlyVariables var(Q);
	constexpr double eps = 1e-8;
	// ExaHyPE workaround:
	// if the input has zeros at weird places, don't do anything
	if(std::abs(pde.gam.det) < eps || var.check()!=magicCheck) {
		printf("Weird FusedSource input (det=%e), skipping. ",pde.gam.det);
		//SYMFOR(i,j) printf("gam(%d,%d)=%e\n", i,j, pde.gam.lo(i,j));
		std::cout << "pos: " << var.pos() << " check: " << var.check() <<  std::endl;
		// set everything to NANs
		NVARS(m) S_[m] = NAN;
		return;
	} else {
		//printf("Good FusedSource input (det=%e)\n",pde.gam.det);
		//NVARS(m) printf("FusedSource input: Q[%d]=%e\n", m, Q[m]);
		
		
		//std::cout << "Good FusedSource, x= " << var.pos() <<  std::endl;
	}
	pde.RightHandSide(gradQ, S);
	S.zero_adm();
	zeroHelpers(S_);
	zero2Din3D(S_);
}
*/
