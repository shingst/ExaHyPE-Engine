#include "GRMHDSolver_FV.h"
#include "GRMHDSolver_FV_Variables.h"

#include <algorithm> // fill_n

#include "peano/utils/Dimensions.h"

#include "PDE/PDE-GRMHD-ExaHyPE.h"
//namespace SVEC::GRMHD::ExaHyPEAdapter = ExaGRMHD;
#define ExaGRMHD SVEC::GRMHD::ExaHyPEAdapter

#include "InitialData/InitialData.h"
#include "BoundaryConditions.h"
#include "DebuggingHelpers.h"

constexpr int nVar = GRMHD::AbstractGRMHDSolver_FV::NumberOfVariables;
constexpr int nDim = DIMENSIONS;
using namespace GRMHD;

tarch::logging::Log GRMHD::GRMHDSolver_FV::_log( "GRMHD::GRMHDSolver_FV" );

void GRMHD::GRMHDSolver_FV::init(const std::vector<std::string>& cmdlineargs, exahype::Parser::ParserView& constants) {
	mexa::mexafile mf = mexa::fromSpecfile(constants.getAllAsOrderedMap(), constants.toString());
	
	std::cout << "Mexa configuration: \n" << mf.toString();
	std::cout << "ID configuration: \n" << mf("initialdata").toString();
	std::cout << "ID NAME: '" << mf.get("initialdata/name").get_string() << "'\n";
	std::cout << "ID subquery NAME: '" << mf("initialdata").get("name").get_string() << "'\n";
	GlobalInitialData::getInstance().setByParameters(mf);
	GlobalBoundaryConditions::getInstance().initializeFV(this).readParameters(mf("boundaries"));
}

void GRMHD::GRMHDSolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* Q) {
	// Set the 9 SRMHD variables (D,S_j,tau,B^j) and the 10 [11] ADM material parameters (N^i,g_ij,[detg])
	if(tarch::la::equals(t,0.0)) {
		InitialData(x,t,Q);
	} else {
		constexpr bool floorAtmosphere = false;
		if(floorAtmosphere) {
			// Try to treat the atmosphere
			SVEC::GRMHD::Shadow qout(Q);
			SVEC::GRMHD::Cons2Prim::Stored prim(Q);
			if(prim.failed || prim.rho < prim.atmo_rho) {
				// floor the primitives
				prim.rho = prim.atmo_rho;
				DFOR(i) prim.vel.up(i) = 0;
				prim.press = prim.atmo_press;
				
				SVEC::GRMHD::Prim2Cons(Q, prim.V);
			}
		}
		
		// Doing this in the FV limiter is probably not a bad idea at all, as
		// initial data might have been broken apart already in the ADER polynomial
		// and then after projection,
		// i.e. we might limit around the NS core
		constexpr bool overwriteADMalways = true;
		if(overwriteADMalways) {
			overwriteADM(x,t,Q);
		}
	}
}

#include "kernels/finitevolumes/riemannsolvers/c/riemannsolvers.h"

double GRMHD::GRMHDSolver_FV::riemannSolver(double* fL, double *fR, const double* qL, const double* qR, int normalNonZero) {

  // Default FV Riemann Solver
  // TODO: We could write an optimized one here which does not ask for the eigenvalues.
  double lambda= kernels::finitevolumes::riemannsolvers::c::rusanov<true, true, GRMHDSolver_FV>(*static_cast<GRMHDSolver_FV*>(this), fL,fR,qL,qR,normalNonZero);
  resetNumericalDiffusionOnADM(fL);
  resetNumericalDiffusionOnADM(fR);

  return lambda;
}



void GRMHD::GRMHDSolver_FV::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
	// Provide eigenvalues for the 9 SRMHD variables (D,S_j,tau,B^j),
	// we split off the 11 ADM material parameters (N^i,g_ij,detg)
	
	// this does a C2P for ones. We don't like that.
	// ExaGRMHD::eigenvalues(Q, dIndex, lambda);
	
	NVARS(i) lambda[i] = 1.0;
	
	
	//NVARS(m) printf("EV[%d]=%f\n", m, lambda[m]);
}


void GRMHD::GRMHDSolver_FV::flux(const double* const Q, double** F) {
	// Provide fluxes for the 9 SRMHD variables (D,S_j,tau,B^j),
	// we split off the 11 ADM material parameters (N^i,g_ij,detg)
	//PDE(Q).flux(F);
	
	// zero detection
	if(isUnphysical(Q)) {
		//printf("WRONG Flux input, all zero!\n");
		//NVARS(m) printf("Q[%d]=%e\n", m, Q[m]);
		//std::abort();
		// Set everything to some neutral value
		DFOR(i) NVARS(m) F[i][m] = 0;
		return;
	}
	
	// as TDIM is 3 and DIMENSIONS is 2, come up with this dirty wrapper:
	#if DIMENSIONS == 2
	double *FT[3], Fz[nVar];
	FT[0] = F[0];
	FT[1] = F[1];
	FT[2] = Fz;
	ExaGRMHD::flux(Q,FT).zeroMaterialFluxes();
	#else
	ExaGRMHD::flux(Q,F).zeroMaterialFluxes();
	#endif
}



void GRMHD::GRMHDSolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateIn,
    double* stateOut) {
	// SET BV for all 9 variables + 11 parameters.
	
	// Let it be managed by user parameters
	GlobalBoundaryConditions::getInstance().apply(FV_BOUNDARY_CALL);
	
	// EXACT FV AlfenWave BC
	// InitialData(x, t, stateOut);
	
	// NEUTRON STAR reflective + outflow BC:
	// setNeutronStarBoundaryConditions(faceIndex, d, stateIn, stateOut);
}

void GRMHD::GRMHDSolver_FV::algebraicSource(const double* const Q,double* S) {
	// set source to zero
	NVARS(i) S[i] = 0;
	return;
	
	if(isUnphysical(Q)) {
		NVARS(i) S[i] = NAN;
		return;
	}
	

	
	// Todo: Should not do a C2P in order to fill out the algebraic source.
	// Or should just return zero here as the Fortran code does.
	
	ExaGRMHD::algebraicSource(Q,S).zero_adm();
}

void GRMHD::GRMHDSolver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
	// set ncp to zero
	//NVARS(i) BgradQ[i] = 0;
	//return;
	
	if(isUnphysical(Q)) {
		NVARS(i) BgradQ[i] = 0;
		return;
	}
	
	// debugging: set BgradQ to zero:
	NVARS(i) BgradQ[i] = 0;
	
	// as we only have DIMENSIONS == 2 but TDIM == 3, construct a zero gradient
	#if DIMENSIONS == 2
	double gradZ[nVar] = {0.};
	
	// deconstruct the gradient matrix because we use more than the 19 variables of GRMHD.
	ExaGRMHD::nonConservativeProduct(Q,gradQ+0,gradQ+nVar,gradZ,BgradQ).zero_adm();
	#else
	ExaGRMHD::nonConservativeProduct(Q,gradQ,BgradQ).zero_adm();
	#endif
}

// Formerly optimized fusedSource
/*
void GRMHD::GRMHDSolver_FV::fusedSource(const double* const Q, const double* const gradQ, double* S_) {
	
	// ExaHyPE workaround:
	// if the input are zeros everywhere, complain
	if(isUnphysical(Q)) {
		//printf("WRONG FusedSource input, all zero!\n");
		//NVARS(m) printf("Q[%d]=%e\n", m, Q[m]);
		//std::abort();
		// Set everything to NAN
		NVARS(m) S_[m] = NAN;
		return;
	}
	
	PDE::Source S(S_);
	PDE pde(Q);
	GRMHD::AbstractGRMHDSolver_FV::ReadOnlyVariables var(Q);

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
}
*/
