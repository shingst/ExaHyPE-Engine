#include "CCZ4Solver_ADERDG.h"

#include "CCZ4Solver_ADERDG_Variables.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/KernelUtils.h" // matrix indexing
#include <cstring> // memset
#include "Fortran/PDE.h"
#include "InitialData/InitialData.h"
#include "Parameters/mexa.h"

// This doesn't seem to be needed anymore
//#include "exahype/disableOptimization.h" // we experience compiler bugs sometimes.

const int nDim = DIMENSIONS;
const int nVar = CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables;
const int order = CCZ4::AbstractCCZ4Solver_ADERDG::Order;
const int basisSize = order+1;

tarch::logging::Log CCZ4::CCZ4Solver_ADERDG::_log( "CCZ4::CCZ4Solver_ADERDG" );

using namespace CCZ4Fortran;


void CCZ4::CCZ4Solver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
	mexa::mexafile mf = mexa::fromSpecfile(constants.getAllAsOrderedMap(), constants.toString());
	
	// Debugging
	std::cout << "Mexa configuration: \n" << mf.toString();
	std::cout << "ID configuration: \n" << mf("initialdata").toString();
	std::cout << "ID NAME: '" << mf.get("initialdata/name").get_string() << "'\n";
	
	// PDE Runtime Parameters
	GlobalPDEParameters::getInstance().setByParameters(mf);
	
	// Initial Data
	GlobalInitialData::getInstance().setByParameters(mf);
	
	// Boundary conditions for CCZ4: Currently always exact.
	//GlobalBoundaryConditions::getInstance().initializeDG(this).readParameters(mf("boundaries"));
	
	// Limiter Runtime Steering
}


bool CCZ4::CCZ4Solver_ADERDG::isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,
      const bool wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t) const {
	// lower left, uppe right radius of cell
	double l = tarch::la::norm2(center - dx/2.0);
	double r = tarch::la::norm2(center + dx/2.0);
	constexpr double radius_shell = 0.5; // limit around this shell
	bool isAdmissible = (l > radius_shell) || (r <= radius_shell);
	//printf("Cell has l=%f,r=%f => isAdmissible=%s\n", l, r, isAdmissible?"true":"false");

	return isAdmissible;
}

void CCZ4::CCZ4Solver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
	AdjustPointSolution(x,t,dt,Q);
}

void CCZ4::CCZ4Solver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
	double Qgp[nVar];

	// Impose exact boundary conditions
	std::memset(stateOut, 0, nVar * sizeof(double));
	std::memset(fluxOut, 0, nVar * sizeof(double));

	//double F[nDim][nVar];
	//kernels::idx2 F_idx(nDim, nVar);

	const int basisSize = order+1;

	// Integrate solution in gauss points (Qgp) in time
	for(int i=0; i < basisSize; i++)  { // i == time
	const double weight = kernels::gaussLegendreWeights[order][i];
	const double xi = kernels::gaussLegendreNodes[order][i];
	double ti = t + xi * dt;

	InitialData(x,ti,Qgp);
	//id->adjustPointSolution(x,ti,Qgp);
	//PDE::flux(Qgp, F);
	for(int m=0; m < nVar; m++) {
		//if(m==checkm) printf("fluxOut[%d] += %.20e\n", m, weight * F[normalNonZero][m]);
		stateOut[m] += weight * Qgp[m];
		//fluxOut[m] += weight * F[d][m];
	}
	}
}


exahype::solvers::Solver::RefinementControl CCZ4::CCZ4Solver_ADERDG::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
	// @todo Please implement/augment if required
	return exahype::solvers::Solver::RefinementControl::Keep;
}


void CCZ4::CCZ4Solver_ADERDG::eigenvalues(const double* const Q,const int d,double* const lambda) {
	PDE::eigenvalues(Q, d, lambda);
}

void CCZ4::CCZ4Solver_ADERDG::algebraicSource(const double* const Q,double* const S) {
	PDE::algebraicSource(Q, S);
}

void  CCZ4::CCZ4Solver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
	PDE::nonConservativeProduct(Q, gradQ, BgradQ);
}

