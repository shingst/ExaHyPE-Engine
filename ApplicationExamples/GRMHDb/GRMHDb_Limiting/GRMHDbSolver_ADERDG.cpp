// This file was generated by the ExaHyPE toolkit.
// It will NOT be regenerated or overwritten.
// Please adapt it to your own needs.
// 
// ========================
//   www.exahype.eu
// ========================

#include "GRMHDbSolver_ADERDG.h"
#include "GRMHDbSolver_ADERDG_Variables.h"

// User defined calls
#include "Tools.h"
#include "PDE.h"
#include "InitialData.h"
#include "tarch/parallel/Node.h"
#include "tarch/la/MatrixVectorOperations.h"

#include "GRMHDbSolver_FV.h"  // for FV patchsize and ghostlayerwidth

#include <algorithm>
#include <cstring> // memset
#include <string>
#include <math.h>
#include "peano/utils/Loop.h"

#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreBasis.h"


#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/multicore/Lock.h"

tarch::logging::Log GRMHDb::GRMHDbSolver_ADERDG::_log( "GRMHDb::GRMHDbSolver_ADERDG" );


void GRMHDb::GRMHDbSolver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // Tip: You find documentation for this method in header file "GRMHDb::GRMHDbSolver_ADERDG.h".
  
  // @todo Please implement/augment if required

	constexpr int order = GRMHDb::AbstractGRMHDbSolver_ADERDG::Order;
	constexpr int basisSize = AbstractGRMHDbSolver_FV::PatchSize;
	constexpr int Ghostlayers = AbstractGRMHDbSolver_FV::GhostLayerWidth;
    int mpirank = tarch::parallel::Node::getInstance().getRank();

	/**************************************************************************/
	static tarch::multicore::BooleanSemaphore initializationSemaphoreDG;
	tarch::multicore::Lock lock(initializationSemaphoreDG);
	/***************************************************/
	// everything in here is thread-safe w.r.t. the lock
	// call Fortran routines
	/***********************/
	pdesetup_(&mpirank);
	printf("\n******************************************************************");
	printf("\n**************<<<    PDESETUP   DONE       >>>****************************");
	printf("\n******************************************************************");
	lock.free();
}

void GRMHDb::GRMHDbSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
  // Tip: You find documentation for this method in header file "GRMHDb::GRMHDbSolver_ADERDG.h".
  // Tip: See header file "GRMHDb::AbstractGRMHDbSolver_ADERDG.h" for toolkit generated compile-time 
  //      constants such as Order, NumberOfVariables, and NumberOfParameters.
  // Dimensions                        = 3
  // Number of variables + parameters  = 19 + 0
  if (tarch::la::equals(t,0.0)) {
	/**************************************************************************/
	static tarch::multicore::BooleanSemaphore initialDataSemaphoreDG;
	tarch::multicore::Lock lock(initialDataSemaphoreDG);
	/***************************************************/
	// everything in here is thread-safe w.r.t. the lock
	// call Fortran routines
	/***********************/
	initialdata_(x, &t, Q);
	lock.free();
  }
}


void GRMHDb::GRMHDbSolver_ADERDG::boundaryValues(const double* const x, const double t, const double dt, const int faceIndex, const int normalNonZero, 
	const double * const fluxIn, const double* const stateIn, const double* const gradStateIn, 
	double *fluxOut, double* stateOut) {
  // Tip: You find documentation for this method in header file "GRMHDb::GRMHDbSolver_ADERDG.h".
  // Tip: See header file "GRMHDb::AbstractGRMHDbSolver_ADERDG.h" for toolkit generated compile-time 
  //      constants such as Order, NumberOfVariables, and NumberOfParameters.
  constexpr int numberOfVariables = AbstractGRMHDbSolver_ADERDG::NumberOfVariables;
  constexpr int numberOfParameters = AbstractGRMHDbSolver_ADERDG::NumberOfParameters;
  constexpr int numberOfData = numberOfVariables + numberOfParameters;
  //constexpr int numberOfData = AbstractGRMHDbSolver_FV::NumberOfVariables;
  const int order = GRMHDb::AbstractGRMHDbSolver_ADERDG::Order;
  const int basisSize = order + 1;
  const int nDim = DIMENSIONS;
  double Qgp[numberOfData],*F[nDim], Fs[nDim][numberOfData];

 //// STANDARD ANALYTICAL BOUNDARY CONDITIONS:
 for(int dd=0; dd<nDim; dd++)
	 F[dd] = Fs[dd];
 for(int i=0; i < basisSize; i++)  { // i == time
	  const double weight = kernels::legendre::weights[order][i];
	  const double xi = kernels::legendre::nodes[order][i];
	  double ti = t + xi * dt;

	  initialdata_(x, &ti, Qgp);
	  //pdeflux_(F[0], F[1], F[2], Qgp);
	  flux(Qgp, F);
          for(int m=0; m < numberOfData; m++) {
		  stateOut[m] += weight * Qgp[m];
		  fluxOut[m] += weight * Fs[normalNonZero][m];
	  }
 }
}

bool isInRefinementZone(const tarch::la::Vector<DIMENSIONS,double>& center){
	double radius = 8.12514;
	// lower left, upper right radius of cell
	double cen = tarch::la::norm2(center);
	double dr = 0.5;
  bool shouldRefine = (cen > (radius -dr) ) && ( cen  <= (radius+dr) ); 
  return shouldRefine;
}

exahype::solvers::Solver::RefinementControl GRMHDb::GRMHDbSolver_ADERDG::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
    if(isInRefinementZone(center))
        return exahype::solvers::Solver::RefinementControl::Refine;
    return exahype::solvers::Solver::RefinementControl::Keep;
}


//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
/// file and its header and rerun the toolkit
//*****************************************************************************
void GRMHDb::GRMHDbSolver_ADERDG::eigenvalues(const double* const Q,const int d,double* const lambda) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 19 + 0
  
	constexpr int numberOfVariables = AbstractGRMHDbSolver_ADERDG::NumberOfVariables;
	constexpr int numberOfParameters = AbstractGRMHDbSolver_ADERDG::NumberOfParameters;
	constexpr int numberOfData = numberOfVariables + numberOfParameters;

	double nv[DIMENSIONS] = {0.};
	nv[d] = 1.0;
	pdeeigenvalues_(lambda, Q, nv);
}



void GRMHDb::GRMHDbSolver_ADERDG::referenceSolution(const double* const x,double t, double* Q) {
	constexpr int numberOfVariables = AbstractGRMHDbSolver_ADERDG::NumberOfVariables;
	constexpr int numberOfParameters = AbstractGRMHDbSolver_ADERDG::NumberOfParameters;
	constexpr int numberOfData = numberOfVariables + numberOfParameters;
	int iErr;
	double Qcons[numberOfData];
	iErr = 0;

	initialdata_(x, &t, &Qcons[0]);
	pdecons2prim_(Q, &Qcons[0], &iErr);
}


void GRMHDb::GRMHDbSolver_ADERDG::flux(const double* const Q,double** const F) {
  // Tip: You find documentation for this method in header file "GRMHDb::GRMHDbSolver_ADERDG.h".
  // Tip: See header file "GRMHDb::AbstractGRMHDbSolver_ADERDG.h" for toolkit generated compile-time 
  //      constants such as Order, NumberOfVariables, and NumberOfParameters.
  // Dimensions                        = 3
  // Number of variables + parameters  = 19 + 0
  pdeflux_(F[0], F[1],F[2], Q);
}



void  GRMHDb::GRMHDbSolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  // Tip: You find documentation for this method in header file "GRMHDb::GRMHDbSolver_ADERDG.h".
  // Tip: See header file "GRMHDb::AbstractGRMHDbSolver_ADERDG.h" for toolkit generated compile-time 
  //      constants such as Order, NumberOfVariables, and NumberOfParameters.
  pdencp_(BgradQ, Q, gradQ);
}

void GRMHDb::GRMHDbSolver_ADERDG::mapDiscreteMaximumPrincipleObservables(
	double* const observables,
	const double* const Q) const {
	if (NumberOfDMPObservables>0) {
		std::copy_n(Q,NumberOfDMPObservables,observables);
	}
}

bool GRMHDb::GRMHDbSolver_ADERDG::vetoDiscreteMaximumPrincipleDecision(
		const double* const                         solution,
		const double* const                         localObservablesMin,
		const double* const                         localObservablesMax,
		const bool                                  wasTroubledInPreviousTimeStep,
		const tarch::la::Vector<DIMENSIONS, double>& center,
		const tarch::la::Vector<DIMENSIONS, double>& dx,
		const double                                timeStamp) const {
	return false; // do not veto = false; veto = true;	
}


bool GRMHDb::GRMHDbSolver_ADERDG::isPhysicallyAdmissible(
	const double* const solution,
	const double* const observablesMin, const double* const observablesMax,
	const bool wasTroubledInPreviousTimeStep,
	const tarch::la::Vector<DIMENSIONS, double>& center,
	const tarch::la::Vector<DIMENSIONS, double>& dx,
	const double t) const {

	double radius = 8.12514;
	double cen = tarch::la::norm2(center);
	double dr = 1.5;

	bool shouldLimit = (cen > (radius -dr) ) && ( cen  <= (radius+dr) ); 
	// return TRUE if the cell does not need limited
	return !shouldLimit;
}

int getGeometricLoadBalancingWeight(
        const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS,double>& cellSize) {
  double cen = tarch::la::norm2(cellCentre);
  double dr = std::max(0.1, tarch::la::norm2(cellSize));
  if(cen < 8.12514+dr)
    return 27;
  return 1;
}

#ifdef OPT_KERNELS
#include "kernels/GRMHDb_GRMHDbSolver_ADERDG/Kernels.h"

using namespace GRMHDb::GRMHDbSolver_ADERDG_kernels::aderdg;

void GRMHDb::GRMHDbSolver_ADERDG::riemannSolver(double* const FL, double* const FR, const double* const QL, const double* const QR, const double t, const double dt, const tarch::la::Vector<DIMENSIONS, double>& lengthScale, const int direction, bool isBoundaryFace, int faceIndex) {
	assertion2(direction >= 0, dt, direction);
	assertion2(direction<DIMENSIONS, dt, direction);
	GRMHDb::GRMHDbSolver_ADERDG_kernels::aderdg::riemannSolver(*static_cast<GRMHDbSolver_ADERDG*>(this), FL, FR, QL, QR, t, dt, direction);
	constexpr int order = GRMHDb::AbstractGRMHDbSolver_ADERDG::Order;
	constexpr int basisSize = order + 1;
	constexpr int numberOfVariables = GRMHDb::GRMHDbSolver_ADERDG_kernels::aderdg::getNumberOfVariablePadded();   //              AbstractGRMHDbSolver_ADERDG::NumberOfVariables;
	// avoid spurious numerical diffusion (ony for Cowling approximation)
	
#ifdef Dim2
	kernels::idx2 idx_FLR(basisSize, numberOfVariables);
        for (int i = 0; i < basisSize; i++) {
                        //resetNumericalDiffusionOnADM(FL + idx_FLR(i, j, 0));
                        //resetNumericalDiffusionOnADM(FR + idx_FLR(i, j, 0));
                        double* FLL = FL + idx_FLR(i, 0);
                        double* FRR = FR + idx_FLR(i, 0);
                        for (int m = 9; m < numberOfVariables; m++) {
                                FLL[m] = 0.0;
                                FRR[m] = 0.0;
                }
        }
#else
	kernels::idx3 idx_FLR(basisSize, basisSize, numberOfVariables);
	for (int i = 0; i < basisSize; i++) {
		for (int j = 0; j < basisSize; j++) {
			//resetNumericalDiffusionOnADM(FL + idx_FLR(i, j, 0));
			//resetNumericalDiffusionOnADM(FR + idx_FLR(i, j, 0));
			double* FLL = FL + idx_FLR(i, j, 0);
			double* FRR = FR + idx_FLR(i, j, 0);
			for (int m = 9; m < numberOfVariables; m++) {
				FLL[m] = 0.0;
				FRR[m] = 0.0;
			}
		}
	}
#endif

}


#else

#include "kernels/aderdg/generic/Kernels.h" 
void GRMHDb::GRMHDbSolver_ADERDG::riemannSolver(double* const FL, double* const FR, const double* const QL, const double* const QR, const double t, const double dt, const tarch::la::Vector<DIMENSIONS, double>& dx, const int direction, bool isBoundaryFace, int faceIndex) {
	assertion2(direction >= 0, dt, direction);
	assertion2(direction<DIMENSIONS, dt, direction);
	kernels::aderdg::generic::c::riemannSolverNonlinear<true, false, GRMHDbSolver_ADERDG>(*static_cast<GRMHDbSolver_ADERDG*>(this), FL, FR, QL, QR, t, dt, dx, direction);
 

	constexpr int order = GRMHDb::AbstractGRMHDbSolver_ADERDG::Order;
	constexpr int basisSize = order + 1;
	constexpr int numberOfVariables = AbstractGRMHDbSolver_ADERDG::NumberOfVariables;
	kernels::idx3 idx_FLR(basisSize, basisSize, numberOfVariables);
	for (int i = 0; i < basisSize; i++) {
		for (int j = 0; j < basisSize; j++) {
			//resetNumericalDiffusionOnADM(FL + idx_FLR(i, j, 0));
			//resetNumericalDiffusionOnADM(FR + idx_FLR(i, j, 0));
			double* FLL = FL + idx_FLR(i, j, 0);
			double* FRR = FR + idx_FLR(i, j, 0);
			for (int m = 9; m < numberOfVariables; m++) {
				FLL[m] = 0.0;
				FRR[m] = 0.0;
			}
		}
	}
}
#endif
