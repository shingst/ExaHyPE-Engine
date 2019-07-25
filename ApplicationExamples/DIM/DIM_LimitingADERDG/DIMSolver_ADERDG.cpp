#include "DIMSolver_ADERDG.h"

#include "DIMSolver_ADERDG_Variables.h"

#include "PDE.h"

#include "InitialData.h"
// Used for the rieman-solver part
#include "peano/utils/Loop.h"
#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreBasis.h"
#include "peano/utils/Loop.h"
tarch::logging::Log DIM::DIMSolver_ADERDG::_log( "DIM::DIMSolver_ADERDG" );


void DIM::DIMSolver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

void DIM::DIMSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 14 + 0
  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
    // Fortran
    initialdata_(x, &t, Q);
  }
}

void DIM::DIMSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut)
  // Dimensions                        = 2
  // Number of variables + parameters  = 14 + 0

  // @todo Please implement/augment if required

  const int nVar = DIM::AbstractDIMSolver_ADERDG::NumberOfVariables;
  const int order = DIM::AbstractDIMSolver_ADERDG::Order;
  const int basisSize = order + 1;
  const int nDim = DIMENSIONS;

  double Qgp[nVar], F[nDim][nVar];

  std::memset(stateOut, 0, nVar * sizeof(double));
  std::memset(fluxOut, 0, nVar * sizeof(double));
  
  for(int i=0; i < basisSize; i++)  { // i == time
     const double weight = kernels::legendre::weights[order][i];
     const double xi = kernels::legendre::nodes[order][i];
     double ti = t + xi * dt;

     initialdata_(x, &ti, Qgp);
     pdeflux_(F[0], F[1], (nDim==3)?F[2]:nullptr, Qgp);
     for(int m=0; m < nVar; m++) {
        stateOut[m] += weight * Qgp[m];
        fluxOut[m] += weight * F[normalNonZero][m];
     }
  }
}

exahype::solvers::Solver::RefinementControl DIM::DIMSolver_ADERDG::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void DIM::DIMSolver_ADERDG::eigenvalues(const double* const Q,const int d,double* const lambda) {
  double nv[3] = {0.};
  nv[d] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}

void DIM::DIMSolver_ADERDG::mapDiscreteMaximumPrincipleObservables(double* const observables, const double* const Q) const {
  assertion(NumberOfDMPObservables==1);
  observables[0]=Q[12]; //extract alpha
}

bool DIM::DIMSolver_ADERDG::isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,
      const bool wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t) const {
  
  // Variant 1 (cheapest, currently works only in 2D)
  //  double outerRadius = 1.25*0.25;
  //  double innerRadius = 0.75*0.25;
  //  double radiusSquared = (center[0])*(center[0])+(center[1])*(center[1]);
  //
  //  if (
  //    radiusSquared<outerRadius*outerRadius &&
  //    radiusSquared>=innerRadius*innerRadius
  //  ) {
  //    return false;
  //  }

  // Variant 2 (comparably cheap, currently works only in 2D)
  //  double largest   = -std::numeric_limits<double>::max();
  //  double smallest  = +std::numeric_limits<double>::max();
  //
  //  kernels::idx3 idx_luh(Order+1,Order+1,NumberOfVariables);
  //  dfor(i,Order+1) {
  //    const double* Q = solution + idx_luh(i(1),i(0),0);
  //
  //    largest  = std::max (largest,   Q[12]);
  //    smallest = std::min (smallest,  Q[12]);
  //  }
  //
  //  if (!tarch::la::equals(largest,smallest,1e-6)) {
  //    return false;
  //  }

  // Variant 3 (expensive since it requires DMP, min and max gets confused by this variant for more than 2 mesh levels. Do not know why yet...)
  if (observablesMin[0] <= 0.0) {
    return false;
  }
  else {
    return true;
  }
}


void  DIM::DIMSolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  pdencp_(BgradQ, Q, gradQ);
}


