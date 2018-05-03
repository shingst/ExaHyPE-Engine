// This file was generated by the ExaHyPE toolkit.
// It will NOT be regenerated or overwritten.
// Please adapt it to your own needs.
// 
// ========================
//   www.exahype.eu
// ========================

#include "MySWESolver_ADERDG.h"
#include "InitialData.h"
#include "MySWESolver_ADERDG_Variables.h"

#include "kernels/KernelUtils.h"

using namespace kernels;

const double grav = 9.81;

const double epsilon = 1e-1;

tarch::logging::Log SWE::MySWESolver_ADERDG::_log( "SWE::MySWESolver_ADERDG" );


void SWE::MySWESolver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required

}

void SWE::MySWESolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 4 + 0

  if (tarch::la::equals(t,0.0)) {
    initialData(x, Q);
  }
}

void SWE::MySWESolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 4 + 0

  //Outflow (not working)
//    for(int i=0; i < NumberOfVariables; i++) {
//        fluxOut[i] = fluxIn[i];
//        stateOut[i] = stateIn[i];
//    }

  //Wall
  std::copy_n(stateIn, NumberOfVariables, stateOut);
  stateOut[1+normalNonZero] =  -stateOut[1+normalNonZero];
  double _F[2][NumberOfVariables]={0.0};
  double* F[2] = {_F[0], _F[1]};
  flux(stateOut,F);
  std::copy_n(F[normalNonZero], NumberOfVariables, fluxOut);
}

exahype::solvers::Solver::RefinementControl SWE::MySWESolver_ADERDG::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void SWE::MySWESolver_ADERDG::eigenvalues(const double* const Q,const int d,double* lambda) {
  /// Dimensions                        = 2
  // Number of variables + parameters  = 4 + 0

  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const double c = std::sqrt(grav*vars.h());
  const double ih = 1./vars.h();
  double u_n = Q[d + 1] * ih;

  eigs.h() = u_n + c;
  eigs.hu() = u_n - c;
  eigs.hv() = u_n;
  eigs.b() = 0.0;
}


void SWE::MySWESolver_ADERDG::flux(const double* const Q,double** F) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 4 + 0

  ReadOnlyVariables vars(Q);

  const double ih = 1./vars.h();

  double* f = F[0];
  double* g = F[1];

  f[0] = vars.hu();
  f[1] = vars.hu()*vars.hu()*ih + 0.5*grav*vars.h()*vars.h();
  f[2] = vars.hu()*vars.hv()*ih;
  f[3] = 0.0;

  g[0] = vars.hv();
  g[1] = vars.hu()*vars.hv()*ih;
  g[2] = vars.hv()*vars.hv()*ih + 0.5*grav*vars.h()*vars.h();
  g[3] = 0.0;
  
}



void  SWE::MySWESolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  idx2 idx_gradQ(DIMENSIONS,NumberOfVariables);

  BgradQ[0] = 0.0;
  BgradQ[1] = grav*Q[0]*gradQ[idx_gradQ(0,3)];
  BgradQ[2] = grav*Q[0]*gradQ[idx_gradQ(1,3)];
  BgradQ[3] = 0.0;
}

void SWE::MySWESolver_ADERDG::mapDiscreteMaximumPrincipleObservables(double* observables,const int numberOfObservables, const double* const Q) const {
    //for (int i=0; i<NumberOfVariables; ++i) {
        observables[0] = Q[0];
    //}
}

bool SWE::MySWESolver_ADERDG::isPhysicallyAdmissible(
        const double* const solution,
        const double* const observablesMin,const double* const observablesMax,
        const int numberOfObservables,
        const tarch::la::Vector<DIMENSIONS,double>& center,
        const tarch::la::Vector<DIMENSIONS,double>& dx,
        const double t, const double dt) const {

    if (observablesMin[0] <= epsilon * 20){
        return false;
    }
    else {
        return true;
    }
}


