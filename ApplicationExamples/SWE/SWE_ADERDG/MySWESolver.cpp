// This file was generated by the ExaHyPE toolkit.
// It will NOT be regenerated or overwritten.
// Please adapt it to your own needs.
// 
// ========================
//   www.exahype.eu
// ========================

#include "MySWESolver.h"
#include "InitialData.h"
#include "MySWESolver_Variables.h"

#include "kernels/KernelUtils.h"

using namespace kernels;

const double grav = 9.81;


tarch::logging::Log SWE::MySWESolver::_log( "SWE::MySWESolver" );


void SWE::MySWESolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required

}

void SWE::MySWESolver::adjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 3 + 1
  if (tarch::la::equals(t,0.0)) {
      initialData(x, Q);
  }
}

void SWE::MySWESolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 3 + 1

    stateOut[0] = stateIn[0];
    stateOut[1] = stateIn[1];
    stateOut[2] = stateIn[2];


    fluxOut[0] = fluxIn[0];
    fluxOut[1] = fluxIn[1];
    fluxOut[2] = fluxIn[2];
}

exahype::solvers::Solver::RefinementControl SWE::MySWESolver::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void SWE::MySWESolver::eigenvalues(const double* const Q,const int d,double* lambda) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 3 + 1

  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const double c = std::sqrt(grav*vars.h());
  const double ih = 1./vars.h();
  double u_n = Q[d + 1] * ih;

  eigs.h() = u_n + c;
  eigs.hu() = u_n - c;
  eigs.hv() = u_n;
}


void SWE::MySWESolver::flux(const double* const Q,double** F) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 3 + 1

  ReadOnlyVariables vars(Q);

  const double ih = 1./vars.h();

  double* f = F[0];
  double* g = F[1];

  f[0] = vars.hu();
  f[1] = vars.hu()*vars.hu()*ih + 0.5*grav*vars.h()*vars.h();
  f[2] = vars.hu()*vars.hv()*ih;

  g[0] = vars.hv();
  g[1] = vars.hu()*vars.hv()*ih;
  g[2] = vars.hv()*vars.hv()*ih + 0.5*grav*vars.h()*vars.h();
  
}



void  SWE::MySWESolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
    idx2 idx_gradQ(DIMENSIONS,NumberOfVariables + NumberOfParameters);

    BgradQ[0] = 0.0;
    BgradQ[1] = grav*Q[0]*gradQ[idx_gradQ(0,3)];
    BgradQ[2] = grav*Q[0]*gradQ[idx_gradQ(1,3)];
}


