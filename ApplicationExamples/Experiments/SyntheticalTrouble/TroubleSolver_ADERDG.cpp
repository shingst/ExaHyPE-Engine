#include "TroubleSolver_ADERDG.h"

#include "TroubleSolver_ADERDG_Variables.h"


tarch::logging::Log Synthetical::TroubleSolver_ADERDG::_log( "Synthetical::TroubleSolver_ADERDG" );

constexpr int nVar = Synthetical::AbstractTroubleSolver_ADERDG::NumberOfVariables;
constexpr int nDim = DIMENSIONS;

// a unique magic number
constexpr double magicProof = 12345;

#define NVARS(i) for(int i=0; i<nVar; i++) 
#define DFOR(i)  for(int i=0; i<nDim; i++)

// Dimensions                        = 2
// Number of variables + parameters  = 10 + 0


void Synthetical::TroubleSolver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

void InitialData(const double* const x,const double t,double* const Q) {
  using namespace Synthetical::TroubleSolver_ADERDG_Variables::shortcuts;
  // initialize everything with rising numbers
  NVARS(i) Q[i] = i;

  // store a magic number to detect whether a position is initialized
  Q[proof] = magicProof;

  // store the coordinate
  DFOR(i) Q[pos+i] = x[i];
  if(nDim<3) Q[pos+2] = 0;
}

bool correctState(const double* const Q) {
  using namespace Synthetical::TroubleSolver_ADERDG_Variables::shortcuts;
  using namespace tarch::la;
  return equals(Q[proof], magicProof);
}

void checkCorrectState(const double* const Q) {
  if(!correctState(Q)) {
     printf("Errnous State:\n");
     NVARS(i) printf("Q[%d] = %e\n", i, Q[i]);
     std::abort();
  }
}

void Synthetical::TroubleSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
  if (tarch::la::equals(t,0.0)) InitialData(x,t,Q);
}

void Synthetical::TroubleSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {

  // trivial BC.

  InitialData(x,t,stateOut);
  NVARS(i) fluxOut[i] = 0;
  /*
  double Fs[nDim][nVar], *F[nDim];
  for(int dd=0; dd<nDim; dd++) F[dd] = Fs[dd];
  F[normalNonZero] = fluxOut;
  flux(stateOut, F);
  */
}

exahype::solvers::Solver::RefinementControl Synthetical::TroubleSolver_ADERDG::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void Synthetical::TroubleSolver_ADERDG::eigenvalues(const double* const Q,const int d,double* const lambda) {
  NVARS(i) lambda[i] = 1;
}


void Synthetical::TroubleSolver_ADERDG::flux(const double* const Q,double** const F) {
  checkCorrectState(Q);

  DFOR(d) NVARS(i) F[d][i] = 0;
}


//You can either implement this method or modify fusedSource
void Synthetical::TroubleSolver_ADERDG::algebraicSource(const double* const Q,double* const S) {
  checkCorrectState(Q);
  NVARS(i) S[i] = 0;
}

void  Synthetical::TroubleSolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  checkCorrectState(Q);
  NVARS(i) BgradQ[i] = 0;
}

