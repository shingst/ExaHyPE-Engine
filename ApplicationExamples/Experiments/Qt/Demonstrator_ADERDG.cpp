#include "Demonstrator_ADERDG.h"

#include "Demonstrator_ADERDG_Variables.h"


tarch::logging::Log Qt::Demonstrator_ADERDG::_log( "Qt::Demonstrator_ADERDG" );


void Qt::Demonstrator_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required

}

void Qt::Demonstrator_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 1 + 0
  const int nVar = Qt::Demonstrator_ADERDG::NumberOfVariables;
  // @todo Please implement/augment if required
  /*if (tarch::la::equals(t,0.0)) {
    for(int i=0; i < nVar; i++)
      Q[i] = sin(M_PI/1000*x[0])*sin(M_PI/1000*x[1])*sin(M_PI/1000*x[2]);
  }*/
}

void Qt::Demonstrator_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 1 + 0
  const int nVar = Qt::Demonstrator_ADERDG::NumberOfVariables;

  // @todo Please implement/augment if required
  for(int i=0; i < nVar; i++)
      stateOut[i] = 0.0;

  for(int i=0; i < nVar; i++)
      fluxOut[i] = 0.0;
}

exahype::solvers::Solver::RefinementControl Qt::Demonstrator_ADERDG::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


// Ackermann function to artificially increase cost of calculation
int ackermannFunction(int m, int n){
    if(m==0)
        return n+1;
    if(n==0)
        return ackermannFunction(m-1,1);
    return ackermannFunction(m-1,ackermannFunction(m,n-1));
}


//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void Qt::Demonstrator_ADERDG::eigenvalues(const double* const Q,const int d,double* lambda) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 1 + 0
  const int nVar = Qt::Demonstrator_ADERDG::NumberOfVariables;

  // @todo Please implement/augment if required
  for(int i=0; i < nVar; i++)
      lambda[i] = 1.0;
}


//You can either implement this method or modify fusedSource
void Qt::Demonstrator_ADERDG::algebraicSource(const double* const Q,double* S) {
  // @todo Please implement/augment if required
  const int nVar = Qt::Demonstrator_ADERDG::NumberOfVariables;
  for(int i=0; i < nVar; i++)
      S[i] = 1.0 + 1e-16*ackermannFunction(3,3);
}

void  Qt::Demonstrator_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  // @todo Please implement/augment if required
  const int nVar = Qt::Demonstrator_ADERDG::NumberOfVariables;
  for(int i=0; i < nVar; i++)
      BgradQ[i] = gradQ[i] + 1e-16*ackermannFunction(3,3);
}


