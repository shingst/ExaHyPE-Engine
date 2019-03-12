#include "DummySolver.h"

#include "DummySolver_Variables.h"


tarch::logging::Log Dummy::DummySolver::_log( "Dummy::DummySolver" );

void Dummy::DummySolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
}


bool isUnphysical(const double* const Q) {
   // check wether state vector holds useful data or garbarage.
  const double eps = 1e-8;
  for(int i=0; i<Dummy::DummySolver::NumberOfVariables; i++) {
    if(Q[i] > 1 + eps || Q[i] < 1 - eps) return true;
  }
  return false;
}

void Dummy::DummySolver::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {
  if(tarch::la::equals(t, 0.0)) {
	for(int i=0; i<NumberOfVariables; i++) Q[i] = 1.0; // Sic, not 0.0
  }
}

void Dummy::DummySolver::eigenvalues(const double* const Q, const int dIndex, double* const lambda) {
  for(int i=0; i<NumberOfVariables; i++) lambda[i] = 1.0;
}

void Dummy::DummySolver::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* const stateOutside) {
	
  for(int i=0; i<NumberOfVariables; i++) stateOutside[i] = stateInside[i];
}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit


void Dummy::DummySolver::flux(const double* const Q,double** const F) {
  if(isUnphysical(Q)) return;

  for(int d=0; d<DIMENSIONS; d++)
  for(int i=0; i<NumberOfVariables; i++)
	F[d][i] = 0.0;
}

//You can either implement this method or modify fusedSource
void Dummy::DummySolver::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
  if(isUnphysical(Q)) return;
  for(int i=0; i<NumberOfVariables; i++) S[i] = 0.0;
}

void  Dummy::DummySolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  if(isUnphysical(Q)) return;
  for(int i=0; i<NumberOfVariables; i++) BgradQ[i] = 0.0;
}

