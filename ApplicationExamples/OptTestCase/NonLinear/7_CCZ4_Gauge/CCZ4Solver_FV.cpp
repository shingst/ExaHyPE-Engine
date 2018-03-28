#include "CCZ4Solver_FV.h"

#include "CCZ4Solver_FV_Variables.h"

#include <cstdlib> // exit()
#include "Fortran/PDE.h"
#include "InitialData/InitialData.h"
#include "Parameters/mexa.h"


tarch::logging::Log CCZ4::CCZ4Solver_FV::_log( "CCZ4::CCZ4Solver_FV" );

// enable nan tracker
//#include <fenv.h>

using namespace CCZ4Fortran;


void CCZ4::CCZ4Solver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
	//feenableexcept(FE_INVALID | FE_OVERFLOW);  // Enable all floating point exceptions but FE_INEXACT
	
	mexa::mexafile mf = mexa::fromSpecfile(constants.getAllAsOrderedMap(), constants.toString());
	
	// PDE Runtime Parameters
	GlobalPDEParameters::getInstance().setByParameters(mf);
	
	// Initial Data
	GlobalInitialData::getInstance().setByParameters(mf);
	
	// Boundary conditions for CCZ4: Currently always exact.
	//GlobalBoundaryConditions::getInstance().initializeDG(this).readParameters(mf("boundaries"));
}

void CCZ4::CCZ4Solver_FV::adjustSolution(const double* const x,const double t,const double dt, double* Q) {
	AdjustPointSolution(x,t,dt,Q);
}

void CCZ4::CCZ4Solver_FV::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
	PDE::eigenvalues(Q, dIndex, lambda);
}

void CCZ4::CCZ4Solver_FV::boundaryValues(const double* const x,  const double t,const double dt, const int faceIndex, const int d, const double* const stateInside, double* stateOutside) {
	// exact BC without integration
	InitialData(x, t, stateOutside);
}

// if you want to use fluxes...
/*
void CCZ4::CCZ4Solver_FV::flux(const double* const Q,double** F) {
  pdeflux_(F[0], F[1], DIMENSIONS==3 ? F[2] : nullptr, Q);
}
*/


void CCZ4::CCZ4Solver_FV::algebraicSource(const double* const Q,double* S) {
	PDE::algebraicSource(Q, S);
}

void  CCZ4::CCZ4Solver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
	PDE::nonConservativeProduct(Q, gradQ, BgradQ);
}

