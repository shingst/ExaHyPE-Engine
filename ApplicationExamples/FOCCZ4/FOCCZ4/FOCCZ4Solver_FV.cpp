#include "FOCCZ4Solver_FV.h"

#include "FOCCZ4Solver_FV_Variables.h"

#include "kernels/finitevolumes/musclhancock/c/musclhancock.h"
#include "InitialData.h"
#include "PDE.h"




tarch::logging::Log FOCCZ4::FOCCZ4Solver_FV::_log( "FOCCZ4::FOCCZ4Solver_FV" );

void FOCCZ4::FOCCZ4Solver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver_FV.h".
  
  // @todo Please implement/augment if required
}

void FOCCZ4::FOCCZ4Solver_FV::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver_FV.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
  
  if (tarch::la::equals(t,0.0)) {
    int md = exahype::solvers::Solver::getMaximumAdaptiveMeshDepth();
    double cms = exahype::solvers::Solver::getCoarsestMeshSize();
    const int order = FOCCZ4::FOCCZ4Solver_FV::PatchSize;
    std::fill_n(Q,96,0.0);

    double x_3[3];
    x_3[2]=0;
    std::copy_n(&x[0],DIMENSIONS,&x_3[0]);
    
    initialdata_(x_3, &t, Q);
  }
  for(int i = 0; i< 96 ; i++){
    assert(std::isfinite(Q[i]));
  }


}

void FOCCZ4::FOCCZ4Solver_FV::eigenvalues(const double* const Q, const int dIndex, double* const lambda) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver_FV.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
  
 
  double nv[3] = {0.};
  nv[dIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);

  for(int i = 0; i< 96 ; i++){
    assert(std::isfinite(lambda[i]));
  }
 
 
}

void FOCCZ4::FOCCZ4Solver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int direction,
    const double* const stateInside,
    double* const stateOutside) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver_FV.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.

    const int nVar = FOCCZ4::FOCCZ4Solver_FV::NumberOfVariables;

	double Qgp[nVar];

	double ti = t + 0.5 * dt;
	// Compute the outer state according to the initial condition
    double x_3[3];
    x_3[2]=0;
    std::copy_n(&x[0],DIMENSIONS,&x_3[0]);
    
    initialdata_(x_3, &t, Qgp);
	// Assign the proper outer state
	for(int m=0; m < nVar; m++) {
        stateOutside[m] = Qgp[m];
	}

}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit


void FOCCZ4::FOCCZ4Solver_FV::flux(const double* const Q,double** const F) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver_FV.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
  
  const int nVar = FOCCZ4::FOCCZ4Solver_FV::NumberOfVariables;
  if(DIMENSIONS == 2){
    double F_3[nVar];
    pdeflux_(F[0], F[1],F_3, Q);
  }else{
    pdeflux_(F[0], F[1],F[2], Q);
  }
  
}




//You can either implement this method or modify fusedSource
void FOCCZ4::FOCCZ4Solver_FV::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver_FV.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.

  pdesource_(S, Q);

}

void  FOCCZ4::FOCCZ4Solver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver_FV.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
  
 pdencp_(BgradQ, Q, gradQ);
 
}

