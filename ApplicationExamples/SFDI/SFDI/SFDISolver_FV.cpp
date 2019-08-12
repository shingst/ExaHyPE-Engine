#include "SFDISolver_FV.h"

#include "SFDISolver_FV_Variables.h"

#include "kernels/finitevolumes/musclhancock/c/musclhancock.h"
#include "InitialData.h"
#include "PDE.h"


tarch::logging::Log SFDI::SFDISolver_FV::_log( "SFDI::SFDISolver_FV" );

void SFDI::SFDISolver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
}

void SFDI::SFDISolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {

  if (tarch::la::equals(t,0.0)) {
    int md = exahype::solvers::Solver::getMaximumAdaptiveMeshDepth();
    double cms = exahype::solvers::Solver::getCoarsestMeshSize();
    const int order = PatchSize;
    std::fill_n(Q,NumberOfVariables,0.0);
    double x_3[3];
    x_3[2]=0;
    std::copy_n(&x[0],DIMENSIONS,&x_3[0]);
    initialdata_(x_3, &t, Q);
  }
  for(int i = 0; i< NumberOfVariables ; i++){
    assert(std::isfinite(Q[i]));
  }
}

void SFDI::SFDISolver_FV::eigenvalues(const double* const Q, const int dIndex, double* const lambda) {
  double nv[3] = {0.};
  nv[dIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
  for(int i = 0; i< NumberOfVariables ; i++){
    assertion2(std::isfinite(lambda[i]),i,lambda[i]);
  }
}

void SFDI::SFDISolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int direction,
    const double* const stateInside,
    double* const stateOutside) {
  const int nVar = NumberOfVariables;
  std::copy_n(stateInside,nVar,stateOutside);
}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit


void SFDI::SFDISolver_FV::flux(const double* const Q,double** const F) {
  const int nVar = SFDI::SFDISolver_FV::NumberOfVariables;
  if(DIMENSIONS == 2){
    double F_3[nVar];
    pdeflux_(F[0], F[1],F_3, Q);
  }else{
    pdeflux_(F[0], F[1],F[2], Q);
  }
  for(int d = 0; d< DIMENSIONS ; d++){
    for(int i = 0; i< NumberOfVariables ; i++){
      assertion3(std::isfinite(F[d][i]),d,i,F[d][i]);
    }
  }
}




//You can either implement this method or modify fusedSource
void SFDI::SFDISolver_FV::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
  pdesource_(S, Q);
  for(int i = 0; i< NumberOfVariables ; i++){
    assertion2(std::isfinite(S[i]),i,S[i]);
  }
}

void  SFDI::SFDISolver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  pdencp_(BgradQ, Q, gradQ);
  for(int i = 0; i< NumberOfVariables ; i++){
    assertion2(std::isfinite(BgradQ[i]),i,BgradQ[i]);
  }
}

double SFDI::SFDISolver_FV::riemannSolver(double* fL, double *fR, const double* qL, const double* qR, const double* gradQL, const double* gradQR, const double* cellSize, int direction) {


    //return kernels::finitevolumes::riemannsolvers::c::rusanov<true, true, false, GRMHDbSolver_FV>(*static_cast<GRMHDbSolver_FV*>(this), fL,fR,qL,qR,gradQL, gradQR, cellSize, direction);
	constexpr int numberOfVariables = AbstractSFDISolver_FV::NumberOfVariables;

	//printf("SONO QUI IN riemannSolver");
	/* HLLEM */
	
	//const int numberOfVariables = GRMHDb::AbstractGRMHDbSolver_FV::NumberOfVariables;
	const int numberOfParameters = SFDI::AbstractSFDISolver_FV::NumberOfParameters;
	const int numberOfData = numberOfVariables + numberOfParameters;
	const int order = 0;  // for finite volume we use one single d.o.f., i.e. the cell average.
	const int basisSize = order + 1;
	// Compute the average variables and parameters from the left and the right
	double QavL[numberOfData] = { 0.0 }; // ~(numberOfVariables+numberOfParameters)
	double QavR[numberOfData] = { 0.0 }; // ~(numberOfVariables+numberOfParameters)
    double lambda = kernels::finitevolumes::riemannsolvers::c::rusanov<true, true, false, SFDISolver_FV>(*static_cast<SFDISolver_FV*>(this), fL, fR, qL, qR, gradQL, gradQR, cellSize, direction);
    //double lambda;
	//hllemfluxfv_(&lambda, fL, fR, qL, qR, &direction);
	
	
	//std::cout << lambda << std::endl;
	//double1 lambda = 10.0;
	return lambda; 
}
