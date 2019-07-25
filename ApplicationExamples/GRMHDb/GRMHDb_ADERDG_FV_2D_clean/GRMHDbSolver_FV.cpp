#include "GRMHDbSolver_FV.h"

#include "GRMHDbSolver_FV_Variables.h"


// User defined calls
//#include "Tools.h"
#include "PDE.h"
#include "InitialData.h"
//#include "TECPLOTinterface.h"
#include "tarch/parallel/Node.h"
#include "tarch/la/MatrixVectorOperations.h"

#include <algorithm>

#include <cstring> // memset

#include <string>

#include <math.h>

#include "peano/utils/Loop.h"

#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreBasis.h"

#include "GRMHDbSolver_ADERDG.h"


#include "tarch/multicore/BooleanSemaphore.h"

#include "tarch/multicore/Lock.h"

tarch::logging::Log GRMHDb::GRMHDbSolver_FV::_log( "GRMHDb::GRMHDbSolver_FV" );

void GRMHDb::GRMHDbSolver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // Tip: You find documentation for this method in header file "GRMHDb::GRMHDbSolver_FV.h".
  
  // @todo Please implement/augment if required

    //const int order = GRMHDb::AbstractGRMHDbSolver_ADERDG::Order;
	//constexpr int basisSize = AbstractGRMHDbSolver_FV::PatchSize;
	//constexpr int Ghostlayers = AbstractGRMHDbSolver_FV::GhostLayerWidth;
    //int mpirank = tarch::parallel::Node::getInstance().getRank();
	//
	//
	///**************************************************************************/
	//static tarch::multicore::BooleanSemaphore initialDataSemaphore;
	//tarch::multicore::Lock lock(initialDataSemaphore);
	///***************************************************/
	//// everything in here is thread-safe w.r.t. the lock
	//// call Fortran routines
	///***********************/
	//
	////printf("\n******************************************************************");
	////printf("\n**************<<<  INIT TECPLOT    >>>****************************");
	////printf("\n******************************************************************");
    ////inittecplot_(&order,&order,&basisSize,&Ghostlayers);
	//printf("\n******************************************************************");
	//printf("\n**************<<<  INIT PDE SETUP  >>>****************************");
	//printf("\n******************************************************************");
    ////pdesetup_(&mpirank);
	//printf("\n******************************************************************");
	//printf("\n**************<<<       DONE       >>>****************************");
	//printf("\n******************************************************************");
    ////fflush(stdout);
	//
	//
	//
	///************/
	//lock.free();
	//// everything afterwards is not thread-safe anymore w.r.t. the lock
	///**************************************************************************/
}

void GRMHDb::GRMHDbSolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {
  // Tip: You find documentation for this method in header file "GRMHDb::GRMHDbSolver_FV.h".
  // Tip: See header file "GRMHDb::AbstractGRMHDbSolver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.

	// Dimensions             = 3
	// Number of variables    = 19 + #parameters

	// @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
	  //const int nVar = GRMHDb::AbstractGRMHDbSolver_ADERDG::NumberOfVariables;
	  
	  constexpr int numberOfVariables = AbstractGRMHDbSolver_FV::NumberOfVariables;
	  constexpr int numberOfParameters = AbstractGRMHDbSolver_FV::NumberOfParameters;
	  constexpr int numberOfData = numberOfVariables + numberOfParameters;
		Q[0] = 0.0;
		Q[1] = 0.0;
		Q[2] = 0.0;
		Q[3] = 0.0;
		Q[4] = 0.0;
		Q[5] = 0.0;
		Q[6] = 0.0;
		Q[7] = 0.0;
		Q[8] = 0.0;
		Q[9] = 0.0;
		Q[10] = 0.0;
		Q[11] = 0.0;
		Q[12] = 0.0;
		Q[13] = 0.0;
		Q[14] = 0.0;
		Q[15] = 0.0;
		Q[16] = 0.0;
		Q[17] = 0.0;
		Q[18] = 0.0;
		/**************************************************************************/
		static tarch::multicore::BooleanSemaphore initialDataSemaphoreFV;
		tarch::multicore::Lock lock(initialDataSemaphoreFV);
		/***************************************************/
		// everything in here is thread-safe w.r.t. the lock
		// call Fortran routines
		/***********************/
/*
        double x3D[3]={0.0};
        for(int i=0;i<DIMENSIONS;i++){
        x3D[i]=x[i];
        }
        //printf("x3d_FV:  %f,  %f",x3D[0],x3D[1]);
        */
        initialdata_(x, &t, Q);
       // printf("\nx,   rho:  %f,  %f",x3D[0],Q[1]);

		/************/
		lock.free();
		// everything afterwards is not thread-safe anymore w.r.t. the lock
		/**************************************************************************/

		/*Q[0] = exp(-(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) / 8.0);
		Q[1] = sin(x[1])*sin(x[0]);
		Q[2] = sin(x[2]);
		for (int i = 3; i < nVar; i++) {
			Q[i] = cos(x[0]);
		}*/
	}
}


void GRMHDb::GRMHDbSolver_FV::eigenvalues(const double* const Q, const int dIndex, double* const lambda) {
  // Tip: You find documentation for this method in header file "GRMHDb::GRMHDbSolver_FV.h".
  // Tip: See header file "GRMHDb::AbstractGRMHDbSolver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.

  // Dimensions             = 3
  // Number of variables    = 19 + #parameters

  //printf("\n******* EIGENVALUES FV *****************");
	constexpr int numberOfVariables = AbstractGRMHDbSolver_FV::NumberOfVariables;
	constexpr int numberOfParameters = AbstractGRMHDbSolver_FV::NumberOfParameters;
	constexpr int numberOfData = numberOfVariables + numberOfParameters;
	for (int i = 0; i < numberOfData; i++) {
		assertion(!std::isnan(Q[i]));
     }
  // @todo Please implement/augment if required
  lambda[0]  = 0.0;
  lambda[1]  = 0.0;
  lambda[2]  = 0.0;
  lambda[3]  = 0.0;
  lambda[4]  = 0.0;
  lambda[5]  = 0.0;
  lambda[6]  = 0.0;
  lambda[7]  = 0.0;
  lambda[8]  = 0.0;
  lambda[9]  = 0.0;
  lambda[10] = 0.0;
  lambda[11] = 0.0;
  lambda[12] = 0.0;
  lambda[13] = 0.0;
  lambda[14] = 0.0;
  lambda[15] = 0.0;
  lambda[16] = 0.0;
  lambda[17] = 0.0;
  lambda[18] = 0.0;
  //return;


  double nv[DIMENSIONS] = {0.};
  nv[dIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}


void GRMHDb::GRMHDbSolver_FV::boundaryValues(
	const double* const x, 
	const double t, const double dt, 
	const int faceIndex, 
	const int d,
	const double* const stateInside,
	double* const stateOutside) {
  // Tip: You find documentation for this method in header file "GRMHDb::GRMHDbSolver_FV.h".
  // Tip: See header file "GRMHDb::AbstractGRMHDbSolver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.

  // Dimensions             = 3
  // Number of variables    = 19 + #parameters 
  // @todo Please implement/augment if required
	constexpr int numberOfVariables = AbstractGRMHDbSolver_FV::NumberOfVariables;
	constexpr int numberOfParameters = AbstractGRMHDbSolver_FV::NumberOfParameters;
	constexpr int numberOfData = numberOfVariables + numberOfParameters;
	double Qtmp[numberOfData];
  stateOutside[0] = stateInside[0];
  stateOutside[1] = stateInside[1];
  stateOutside[2] = stateInside[2];
  stateOutside[3] = stateInside[3];
  stateOutside[4] = stateInside[4];
  stateOutside[5] = stateInside[5];
  stateOutside[6] = stateInside[6];
  stateOutside[7] = stateInside[7];
  stateOutside[8] = stateInside[8];
  stateOutside[9] = stateInside[9];
  stateOutside[10] = stateInside[10];
  stateOutside[11] = stateInside[11];
  stateOutside[12] = stateInside[12];
  stateOutside[13] = stateInside[13];
  stateOutside[14] = stateInside[14];
  stateOutside[15] = stateInside[15];
  stateOutside[16] = stateInside[16];
  stateOutside[17] = stateInside[17];
  stateOutside[18] = stateInside[18];

  //return;
  // d is one between 0,1,2
  // THIS IS FOR 1D Riemann problems. (inviscid reflection at the y boundaries)
  if (d==0) {
	  double ti = t + 0.5 * dt;
	  initialdata_(x, &ti, &Qtmp[0]);
	  for (int m = 0; m < numberOfData; m++) {
		  stateOutside[m] = Qtmp[m];
	  }
  }
  else {
	  int iErr;
	  double VtmpIn[numberOfData];
	  double VtmpOut[numberOfData];
	  for (int i = 0; i < numberOfData; i++) {
		  Qtmp[i] = stateInside[i];
	  }
	  pdecons2prim_(&VtmpIn[0], &Qtmp[0], &iErr);
	  for (int i = 0; i < numberOfData; i++) {
		  VtmpOut[i] = VtmpIn[i];
	  }
	  VtmpOut[1 + d] = -VtmpIn[1 + d];
	  pdeprim2cons_(&stateOutside[0], &VtmpOut[0]);
	  //stateOutside[1 + d] = -stateInside[1 + d];
  }

// THIS IS FOR ANALYTICAL BOUNDARY CONDITIONS:
 // double ti = t + 0.5 * dt;
 ///* 
 //       double x3D[3]={0.0};
 //       for(int i=0;i<DIMENSIONS;i++){
 //       x3D[i]=x[i];
 //       }*/

 // initialdata_(x, &ti, &Qgp[0]);
 // for(int m=0; m < numberOfData; m++) {
 //       stateOutside[m] = Qgp[m];
 // }


  /*stateOutside[0] = exp(-(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) / 8.0);
  stateOutside[1] = sin(x[1])*sin(x[0]);
  stateOutside[2] = sin(x[2]);
  for (int i = 3; i < nVar; i++) {
	  stateOutside[i] = cos(x[0]);
  }*/

}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit


void GRMHDb::GRMHDbSolver_FV::flux(const double* const Q,double** const F) {
  // Tip: You find documentation for this method in header file "GRMHDb::GRMHDbSolver_FV.h".
  // Tip: See header file "GRMHDb::AbstractGRMHDbSolver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.

  // Dimensions                        = 3
  // Number of variables + parameters  = 19 + 0
	constexpr int numberOfVariables = AbstractGRMHDbSolver_FV::NumberOfVariables;
	constexpr int numberOfParameters = AbstractGRMHDbSolver_FV::NumberOfParameters;
	constexpr int numberOfData = numberOfVariables + numberOfParameters;
  
  
       // printf("\n******* FLUXES *****************");

  // @todo Please implement/augment if required
  F[0][0] = 0.0;
  F[0][1] = 0.0;
  F[0][2] = 0.0;
  F[0][3] = 0.0;
  F[0][4] = 0.0;
  F[0][5] = 0.0;
  F[0][6] = 0.0;
  F[0][7] = 0.0;
  F[0][8] = 0.0;
  F[0][9] = 0.0;
  F[0][10] = 0.0;
  F[0][11] = 0.0;
  F[0][12] = 0.0;
  F[0][13] = 0.0;
  F[0][14] = 0.0;
  F[0][15] = 0.0;
  F[0][16] = 0.0;
  F[0][17] = 0.0;
  F[0][18] = 0.0;

  F[1][0] =  0.0;
  F[1][1] = 0.0;
  F[1][2] = 0.0;
  F[1][3] = 0.0;
  F[1][4] = 0.0;
  F[1][5] = 0.0;
  F[1][6] = 0.0;
  F[1][7] = 0.0;
  F[1][8] = 0.0;
  F[1][9] = 0.0;
  F[1][10] = 0.0;
  F[1][11] = 0.0;
  F[1][12] = 0.0;
  F[1][13] = 0.0;
  F[1][14] = 0.0;
  F[1][15] = 0.0;
  F[1][16] = 0.0;
  F[1][17] = 0.0;
  F[1][18] = 0.0;

#ifdef Dim3
  F[2][0] = 0.0;
  F[2][1] = 0.0;
  F[2][2] = 0.0;
  F[2][3] = 0.0;
  F[2][4] = 0.0;
  F[2][5] = 0.0;
  F[2][6] = 0.0;
  F[2][7] = 0.0;
  F[2][8] = 0.0;
  F[2][9] = 0.0;
  F[2][10] = 0.0;
  F[2][11] = 0.0;
  F[2][12] = 0.0;
  F[2][13] = 0.0;
  F[2][14] = 0.0;
  F[2][15] = 0.0;
  F[2][16] = 0.0;
  F[2][17] = 0.0;
  F[2][18] = 0.0; 
#endif
   


    if(DIMENSIONS == 2){
        //const int nVar = GRMHDb::AbstractGRMHDbSolver_FV::NumberOfVariables;
		double F_3[numberOfData];
		pdeflux_(F[0], F[1],F_3, Q);
	}else{
		pdeflux_(F[0], F[1],F[2], Q);
	}


  //constexpr int nFF = numberOfVariables * DIMENSIONS;
  //double FF[nFF];
  //int count;
  ////count = 0;
  ////for (int dloc = 0; dloc < DIMENSIONS; dloc++) {
	 //// for (int i = 0; i < numberOfVariables; i++) {
		////  FF[count]= F[dloc][i];
		////  printf("FF(count) is  %d,    %f", count, FF[count]);
		////  count++;
	 //// }
  ////}
  //pdeflux_(&FF[0], Q);


  //count = 0;
  //for (int dloc = 0; dloc < DIMENSIONS; dloc++) {
  //for (int i = 0; i < numberOfVariables; i++) {
	 // F[dloc][i] = FF[count];
	 // count++;
	 // }
  //}

  //F[1][0] = 0.0;
  //F[1][1] = 0.0;
  //F[1][2] = 0.0;
  //F[1][3] = 0.0;
  //F[1][4] = 0.0;
  //F[1][5] = 0.0;
  //F[1][6] = 0.0;
  //F[1][7] = 0.0;
  //F[1][8] = 0.0;
  //F[1][9] = 0.0;
  //F[1][10] = 0.0;
  //F[1][11] = 0.0;
  //F[1][12] = 0.0;
  //F[1][13] = 0.0;
  //F[1][14] = 0.0;
  //F[1][15] = 0.0;
  //F[1][16] = 0.0;
  //F[1][17] = 0.0;
  //F[1][18] = 0.0;




 //   if(DIMENSIONS == 2){
 //       //const int nVar = GRMHDb::AbstractGRMHDbSolver_FV::NumberOfVariables;
	//	double F_3[numberOfData];
	//	pdeflux_(F[0], F[1],F_3, Q);
	//}else{
	//	pdeflux_(F[0], F[1],F[2], Q);
	//}
  
}




void  GRMHDb::GRMHDbSolver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  // Tip: You find documentation for this method in header file "GRMHDb::GRMHDbSolver_FV.h".
  // Tip: See header file "GRMHDb::AbstractGRMHDbSolver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.

       // printf("\n*******nonConservativeProduct *****************");
  
  // @todo Please implement/augment if required
  BgradQ[0] = 0.0;
  BgradQ[1] = 0.0;
  BgradQ[2] = 0.0;
  BgradQ[3] = 0.0;
  BgradQ[4] = 0.0;
  BgradQ[5] = 0.0;
  BgradQ[6] = 0.0;
  BgradQ[7] = 0.0;
  BgradQ[8] = 0.0;
  BgradQ[9] = 0.0;
  BgradQ[10] = 0.0;
  BgradQ[11] = 0.0;
  BgradQ[12] = 0.0;
  BgradQ[13] = 0.0;
  BgradQ[14] = 0.0;
  BgradQ[15] = 0.0;
  BgradQ[16] = 0.0;
  BgradQ[17] = 0.0;
  BgradQ[18] = 0.0;
  
  //return;
  //
  ////printf("\n QUI QUO QUA");
  //const int numberOfVariables = AbstractGRMHDbSolver_ADERDG::NumberOfVariables;
  //constexpr int tot3Dvar = 3*numberOfVariables;
  ////printf("\nTOT 3D var = %d",tot3Dvar);
  //double gradQ3D[tot3Dvar]={0.0};
  //int count;
  //count = 0;
  //for(int i=0;i<DIMENSIONS;i++){
  //for(int m=0;m<numberOfVariables;m++){
  //      //printf("\n NCP count=%d    inizio",count);
  //      //fflush(stdout);
  //      gradQ3D[count]=1.0;
  //      //printf("\n NCP count=%d  mezzo",count); 
  //      gradQ3D[count]=1.0;
  //      //fflush(stdout);

  //      double myself = gradQ[count];      
  //      //printf("\n NCP count=%d  fine ",count);
  //      //fflush(stdout);
  //      gradQ3D[count]=gradQ[count];
  //      count++;
  //  }
  //  }
   //printf("\n ARRIVATO QUI");   
  pdencp_(BgradQ, Q, gradQ);


      //  printf("\n******* FV *****************");

}


void GRMHDb::GRMHDbSolver_FV::referenceSolution(const double* const x, double t, double* Q) {
	constexpr int numberOfVariables = AbstractGRMHDbSolver_FV::NumberOfVariables;
	constexpr int numberOfParameters = AbstractGRMHDbSolver_FV::NumberOfParameters;
	constexpr int numberOfData = numberOfVariables + numberOfParameters;
	
	
       // printf("\n*******referenceSolution*****************");

	
	int iErr;
	double Qcons[numberOfData];
	iErr = 0;
        /*
        double x3D[3]={0.0};
        for(int i=0;i<DIMENSIONS;i++){
        x3D[i]=x[i];
        }*/
	//printf("\nSONO QUI IN REFERENCE SOLUTION");
	initialdata_(x, &t, &Qcons[0]);

	//// test:
	//Q[0] = exp(-(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) / 8.0);
	//	Q[1] = sin(x[1])*sin(x[0]);
	//	Q[2] = sin(x[2]);
	//for (int i = 3; i < nVar; i++) {
	//	Q[i] = cos(x[0]);
	//}
	pdecons2prim_(Q, &Qcons[0], &iErr);
}


#include "kernels/finitevolumes/riemannsolvers/c/riemannsolvers.h"
double GRMHDb::GRMHDbSolver_FV::riemannSolver(double* fL, double *fR, const double* qL, const double* qR, const double* gradQL, const double* gradQR, const double* cellSize, int direction) {
//double GRMHDb::GRMHDbSolver_FV::riemannSolver(double* const fL, double* const fR, const double* const qL, const double* const qR, const double* gradQL, const double* gradQR, int direction) {
	//// Default FV Riemann Solver
        // printf("\n*******riemannSolver(*****************");

    //return kernels::finitevolumes::riemannsolvers::c::rusanov<true, true, false, GRMHDbSolver_FV>(*static_cast<GRMHDbSolver_FV*>(this), fL,fR,qL,qR,gradQL, gradQR, cellSize, direction);
	constexpr int numberOfVariables = AbstractGRMHDbSolver_FV::NumberOfVariables;

	//printf("SONO QUI IN riemannSolver");
	/* HLLEM */
	
	//const int numberOfVariables = GRMHDb::AbstractGRMHDbSolver_FV::NumberOfVariables;
	const int numberOfParameters = GRMHDb::AbstractGRMHDbSolver_FV::NumberOfParameters;
	const int numberOfData = numberOfVariables + numberOfParameters;
	const int order = 0;  // for finite volume we use one single d.o.f., i.e. the cell average.
	const int basisSize = order + 1;
	// Compute the average variables and parameters from the left and the right
	double QavL[numberOfData] = { 0.0 }; // ~(numberOfVariables+numberOfParameters)
	double QavR[numberOfData] = { 0.0 }; // ~(numberOfVariables+numberOfParameters)

										 // std::cout << "opened ---------------------"<< std::endl;
         
        // printf("\n******* RIEMANN SOLVER FV*****************");

	kernels::idx2 idx_QLR(basisSize, numberOfData);
	for (int j = 0; j < basisSize; j++) {
		const double weight = kernels::legendre::weights[order][j];

		for (int k = 0; k < numberOfData; k++) {
			QavL[k] += weight * qL[idx_QLR(j, k)];
			QavR[k] += weight * qR[idx_QLR(j, k)];
		}
	}
	
        // printf("\n***DONE*****");

	double lambda = 2.0;
	hllemfluxfv_(fL, fR, qL, qR, QavL, QavR, &direction);

	/* OSHER */
	//double lambda = kernels::finitevolumes::riemannsolvers::c::generalisedOsherSolomon<false, true, false, 3, EulerSolver_FV>(*static_cast<EulerSolver_FV*>(this), fL, fR, qL, qR, direction);
	/* RUSANOV */
	//double lambda = kernels::finitevolumes::riemannsolvers::c::rusanov<true, true, false, GRMHDbSolver_FV>(*static_cast<GRMHDbSolver_FV*>(this), fL, fR, qL, qR, gradQL, gradQR, cellSize, direction);
	// avoid spurious numerical diffusion (ony for Cowling approximation)
	for (int m = 9; m < numberOfVariables; m++) {
		fL[m] = 0.0;
		fR[m] = 0.0;
	}
	return lambda; 
}


