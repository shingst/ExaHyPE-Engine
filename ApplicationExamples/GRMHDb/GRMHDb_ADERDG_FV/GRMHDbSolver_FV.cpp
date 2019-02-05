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
#include "kernels/GaussLegendreQuadrature.h"

#include "GRMHDbSolver_ADERDG.h"


#include "tarch/multicore/BooleanSemaphore.h"

#include "tarch/multicore/Lock.h"

tarch::logging::Log GRMHDb::GRMHDbSolver_FV::_log( "GRMHDb::GRMHDbSolver_FV" );

void GRMHDb::GRMHDbSolver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
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

void GRMHDb::GRMHDbSolver_FV::adjustSolution(const double* const x, const double t, const double dt, double* Q) {
	// Dimensions             = 3
	// Number of variables    = 19 + #parameters

	// @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
	  const int nVar = GRMHDb::AbstractGRMHDbSolver_ADERDG::NumberOfVariables;
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
		static tarch::multicore::BooleanSemaphore initialDataSemaphore;
		tarch::multicore::Lock lock(initialDataSemaphore);
		/***************************************************/
		// everything in here is thread-safe w.r.t. the lock
		// call Fortran routines
		/***********************/

                initialdata_(x, &t, Q);

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


void GRMHDb::GRMHDbSolver_FV::referenceSolution(const double* const x, double t, double* Q) {
	const int nVar = GRMHDb::AbstractGRMHDbSolver_ADERDG::NumberOfVariables;
	int iErr = 0;
	double* Qcons;
	initialdata_(x, &t, Qcons);
	/*Q[0] = exp(-(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) / 8.0);
	Q[1] = sin(x[1])*sin(x[0]);
	Q[2] = sin(x[2]);
	for (int i = 3; i < nVar; i++) {
		Q[i] = cos(x[0]);
	}*/
	pdecons2prim_(Q, Qcons, &iErr);
}


void GRMHDb::GRMHDbSolver_FV::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
  // Dimensions             = 3
  // Number of variables    = 19 + #parameters
  
  // @todo Please implement/augment if required
  lambda[0] = 1.0;
  lambda[1] = 1.0;
  lambda[2] = 1.0;
  lambda[3] = 1.0;
  lambda[4] = 1.0;
  lambda[5] = 1.0;
  lambda[6] = 1.0;
  lambda[7] = 1.0;
  lambda[8] = 1.0;
  lambda[9] = 1.0;
  lambda[10] = 1.0;
  lambda[11] = 1.0;
  lambda[12] = 1.0;
  lambda[13] = 1.0;
  lambda[14] = 1.0;
  lambda[15] = 1.0;
  lambda[16] = 1.0;
  lambda[17] = 1.0;
  lambda[18] = 1.0;
  double nv[3] = {0.};
  nv[dIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}

void GRMHDb::GRMHDbSolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* stateOutside) {
  // Dimensions             = 3
  // Number of variables    = 19 + #parameters 
  // @todo Please implement/augment if required
	const int nVar = GRMHDb::AbstractGRMHDbSolver_FV::NumberOfVariables;
	double Qgp[nVar];
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

  double ti = t + 0.5 * dt;
  initialdata_(x, &ti, &Qgp[0]);
  for(int m=0; m < nVar; m++) {
        stateOutside[m] = Qgp[m];
  }


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


void GRMHDb::GRMHDbSolver_FV::flux(const double* const Q,double** F) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 19 + 0
  
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
  
  F[1][0] = 0.0;
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
  
    if(DIMENSIONS == 2){
        const int nVar = GRMHDb::AbstractGRMHDbSolver_FV::NumberOfVariables;
		double F_3[nVar];
		pdeflux_(F[0], F[1],F_3, Q);
	}else{
		pdeflux_(F[0], F[1],F[2], Q);
	}

}



void  GRMHDb::GRMHDbSolver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
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
  pdencp_(BgradQ, Q, gradQ);
}

