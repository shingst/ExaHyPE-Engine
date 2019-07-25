#include "GRMHDSolver_FV.h"

#include "GRMHDSolver_FV_Variables.h"

#include "AbstractGRMHDSolver_ADERDG.h"
#include "InitialData/InitialData.h"
#include "Fortran/PDE.h"

#include <stdio.h>
#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreBasis.h"

#include "exahype/disableOptimization.h" // bugs when limiting is on. whatevers

const double excision_radius = 1.0;

tarch::logging::Log GRMHD::GRMHDSolver_FV::_log("GRMHDSolver_FV");
constexpr int nVar = GRMHD::AbstractGRMHDSolver_FV::NumberOfVariables;

// enable nan tracker
#include <fenv.h>

void GRMHD::GRMHDSolver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
    prepare_id("TOVSolver");
}

void __attribute__((optimize("O0"))) initialData_FV(const double* const x,const double t,const double dt,double* const Q) {
  id->Interpolate(x, t, Q);
  //Debugging
  /*for(int i=0; i<nVar; i++) {
    if(!std::isfinite(Q[i])) {
      printf("NAN in i=%d at t=%f, x=[%f,%f,%f], Q[%d]=%f\n", i, t, x[0],x[1],x[2], i, Q[i]);
    }
  }*/

}

bool doneOnce=false;
void GRMHD::GRMHDSolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {
  using namespace tarch::la;
  // Do excision only in 3D.

  if(equals(t,0.0)) {
    if(!doneOnce) { printf("Calling ID on FV grid\n"); doneOnce=true; }
    initialData_FV(x, t, dt, Q);

    if( (x[1] > -1.e-10) && (x[1]  < 1.e-10)) {
      if( (x[2] > -1.e-10) && (x[2]  < 1.e-10)) {
//        std::cout << "adjusting FV at : " <<   x[0]  << " t = " << t <<  " Q[0] = " << Q[0] << std::endl;
        printf("adjusting FV at : %.5e , %.5e\n",x[0],Q[0]);
      }
    }
  }
}

void GRMHD::GRMHDSolver_FV::eigenvalues(const double* const Q, const int dIndex, double* const lambda) {
  double nv[3] = {0.};
  nv[dIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}


// Detection of unphysical states. In these cases, the user PDE functions shall never be called.
// We workaround by returning some kind of "neutral" values which go well with the scheme.
inline bool isAllZero(const double* const Q) {
	// TODO: Check only the metric, since we see 1e-14 values in Q despite
	// the state vector is not coming from the ID.
	for(int i=9; i<GRMHD::GRMHDSolver_FV::NumberOfVariables; i++){
        if(Q[i]>0)
            return false; 
    }
	return true;
}


void GRMHD::GRMHDSolver_FV::flux(const double* const Q, double** const F) {
  if(!isAllZero(Q))
    pdeflux_(F[0], F[1], (DIMENSIONS==3)?F[2]:nullptr, Q);
}

void GRMHD::GRMHDSolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateIn,
    double* const stateOut) {
  double Qgp[nVar];

  for(int m=0;m<nVar;m++) {
    stateOut[m] = stateIn[m];
  }

  double ti = t + 0.5 * dt;
  initialData_FV(x, ti, dt, Qgp);
  for(int m=0; m < nVar; m++) {
    stateOut[m] = Qgp[m];
  } 
}



void GRMHD::GRMHDSolver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  if(!isAllZero(Q))
    pdencp_(BgradQ, Q, gradQ);
}

/*
void GRMHD::GRMHDSolver_FV::coefficientMatrix(const double* const Q,const int d,double* const Bn) {
  // new FV scheme has no coefficient matrix
  static tarch::logging::Log _log("GRMHDSolver");
  logError("coefficientMatrix()", "Coefficient Matrix invoked");
  exit(-2);

  double nv[3] = {0.};
  nv[d] = 1;
  pdematrixb_(Bn, Q, nv);
}
*/
