
#include <algorithm>
#include <cstring>

#include "kernels/Euler_EulerSolver_ADERDG/Kernels.h"
#include "kernels/Euler_EulerSolver_ADERDG/Quadrature.h"

#include "kernels/Euler_EulerSolver_ADERDG/gemmsCPP.h"

#include "EulerSolver_ADERDG.h"

//Fortran (Limiter.f90): GetSubcellData
void Euler::EulerSolver_ADERDG_kernels::aderdg::projectOnFVLimiterSpace(const double* const luh, double* const lim) {
  
  //compact projection without ghostlayer
  // x
  double tmpX[224] __attribute__((aligned(ALIGNMENT)));
  for (int zy = 0; zy < 4; zy++) {
    // will overwrite tmpX, no need to set to 0
    #ifdef USE_IPO
    #pragma forceinline
    #endif
    gemm_5_7_4_dg2fv_x(luh+zy*20, dg2fv, tmpX+zy*56);
  }
  
  // y
  for (int x = 0; x < 7; x++) {
    // will overwrite lim, no need to set to 0
    #ifdef USE_IPO
    #pragma forceinline
    #endif
    gemm_5_7_4_dg2fv_y(tmpX+x*8, dg2fv, lim+x*5);
  }
    
  //space out for ghostlayer (in place)
  for (int z = 0; z >= 0; z--) {
    for (int y = 6; y >= 0; y--) {
      std::copy_n(lim+(z*7+y)*35, 35, lim+(((z+0)*7+y+0)*7+0)*5); // no overlap since we have ghostlayers at the beginning
      std::memset(lim+(z*7+y)*35, 0, sizeof(double)*35); //delete the memory block that was moved to ensure ghostlayers are set at 0.
    }
  }
  
}

//Fortran (Limiter.f90): PutSubcellData
void Euler::EulerSolver_ADERDG_kernels::aderdg::projectOnDGSpace(const double* const lim, double* const luh) {
  
  // x
  // ignore and remove ghostlayers
  double tmpX[224] __attribute__((aligned(ALIGNMENT)));
  for (int z = 0; z < 1; z++) {
    for (int y = 0; y < 7; y++) {
      // will overwrite tmpX, no need to set to 0
      #ifdef USE_IPO
      #pragma forceinline
      #endif
      gemm_5_4_7_fv2dg_x(lim+((z+0)*7+y+0)*35+0, fv2dg, tmpX+(z*7+y)*32);
    }
  }
  
  // y
  for (int x = 0; x < 4; x++) {
    // will overwrite luh, no need to set to 0
    #ifdef USE_IPO
    #pragma forceinline
    #endif
    gemm_5_4_7_fv2dg_y(tmpX+x*8, fv2dg, luh+x*5);
  }
    
}


bool Euler::EulerSolver_ADERDG_kernels::aderdg::discreteMaximumPrincipleAndMinAndMaxSearch(
    const double* const luh,
    const exahype::solvers::ADERDGSolver* solver,
    const double relaxationParameter,const double differenceScaling,
    double* boundaryMinPerObservable, double* boundaryMaxPerObservable
) {

  double localMinPerObservable[0] __attribute__((aligned(ALIGNMENT)));
  double localMaxPerObservable[0] __attribute__((aligned(ALIGNMENT)));

  // 1. Determine the new cell-local -minimum and maximummin and max
  findCellLocalMinAndMax(luh,solver,localMinPerObservable,localMaxPerObservable);
  
  // 2. Compare to the boundary minimum and maximum
  bool discreteMaximumPrincipleSatisfied=true;
  for(int v = 0; v < 0; v++) {
    double boundaryMin = boundaryMinPerObservable[v];
    for (int i=1; i<4; i++) {
      boundaryMin = std::min( boundaryMin, boundaryMinPerObservable[i*0+v] );
    }
    double boundaryMax = boundaryMaxPerObservable[v];
    for (int i=1; i<4; i++) {
      boundaryMax = std::max( boundaryMax, boundaryMaxPerObservable[i*0+v] );
    }
    
    const double scaledRelaxationParameter =
        solver->getDiscreteMaximumPrincipleRelaxationParameter(
            relaxationParameter, v,
            localMinPerObservable[v],localMaxPerObservable[v],
            boundaryMin,boundaryMax);
    double scaledDifference = (boundaryMax - boundaryMin) * differenceScaling;
    scaledDifference = std::max( scaledDifference, scaledRelaxationParameter );

    if((localMinPerObservable[v] < (boundaryMin - scaledDifference)) ||
       (localMaxPerObservable[v] > (boundaryMax + scaledDifference))) {
      discreteMaximumPrincipleSatisfied=false;
    }

    // check for nans and infinity values
    discreteMaximumPrincipleSatisfied &= std::isfinite(localMinPerObservable[v]) &&
                                         std::isfinite(localMaxPerObservable[v]);

    // TODO(Dominic): A little hacky

    // We have the new min and max directly available now and
    // overwrite the block for face 0 with it
    boundaryMinPerObservable[v] = localMinPerObservable[v];
    boundaryMaxPerObservable[v] = localMaxPerObservable[v];

    // In the block for face 1, we write the boundary min and max
    boundaryMinPerObservable[v+0] = boundaryMin;
    boundaryMaxPerObservable[v+0] = boundaryMax;
  }

  return discreteMaximumPrincipleSatisfied;
}


void Euler::EulerSolver_ADERDG_kernels::aderdg::findCellLocalLimiterMinAndMax(
    const double* const lim,
    const exahype::solvers::ADERDGSolver* solver,
    double* const localMinPerVariables, double* const localMaxPerVariables
) {
  
  std::fill_n(localMinPerVariables,0,+std::numeric_limits<double>::infinity());
  std::fill_n(localMaxPerVariables,0,-std::numeric_limits<double>::infinity());

  double observables[0] __attribute__((aligned(ALIGNMENT)));
  
  for (int z = 0; z < 1; z++) { // skip the last element
    for (int y = 0; y < 7; y++) {
      for (int x = 0; x < 7; x++) {
        solver->mapDiscreteMaximumPrincipleObservables(observables,lim+(((z*7+y)*7+x)*5));

        for (int v = 0; v < 0; v++) {
          localMinPerVariables[v] = std::min ( localMinPerVariables[v], observables[v] );
          localMaxPerVariables[v] = std::max ( localMaxPerVariables[v], observables[v] );
        }
      }
    }
  }
  
}

/**
 * localMinPerVariables, localMaxPerVariables are double[numberOfVariables]
 */
void Euler::EulerSolver_ADERDG_kernels::aderdg::findCellLocalMinAndMax(
    const double* const luh,
    const exahype::solvers::ADERDGSolver* solver,
    double* const localMinPerVariables, double* const localMaxPerVariables
) {
  
  std::fill_n(localMinPerVariables,0,+std::numeric_limits<double>::infinity());
  std::fill_n(localMaxPerVariables,0,-std::numeric_limits<double>::infinity());

  double observables[0] __attribute__((aligned(ALIGNMENT)));

  for (int zyx = 0; zyx < 16; zyx++) {
    solver->mapDiscreteMaximumPrincipleObservables(observables,luh+(zyx*5));

    for (int v = 0; v < 0; v++) {
      localMinPerVariables[v] = std::min ( localMinPerVariables[v], observables[v] );
      localMaxPerVariables[v] = std::max ( localMaxPerVariables[v], observables[v] );
    }
  }
  compareWithADERDGSolutionAtGaussLobattoNodes(luh, solver, localMinPerVariables, localMaxPerVariables);
  compareWithADERDGSolutionAtFVSubcellCenters (luh, solver, localMinPerVariables, localMaxPerVariables);

}

//*************************
//*** Private functions ***
//*************************
/**
 * Auxilliary function to findMinMax
 * Project to GaussLobatto and modify the min/max if required
 */
void Euler::EulerSolver_ADERDG_kernels::aderdg::compareWithADERDGSolutionAtGaussLobattoNodes(
    const double* const luh,
    const exahype::solvers::ADERDGSolver* solver,
    double* min, 
    double* max
) {

  // x
  double tmpX[128] __attribute__((aligned(ALIGNMENT)));
  for (int zy = 0; zy < 4; zy++) {
    // will overwrite tmpX, no need to set to 0
    #ifdef USE_IPO
    #pragma forceinline
    #endif
    gemm_5_4_4_uh2lob_x(luh+zy*20, uh2lob, tmpX+zy*32);
  }
  
  // y
  double observables[0] __attribute__((aligned(ALIGNMENT)));
  double lob[32] __attribute__((aligned(ALIGNMENT)));
 // constant x slice of projected solution
  for (int x = 0; x < 4; x++) {
    // will overwrite lop, no need to set to 0
    #ifdef USE_IPO
    #pragma forceinline
    #endif
    gemm_8_4_4_uh2lob_y_slice(tmpX+x*8, uh2lob, lob);
 
    for(int p = 0; p < 4; p++) {
      solver->mapDiscreteMaximumPrincipleObservables(observables, lob+p*8);
      for (int n = 0; n < 0; n++) {
        min[n] = std::min( min[n], observables[n] );
        max[n] = std::max( max[n], observables[n] );
      }
    }
  }
    
}

/**
 * Auxilliary function to findMinMax
 * Project onto FV subcell nodes and modify the min/max if required
 */
void Euler::EulerSolver_ADERDG_kernels::aderdg::compareWithADERDGSolutionAtFVSubcellCenters(
    const double* const luh,
    const exahype::solvers::ADERDGSolver* solver,
    double* min, 
    double* max
) {
  
  // x
  double tmpX[224] __attribute__((aligned(ALIGNMENT)));
  for (int zy = 0; zy < 4; zy++) {
    // will overwrite tmpX, no need to set to 0
    #ifdef USE_IPO
    #pragma forceinline
    #endif
    gemm_5_7_4_dg2fv_x(luh+zy*20, dg2fv, tmpX+zy*56);
  }
  
  // y
  double observables[0] __attribute__((aligned(ALIGNMENT)));
  double lim[56] __attribute__((aligned(ALIGNMENT)));
 // constant x slice of projected solution
  for (int x = 0; x < 7; x++) {
    // will overwrite lim, no need to set to 0
    #ifdef USE_IPO
    #pragma forceinline
    #endif
    gemm_8_7_4_dg2fv_y_slice(tmpX+x*8, dg2fv, lim);
    for(int p = 0; p < 7; p++) {
      solver->mapDiscreteMaximumPrincipleObservables(observables, lim+p*8);
      for (int n = 0; n < 0; n++) {
        min[n] = std::min( min[n], observables[n] );
        max[n] = std::max( max[n], observables[n] );
      }
    }
  }
    
}