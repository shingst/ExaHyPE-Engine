
#ifndef _EXAHYPE_SWE_MYSWESOLVER_ADERDG_KERNELS_ADERDG_KERNELS_H_
#define _EXAHYPE_SWE_MYSWESOLVER_ADERDG_KERNELS_ADERDG_KERNELS_H_

#ifndef __INTEL_COMPILER
#include <mm_malloc.h>
#endif
#include <vector>

//forward declaration of the user solver
namespace SWE {
  class MySWESolver_ADERDG;
}

//forward declaration of the ADERDGSolver for limiter functions
namespace exahype {
namespace solvers {
 class ADERDGSolver;
}
}


#define NDEBUG
namespace SWE {
namespace MySWESolver_ADERDG_kernels {
namespace aderdg {
  int fusedSpaceTimePredictorVolumeIntegral(
    SWE::MySWESolver_ADERDG& solver,
    double* restrict lduh,
    double* restrict lQhbnd,
    double* restrict lFhbnd,
    double* restrict lQi,
    double* restrict rhs,
    double* restrict lFi,
    double* restrict lSi,   // for NCP or Source
    double* restrict lQhi,
    double* restrict lFhi,
    double* restrict lShi,  // for NCP or Source
    double* restrict gradQ, // for NCP or Source
    const double* const restrict luh,
    const double inverseDx, //Assume dx[0] == dx[1] == dx[2]
    const double dt
  );

  void solutionUpdate( 
    double* restrict luh,
    const double* restrict const luhOld, 
    const double* restrict const lduh, 
    const double dt
  );
  
  void surfaceIntegral( 
    double* restrict lduh, 
    const double* restrict const lFhbnd, 
    const double inverseDx //Assume dx[0] == dx[1] == dx[2]
  );
  
  void faceIntegral(
    double *lduh, 
    const double *const lFhbnd,
    const int direction, 
    const int orientation,
    const double inverseDxDirection
  );

  void solutionAdjustment(
    SWE::MySWESolver_ADERDG& solver,
    double* luh,
    const double* const center,
    const double dx, //Assume dx[0] == dx[1] == dx[2]
    const double t,
    const double dt
  );

  void riemannSolver( 
    SWE::MySWESolver_ADERDG& solver,
    double* restrict FL,
    double* restrict FR,
    const double* restrict const QL,
    const double* restrict const QR,
    const double t,
    const double dt,
    const int direction
  );

  double stableTimeStepSize(
    SWE::MySWESolver_ADERDG& solver,
    const double* restrict const luh,
    const double inverseDx //Assume dx[0] == dx[1] == dx[2]
  );

  void boundaryConditions(
    SWE::MySWESolver_ADERDG& solver,
    double* fluxOut, 
    double* stateOut, 
    const double* const fluxIn, 
    const double* const stateIn, 
    const double* const cellCentre, 
    const double* const cellSize, 
    const double t,const double dt, 
    const int faceIndex, 
    const int normalNonZero 
  );

  
//AMR Routines
//------------

  void faceUnknownsProlongation(
    double* restrict lQhbndFine,
    double* restrict lFhbndFine,
    const double* const restrict lQhbndCoarse,
    const double* const restrict lFhbndCoarse,
    const int coarseGridLevel,
    const int fineGridLevel,
    const int* const subfaceIndex
  );
  
  // used by faceIntegral, only restrict the flux on the face
  void faceFluxRestriction(
    double* restrict lFhbndCoarse,
    const double* const restrict lFhbndFine,
    const int* const subfaceIndex,
    const int levelDelta
  );

  void volumeUnknownsProlongation(
    double* restrict luhFine,
    const double* const restrict luhCoarse,
    const int coarseGridLevel,
    const int fineGridLevel,
    const int* const subcellIndex
  );
  
  void volumeUnknownsRestriction(
    double* restrict luhCoarse,
    const double* const restrict luhFine,
    const int coarseGridLevel,
    const int fineGridLevel,
    const int* const subcellIndex
  );

//Limiter
//-------

  void projectOnFVLimiterSpace(const double* const luh, double* const lim);

  void projectOnDGSpace(const double* const lim, double* const luh);

  void findCellLocalMinAndMax(
    const double* const luh,
    const exahype::solvers::ADERDGSolver* solver,
    double* const localMinPerVariables, 
    double* const localMaxPerVariable
  );

  void findCellLocalLimiterMinAndMax(
    const double* const lim,
    const exahype::solvers::ADERDGSolver* solver,
    double* const localMinPerObservable, 
    double* const localMaxPerObservable
  );

  bool discreteMaximumPrincipleAndMinAndMaxSearch(
    const double* const luh,
    const exahype::solvers::ADERDGSolver* solver,
    const double relaxationParameter,
    const double differenceScaling,
    double* boundaryMinPerVariables, 
    double* boundaryMaxPerVariables
  );

  //private
  void compareWithADERDGSolutionAtGaussLobattoNodes(
    const double* const luh,
    const exahype::solvers::ADERDGSolver* solver,
    double* min, 
    double* max
  );

  //private
  void compareWithADERDGSolutionAtFVSubcellCenters(
    const double* const luh,
    const exahype::solvers::ADERDGSolver* solver,
    double* min, 
    double* max
  );
  
}
}
}

#include "kernels/SWE_MySWESolver_ADERDG/ConfigurationParameters.cpph"
#include "kernels/SWE_MySWESolver_ADERDG/matrixUtils.cpph"

#endif // _EXAHYPE_SWE_MYSWESOLVER_ADERDG_KERNELS_ADERDG_KERNELS_H_