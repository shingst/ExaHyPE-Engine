#ifndef __NavierStokesSolver_CLASS_HEADER__
#define __NavierStokesSolver_CLASS_HEADER__

// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================

#include <string>

#include "exahype/solvers/LimitingADERDGSolver.h"
#include "NavierStokesSolver_ADERDG.h"
#include "NavierStokesSolver_FV.h"

namespace NavierStokes{
  class NavierStokesSolver;
}

class NavierStokes::NavierStokesSolver: public exahype::solvers::LimitingADERDGSolver {  
  public:
    static constexpr int NumberOfVariables         = NavierStokes::AbstractNavierStokesSolver_ADERDG::NumberOfVariables;
    static constexpr int NumberOfParameters        = NavierStokes::AbstractNavierStokesSolver_ADERDG::NumberOfParameters;
    static constexpr int Order                     = NavierStokes::AbstractNavierStokesSolver_ADERDG::Order;
    static constexpr int NumberOfGlobalObservables = NavierStokes::AbstractNavierStokesSolver_ADERDG::NumberOfGlobalObservables;
    static constexpr int NumberOfDMPObservables    = NavierStokes::AbstractNavierStokesSolver_ADERDG::NumberOfDMPObservables;
    static constexpr int GhostLayerWidth           = NavierStokes::AbstractNavierStokesSolver_FV::GhostLayerWidth;
      
    NavierStokesSolver(
        const double maximumMeshSize,
        const int maximumMeshDepth,
        const int haloCells,
        const int regularisedFineGridLevels,
        const exahype::solvers::Solver::TimeStepping timeStepping,
        const int DMPObservables,
        const double DMPRelaxationParameter,
        const double DMPDifferenceScaling
);
    
    void projectOnFVLimiterSpace(const double* const luh, double* const lim) const override;
    void projectOnDGSpace(const double* const lim, double* const luh) const override;
    bool discreteMaximumPrincipleAndMinAndMaxSearch(const double* const luh, double* const boundaryMinPerVariables, double* const boundaryMaxPerVariables) override;
    void findCellLocalMinAndMax(const double* const luh, double* const localMinPerVariables, double* const localMaxPerVariable) override;
    void findCellLocalLimiterMinAndMax(const double* const lim, double* const localMinPerObservable, double* const localMaxPerObservable) override;
};

#endif // __NavierStokesSolver_CLASS_HEADER__