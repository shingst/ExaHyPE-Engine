// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#ifndef __AbstractMyEulerSolver_ADERDG_CLASS_HEADER__
#define __AbstractMyEulerSolver_ADERDG_CLASS_HEADER__

#include <ostream>
#include <algorithm>

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/UserSolverInterface.h"

/**
 * We include Peano's assertion collection here.
 */
#include "tarch/Assertions.h"
#include "tarch/la/Vector.h"

namespace EulerADERDG {
  class MyEulerSolver_ADERDG;
  class AbstractMyEulerSolver_ADERDG;
}

class EulerADERDG::AbstractMyEulerSolver_ADERDG: public exahype::solvers::ADERDGSolver, public exahype::solvers::UserADERDGSolverInterface {
  public:
    static constexpr int Order                     = 5;
    static constexpr int NumberOfVariables         = 5;
    static constexpr int NumberOfParameters        = 0;
    static constexpr int NumberOfGlobalObservables = 0;
    static constexpr int NumberOfDMPObservables    = 2; // only of interest if this ADERDGSolver is a component of a LimitingADERDSolver 
    static constexpr int MaxPicardIterations       = -1;
    static constexpr bool UseMaxPicardIterations   = false;
    static constexpr double CFL                    = 0.9;
    
    // virtual getters for the constexpr's
    int constexpr_getNumberOfVariables()  const override { return NumberOfVariables; };
    int constexpr_getNumberOfParameters() const override { return NumberOfParameters; };
    int constexpr_getOrder()              const override { return Order; };
    double constexpr_getCFLNumber()       const override { return CFL; };  
  
    class VariableMetrics;
    class Variables;
    class ReadOnlyVariables;
    class Fluxes;
    class VariableShortcuts;
    class VariableMultiplicities;
    class VariableNames;
    
    AbstractMyEulerSolver_ADERDG(
      const double maximumMeshSize,
      const int maximumMeshDepth,
      const int haloCells,
      const int regularisedFineGridLevels,
      const exahype::solvers::Solver::TimeStepping timeStepping,const int DMPObservables
);

    /**
     * This operation should be overwritten in your application-specific 
     * solver. Alternatively, make your own subclass useConservativeFlux()
     * return false.
     */
    void flux(const double* const Q,double** const F) override;

    /**
     * This operation should be overwritten in your application-specific 
     * solver. Alternatively, make your own subclass useConservativeFlux()
     * return false.
     */
    void viscousFlux(const double* const Q,const double* const gradQ,double** const F) override;

    /**
     * This operation should be overwritten in your application-specific 
     * solver. Alternatively, make your own subclass useConservativeFlux()
     * return false.
     */
    void viscousEigenvalues(const double* const Q,const int d,double* const lambda) override;

    
    /**
     * Default implementation. Has to be be overwritten by user's solver if you 
     * make the corresponding use operation activate the feature.
     *
     * See superclass for function's semantics.
     */
    void algebraicSource(const double* const Q,double* const S) override;
        
    /**
     * Default implementation. Has to be be overwritten by user's solver if you 
     * make the corresponding use operation activate the feature.
     *
     * See superclass for function's semantics.
     */
    void nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) override;
        
    /**
     * Default implementation. Has to be be overwritten by user's solver if you 
     * make the corresponding use operation activate the feature.
     *
     * See superclass for function's semantics.
     */
        
    /**
     * Default implementation. Has to be be overwritten by user's solver if you 
     * make the corresponding use operation activate the feature.
     *
     * See superclass for function's semantics.
     */
    void pointSource(const double* const Q,const double* const x,const double t,const double dt, double* const forceVector,int n) override;
       
    /**
     * Default implementation. Please overwrite.
     *
     * See superclass for function's semantics.
     */
    bool isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,
      const bool wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t, const double dt) const override { return true; }
    
    /**
     * Default implementation. Please overwrite.
     *
	 * See superclass for function's semantics.
	 */
    void mapDiscreteMaximumPrincipleObservables(double* const observables, const double* const Q) const override {
      if (NumberOfDMPObservables>0) {
      	std::copy_n(Q,NumberOfDMPObservables,observables);
      }
  	}
     
    int fusedSpaceTimePredictorVolumeIntegral(double* const lduh, double* const lQhbnd, double* const lFhbnd, const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, const double t, const double dt) override;
    void solutionUpdate(double* const luh,const double* const lduh,const double dt) override;
    void riemannSolver(double* const FL,double* const FR,const double* const QL,const double* const QR,const double dt,const int direction, bool isBoundaryFace, int faceIndex) override;
    void boundaryConditions(double* const update, double* const fluxIn,const double* const stateIn,const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& cellCentre,const tarch::la::Vector<DIMENSIONS,double>& cellSize,const double t,const double dt,const int direction,const int orientation) override;
    void faceIntegral(double* const lduh,const double* const lFhbnd,const int direction, const int orientation,const tarch::la::Vector<DIMENSIONS-1,int>& subfaceIndex,const int levelDelta,const tarch::la::Vector<DIMENSIONS, double>& cellSize) override;    
    double stableTimeStepSize(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& dx) override;
    void adjustSolution(double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) override; 
    void faceUnknownsProlongation(double* const lQhbndFine,double* const lFhbndFine,const double* const lQhbndCoarse,const double* const lFhbndCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1,int>& subfaceIndex) override;
    void volumeUnknownsProlongation(double* const luhFine,const double* const luhCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) override;
    void volumeUnknownsRestriction(double* const luhCoarse,const double* const luhFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) override;

    static void constantsToString(std::ostream& os);
    static void abortWithMsg(const char* const msg);
    
    //override the size of unused data storage to -1 to not allocate it
    int getTempPointForceSourcesSize() const {return -1;} //pointSource not required
    int getTempSpaceTimeFluxUnknowns1Size() const {return -1;} //gradQ not required
    
    //not used PDE terms
    void multiplyMaterialParameterMatrix(const double* const Q, double* const rhs) final {}
    
};

#endif // __AbstractMyEulerSolver_ADERDG_CLASS_HEADER__