#ifndef __AbstractCCZ4Solver_ADERDG_CLASS_HEADER__
#define __AbstractCCZ4Solver_ADERDG_CLASS_HEADER__

// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include <ostream>
#include <algorithm>

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/UserSolverInterface.h"

/**
 * We include Peano's assertion collection here.
 */
#include "tarch/Assertions.h"
#include "tarch/la/Vector.h"

namespace CCZ4{
  class CCZ4Solver_ADERDG;
  class AbstractCCZ4Solver_ADERDG;
}

class CCZ4::AbstractCCZ4Solver_ADERDG: public exahype::solvers::ADERDGSolver, public exahype::solvers::UserADERDGSolverInterface {
  public:
    static constexpr int NumberOfVariables       = 59;
    static constexpr int NumberOfParameters      = 0;
    static constexpr int Order                   = 3; 
    static constexpr int MaxPicardIterations     = -1;
    static constexpr bool UseMaxPicardIterations = false;
    static constexpr double CFL                  = 0.9;
    
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
    
    AbstractCCZ4Solver_ADERDG(double maximumMeshSize,int maximumAdaptiveMeshDepth,int DMPObservables,int limiterHelperLayers,exahype::solvers::Solver::TimeStepping timeStepping);

    /**
     * This operation should be overwritten in your application-specific 
     * solver. Alternatively, make your own subclass useConservativeFlux()
     * return false.
     */
    void flux(const double* const Q,double** F) override;

    /**
     * Default implementation. Has to be be overwritten by user's solver if you 
     * make the corresponding use operation activate the feature.
     *
     * See superclass for function's semantics.
     */
    void algebraicSource(const double* const Q,double* S) override;
        
    /**
     * Default implementation. Has to be be overwritten by user's solver if you 
     * make the corresponding use operation activate the feature.
     *
     * See superclass for function's semantics.
     */
    void nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) override;
        
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
    void pointSource(const double* const Q,const double* const x,const double t,const double dt, double* forceVector,int n) override;
       
    /**
     * Default implementation. Please overwrite.
     *
     * See superclass for function's semantics.
     */
    bool isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,const int numberOfObservables,
      const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t, const double dt) const override { return true; }
    
    /**
     * Default implementation. Please overwrite.
     *
	 * See superclass for function's semantics.
	 */
    void mapDiscreteMaximumPrincipleObservables(double* observables,const int numberOfObservables,const double* const Q) const override {
      if (numberOfObservables>0) {
      	std::copy_n(Q,numberOfObservables,observables);
      }
  	}
     
    int fusedSpaceTimePredictorVolumeIntegral(double* lduh, double* lQhbnd, double* lFhbnd, const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, const double t, const double dt) override;
    void solutionUpdate(double* luh,const double* const lduh,const double dt) override;
    void surfaceIntegral(double* lduh,const double* const lFhbnd,const tarch::la::Vector<DIMENSIONS,double>& dx) override;
    void riemannSolver(double* FL,double* FR,const double* const QL,const double* const QR,const double dt,const int direction, bool isBoundaryFace, int faceIndex) override;
    void boundaryConditions(double* fluxIn,const double* const stateIn,const tarch::la::Vector<DIMENSIONS,double>& cellCentre,const tarch::la::Vector<DIMENSIONS,double>& cellSize,const double t,const double dt,const int faceIndex,const int direction) override;
    double stableTimeStepSize(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& dx) override;
    void adjustSolution(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) override; 
    void faceUnknownsProlongation(double* lQhbndFine,double* lFhbndFine,const double* lQhbndCoarse,const double* lFhbndCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1,int>& subfaceIndex) override;
    void faceUnknownsRestriction(double* lQhbndCoarse,double* lFhbndCoarse,const double* lQhbndFine,const double* lFhbndFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1,int>& subfaceIndex) override;
    void volumeUnknownsProlongation(double* luhFine,const double* luhCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) override;
    void volumeUnknownsRestriction(double* luhCoarse,const double* luhFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) override;

    static void constantsToString(std::ostream& os);
    static void abortWithMsg(const char* const msg);
    
    //override the size of unused data storage to -1 to not allocate it
    int getTempPointForceSourcesSize() const {return -1;} //pointSource not required
    
    //not used PDE
    void multiplyMaterialParameterMatrix(const double* const Q, double* rhs) final {}
    
};

#endif // __AbstractCCZ4Solver_ADERDG_CLASS_HEADER__
