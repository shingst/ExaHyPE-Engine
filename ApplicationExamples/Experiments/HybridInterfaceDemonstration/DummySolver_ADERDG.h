#ifndef __DummySolver_ADERDG_CLASS_HEADER__
#define __DummySolver_ADERDG_CLASS_HEADER__

// This file was initially generated by the ExaHyPE toolkit.
// You can modify it in order to extend your solver with features.
// Whenever this file is present, a re-run of the ExaHyPE toolkit will
// not overwrite it. Delete it to get it regenerated.
//
// ========================
//   www.exahype.eu
// ========================

#include <ostream>

#include "AbstractDummySolver_ADERDG.h"
#include "exahype/parser/ParserView.h"

/**
 * We use Peano's logging
 */
#include "tarch/logging/Log.h"

namespace HybridInterfaceDemonstration{
  class DummySolver_ADERDG;
}

class HybridInterfaceDemonstration::DummySolver_ADERDG : public HybridInterfaceDemonstration::AbstractDummySolver_ADERDG {
  private:
    /**
     * Log device
     */
    static tarch::logging::Log _log;
  public:
    DummySolver_ADERDG(
        const double maximumMeshSize,
        const int maximumMeshDepth,
        const int haloCells,
        const int haloBufferCells,
        const int limiterBufferCells,
        const int regularisedFineGridLevels,
        const exahype::solvers::Solver::TimeStepping timeStepping,
        const int DMPObservables
);

    /**
     * Initialise the solver.
     *
     * @param[in] cmdlineargs the command line arguments.
     * @param[in] constants   access to the constants specified for the solver.
     */
    void init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) final override;
    
    /**
     * Adjust the conserved variables and parameters (together: Q) at a given time t at the (quadrature) point x.
     *
     * @note Please overwrite function adjustSolution(...) if you want to
     * adjust the solution degrees of freedom in a cellwise manner.
     *
     * @param[in]    x   physical coordinate on the face.
     * @param[in]    t   start of the time interval.
     * @param[in]    dt  width of the time interval.
     * @param[in]    Q   vector of state variables (plus parameters); 
     *                   range: [0,nVar+nPar-1], already allocated.
     */
    void adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) final override;

    /**
     * Compute the eigenvalues of the flux tensor per coordinate  @p direction.
     *
     * @param[in]    Q          vector of state variables (plus parameters); 
     *                          range: [0,nVar+nPar-1], already allocated.
     * @param[in]    direction  normal direction of the face / column of the flux vector (range: [0,nDims-1]).
     * @param[inout] lambda     eigenvalues as C array;
     *                          range: [0,nVar-1], already allocated.
     */
    void eigenvalues(const double* const Q,const int d,double* const lambda) final override;
    

    /**
     * Impose boundary conditions at a point on a boundary face
     * within the time interval [t,t+dt].
     *
     * @param[in]    x         physical coordinate on the face; 
     *                         range: [0,nDims-1].
     * @param[in]    t         start of the time interval.
     * @param[in]    dt        width of the time interval.
     * @param[in]    faceIndex indexing of the face (0 -- {x[0]=xmin}, 1 -- {x[1]=xmax}, 2 -- {x[1]=ymin}, 3 -- {x[2]=ymax}, and so on,
     *                         where xmin,xmax,ymin,ymax are the bounds of the cell containing point x.
     * @param[in]    direction coordinate direction the face normal is pointing to.
     * 
     * @param[in]    QIn       the conserved variables (plus parameters) at point x from inside of the domain
     *                         and time-averaged (over [t,t+dt]);
     *                         range: [0,nVar+nPar-1], already allocated.
     * @param[in]    FIn       the normal fluxes at point x from inside of the domain
     *                         and time-averaged (over [t,t+dt]);
     *                         range: [0,nVar-1], already allocated.
     *                         
     * @param[inout] QOut      the conserved variables at point x from outside of the domain
     *                         and time-averaged (over [t,t+dt]);
     *                         range: [0,nVar+nPar-1], already allocated.
     * @param[inout] FOut      the normal fluxes at point x from outside of the domain;
     *                         range: [0,nVar-1], already allocated.
     */
    void boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) final override;
    
    /**
     * Evaluate the refinement criterion within a cell.
     *
     * @note Instead of a variables array at a single quadrature point we give
     * you all NumberOfVariables*(Order+1)^DIMENSIONS solution degrees of freedom.
     *
     * @note Use this function and ::adjustSolution to set initial conditions.
     *
     * @param[in]    luh    all state variables (plus parameters) of a cell; 
     *                      range[outer->inner]: [0,(p+1)^nDim-1]x[0,nVar+nPar-1], 
     *                      already allocated.
     *                      
     * @param[in]    centre centre of the cell.
     * @param[in]    dx     extent of the cell.
     * @param[in]    t      start of the time interval.
     * @param[in]    dt     width of the time interval.
     * 
     * @return One of exahype::solvers::Solver::RefinementControl::{Erase,Keep,Refine}.
     */
    exahype::solvers::Solver::RefinementControl refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& centre,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) override;
    
    //PDE
    
    /**
     * Compute the flux tensor.
     *
     * @param[in]    Q vector of state variables (plus parameters); 
     *                 range: [0,nVar+nPar-1], already allocated.
     *                 
     * @param[inout] F flux at that point;
     *                 range[outer->inner]: [0,nDim-1]x[0,nVar-1], 
     *                 already allocated.
     */
    void flux(const double* const Q,double** const F) final override;

/* viscousFlux() function not included, as requested in the specification file */


    
     /**
     * Compute the Algebraic Sourceterms.
     * 
     * You may want to overwrite this with your PDE Source (algebraic RHS contributions).
     * However, in all schemes we have so far, the source-type contributions are
     * collected with non-conservative contributions into a fusedSource, see the
     * fusedSource method. From the kernels given with ExaHyPE, only the fusedSource
     * is called and there is a default implementation for the fusedSource calling
     * again seperately the nonConservativeProduct function and the algebraicSource
     * function.
     *
     * @param[in]    Q vector of state variables (plus material 
     *                 parameters); range: [0,nVar+nPar-1], already allocated.
     * @param[inout] S source term; range: [0,nVar-1], already allocated.
     */
    void algebraicSource(const double* const Q,double* const S) final override;

    /**
     * Compute the nonconservative term $B(Q) \nabla Q$.
     * 
     * This function shall return a vector BgradQ which holds the result
     * of the full term. To do so, it gets the vector Q and the matrix
     * gradQ which holds the derivative of Q in each spatial direction.
     * Currently, the gradQ is a continous storage and users can use the
     * kernels::idx2 class in order to compute the positions inside gradQ.
     *
     * @TODO: Check if the following is still right:
     * 
     * !!! Warning: BgradQ is a vector of size NumberOfVariables if you
     * use the ADER-DG kernels for nonlinear PDEs. If you use
     * the kernels for linear PDEs, it is a tensor with dimensions
     * Dim x NumberOfVariables.
     * 
     * @param[in]  Q        vector of state variables (plus material 
     *                      parameters); range: [0,nVar+nPar-1], already allocated.
     * @param[in]    gradQ  the gradients of the vector of unknowns, stored 
     *                      in a linearized array. (range: [0,dim*(nVar+nPar)-1].
     * @param[inout] BgradQ the nonconservative product (extends nVar), 
     *                      already allocated. 
     */
    void nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) final override;

/* pointSource() function not included, as requested in the specification file */

/* multiplyMaterialParameterMatrix() not included, as requested in the specification file */


    bool isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,
      const bool wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t) const override;
      
      
   void beginTimeStep(const double minTimeStamp,const bool isFirstTimeStepOfBatchOrNoBatch) override;
   void endTimeStep(const double minTimeStamp,const bool isLastTimeStepOfBatchOrNoBatch) override;
};

#endif // __DummySolver_ADERDG_CLASS_HEADER__
