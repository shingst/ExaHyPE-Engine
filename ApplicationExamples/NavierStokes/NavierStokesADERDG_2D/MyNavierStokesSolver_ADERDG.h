#ifndef __MyNavierStokesSolver_ADERDG_CLASS_HEADER__
#define __MyNavierStokesSolver_ADERDG_CLASS_HEADER__

// This file was initially generated by the ExaHyPE toolkit.
// You can modify it in order to extend your solver with features.
// Whenever this file is present, a re-run of the ExaHyPE toolkit will
// not overwrite it. Delete it to get it regenerated.
//
// ========================
//   www.exahype.eu
// ========================

#include <ostream>

#include "AbstractMyNavierStokesSolver_ADERDG.h"
#include "exahype/parser/ParserView.h"

/**
 * We use Peano's logging
 */
#include "tarch/logging/Log.h"

namespace NavierStokesADERDG{
  class MyNavierStokesSolver_ADERDG;
}

class NavierStokesADERDG::MyNavierStokesSolver_ADERDG : public NavierStokesADERDG::AbstractMyNavierStokesSolver_ADERDG {
  private:
    /**
     * Log device
     */
    static tarch::logging::Log _log;
  public:
    MyNavierStokesSolver_ADERDG(
        const double maximumMeshSize,
        const int maximumMeshDepth,
        const int haloCells,
        const int regularisedFineGridLevels,
        const exahype::solvers::Solver::TimeStepping timeStepping,
        const int limiterHelperLayers,
        const int DMPObservables
        );

    /**
     * Initialise the solver.
     *
     * \param[in] cmdlineargs the command line arguments.
     * \param[in] constants   access to the constants specified for the solver.
     */
    void init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) final override;

    /**
     * Adjust the conserved variables and parameters (together: Q) at a given time t at the (quadrature) point x.
     *
     * \note Please overwrite function adjustSolution(...) if you want to
     * adjust the solution degrees of freedom in a cellwise manner.
     *
     * \param[in]    x         the physical coordinate on the face.
     * \param[in]    t         the start of the time interval.
     * \param[in]    dt        the width of the time interval.
     * \param[inout] Q         the conserved variables (and parameters) associated with a quadrature point
     *                         as C array (already allocated).
     */
    void adjustPointSolution(const double* const x,const double t,const double dt,double* Q) final override;

    /**
     * Compute the eigenvalues of the flux tensor per coordinate direction \p d.
     *
     * \param[in] Q  the conserved variables associated with a quadrature node
     *               as C array (already allocated).
     * \param[in] d  the column of the flux vector (d=0,1,...,DIMENSIONS).
     * \param[inout] lambda the eigenvalues as C array (already allocated).
     */
    void eigenvalues(const double* const Q,const int d,double* lambda) final override;
    void viscousEigenvalues(const double* const Q,const int d,double* lambda) final override;
    
    /**
     * Impose boundary conditions at a point on a boundary face
     * within the time interval [t,t+dt].
     *
     * \param[in]    x         the physical coordinate on the face.
     * \param[in]    t         the start of the time interval.
     * \param[in]    dt        the width of the time interval.
     * \param[in]    faceIndex indexing of the face (0 -- {x[0]=xmin}, 1 -- {x[1]=xmax}, 2 -- {x[1]=ymin}, 3 -- {x[2]=ymax}, and so on,
     *                         where xmin,xmax,ymin,ymax are the bounds of the cell containing point x.
     * \param[in]    d         the coordinate direction the face normal is pointing to.
     * \param[in]    QIn       the conserved variables at point x from inside of the domain
     *                         and time-averaged (over [t,t+dt]) as C array (already allocated).
     * \param[in]    FIn       the normal fluxes at point x from inside of the domain
     *                         and time-averaged (over [t,t+dt]) as C array (already allocated).
     * \param[inout] QOut      the conserved variables at point x from outside of the domain
     *                         and time-averaged (over [t,t+dt]) as C array (already allocated).
     * \param[inout] FOut      the normal fluxes at point x from outside of the domain
     *                         and time-averaged (over [t,t+dt]) as C array (already allocated).
     */
    void boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double * const fluxIn,const double* const stateIn,double *fluxOut,double* stateOut) final override;
    
    /**
     * Evaluate the refinement criterion within a cell.
     *
     * \note Instead of a variables array at a single quadrature point we give
     * you all NumberOfVariables*(Order+1)^DIMENSIONS solution degrees of freedom.
     *
     * \note Use this function and ::adjustSolution to set initial conditions.
     *
     * \param[in]    centre    The centre of the cell.
     * \param[in]    dx        The extent of the cell.
     * \param[in]    t         the start of the time interval.
     * \param[in]    dt        the width of the time interval.
     * \return One of exahype::solvers::Solver::RefinementControl::{Erase,Keep,Refine}.
     */
    exahype::solvers::Solver::RefinementControl refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& centre,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) override;
    
    //PDE

/* flux() function not included, as requested in the specification file */

    /**
     * Compute the flux tensor with diffusive components.
     *
     * \param[in]    Q the conserved variables (and parameters) associated with a quadrature point
     *                 as C array (already allocated).
     *               gradQ the gradient of the conserved variables (and parameters) associated with a quadrature point
     *                     as C array (already allocated).
     * \param[inout] F the fluxes at that point as C array (already allocated).
     */
    void viscousFlux(const double* const Q,const double* gradQ, double** F) final override;

/* algebraicSource() function not included, as requested by the specification file */

/* nonConservativeProduct() function is not included, as requested in the specification file */

/* pointSource() function not included, as requested in the specification file */

/* multiplyMaterialParameterMatrix() not included, as requested in the specification file */

    void mapDiscreteMaximumPrincipleObservables(
            double* observables,const int numberOfObservables,
            const double* const Q) const override;
    bool isPhysicallyAdmissible(
            const double* const solution,
            const double* const observablesMin,const double* const observablesMax,
            const int numberOfObservables,
            const tarch::la::Vector<DIMENSIONS,double>& center,
            const tarch::la::Vector<DIMENSIONS,double>& dx,
            const double t, const double dt) const override;

};

#endif // __MyNavierStokesSolver_ADERDG_CLASS_HEADER__
