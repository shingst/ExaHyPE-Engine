/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
 
#ifndef _EXAHYPE_SOLVERS_BASISSOLVER_H
#define _EXAHYPE_SOLVERS_BASISSOLVER_H

namespace exahype {
namespace solvers {

class UserSolverInterface;
class UserADERDGSolverInterface;
class UserFiniteVolumesSolverInterface;

} /* namespace solvers */
} /* namespace exahype */



/**
 * The Basis API for User solvers, purely virtual. New from 2017-05-14.
 * Cf. https://gitlab.lrz.de/exahype/ExaHyPE-Engine/issues/143
 * 
 * This is for a template-free glue code (Abstract*Solver).
 * 
 * Direct classes which inherit UserSolverInterface:
 *   1) ADERDGSolver
 *   3) FVSolver
 *
 * TODO: The UseAdjustSolution() user functions should be unified accross
 *   FV/ADERDG solvers and also be added here.
 **/
class exahype::solvers::UserSolverInterface {
public:
  virtual ~UserSolverInterface() {};

 /**
  * @defgroup Theoretically-Constexpr-Getters
  */
  ///@{
  // Read off the constexpr's in the Abstract*Solver
  virtual int constexpr_getNumberOfVariables()  const = 0;
  virtual int constexpr_getNumberOfParameters() const = 0;
  virtual double constexpr_getCFLNumber()       const = 0;
  ///@}

  
  /**
     * Compute the eigenvalues of the flux tensor per coordinate direction \p d.
     *
     * \param[in] Q  the conserved variables associated with a quadrature node
     *               as C array (already allocated).
     * \param[in] d  the column of the flux vector (d=0,1,...,DIMENSIONS).
     * \param[inout] lambda the eigenvalues as C array (already allocated).
     */
  virtual void eigenvalues(const double* const Q,const int d,double* lambda) = 0;

 /**
  * @defgroup User PDE
  */
  ///@{
  /**
   * Compute a pointSource contribution.
   * 
   * @TODO: Document me, please.
   **/
   virtual void pointSource(const double* const Q,const double* const x,const double t,const double dt, double* forceVector,int n) = 0;

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
   * \param[in]    Q the conserved variables (and parameters) associated with a quadrature point
   *                 as C array (already allocated).
   * \param[inout] S the source point as C array (already allocated).
   */
  virtual void algebraicSource(const double* const Q,double* S) = 0;
  
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
   * \param[in]   Q   the vector of unknowns at the given position
   * \param[in]   gradQ   the gradients of the vector of unknowns,
   *                  stored in a linearized array.
   * \param[inout]  The vector BgradQ (extends nVar), already allocated. 
   *
   **/
  virtual void nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) = 0;

  virtual void multiplyMaterialParameterMatrix(const double* const Q, double* rhs) = 0;
  
  /**
   * Compute the conserved flux.
   * 
   * \param[in]  Q the conserved variabels (and parameters) associated with a
   *               quadrature point as C array.
   * \param[inout] F a C array with shape [nDim][nVars]. That is, this is an C list
   *               holding pointers to actual lists. Thus, the storage may be noncontinous.
   *               In any case, the storage has already been allocated.
   **/
  virtual void flux(const double* const Q,double** F) = 0;
  ///@}
};
 // UserSolverInterface

class exahype::solvers::UserADERDGSolverInterface : public exahype::solvers::UserSolverInterface {
public:
  virtual ~UserADERDGSolverInterface() {};

  virtual int constexpr_getOrder()  const  = 0;

  /**
   * Adjust the conserved variables and parameters (together: Q) at a given time t at the (quadrature) point x.
   *
   * \note Use this function and ::useAdjustSolution to set initial conditions.
   *
   * \param[in]    x         the physical coordinate on the face.
   * \param[in]    t         the start of the time interval.
   * \param[in]    dt        the width of the time interval.
   * \param[inout] Q         the conserved variables (and parameters) associated with a quadrature point
   *                         as C array (already allocated).
   */
  virtual void adjustPointSolution(const double* const x,const double t,const double dt,double* Q) = 0;
  ///@}

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
  virtual void boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double * const fluxIn,const double* const stateIn,double *fluxOut,double* stateOut) = 0;
};

class exahype::solvers::UserFiniteVolumesSolverInterface : public exahype::solvers::UserSolverInterface {
public:
  virtual ~UserFiniteVolumesSolverInterface() {};

  virtual int constexpr_getPatchSize()  const  = 0;
  virtual int constexpr_getGhostLayerWidth() const  = 0;



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
  virtual void boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const stateIn,double* stateOut) = 0;
};

#endif /* _EXAHYPE_SOLVERS_BASISSOLVER_H */
