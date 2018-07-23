/*
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
 *
 * \author Dominic E. Charrier, Tobias Weinzierl
 **/
#ifndef LIMITEDADERDGSOLVER_H_
#define LIMITEDADERDGSOLVER_H_

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

#include "tarch/multicore/BooleanSemaphore.h"

#include "exahype/plotters/Plotter.h"

#include "exahype/profilers/simple/NoOpProfiler.h"

namespace exahype {
namespace solvers {

/**
 * A solver that combines high-order ADER-DG
 * with a more robust Finite Volumes solver in areas
 * where shocks are present.
 *
 * The ADER-DG solver is regarded as the main solver
 * that deals with time step sizes, mesh update requests etc.
 *
 * <h1>Algorithm sections</h1>
 * The solver might be active in one of the
 * following algorithm sections.
 *
 * <h2>Mesh refinement</h2>
 * If a solver requests regular mesh refinement, no limiter status spreading must be
 * performed! We then simply adjust the mesh.
 * No limiter patches are allocated!
 *
 * <h2>Local recomputation</h2>
 * If a solver requests local recomputation, we perform limiter status spreading.
 * No limiter patches are allocated during these spreading iterations.
 * Only if it is decided that no global recomputation is performed,
 * limiter patches are allocated.
 *
 * <h2>Global recomputation</h2>
 * The solver is redoing the last ADER-DG time
 * step completely but performs some mesh
 * refinement beforehand.
 *
 * More precisely, the solver will fist perform a rollback to
 * the previous time step> It will then perform mesh refinement
 * until a troubled compute cell and all its compute cell
 * neighours are placed on the finest mesh level.
 * Next it will perform the computation of a new
 * time step size and then, the computation of
 * a new space-time predictor.
 * Then, the last time step is redone.
 *
 * The following scenario causes the solver
 * to switch to this algorithmic section:
 *
 * Scenario 1:
 * A compute cell was marked as troubled on a
 * mesh level coarser than the finest one.
 */
class LimitingADERDGSolver;

} /* namespace solvers */
} /* namespace exahype */

class exahype::solvers::LimitingADERDGSolver : public exahype::solvers::Solver {

protected:
  /**
   * The ADERDG solver.
   */
  std::unique_ptr<exahype::solvers::ADERDGSolver> _solver;

  /**
   * The finite volumes solver used for the a posteriori subcell limiting.
   */
  std::unique_ptr<exahype::solvers::FiniteVolumesSolver> _limiter;
  /**
   * The maximum relaxation parameter
   * used for the discrete maximum principle.
   */
  const double _DMPMaximumRelaxationParameter;

  /**
   * The difference scaling
   * used for the discrete maximum principle.
   */
  const double _DMPDifferenceScaling;

  /**
   * A counter holding the number of iterations to
   * cure a troubled cell.
   * This counter will be initialised to a certain
   * (user-dependent?) value if a cell is flagged as troubled.
   *
   * If the cell is not troubled for one iteration, the counter is
   * decreased until it reaches 0. Then, the
   * cell is considered as cured.
   * Note that the counter can be reset to the maximum value
   * in the meantime if the cell is marked again as troubled.
   *
   * This counter prevents that a cell is toggling between
   * troubled and Ok (cured).
   */
  int _iterationsToCureTroubledCell;

private:
  typedef exahype::records::ADERDGCellDescription SolverPatch;

  typedef exahype::records::FiniteVolumesCellDescription LimiterPatch;

  /**
   * Log device.
   */
  static tarch::logging::Log _log;


  #ifdef Parallel
  std::vector<double> _receivedMax;
  std::vector<double> _receivedMin;
  /**
   * TODO(WORKAROUND): We store these fields in order
   * to use the symmetric boundary exchanger of Peano
   * which does not yet support asymmetric send buffers.
   */
  std::vector<double> _invalidObservables;
  #endif

  /**
   * TODO(Dominc): Remove after docu is recycled.
   *
   * This operation sets the solutions' minimum and maximum value on a cell.
   * The routine is to be invoked after the code has determined the new minimum
   * and maximum value within a cell. In turn, it evaluates whether the new
   * minimum and maximum value have decreased or grown, respectively.
   *
   * If the new min/max values indicate that the new solution comprises
   * oscillations, the routine returns false. This is an indicator that the
   * solution should be limited.
   *
   * If the new min/max values fit, the routine returns true.
   *
   * <h2>Implementation</h2>
   * We hold the min/max information exclusively on the faces. The first thing
   * the routine does is to project the min/max values into the cell. For this
   * it evaluates the 2d faces. The projected value then is compared to the
   * arguments. Once the results of the operation is determined, the routine
   * writes the new arguments onto the 2d face entries. This, on the one hand,
   * stores the data for the subsequent time step, but it also propagates the
   * min/max information into the face-connected neighbours.
   *
   * @param  min          New minimum values within the cell. Array of length
   *                      _numberOfUnknowns.
   * @param  max          New maximum values within the cell
   * @param  solverIndex  Number of the solver within the cell. Please ensure
   *                      that solverIndex refers to an ADER-DG solver.
   * @return True if the new min and max values fit into the restricted min
   *   max solutions. Return false if we seem to run into oscillations.
   */
  //  void setSolutionMinMax(double* min, double* max) const;

  /**
   * Merge the solution min and max values on a face between two cell
   * descriptions. Signature is similar to that of the solver of a Riemann problem.
   */
  void mergeSolutionMinMaxOnFace(
      SolverPatch& pLeft,
      SolverPatch& pRight,
      const int faceIndexLeft,
      const int faceIndexRight) const;

  /**
   * Checks if the updated solution
   * of the ADER-DG solver contains unphysical oscillations (false)
   * or not (true).
   *
   * Compute the new min and max at the same time.
   */
  bool evaluateDiscreteMaximumPrincipleAndDetermineMinAndMax(SolverPatch& solverPatch);

  /**
   * Checks if the updated solution
   * of the ADER-DG solver contains unphysical oscillations (false)
   * or not (true).
   */
  bool evaluateDiscreteMaximumPrinciple(SolverPatch& solverPatch);

  /**
   * Checks if the updated solution
   * of the ADER-DG solver is
   * a physically admissible one (true).
   *
   * \note This method assumes the (previous) refinement status
   * was not modified yet by another routine.
   */
  bool evaluatePhysicalAdmissibilityCriterion(SolverPatch& solverPatch);

  /**
   * Computes the cell-local minimum and maximum
   * of the solution (per variable)
   * and makes them accessible per face.
   */
  void determineSolverMinAndMax(SolverPatch& solverPatch);

  /**
   * Computes the cell-local minimum and maximum
   * of the limiter's solution (per variable)
   * and makes them accessible per face.
   *
   * This method is used for troubled cells that
   * do not hold a valid ADER-DG solution,
   * as well as their neighbours.
   */
  void determineLimiterMinAndMax(SolverPatch& solverPatch,LimiterPatch& limiterPatch);

  /**
   * Updates the merged limiter status based on the cell-local ADER-DG solution
   * values,
   */
  void updateLimiterStatusAfterSolutionUpdate(SolverPatch& solverPatch,const bool isTroubled);

  /**
   * Deallocates a limiter patch.
   *
   * \note Thread-safe.
   */
  void deallocateLimiterPatch(
      const SolverPatch& solverPatch,
      const int cellDescriptionsIndex) const;

  /**
   * Allocates a new limiter patch,
   * copies geometry information from the solver
   * patch, and projects the solver patch's DG solution
   * onto the FV limiter patch.
   *
   * Adjusts the limiter patch's solution
   * as long if that is requested.
   *
   * \return The index of the patch in the heap
   * vector at address \p cellDescriptionsIndex.
   *
   * \note Thread-safe.
   *
   * TODO(Dominic): A few lines are candidates for
   * being run as background job.
   */
  int allocateLimiterPatch(
      const SolverPatch&  solverPatch,
      const int cellDescriptionsIndex) const;

  /**
   * Project the DG solution onto the FV
   * limiter space and further adjust
   * the so-obtained FV solution if requested
   * by the user.
   */
  void adjustLimiterSolution(
      SolverPatch&  solverPatch,
      LimiterPatch& limiterPatch) const;

  /**
   * Deallocates the limiter patch for solver patches
   * of type which have a limiterStatus and
   * previousLimiterStatus with value 0.
   *
   * Otherwise, a limiter patch holding a valid FV solution
   * might be removed in one of the first iterations but is
   * then required after the limiter status spreading has converged to
   * perform a troubled cell recomputation or is required
   * for performing a global rollback.
   */
  void ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(
      const SolverPatch&  solverPatch,
      const int cellDescriptionsIndex) const;

  /**
   * Update the limiter status based on the cell-local solution values.
   *
   * If the new limiter status is changed to or remains troubled,
   * set the iterationsToCureTroubledCell counter to a certain
   * maximum value.
   * If the limiter status changes from troubled to something else,
   * decrease the iterationsToCureTroubledCell counter.
   * If the counter is set to zero, change a troubled cell
   * to NeighbourOfCellIsTroubled1.
   *
   * \return True if the limiter domain changes irregularly in the cell, i.e.,
   * if a patch with status Ok, NeighbourOfTroubled3, NeighbourOfTroubled4
   * changes its status to Troubled.
   *
   * If the limiter status changes regularly, i.e., from NeighbourOfTroubled1
   * to Troubled or from Troubled to NeighbourOfTroubled3, NeighbourOfTroubled4, this
   * methods returns false.
   */
  MeshUpdateEvent determineRefinementStatusAfterSolutionUpdate(
      SolverPatch& solverPatch,
      const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed);

  /**
   * Takes the FV solution from the limiter patch and projects it on the
   * DG space, overwrites the DG solution on the solver patch with the projected values.
   */
  void projectFVSolutionOnDGSpace(SolverPatch& solverPatch,LimiterPatch& limiterPatch) const;

  /**
   * Takes the DG solution from the solver patch and projects it on the
   * FV space, overwrites the FV solution on the limiter patch with the projected values.
   */
  void projectDGSolutionOnFVSpace(SolverPatch& solverPatch,LimiterPatch& limiterPatch) const;

  /**
   * Ensure that the time step data of the limiter is
   * consistent with the one of the solver.
   */
  void ensureLimiterTimeStepDataIsConsistent() const;

  /**
   * Returns the index of the limiter patch corresponding to
   * the solver patch with index \p solverElement.
   * Both patches link to the same ::LimitingADERDGSolver.
   * If no limiter patch is found, this function
   * returns exahype::solvers::Solver::NotFound.
   */
  int tryGetLimiterElementFromSolverElement(
      const int cellDescriptionsIndex,
      const int solverElement) const {
    SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
    return _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  }

  /**
   * If a limiter patch is allocated for the solver patch,
   * ensure that it's time step data is consistent
   * with the solver patch's time step data.
   */
  void ensureLimiterPatchTimeStepDataIsConsistent(
        const int cellDescriptionsIndex,
        const int solverElement) const;

  void ensureLimiterPatchTimeStepDataIsConsistent(
      const SolverPatch& solverPatch,
      const int cellDescriptionsIndex) const;

  /**
   * Copies the time stamp and the time step sizes from the solver patch
   * to the limiter patch.
   */
  static void copyTimeStepDataFromSolverPatch(
      const SolverPatch& solverPatch, const int cellDescriptionsIndex, const int limiterElement);

  /**
   * Copies the time stamp and the time step sizes from the solver patch
   * to the limiter patch.
   */
  static void copyTimeStepDataFromSolverPatch(
      const SolverPatch& solverPatch, LimiterPatch& limiterPatch);

  /**
   * Body of LimitingADERDGSolver::fusedTimeStep(...).
   */
  UpdateResult fusedTimeStepBody(
      const int   cellDescriptionsIndex,
      const int   element,
      const bool  isFirstIterationOfBatch,
      const bool  isLastIterationOfBatch,
      const bool  isSkeletonJob,
      const bool  mustBeDoneImmediately,
      const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed);

  /**
   * Body of LimitingADERDGSolver::adjustSolutionDuringMeshRefinement(int,int).
   */
  void adjustSolutionDuringMeshRefinementBody(
      SolverPatch& solverPatch,
      const int cellDescriptionsIndex,
      const bool isInitialMeshRefinement);

#ifdef Parallel

  /**
   * Send the solution minimum and maximum values per variable
   * and further the merged limiter status of the solver patch
   * \p element in heap array \p cellDescriptionsIndex to the
   * neighbour.
   *
   * We send the min and max values to all neighbours if the current
   * solver patch holds face data.
   *
   * \see exahype::solvers::Solver::sendDataToNeighbour
   * for a description of the parameters.
   */
  void sendMinAndMaxToNeighbour(
      const int                                    toRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const;

  /**
   * Receive the solution minimum and maximum values per variable
   * and further the merged limiter status for the solver patch
   * \p element in heap array \p cellDescriptionsIndex from the
   * neighbour.
   *
   * We only merge the received values if we consider
   * a solver patch of type Cell. Otherwise,
   * we just drop them.
   *
   * \see exahype::solvers::Solver::sendDataToNeighbour
   * for a description of the parameters.
   */
  void mergeWithNeighbourMinAndMax(
      const int                                    fromRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const;

  /**
   * Single-sided variant of mergeSolutionMinMaxOnFace() that is required
   * for MPI where min and max value are explicitly exchanged through messages.
   */
  void mergeSolutionMinMaxOnFace(
      SolverPatch&  solverPatch,
      const int     faceIndex,
      const double* const min, const double* const  max) const;

#endif

  /**
   * A job which performs a fused limiting ADER-DG time step, i.e., it performs the solution update,
   * updates the local time stamp, performs the space-time predictor commputation,
   * and finally evaluates the a-posteriori limiting criterion.
   *
   * \note Spawning these operations as background job makes only sense if you
   * do not plan to reduce the admissible time step size, refinement requests,
   * or limiter requests within a consequent reduction step.
   *
   * \note We have to copy the neighbourMergePerformed flags of a cell description
   * as they are requrired when determining a new limiter status.
   * We exclude here faces where no merge has been performed, boundary faces e.g.
   *
   * TODO(Dominic): Minimise time step sizes and refinement requests per patch
   * (->transpose the typical minimisation order)
   */
  class FusedTimeStepJob {
  private:
    LimitingADERDGSolver&             _solver;
    const int                         _cellDescriptionsIndex;
    const int                         _element;
    std::bitset<DIMENSIONS_TIMES_TWO> _neighbourMergePerformed;
    const bool                        _isSkeletonJob;
  public:
    FusedTimeStepJob(
        LimitingADERDGSolver&                    solver,
        const int                                cellDescriptionsIndex,
        const int                                element,
        const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed,
        const bool                               isSkeletonJob);

    bool operator()();
  };

  /**
   * A job that calls Solver::adjustSolutionDuringMeshRefinementBody(...).
   */
  class AdjustSolutionDuringMeshRefinementJob {
  private:
    LimitingADERDGSolver& _solver;
    SolverPatch&          _solverPatch;
    const int             _cellDescriptionsIndex;
    const bool            _isInitialMeshRefinement;
  public:
    AdjustSolutionDuringMeshRefinementJob(
        LimitingADERDGSolver& solver,
        SolverPatch&          solverPatch,
        const int             cellDescriptionsIndex,
        const bool            isInitialMeshRefinement);

    bool operator()();
  };

public:

  //Virtual methods to bind the limiter kernels in the abstractLimiterSolver (generated by the Toolkit and application specific)
  virtual void projectOnFVLimiterSpace(const double* const luh, double* const lim) const = 0;
  virtual void projectOnDGSpace(const double* const lim, double* const luh) const = 0;
  virtual bool discreteMaximumPrincipleAndMinAndMaxSearch(const double* const luh, double* boundaryMinPerVariables, double* boundaryMaxPerVariables) = 0;
  virtual void findCellLocalMinAndMax(const double* const luh, double* const localMinPerVariables, double* const localMaxPerVariable) = 0;
  virtual void findCellLocalLimiterMinAndMax(const double* const lim, double* const localMinPerObservable, double* const localMaxPerObservable) = 0;

  /**
   * Create a limiting ADER-DG solver.
   *
   * <h2>Discrete maximum principle</h2>
   * By default this constructor initialises the maximum relaxation
   * parameter to the value to \f$ \delta_0 = 1\cdot 10^{-4} \f$
   * and the difference scaling parameter to \f$ \epsilon = 1\cdot 10^{-3} \f$.
   * See Dumbser et al., 2014. doi:10.1016/j.jcp.2014.08.009 for more details on
   * the notation.
   */
  LimitingADERDGSolver(
      const std::string& identifier,
      exahype::solvers::ADERDGSolver* solver,
      exahype::solvers::FiniteVolumesSolver* limiter,
      const double DMPRelaxationParameter=1e-4,
      const double DMPDifferenceScaling=1e-3,
      const int iterationsToCureTroubledCell=2);

  virtual ~LimitingADERDGSolver() {
    _solver.reset();
    _limiter.reset();
  }

  // Disallow copy and assignment
  LimitingADERDGSolver(const ADERDGSolver& other) = delete;
  LimitingADERDGSolver& operator=(const ADERDGSolver& other) = delete;
  
  /** Wire through to ADER-DG Solver */
  MeshUpdateEvent getNextMeshUpdateEvent() const final override;
  void updateNextMeshUpdateEvent(MeshUpdateEvent meshUpdateEvent) final override;
  void setNextMeshUpdateEvent() final override;
  MeshUpdateEvent getMeshUpdateEvent() const final override;
  void overwriteMeshUpdateEvent(MeshUpdateEvent newMeshUpdateEvent) final override;

  /*
   * A time stamp minimised over all the ADERDG and FV solver
   * patches.
   */
  double getMinTimeStamp() const final override;

  /**
   * Run over all solvers and identify the minimal time step size.
   */
  double getMinTimeStepSize() const final override;

  double getMinNextTimeStepSize() const final override;

  void updateMinNextTimeStepSize(double value) final override;

  /**
   * \copydoc ::exahype::solvers::Solver::initSolver
   *
   * Additionally, set the
   * ::_limiterDomainChangedIrregularly flag to true
   * since the limiter mappings all check this flag
   * in order to distinguish between solvers in a multisolver
   * run.
   */
  void initSolver(
      const double timeStamp,
      const tarch::la::Vector<DIMENSIONS,double>& domainOffset,
      const tarch::la::Vector<DIMENSIONS,double>& domainSize,
      const tarch::la::Vector<DIMENSIONS,double>& boundingBoxSize,
      const std::vector<std::string>& cmdlineargs,
      const exahype::parser::ParserView& parserView) override;

  bool isPerformingPrediction(const exahype::State::AlgorithmSection& section) const final override;
  bool isMergingMetadata(const exahype::State::AlgorithmSection& section) const final override;

  void synchroniseTimeStepping(
          const int cellDescriptionsIndex,
          const int element) const final override;

  void synchroniseTimeStepping(
      SolverPatch& solverPatch,
      const int cellDescriptionsIndex) const;

  /**
   * We always override the limiter time step
   * data by the ADER-DG one before a solution update.
   */
  void startNewTimeStep() final override;

  void startNewTimeStepFused(
      const bool isFirstIterationOfBatch,
      const bool isLastIterationOfBatch) final override;

  void updateTimeStepSizes() final override;

  void updateTimeStepSizesFused() final override;

  void zeroTimeStepSizes() final override;

  /**
   * Roll back the time step data to the
   * ones of the previous time step.
   */
  void rollbackToPreviousTimeStep() final override;

  /**
   * Same as LimitingADERDGSolver::rollbackToPreviousTimeStep
   * but for the fused time stepping scheme.
   */
  void rollbackToPreviousTimeStepFused() final override;

  void updateNextMaxLevel(int maxLevel) final override;
  int getNextMaxLevel() const final override;
  int getMaxLevel() const final override;

  /**
   * Returns the index of the solver patch registered for the solver with
   * index \p solverNumber in exahype::solvers::RegisteredSolvers.
   * If no patch is registered for the solver, this function
   * returns exahype::solvers::Solver::NotFound.
   */
  int tryGetElement(
      const int cellDescriptionsIndex,
      const int solverNumber) const final override {
    return _solver->tryGetElement(cellDescriptionsIndex,solverNumber);
  }

  /**
   * Returns the index of the limiter patch registered for the solver with
   * index \p solverNumber in exahype::solvers::RegisteredSolvers.
   * If no patch is registered for the solver, this function
   * returns exahype::solvers::Solver::NotFound.
   */
  int tryGetLimiterElement(
      const int cellDescriptionsIndex,
      const int solverNumber) const {
    return _limiter->tryGetElement(cellDescriptionsIndex,solverNumber);
  }

  /**
   * Looks up the limiter patch for the given solver patch.
   *
   * Further copies the time step sizes from the solver patch
   * to the limiter patch.
   *
   * Assumes that \p solverPatch has a limiter status > 0 and thus
   * has a limiter patch assigned.
   */
  LimiterPatch& getLimiterPatchForSolverPatch(const int cellDescriptionsIndex, const SolverPatch& solverPatch) const;

  /**
   * Similar to ::getLimiterPatchforSolverPatch but does not lookup the
   * \p limiterElement by itself.
   *
   * Assumes that \p limiterElement is valid.
   */
  LimiterPatch& getLimiterPatchForSolverPatch(
      const SolverPatch& solverPatch, const int cellDescriptionsIndex, const int limiterElement) const;

  ///////////////////////////////////
  // MODIFY CELL DESCRIPTION
  ///////////////////////////////////

  /**
   * TODO(Dominic): Update docu.
   *
   * Based on the limiter status of a solver patch
   * and the solver patch's type, we perform the
   * following actions:
   *
   * | New Status | Type                        | Action                                                                                           |
   * ----------------------------------------------------------------------------------------------------------------------------------------------|
   * | O/NNT      | Any                         | Do nothing.                                                                                      |
   * | T/NT       | Cell                        | Set RefinementRequested event on parent cell if its current event is None or AugmentingRequested |
   * | T/NT       | Descendant                  | Set RefinementRequested event if current event is None or AugmentingRequested                    |
   * | T/NT       | Else                        | Do nothing                                                                                       |
   *
   * \note Currently we assume that the problem and load-balancing is so well-behaved that
   * we always find a Cell as parent of a Descendant on the same MPI rank. We further do not
   * consider Master-Worker boundaries in the lookup of the parent.
   *
   * Legend: O: Ok, T: Troubled, NT: NeighbourOfTroubled1..2, NNT: NeighbourOfTroubled3..4
   */
  void markForRefinementBasedOnLimiterStatus(SolverPatch& solverPatch) const;

  /**
   * Update the cellwise limiter status using the facewise limiter status
   * values.
   *
   * Ensure a limiter patch is allocated on all cells
   * on the finest level of the grid which have
   * a limiter status greater than 0.
   *
   * Deallocate the limiter patch on helper cells.
   * Deallocate the limiter patch also on compute cells if
   * both, the current and previous, limiter status flags
   * are equalling Ok (0).
   * Otherwise, the limiter patch is kept since
   * it might turn out that it is still required in
   * one of the next limiter status spreading or
   * mesh refinement iterations.
   *
   * \note We overwrite the facewise limiter status values with the new value
   * in order to use the updateLimiterStatusAfterSetInitialConditions function
   * afterwards which calls determineLimiterStatus(...) again.
   *
   * \note calls synchroniseTimeStepping
   *
   * returns true if a new limiter patch was allocated.
   */
  MeshUpdateEvent  updateRefinementStatusDuringRefinementStatusSpreading(SolverPatch& solverPatch) const;

  /**\copydoc exahype::solvers::Solver::progressMeshRefinementInEnterCell
   *
   * <h2>Limiter status-based refinement</h2>
   *
   * We do not mark for refinement based on the limiter status
   * as long as the refinement request of the cell description
   * is still in the pending state. In this case,
   * the physical admissibility detection and the user's refinement
   * criterion have not been evaluated.
   */
  bool progressMeshRefinementInEnterCell(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const int  solverNumber,
      const bool stillInRefiningMode) override;

 bool progressMeshRefinementInLeaveCell(
     exahype::Cell& fineGridCell,
     exahype::Vertex* const fineGridVertices,
     const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
     exahype::Cell& coarseGridCell,
     const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
     const int solverNumber) override;

 exahype::solvers::Solver::RefinementControl eraseOrRefineAdjacentVertices(
     const int cellDescriptionsIndex,
     const int solverNumber,
     const tarch::la::Vector<DIMENSIONS, double>& cellOffset,
     const tarch::la::Vector<DIMENSIONS, double>& cellSize,
     const bool checkThoroughly) const final override;

  bool attainedStableState(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      const int solverNumber) const final override;

  void finaliseStateUpdates(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) final override;

  ///////////////////////////////////
  // CELL-LOCAL
  //////////////////////////////////
  double startNewTimeStep(
      const int cellDescriptionsIndex,
      const int element) final override;

  double startNewTimeStepFused(
        const int cellDescriptionsIndex,
        const int element,
        const bool isFirstIterationOfBatch,
        const bool isLastIterationOfBatch) final override;

  double updateTimeStepSizesFused(
      const int cellDescriptionsIndex,
      const int element) final override;

  double updateTimeStepSizes(
        const int cellDescriptionsIndex,
        const int element) final override;

  void zeroTimeStepSizes(
      SolverPatch& solverPatch,
      const int cellDescriptionsIndex) const;

 /**
   * Rollback to the previous time step, i.e,
   * overwrite the time step size and time stamp
   * fields of the solver and limiter patches
   * by the values used in the previous iteration.
   *
   * TODO(Dominic): Move into solver
   */
  void rollbackToPreviousTimeStep(
      const int cellDescriptionsIndex,
      const int solverElement) const;

  /*
   * Same as LimitingADERDGSolver::rollbackToPreviousTimeStep
   * but for the fused time stepping scheme.
   */
  void rollbackToPreviousTimeStepFused(
      const int cellDescriptionsIndex,
      const int solverElement) const;

  UpdateResult fusedTimeStep(
      const int cellDescriptionsIndex,
      const int element,
      const bool isFirstIterationOfBatch,
      const bool isLastIterationOfBatch,
      const bool isAtRemoteBoundary) final override;

  UpdateResult update(
        const int cellDescriptionsIndex,
        const int element,
        const bool isAtRemoteBoundary) final override;

  void compress(
      const int cellDescriptionsIndex,
      const int element,
      const bool isAtRemoteBoundary) const final override;

  void adjustSolutionDuringMeshRefinement(
      const int cellDescriptionsIndex,const int element) final override;

  /**
   * Update the solution of a solver patch and or
   * its associated limiter patch
   *
   * This method assumes the ADERDG solver's cell-local limiter status has
   * already been determined.
   *
   * Before performing an update with the limiter,
   * set the ADER-DG time step sizes for the limiter patch.
   * (ADER-DG is always dictating the time step sizes.)
   *
   * \see determineLimiterStatusAfterRefinementStatusSpreading(...)
   *
   * \note Make sure to reset neighbour merge
   * helper variables in this method call.
   *
   * \note Has no const modifier since kernels are not const functions yet.
   *
   * \param[in] backupPreviousSolution Set to true if the solution should be backed up before
   *                                   we overwrite it by the updated solution.
   *
   */
  void updateSolution(
      const int cellDescriptionsIndex,
      const int element,
      const bool backupPreviousSolution);

  /**
   * Determine the new cell-local min max values.
   *
   * Must be invoked after ::determineLimiterStatusAfterSolutionUpdate.
   *
   * TODO(Dominic): Tobias's integer
   * flagging idea might reduce complexity here
   */
  void determineMinAndMax(
      const int cellDescriptionsIndex,
      const int element);

  /**
   * Evaluates a discrete maximum principle (DMP) and
   * the physical admissibility detection (PAD) criterion for
   * the solution values stored for any solver patch
   * that is of type Cell independent of the mesh level
   * it is located at.
   * This method then invokes
   * ::determinLimiterStatusAfterSolutionUpdate(SolverPatch&,const bool)
   * with the result of these checks.
   *
   * For solver patches of a type other than Cell,
   * we simply update the limiter status using
   * the information taken from the neighbour
   * merging.
   *
   * \note Must be called after starting a new time step for the patch.
   */
  MeshUpdateEvent
  updateRefinementStatusAndMinAndMaxAfterSolutionUpdate(
      const int cellDescriptionsIndex,
      const int element,
      const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed);

  /**
   * Similar to ::determineLimiterStatusAfterSolutionUpdate(const int,const int)
   * Does only evaluate the physical admissibility detection (PAD) but not the
   * discrete maximum principle (DMP).
   *
   * \note We overwrite the facewise limiter status values with the new value
   * in order to reuse the determineLimiterStatusAfterSolutionUpdate function
   * which calls determineLimiterStatus(...) again.
   */
  void updateLimiterStatusAndMinAndMaxAfterAdjustSolution(
      const int cellDescriptionsIndex,
      const int element);

  /*
   * Deallocate the limiter patch on all AMR related
   * helper cells.
   *
   * It is safe to use this method during
   * the mesh refinement iterations.
   *
   * \note Thread-safe.
   */
   void ensureNoLimiterPatchIsAllocatedOnHelperCell(
       const SolverPatch& solverPatch,
       const int cellDescriptionsIndex) const;

   /*
    * Ensures that a limiter patch is allocated
    * on all compute cells (Cell) on the finest mesh
    * level that are flagged
    * with a limiter status other than Ok.
    *
    * It is safe to use this method during
    * the mesh refinement iterations.
    */
   bool ensureRequiredLimiterPatchIsAllocated(
           const SolverPatch& solverPatch,
           const int cellDescriptionsIndex,
           const int limiterStatus) const;


   /**
    * Go back to previous time step with
    * time step data and solution.
    *
    * Keep the new refinement status.
    *
    * Allocate necessary new limiter patches.
    */
   void rollbackSolutionGlobally(
          const int cellDescriptionsIndex,
          const int element,
          const bool fusedTimeStepping) const final override;

  /**
   * Reinitialises cells that have been subject to a limiter status change.
   * This method is invoked (during and??) after the limiter status spreading.
   *
   * The method has to take into account which solution, the solver's
   * or the limiter's, was populated with valid solution values
   * in the last iteration. The action of this method is
   * thus based on the new and old limiter status.
   *
   * We perform the following actions based on the
   * old and new limiter status:
   *
   * | New Status | Old Status | Action                                                                                                                                        |
   * ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
   * | O          | O          | Do nothing.                                                                                                                                   |
   * | ~          | T/NT/NNT   | Remove the limiter patch.                                                                                                                     |
   * | T/NT/NNT   | T/NT       | Roll back the limiter solution.                                                                                                               |
   * | ~          | O/NNT      | Roll back the solver solution. Initialise limiter patch if necessary. Project (old, valid) solver solution onto the limiter's solution space. |
   *
   * Legend: O: Ok, T: Troubled, NT: NeighbourOfTroubled1..2, NNT: NeighbourOfTroubled3..4
   *
   * We do not overwrite the old limiter status set in this method.
   * We compute the new limiter status based on the merged limiter statuses associated
   * with the faces.
   *
   * <h2>A-posteriori refinement</h2>
   * In case of a-posteriori refinement, we also perform a rollback
   * in the Ok cells. Then, the global time step size used by the predictor
   * is not valid anymore (assumption: global time stepping)
   *  and the last solution update has to be redone.
   *
   * <h2>Compute cell limiter patch deallocation</h2>
   * It is only safe to deallocate unrequired compute cell limiter patches after
   * the mesh refinement iterations since we might throw away valid
   * FV values during the first iterations. However, then find out later
   * that we need them after the limiter status diffusion
   * has converged.
   * Helper cell limiter patches can be deallocated during
   * the mesh refinement iterations.
   */
  void rollbackSolutionLocally(
      const int cellDescriptionsIndex,
      const int element,
      const bool fusedTimeStepping) const;

  /**
   * Recompute the solution in cells that have been subject to a limiter status change
   * This method is invoked after the solver reinitialisation
   * (see exahype::solvers::LimitingADERDGSolver::reinitialiseSolvers).
   *
   * It evolves the solution of the solver and limiter in the reinitialised cells to the
   * correct time stamp.
   *
   * We perform the following actions based on the
   * new limiter status:
   *
   * |New Status | Action                                                                                                                                      |
   * ----------------------------------------------------------------------------------------------------------------------------------------------------------|
   * |O          | Do nothing. Solver solution has been evolved correctly before. DG solution is correct.                                                                             |
   * |T/NT       | Evolve FV solver. Project result onto the ADER-DG space. Recompute the space-time predictor if not initial recomputation.                                                                                  |
   * |NNT        | Evolve solver and project its solution onto the limiter solution space. (We had to do a rollback beforehand in the reinitialisation phase.) |
   *
   * Legend: O: Ok, T: Troubled, NT: NeighbourOfTroubled1..2, NNT: NeighbourOfTroubled3..4
   */
  void recomputeSolutionLocally(
      const int cellDescriptionsIndex,
      const int element);

  /**
   * !!! Only for fused time stepping !!!
   *
   * Recompute the predictor in particular cells.
   *
   * We perform the following actions based on the
   * new and old limiter status:
   *
   * |New Status | Old Status | Action                                                                                                                                      |
   * ----------------------------------------------------------------------------------------------------------------------------------------------------------|
   * |O/NNTT     | O/NT/NNT   | Do nothing. Underlying DG solution has not changed.
   * |           | T          | Recompute predictor. Cell was skipped before in predictor computation
   * |           |            | since it is marked T. See mapping Prediction::enterCell.
   * |NT         | *          | Recompute predictor. DG solution has been recomputed.
   * |T          | *          | Not necesssary to compute predictor
   *
   * Legend: O: Ok, T: Troubled, NT: NeighbourIsTroubledCell, NNT: NeighbourIsNeighbourOfTroubledCell
   *
   * We do not overwrite the old limiter status set in this method.
   * We compute the new limiter status based on the merged limiter statuses associated
   * with the faces.
   *
   * \param[in] isAtRemoteBoundary Flag indicating that the cell hosting the
   *                                   cell description is adjacent to a remote rank.
   */
  void recomputePredictorLocally(
      const int cellDescriptionsIndex,
      const int element,
      const bool isAtRemoteBoundary);

  void prolongateFaceData(
      const int cellDescriptionsIndex,
      const int element) final override;

  void restriction(
        const int cellDescriptionsIndex,
        const int element) final override;

  ///////////////////////////////////
  // NEIGHBOUR
  ///////////////////////////////////
  void mergeNeighboursMetadata(
      const int                                 cellDescriptionsIndex1,
      const int                                 element1,
      const int                                 cellDescriptionsIndex2,
      const int                                 element2,
      const tarch::la::Vector<DIMENSIONS, int>& pos1,
      const tarch::la::Vector<DIMENSIONS, int>& pos2) const final override;

  void mergeNeighbours(
      const int                                 cellDescriptionsIndex1,
      const int                                 element1,
      const int                                 cellDescriptionsIndex2,
      const int                                 element2,
      const tarch::la::Vector<DIMENSIONS, int>& pos1,
      const tarch::la::Vector<DIMENSIONS, int>& pos2) final override;

  /**
   * Merge solver boundary data (and other values) of two adjacent
   * cells based on their limiter status.
   *
   * The solver involved in the neighbour merge
   * is selected according to the following scheme:
   *
   * | Status 1 | Status 2 | Solver to Merge
   * ---------------------------------------
   * | O        | O        | ADER-DG       |
   * | O        | NNT      | ADER-DG       |// O|NNT x O|NNT
   * | NNT      | O        | ADER-DG       |
   * | NNT      | NNT      | ADER-DG       |
   *
   * | NNT      | NT       | FV            |
   * | NT       | NNT      | FV            | // NT&NNT | N&NNT
   *
   * | NT       | NT       | FV            |
   * | NT       | T        | FV            |
   * | T        | NT       | FV            | // T|NT x T|NT
   * | T        | T        | FV            |
   *
   * Legend: O: Ok, T: Troubled, NT: NeighbourIsTroubledCell, NNT: NeighbourIsNeighbourOfTroubledCell
   *
   * <h2>Solution Recomputation</h2>
   * If we perform a solution recomputation, we do not need to perform
   * a solution update in the non-troubled solver patches and the patches
   * with status neighbourIsNeighbourOfTroubledCell. Instead we
   * can simply reuse the already computed update to advance in
   * time to the desired time stamp.
   *
   * We thus do not need to merge these patches with neighbour data
   * in the recomputation phase.
   *
   * \note Limiting is only performed on the finest level
   * of the mesh. The other levels work only with the ADER-DG
   * solver.
   *
   * TODO(Dominic): Remove limiterstatus1 and limiterStatus2 argument.
   * They depend on the isRecomputation value
   */
  void mergeNeighboursBasedOnLimiterStatus(
      const int                                 cellDescriptionsIndex1,
      const int                                 element1,
      const int                                 cellDescriptionsIndex2,
      const int                                 element2,
      const tarch::la::Vector<DIMENSIONS, int>& pos1,
      const tarch::la::Vector<DIMENSIONS, int>& pos2,
      const bool                                isRecomputation) const;

  /**
   * Merges the min max of two neighbours sharing a face.
   *
   * This method is used to detect cells that are
   * troubled after the imposition of initial conditions.
   */
  void mergeSolutionMinMaxOnFace(
      const int                                 cellDescriptionsIndex1,
      const int                                 element1,
      const int                                 cellDescriptionsIndex2,
      const int                                 element2,
      const tarch::la::Vector<DIMENSIONS, int>& pos1,
      const tarch::la::Vector<DIMENSIONS, int>& pos2) const;

  /**
   * Depending on the limiter status, we impose boundary conditions
   * onto the solution of the solver or of the limiter.
   */
  void mergeWithBoundaryData(
      const int                                 cellDescriptionsIndex,
      const int                                 element,
      const tarch::la::Vector<DIMENSIONS, int>& posCell,
      const tarch::la::Vector<DIMENSIONS, int>& posBoundary) final override;

  /**
   * Merge solver boundary data (and other values) of a
   * cell with the boundary conditions based on the cell's
   * limiter status.
   *
   * The solver involved in the merge
   * is selected according to the following scheme:
   *
   * | Status   | Solver to Merge |
   * ------------------------------
   * | O        | ADER-DG         |
   * | NNT      | ADER-DG         |
   *
   * | NT       | FV              |
   * | T        | FV              |
   *
   * Legend: O: Ok, T: Troubled, NT: NeighbourIsTroubledCell, NNT: NeighbourIsNeighbourOfTroubledCell
   *
   * <h2>Solution Recomputation</h2>
   * If we perform a solution recomputation, we do not need to perform
   * a solution update in the non-troubled solver patches and the patches
   * with status neighbourIsNeighbourOfTroubledCell. Instead we
   * can simply reuse the already computed update to advance in
   * time to the desired time stamp.
   *
   * We thus do not need to merge these patches with boundary data
   * in the recomputation phase.
   *
   * \param[in] isRecomputation Flag indicating if this merge is part of a solution recomputation phase.
   */
  void mergeWithBoundaryDataBasedOnLimiterStatus(
        const int                                 cellDescriptionsIndex,
        const int                                 element,
        const int                                 limiterStatusAsInt,
        const tarch::la::Vector<DIMENSIONS, int>& posCell,
        const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
        const bool                                isRecomputation);

#ifdef Parallel
  ///////////////////////////////////
  // NEIGHBOUR
  ///////////////////////////////////
  void appendNeighbourCommunicationMetadata(
      exahype::MetadataHeap::HeapEntries& metadata,
      const tarch::la::Vector<DIMENSIONS,int>& src,
      const tarch::la::Vector<DIMENSIONS,int>& dest,
      const int cellDescriptionsIndex,
      const int solverNumber) const final override;

  void mergeWithNeighbourMetadata(
      const exahype::MetadataHeap::HeapEntries& neighbourMetadata,
      const tarch::la::Vector<DIMENSIONS, int>& src,
      const tarch::la::Vector<DIMENSIONS, int>& dest,
      const int                                 cellDescriptionsIndex,
      const int                                 element) const final override;

  void sendDataToNeighbour(
      const int                                     toRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) final override;

  /**
   * Send data or empty data to the neighbour data based
   * on the limiter status.
   *
   * \param[in] isRecomputation Indicates if this called within a local recomputation
   *                            process.
   * \param[in] limiterStatus   This method is used in two different contexts (see \p isRecomputation)
   *                            which either make use of the unified face-wise limiter status (isRecomputation)
   *                            or make use of the cell-wise limiter status (!isRecomputation).
   *
   * <h2>Local Recomputation</h2>
   * In case this method is called during a local recomputation,
   * the ADER-DG solver does only send empty messages to the neighbour.
   * Otherwise it merges the received data and adds it to the update.
   *
   * <h2>Fused Time Stepping</h2>
   * If the limiter status of a cell changes dramatically, a cell might have
   * not allocated a limiter patch yet when fused time stepping is used.
   * In this case, we send out an empty FV boundary data message.
   */
  void sendDataToNeighbourBasedOnLimiterStatus(
        const int                                    toRank,
        const int                                    cellDescriptionsIndex,
        const int                                    element,
        const tarch::la::Vector<DIMENSIONS, int>&    src,
        const tarch::la::Vector<DIMENSIONS, int>&    dest,
        const bool                                   isRecomputation,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const int                                    level) const;

  void sendEmptyDataToNeighbour(
      const int                                     toRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) const final override;

  void mergeWithNeighbourData(
      const int                                    fromRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) final override;

  /**
   * Merge or drop received neighbour data based
   * on the limiter status.
   *
   * \param isRecomputation Indicates if this called within a solution recomputation
   *                        process.
   * \param limiterStatus   This method is used in two different contexts (see \p isRecomputation)
   *                        which either make use of the unified face-wise limiter status (isRecomputation)
   *                        or make use of the cell-wise limiter status (!isRecomputation).
   *
   * <h2>SolutionRecomputation</h2>
   * In case this method is called within a solution recomputation,
   * the ADER-DG solver drops the received boundary data.
   * Otherwise it merges the received data and adds it to the update.
   *
   *  \note This method assumes that there has been a unified face-wise limiter status value
   *  determined and written back to the faces.
   */
  void mergeWithNeighbourDataBasedOnLimiterStatus(
      const int                                    fromRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const bool                                   isRecomputation,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const;

  void dropNeighbourData(
      const int                                     fromRank,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) const final override;


  ///////////////////////////////////////
  // NEIGHBOUR - Solution Recomputation
  ///////////////////////////////////////
  /**
   * TODO(Dominic):
   * Add more docu.
   *
   * We do not send the minimum and maximum solution values
   * during the solution recomputation process.
   */
  void sendEmptySolverAndLimiterDataToNeighbour(
      const int                                     toRank,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) const;

  /**
   * TODO(Dominic):
   * Add more docu.
   *
   * We do not send the minimum and maximum solution values
   * during the solution recomputation process.
   */
  void dropNeighbourSolverAndLimiterData(
        const int                                     fromRank,
        const tarch::la::Vector<DIMENSIONS, int>&     src,
        const tarch::la::Vector<DIMENSIONS, int>&     dest,
        const tarch::la::Vector<DIMENSIONS, double>&  x,
        const int                                     level) const;

  /////////////////////////////////////
  // MASTER<=>WORKER
  /////////////////////////////////////
  /**
   * Kind of similar to progressMeshRefinementInPrepareSendToWorker
   * but performs a few additional operations in order to
   * notify the worker about some coarse grid operations only
   * the master knows.
   *
   * \note This function sends out MPI messages.
   */
  void progressMeshRefinementInPrepareSendToWorker(
      const int workerRank,
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const int solverNumber) final override;

  /**
   * Forward to ADERDGSolver
   */
  void sendDataToWorkerIfProlongating(
      const int                                     toRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) const final override;

  /**
   * Forward to ADERDGSolver
   */
  void receiveDataFromMasterIfProlongating(
      const int masterRank,
      const int receivedCellDescriptionsIndex,
      const int receivedElement,
      const tarch::la::Vector<DIMENSIONS,double>& x,
      const int level) const final override;

  /**
   * Finish prolongation operations started on the master.
   *
   * TODO(Dominic): No const modifier const as kernels are not const yet
   */
  bool progressMeshRefinementInMergeWithWorker(
      const int localCellDescriptionsIndex,
      const int receivedCellDescriptionsIndex, const int receivedElement) final override;

  /**
   * Finish erasing operations on the worker side and
   * send data up to the master if necessary.
   * This data is then picked up to finish restriction
   * operations.
   */
  void progressMeshRefinementInPrepareSendToMaster(
      const int masterRank,
      const int cellDescriptionsIndex, const int element,
      const tarch::la::Vector<DIMENSIONS,double>& x,
      const int level) const final override;

  /**
   * Finish prolongation operations started on the master.
   *
   * \return If we the solver requires master worker communication
   * at this cell
   *
   * TODO(Dominic): No const modifier const as kernels are not const yet
   */
  bool progressMeshRefinementInMergeWithMaster(
      const int worker,
      const int localCellDescriptionsIndex,
      const int localElement,
      const int coarseGridCellDescriptionsIndex,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int  level,
      const bool stillInRefiningMode) final override;

  void appendMasterWorkerCommunicationMetadata(
      exahype::MetadataHeap::HeapEntries& metadata,
      const int cellDescriptionsIndex,
      const int solverNumber) const final override;

  void sendDataToWorkerOrMasterDueToForkOrJoin(
      const int                                     toRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const peano::heap::MessageType&               messageType,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) const override;

  void mergeWithWorkerOrMasterDataDueToForkOrJoin(
      const int                                     fromRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const peano::heap::MessageType&               messageType,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) const override;

  ///////////////////////////////////
  // WORKER->MASTER
  ///////////////////////////////////
  void sendDataToMaster(
      const int                                    masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const final override;

  void mergeWithWorkerData(
      const int                                    workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) final override;

  ///////////////////////////////////
  // MASTER->WORKER
  ///////////////////////////////////
  void sendDataToWorker(
      const                                        int workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const final override;

  void mergeWithMasterData(
      const                                        int masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) final override;

#endif

  std::string toString() const final override;

  void toString (std::ostream& out) const final override;

  const std::unique_ptr<exahype::solvers::FiniteVolumesSolver>&
  getLimiter () const {
    return _limiter;
  }

  const std::unique_ptr<exahype::solvers::ADERDGSolver>&
  getSolver () const {
    return _solver;
  }
};


#endif /* LIMITEDADERDGSOLVER_H_ */
