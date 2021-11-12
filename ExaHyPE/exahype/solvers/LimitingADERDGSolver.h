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

#include "exahype/solvers/Solver.h"
#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

#include "tarch/multicore/BooleanSemaphore.h"

#include "exahype/plotters/Plotter.h"

#include "exahype/profilers/simple/NoOpProfiler.h"

#ifdef USE_ITAC
#include "VT.h"
#endif

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

public:
  #ifdef USE_ITAC
  /**
   * These handles are used to trace solver events with Intel Trace Analyzer and Collector.
   */
  static int adjustSolutionHandle;
  static int fusedTimeStepBodyHandle;
  static int fusedTimeStepBodyHandleSkeleton;
  static int updateBodyHandle;
  static int mergeNeighboursHandle;
  #endif
  
  /**
   * The ADERDG solver.
   */
  std::unique_ptr<exahype::solvers::ADERDGSolver> _solver;

protected:

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
  //  void setSolutionMinMax(double* const min, double* const max) const;

  /**
   * Merge the solution min and max values on a face between two cell
   * descriptions. Signature is similar to that of the solver of a Riemann problem.
   */
  void mergeSolutionMinMaxOnFace(
      SolverPatch& solverPatch1,
      SolverPatch& solverPatch2,
      Solver::InterfaceInfo& face) const;

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
   * @note This method assumes the (previous) refinement status
   * was not modified yet by another routine.
   */
  bool evaluatePhysicalAdmissibilityCriterion(SolverPatch& solverPatch,const double t);

  /**
   * Computes the cell-local minimum and maximum
   * of the solution (per variable)
   * and makes them accessible per face.
   */
  void determineSolverMinAndMax(SolverPatch& solverPatch,const bool validate);

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
   * Deallocates a limiter patch.
   *
   * @note Thread-safe.
   */
  void deallocateLimiterPatch(const SolverPatch& solverPatch,CellInfo& cellInfo) const;

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
   * @note Thread-safe.
   *
   * TODO(Dominic): A few lines are candidates for
   * being run as background job.
   */
  void allocateLimiterPatch(const SolverPatch&  solverPatch,CellInfo& cellInfo) const;

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
      const SolverPatch&  solverPatch,CellInfo& cellInfo) const;

  /**
   * Update the limiter status based on the cell-local solution values.
   *
   * \return MeshUpdateEvent::IrregularLimiterDomainChange if the limiter domain changes irregularly on the finest grid level, i.e.,
   * if a cell outside of the buffer layers is set to troubled.
   * Returns MeshUpdateEvent::RefinementRequested if a coarse grid cell is flagged as a troubled or
   * if feature-based refinement was requested by the user.
   * Returns MeshUpdateEvent::None in all other cases.
   */
  MeshUpdateEvent updateRefinementStatusAfterSolutionUpdate(
      SolverPatch& solverPatch,
      CellInfo&    cellInfo,
      const bool   isTroubled);
  
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
   * Compute a new admissible time step size.
   * Depending on the refinement status, either the
   * limiter or solver solution's eigenvalues are
   * used.
   *
   * @param solverPatch solver patch linking to the data
   * @param cellInfo    refers to all cell descriptions/patches found in the cell holding the solver patch
   *
   * @note Doesn't have const modifier as kernels are not const yet.
   *
   * @return an admissible time step size.
   */
  double computeTimeStepSize(SolverPatch& solverPatch,CellInfo& cellInfo);

  /**
   * If a limiter patch is allocated for the solver patch,
   * ensure that it's time step data is consistent
   * with the solver patch's time step data.
   */
  void ensureLimiterPatchTimeStepDataIsConsistent(
      const SolverPatch& solverPatch,CellInfo& cellInfo) const;

  /**
   * Uncompress solver and limiter degrees of freedom.
   */
  void uncompress(
      SolverPatch& solverPatch,
      CellInfo& cellInfo) const;

  /**
   * Compress solver and limiter degrees of freedom.
   */
  void compress(
      SolverPatch& solverPatch,
      CellInfo& cellInfo,
      const bool isAtRemoteBoundary) const;

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
   * Calls the updateGlobalObservables on the solver if the
   * cell is not troubled. Otherwise, calls the routine
   * on the limiter.
   *
   * @param solverPatch a solver patch of type Cell.
   * @param cellInfo    struct referring to all cell descriptions registered for a cell
   */
  void reduce(
      const SolverPatch&  solverPatch,
      CellInfo&           cellInfo,
      const UpdateResult& updateResult);

  /**
   * Body of LimitingADERDGSolver::fusedTimeStepOrRestrict(...).
   *
   * <h2> Order of operations</h2>
   * Data stored on a patch must be compressed by the last operation touching
   * the patch. If we spawn the prediction as background job, it is very likely
   * that it is executed last. In order to have a deterministic order of
   * operations, we thus always run the prediction last.
   *
   * This decision implies that the time step data is updated before running the prediction.
   * We thus need to memorise the prediction time stamp and time step size before performing
   * the time step update. Fortunately, it is already memorised as it is copied
   * into the correction time step data fields of the patch
   * after the time step data update.
   *
   * @param solverPatch             an ADER-DG cell description of type Cell
   * @param cellInfo                struct referring to all cell descriptions registered for a cell
   * @param isFirstTimeStepOfBatch  if this the first time step in a batch (at spawn time if run by job)
   * @param isLastTimeStepOfBatch   if this the last time step in a batch  (at spawn time if run by job)
   * @param predictionTimeStamp     the time stamp which should be used for the prediction (at spawn time if run by job)
   * @param predictionTimeStepSize  the time step size which should be used for the prediction (at spawn time if run by job)
   * @param isSkeletonCell          if this cell description belongs to the MPI or AMR skeleton.
   * @param mustBeDoneImmediately   if the prediction has to be performed immediately and cannot be spawned as background job
   *
   * @note Might be called by background task. Do not synchronise time step data here.
   */
  void fusedTimeStepBody(
      SolverPatch&                                       solverPatch,
      CellInfo&                                          cellInfo,
      const double                                       predictionTimeStamp,
      const double                                       predictionTimeStepSize,
      const bool                                         isFirstTimeStepOfBatch,
      const bool                                         isLastTimeStepOfBatch,
      const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers,
      const bool                                         isSkeletonCell,
      const bool                                         mustBeDoneImmediately);

  /**
   * @return Computes a merged limiter status as a maximum of
   * the current cell value and the neighbours' values decremented by 1.
   *
   * @note This method is required as we cure troubled cells
   * on the fly by lowering their limiter status by 2.
   * If the
   *
   * @param[in] solverPatch A solver patch.
   */
  static int computeMergedRefinementStatus(const SolverPatch& solverPatch);

  /**
   * Body of LimitingADERDGSolver::updateOrRestrict(...).
   *
   * @param solverPatch        a solver patch which may or may not have an associated limiter patchs
   * @param cellInfo           links to all solver and limiter patches registered for a cell
   *
   *
   * @return an admissible time step size and a mesh update event for the solver patch
   */
  void updateBody(
      SolverPatch&                                       solverPatch,
      CellInfo&                                          cellInfo,
      const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers);

 /**
   * Rollback to the previous time step, i.e,
   * overwrite the time step size and time stamp
   * fields of the solver and limiter patches
   * by the values used in the previous iteration.
   */
  void rollbackToPreviousTimeStep(SolverPatch& solverPatch,CellInfo& cellInfo) const;

  /**
   * Body of LimitingADERDGSolver::adjustSolutionDuringMeshRefinement(int,int).
   *
   * @note May be called from background task. Do not synchronise time step data here.
   */
  void adjustSolutionDuringMeshRefinementBody(
      SolverPatch& solverPatch,
      CellInfo& cellInfo,
      const bool isInitialMeshRefinement);

  /**
   * @return if the given solver patch is involved in a local
   * recomputation. (Makes only sense to call this function if
   * one is currently ongoing.)
   *
   * @param solverPatch a solver patch.
   */
  bool isInvolvedInLocalRecomputation(const SolverPatch& solverPatch) const {
    return
        solverPatch.getType()==SolverPatch::Type::Leaf &&
        solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
        solverPatch.getRefinementStatus()>=_solver->_minRefinementStatusForTroubledCell-2; 
        //solverPatch.getCorruptionStatus()!=ADERDGSolver::CorruptedAndCorrected &&
        //solverPatch.getCorruptionStatus()!=ADERDGSolver::CorruptedNeighbour;
  }

  /**
   * Recompute the FV solution in cells that have been subject to a limiter status change
   *
   * We perform the following actions based on the new limiter status:
   *
   * |New Status | Action                                                                                                                                      |
   * ----------------------------------------------------------------------------------------------------------------------------------------------------------|
   * |O          | Do nothing. Solver solution has been evolved correctly before.
   * |T/NT       | Evolve limiter. Project result onto the ADER-DG space. Recompute the space-time predictor in NT cells if fused scheme is used.
   * |NNT        | Rollback solver and project result onto limiter solution
   *
   * Legend: O: Ok (ADER-DG cells), T: Troubled (FV cells), NT: FV->DG cells, NNT: DG->FV cells
   *
   * @param solverPatch a solver patch
   * @param cellInfo    struct referring to all cell description associated with a cell
   */
  void localRecomputation(
      SolverPatch&                                       solverPatch,
      CellInfo&                                          cellInfo,
      const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers);

  /**
   * Function Body for public recomputeSolutionLocally function.
   */
  double localRecomputationBody(
      SolverPatch&                                       solverPatch,
      Solver::CellInfo&                                  cellInfo,
      const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers);

  /**
   * Veto a coarsening attempt of a Parent cell description if
   * the restricted solution is troubled.
   *
   * @param cellDescriptionsIndex heap index associated with a mesh cell
   * @param solverNumber          identification number of a solver.
   */
  void vetoCoarseningIfRestrictedSolutionIsTroubled(
      const int cellDescriptionsIndex,
      const int solverNumber);

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
      const SolverPatch&                           solverPatch,
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
      SolverPatch&                                 solverPatch,
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
      Solver::BoundaryFaceInfo& face,
      const double* const min, const double* const  max) const;

#endif

  /**
   * A job which performs a fused limiting ADER-DG time step, i.e., it performs the solution update,
   * updates the local time stamp, performs the space-time predictor commputation,
   * and finally evaluates the a-posteriori limiting criterion.
   *
   * TODO(Dominic): Minimise time step sizes and refinement requests per patch
   * (->transpose the typical minimisation order)
   */
  class FusedTimeStepJob: public tarch::multicore::jobs::Job {
  private:
    LimitingADERDGSolver&                             _solver;
    SolverPatch&                                      _solverPatch;
    CellInfo                                          _cellInfo;               // copy
    const double                                      _predictionTimeStamp;    // copy
    const double                                      _predictionTimeStepSize; // copy
    const bool                                        _isFirstTimeStepOfBatch; // copy
    const bool                                        _isLastTimeStepOfBatch;  // copy
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int> _boundaryMarkers;        // copy
    const bool                                        _isSkeletonJob;
  public:

  /**
   * @note Job is spawned as high priority job if spawned in the last time step.
   * It further spawns a prediction job in this case in order
   * to overlap work with the reduction of time step size
   * and mesh update events.
   *
   * @param solver                 the solver the patch is associated with
   * @param solverPatch            solver patch linking to the data
   * @param cellInfo               refers to all cell descriptions/patches found in the cell holding the solver patch
   * @param predictionTimeStamp    the time stamp used for the prediction at time of the job spawning.
   * @param predictionTimeStepSize the time step size used for the prediction at time of the job spawning.
   * @param isFirstTimeStepOfBatch if this is the first time step in a batch
   * @param isLastTimeStepOfBatch  if this is the last time step in a batch
   * @param isSkeletonJob          if this job was spawned in a cell belonging to the MPI or AMR skeleton
     */
    FusedTimeStepJob(
        LimitingADERDGSolver&                              solver,
        SolverPatch&                                       solverPatch,
        CellInfo&                                          cellInfo,
        const double                                       predictionTimeStamp,
        const double                                       predictionTimeStepSize,
        const bool                                         isFirstTimeStepOfBatch,
        const bool                                         isLastTimeStepOfBatch,
        const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers,
        const bool                                         isSkeletonJob);

    bool run(bool runOnMasterThread) override;
  };

  /**
   * A job which performs the solution update and computes a new time step size.
   *
   * @note Spawning these operations as background job makes only sense if you
   * wait in endIteration(...) on the completion of the job.
   * It further important to flag this job as high priority job to
   * ensure completion before the next reduction.
   *
   * @note All update jobs have the same priority as their outcome is not directly
   * piped into an MPI send task. This is different to the result of the FusedTimeStepJob.
   */
  class UpdateJob: public tarch::multicore::jobs::Job {
    private:
      LimitingADERDGSolver&                             _solver; // TODO not const because of kernels
      SolverPatch&                                      _solverPatch;
      CellInfo                                          _cellInfo;
      const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int> _boundaryMarkers; // copy
    public:
      /**
       * Construct an UpdateJob.
       *
       * @note Job is always spawned as high priority job.
       *
       * @param     solver          the spawning solver
       * @param     solverPatch     a cell description
       * @param     cellInfo        links to all cell descriptions associated with the cell
       * @param[in] boundaryMarkers per face, a flag indicating if the cell description is adjacent to a remote or domain boundary.
       */
      UpdateJob(
          LimitingADERDGSolver&                              solver,
          SolverPatch&                                       solverPatch,
          CellInfo&                                          cellInfo,
          const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers);

      bool run(bool runOnMasterThread) override;
      void prefetchData() override;
  };

  /**
   * A job which performs the local recomputation and computes a new time step size.
   *
   * @note Spawning these operations as background job makes only sense if you
   * wait in endIteration(...) on the completion of the job.
   * It further important to flag this job as high priority job to
   * ensure completion before the next reduction.
   */
  class LocalRecomputationJob: public tarch::multicore::jobs::Job {
  private:
    LimitingADERDGSolver&                             _solver; // TODO not const because of kernels
    SolverPatch&                                      _solverPatch;
    CellInfo                                          _cellInfo;
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int> _boundaryMarkers; // copy
  public:
    /**
     * Construct an LocalRecomputationJob.
     *
     * @note Job is always spawned as high priority job.
     *
     * @param solver             the spawning solver
     * @param solverPatch        a cell description
     * @param cellInfo           links to all cell descriptions associated with the cell
     * @param isAtRemoteBoundary if the cell is at boundary to a remote rank
     */
    LocalRecomputationJob(
        LimitingADERDGSolver&                              solver,
        SolverPatch&                                       solverPatch,
        CellInfo&                                          cellInfo,
        const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers);

    bool run(bool runOnMasterThread) override;
    void prefetchData() override;
  };

  /**
   * A job that calls Solver::adjustSolutionDuringMeshRefinementBody(...).
   */
  class AdjustSolutionDuringMeshRefinementJob: public tarch::multicore::jobs::Job {
  private:
    LimitingADERDGSolver& _solver;
    SolverPatch&          _solverPatch;
    CellInfo              _cellInfo; // copy
    const bool            _isInitialMeshRefinement;
  public:
    AdjustSolutionDuringMeshRefinementJob(
        LimitingADERDGSolver& solver,
        SolverPatch&          solverPatch,
        CellInfo&             cellInfo,
        const bool            isInitialMeshRefinement);

    bool run(bool runOnMasterThread) override;
  };

#if defined(Parallel)
  class CheckAndCorrectSolutionJob : public tarch::multicore::jobs::Job {
    enum class SDCCheckResult {NoCorruption, OutcomeSaneAsLimiterNotActive, UncorrectableSoftError};

    private:
      LimitingADERDGSolver&               _solver;
      SolverPatch&                _solverPatch;
      CellInfo                    _cellInfo;
      const double                _predictorTimeStamp;
      const double                _predictorTimeStepSize;

      // actual execution of a STP job
      bool run(bool runOnMasterThread) override;
      SDCCheckResult checkAgainstOutcome(ADERDGSolver::MigratablePredictionJobData *outcome);
      void correctWithOutcomeAndDeleteLimiterStatus(ADERDGSolver::MigratablePredictionJobData *outcome);

    public:
      // constructor for local jobs that can be stolen
      CheckAndCorrectSolutionJob(
        LimitingADERDGSolver&     solver,
        SolverPatch& solverPatch,
        CellInfo& cellInfo,
        const double predictorTimeStamp,
        const double predictorTimeStepSize
      );
      ~CheckAndCorrectSolutionJob();

      CheckAndCorrectSolutionJob(const CheckAndCorrectSolutionJob& stp) = delete;
  };

#endif

public:
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
      const double DMPDifferenceScaling=1e-3);

  virtual ~LimitingADERDGSolver() {
    _solver.reset();
    _limiter.reset();
  }

  // Disallow copy and assignment
  LimitingADERDGSolver(const ADERDGSolver& other) = delete;
  LimitingADERDGSolver& operator=(const ADERDGSolver& other) = delete;
  
  /** Wire through to ADER-DG Solver */
  void updateMeshUpdateEvent(MeshUpdateEvent meshUpdateEvent) final override;
  void resetMeshUpdateEvent() final override;
  MeshUpdateEvent getMeshUpdateEvent() const final override;

  double getMinTimeStamp() const final override;
  double getMinTimeStepSize() const final override;
  double getAdmissibleTimeStepSize() const override;
  void updateAdmissibleTimeStepSize(double value) override;
  void resetAdmissibleTimeStepSize() final override;

  void initSolver(
      const double                                timeStamp,
      const tarch::la::Vector<DIMENSIONS,double>& domainOffset,
      const tarch::la::Vector<DIMENSIONS,double>& domainSize,
      const double                                boundingBoxSize,
      const double                                boundingBoxMeshSize,
      const std::vector<std::string>&             cmdlineargs,
      const exahype::parser::ParserView&          parserView) final override;

  /** @copydoc Solver::kickOffTimeStep
   *
   * Kick off time step on solver and limiter.
   */
  void kickOffTimeStep(const bool isFirstTimeStepOfBatchOrNoBatch) final override;

  /** @copydoc Solver::kickOffTimeStep
   *
   * Wrap up solver and limiter time step.
   * Compare the global observables of both and
   * copy the result to both solver.s
   */
  void wrapUpTimeStep(const bool isFirstTimeStepOfBatchOrNoBatch,const bool isLastTimeStepOfBatchOrNoBatch) final override;

  bool isPerformingPrediction(const exahype::State::AlgorithmSection& section) const final override;
  bool isMergingMetadata(const exahype::State::AlgorithmSection& section) const final override;

  /**
   * Synchronies the solver patch with the global time step data and
   * further copies this data to the FV patch if one is allocated
   * for this solver patch.
   *
   * @param solverPatch a solution patch of the main solver
   * @param limiterPatches list of limiter patches registered at the cell
   * @param limiterElement element index of the limiter patch allocated for
   *                       the solver patch, or NotFound if none has been allocated.
   */
  void synchroniseTimeStepping(SolverPatch& solverPatch,Solver::CellInfo& cellInfo) const;

  void updateTimeStepSize() final override;

  void updateGlobalObservables() final override;

  /**
   * Roll back the time step data to the
   * ones of the previous time step.
   */
  void rollbackToPreviousTimeStep() final override;

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
   * @return  the limiter patch matching the solver patch.
   *
   * @note Copies the solver patch time step data onto the limiter patch
   *
   * @param solverPatch a solver patch
   * @param cellInfo    holds references to the cell descriptions associated with a mesh cell
   *
   * @note Only use this function if it sure that there is a limiter patch for the
   * given solver patch.
   */
  LimiterPatch& getLimiterPatch(const SolverPatch& solverPatch,CellInfo& cellInfo) const;

  /**
   * Similar to @see getLimiterPatch(const SolverPatch& solverPatch,CellInfo& cellInfo)
   * but here the limiterElement was already obtained.
   *
   * @note Copies the solver patch time step data onto the limiter patch.
   */
  LimiterPatch& getLimiterPatch(const SolverPatch& solverPatch,CellInfo& cellInfo,const int limiterElement) const;

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
   * @note Currently we assume that the problem and load-balancing is so well-behaved that
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
   * @note We overwrite the facewise limiter status values with the new value
   * in order to use the updateLimiterStatusAfterSetInitialConditions function
   * afterwards which calls determineLimiterStatus(...) again.
   *
   * @note calls synchroniseTimeStepping
   *
   * returns true if a new limiter patch was allocated.
   */
  MeshUpdateEvent  updateRefinementStatusDuringRefinementStatusSpreading(SolverPatch& solverPatch) const;

  void updateCorruptionStatusDuringRefinementStatusSpreading( SolverPatch& solverPatch) const;


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
      const int  solverNumber) override;

 void progressMeshRefinementInLeaveCell(
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
     const int level,
     const bool checkThoroughly,
     bool& checkSuccessful) const final override;

  void finaliseStateUpdates(
      const int solverNumber,
      CellInfo& cellInfo) final override;

  ///////////////////////////////////
  // CELL-LOCAL
  //////////////////////////////////

  double startNewTimeStep(
      SolverPatch& solverPatch,Solver::CellInfo& cellInfo,
      const bool isFirstTimeStepOfBatch);

  void updateTimeStepSize(const int solverNumber,CellInfo& cellInfo) final override;

  void updateGlobalObservables(const int solverNumber,CellInfo& cellInfo) final override;

  /**
   * Mostly identical but slightly different to the ADER-DG equivalent
   * as only non-troubled cells compute a predictor.
   *
   * Backs up solver and limiter solution if not all ADER-DG phases are fused.
   * Adds volume integral result directly to solution if not all ADER-DG phases are fused.
   */
  void predictionAndVolumeIntegral(
      const int solverNumber,
      CellInfo& cellInfo,
      const bool isAtRemoteBoundary);

  void fusedTimeStepOrRestrict(
      const int                                          solverNumber,
      CellInfo&                                          cellInfo,
      const bool                                         isFirstTimeStepOfBatch,
      const bool                                         isLastTimeStepOfBatch,
      const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers) final override;

  void updateOrRestrict(
      const int                                          solverNumber,
      CellInfo&                                          cellInfo,
      const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers) final override;

  void compress(
      const int  solverNumber,
      CellInfo&  cellInfo,
      const bool isAtRemoteBoundary) const final override;

  void adjustSolutionDuringMeshRefinement(const int solverNumber,CellInfo& cellInfo) final override;

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
   * @note Has no const modifier since kernels are not const functions yet.
   *
   * @param[in] backupPreviousSolution            Set to true if the solution should be backed up before
   *                                              we overwrite it by the updated solution.
   * @param[in] addSurfaceIntegralResultToUpdate  set to true if the surface integral result should be added to the update. Otherwise, is added directly to the solution.
   *                                              (Fused time stepping for nonlinear PDEs is the only time stepping variant where we need to use an update vector.)
   */
  void updateSolution(
      SolverPatch&                                       solverPatch,
      CellInfo&                                          cellInfo,
      const bool                                         isFirstTimeStep,
      const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers,
      const bool                                         addSurfaceIntegralResultToSolution);

  /**
   * Evaluate the discrete maximum principle and the physically admissibility detection criterion.
   *
   * Furthermore, computes the new cell local minimum and maximum of the DMP observables.
   * These are written to the solutionMin and solutionMax arrays of the cell solver patch.
   * If a cell is troubled, the limiter patch is used to compute the minimum and maximum.
   *
   * @param solverPatch a solver patch
   * @param cellInfo    a struct referring to all cell descriptions
   *
   * @return true if the cell is troubled.
   */
  bool checkIfCellIsTroubledAndDetermineMinAndMax(SolverPatch& solverPatch,CellInfo& cellInfo);

  /**
   * Determine the new cell-local min max values.
   *
   * Must be invoked after ::determineLimiterStatusAfterSolutionUpdate.
   *
   * TODO(Dominic): Tobias's integer
   * flagging idea might reduce complexity here
   */
  void determineMinAndMax(const int solverNumber,Solver::CellInfo& cellInfo);

  /**
   * Deallocate the limiter patch on all AMR related
   * helper cells.
   *
   * It is safe to use this method during
   * the mesh refinement iterations.
   *
   * @note Thread-safe.
   */
   void ensureNoLimiterPatchIsAllocatedOnHelperCell(
       const SolverPatch& solverPatch,CellInfo& cellInfo) const;

   /**
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
       CellInfo&          cellInfo,
       const int          limiterStatus) const;


   /**
    * Go back to previous time step with
    * time step data and solution.
    *
    * Reset the refinement status to the old value.
    * (Nothing has happened but we have a new mesh.)
    *
    * Allocate necessary new limiter patches.
    */
   void rollbackSolutionGlobally(const int  solverNumber,CellInfo&  cellInfo) const final override;

  /**
   * Reinitialises cells that have been subject to a limiter status change.
   * This method is invoked after the limiter status spreading.
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
   * and the last solution update has to be redone.
   *
   * In this case, we perform a complete rollback, i.e. we redo the last time step (with a new mesh).
   *
   * @note Solely the solution arrays are swapped and projections are performed where necessary. The time step
   * data on the patches is not touched. This has to be taken into account when the local recomputation is performed.
   *
   * @see recomputeSolutionLocally (which is run on iteration after).
   */
  void rollbackSolutionLocally(
      const int  solverNumber,
      CellInfo&  cellInfo,
      const bool fusedTimeStepping) const;

  /**
   * Invoke ::localRecomputationBody(...)
   *
   * Fused Time Stepping
   * -------------------
   *
   * Recomputes a new predictor as well if necessary.
   * Further see ::fusedTimeBody regarding order of operations.
   */
  void localRecomputation(
      const int                                          solverNumber,
      Solver::CellInfo&                                  cellInfo,
      const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers);

  ///////////////////////////////////
  // NEIGHBOUR
  ///////////////////////////////////

  /**
   * Merge two neighbouring cell descriptions during the local recomputation phase of the limiting ADER-DG
   * algorithm.
   *
   * As it is performed during the local recomputation phase, the merge only considers
   * cells with troubled status - i, where i=0,1,2.
   * All cells with one of these statuses, perform a merge using the limiter.
   *
   * @param solverNumber identifier for this solver
   * @param cellInfo1    cell descriptions associated with one cell
   * @param cellInfo2    cell descriptions associated with the other cell
   * @param pos1         relative position of the first cell to the vertex issuing the neighbour merge
   * @param pos2         relative position of the first cell to the vertex issuing the neighbour merge
   */
  void mergeNeighboursDataDuringLocalRecomputation(
      const int                                  solverNumber,
      Solver::CellInfo&                          cellInfo1,
      Solver::CellInfo&                          cellInfo2,
      const tarch::la::Vector<DIMENSIONS, int>&  pos1,
      const tarch::la::Vector<DIMENSIONS, int>&  pos2);


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
   * | O        | Tm2      | ADER-DG       |// O|Tm2 x O|Tm2
   * | Tm2      | O        | ADER-DG       |
   * | Tm2      | Tm2      | ADER-DG       |
   *
   * | Tm2      | Tm1       | FV            |
   * | Tm1      | Tm2       | FV            | // Tm1&Tm2 | N&Tm2
   *
   * | Tm1       | Tm1      | FV            |
   * | Tm1       | T        | FV            |
   * | T         | Tm1      | FV            | // T|Tm1 x T|Tm1
   * | T         | T        | FV            |
   *
   * Legend: O: Ok, T: Troubled, Tm1: Troubled (status) minus 1, Tm2: Troubled (status) minus 2
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
   * @note Limiting is only performed on the finest level
   * of the mesh. The other levels work only with the ADER-DG
   * solver.
   *
   * TODO(Dominic): Remove limiterstatus1 and limiterStatus2 argument.
   * They depend on the isRecomputation value
   *
   * @param solverNumber
   * @param cellInfo1
   * @param cellInfo2
   * @param pos1
   * @param pos2
   * @param isRecomputation
   */
  void mergeNeighboursData(
      const int                                  solverNumber,
      Solver::CellInfo&                          cellInfo1,
      Solver::CellInfo&                          cellInfo2,
      const tarch::la::Vector<DIMENSIONS, int>&  pos1,
      const tarch::la::Vector<DIMENSIONS, int>&  pos2);

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

  /**
   * Send data to a neighbouring rank during the local recomputation phase.
   *
   * Merge it if the limiter status of the solver patch is troubled or one or two less than
   * troubled. All other cells send empty messages.
   *
   * @param toRank       the rank we send data to
   * @param solverNumber identifier for this solver
   * @param cellInfo     a struct linking to the cell descriptions associated with a cell
   * @param src          position of the src
   * @param dest         position of the dest
   * @param x            barycentre of the face
   * @param level        mesh level of the face
   */
  void sendDataToNeighbourDuringLocalRecomputation(
      const int                                    toRank,
      const int                                    solverNumber,
      CellInfo&                                    cellInfo,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level);

  /**
   * Receive limiter data from a neighbouring rank durin the local
   * recomputation phase.
   *
   * Merge it if the limiter status of the solver patch is troubled or one less than
   * troubled. All other cells drop it.
   *
   * @param fromRank     the rank we receive data from
   * @param solverNumber identifier for this solver
   * @param cellInfo     a struct linking to the cell descriptions associated with a cell
   * @param src          position of the src
   * @param dest         position of the dest
   * @param x            barycentre of the face
   * @param level        mesh level of the face
   */
  void mergeWithNeighbourDataDuringLocalRecomputation(
      const int                                    fromRank,
      const int                                    solverNumber,
      CellInfo&                                    cellInfo,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level);

  /**
   * Send data to a neighbouring rank.
   *
   * send order:   minAndMax,solver,limiter
   * receive order limiter,solver,minAndMax
   *
   * @note The limiter is only active on the finest mesh level
   * in ExaHyPE. All other levels use the ADER-DG solver.
   * All level send the minimum and maximum DMP
   * observables around, which are used for
   * the shock detection.
   *
   * @note If we send out data during a fused time step,
   * a cell might be newly marked as troubled but a new limiter patch
   * is not allocated yet. This will be
   * done in the following local recomputation
   * step. We send an empty limiter message in this case.
   *
   * @param toRank       the rank we send data to
   * @param solverNumber identification number of this solver
   * @param cellInfo     links to the data assocated with the source cell
   * @param src          relative position of message source to vertex
   * @param dest         relative position of message destination to vertex
   * @param x            face barycentre
   * @param level        mesh level
   */
  void sendDataToNeighbour(
      const int                                     toRank,
      const int                                     solverNumber,
      Solver::CellInfo&                             cellInfo,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level);

  void mergeWithNeighbourData(
      const int                                    fromRank,
      const int                                    solverNumber,
      Solver::CellInfo&                            cellInfo,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level);

  /**
   * @note only used when no recomputation.
   *
   * @param fromRank the rank we expect data from
   * @param x
   * @param level
   */
  void dropNeighbourData(
      const int                                    fromRank,
      const int                                    solverNumber,
      Solver::CellInfo&                            cellInfo,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const;

  /////////////////////////////////////
  // MASTER<=>WORKER
  /////////////////////////////////////
  /**
   * Kind of similar to progressMeshRefinementInPrepareSendToWorker
   * but performs a few additional operations in order to
   * notify the worker about some coarse grid operations only
   * the master knows.
   *
   * @note This function sends out MPI messages.
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
  void progressMeshRefinementInMergeWithMaster(
      const int worker,
      const int localCellDescriptionsIndex,
      const int localElement,
      const int coarseGridCellDescriptionsIndex,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int  level) final override;

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

  /** @copydoc Solver::mergeWithMasterData
   *
   * Send the solver's global observables to the
   * the master rank.
   *
   * @note Solver and limiter observables have been merged in
   * wrapUpTimeStep.
   *
   * @see wrapUpTimeStep
   */
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

  /** @copydoc Solver::mergeWithMasterData
   *
   * Overwrites the limiter's global observables with global observables
   * received by the solver.
   *
   * @see wrapUpTimeStep
   */
  void mergeWithMasterData(
      const                                        int masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) final override;

   void initOffloadingManager();
   void stopOffloadingManager();
#ifndef OffloadingUseProgressThread
   void resumeOffloadingManager();
   void pauseOffloadingManager();
#endif
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

  ///////////////////////
  // PROFILING
  ///////////////////////

  CellProcessingTimes measureCellProcessingTimes(const int numberOfRuns=100) override;

  /////////////////////////////////
  // USER HOOKS - Make inaccessible
  /////////////////////////////////

  int getGeometricLoadBalancingWeight(
      const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>& cellSize) final override {
    return _solver->getGeometricLoadBalancingWeight(cellCentre,cellSize);
  }

  void resetGlobalObservables(double* const globalObservables) final override;

  void updateGlobalObservables(
      double* const                               globalObservables,
      const double* const                         luh,
      const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>& cellSize,
      const double t,
      const double dt) final override;

  void mergeGlobalObservables(
      double* const       globalObservables,
      const double* const otherObservables) final override;

  void wrapUpGlobalObservables(
      double* const globalObservables) final override;

protected:

  /** @name Plugin points for derived solvers.
   *
   *  These are the macro kernels solvers derived from
   *  ADERDGSolver need to implement.
   */
  ///@{
  //Virtual methods to bind the limiter kernels in the abstractLimiterSolver (generated by the Toolkit and application specific)
  virtual void projectOnFVLimiterSpace(const double* const luh, double* const lim) const = 0;
  virtual void projectOnDGSpace(const double* const lim, double* const luh) const = 0;
  virtual bool discreteMaximumPrincipleAndMinAndMaxSearch(const double* const luh, double* const boundaryMinPerVariables, double* const boundaryMaxPerVariables) = 0;
  virtual void findCellLocalMinAndMax(const double* const luh, double* const localMinPerVariables, double* const localMaxPerVariable) = 0;
  virtual void findCellLocalLimiterMinAndMax(const double* const lim, double* const localMinPerObservable, double* const localMaxPerObservable) = 0;
  ///@}
};


#endif /* LIMITEDADERDGSOLVER_H_ */
