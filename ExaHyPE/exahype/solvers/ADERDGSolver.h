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
 *
 * \author Dominic E. Charrier, Tobias Weinzierl, Jean-Matthieu Gallard, Fabian Güra
 **/

#ifndef _EXAHYPE_SOLVERS_ADERDG_SOLVER_H_
#define _EXAHYPE_SOLVERS_ADERDG_SOLVER_H_

#include <iostream>
#include <string>
#include <vector>
#include <tuple>

#include "exahype/solvers/Solver.h"

#include "peano/heap/Heap.h"
#include "peano/utils/Globals.h"

#include "tarch/Assertions.h"
#include "tarch/la/Vector.h"

#include "exahype/profilers/simple/NoOpProfiler.h"
#include "exahype/records/ADERDGCellDescription.h"

namespace exahype {
  namespace parser {
    class ParserView;
  }
  namespace solvers {
    class ADERDGSolver;
  }
}

/**
 * ADER-DG solver base class.
 *
 *
 * Mesh Refinement
 * ===============
 *
 * Some notes on the solver which will hopefully grow into a full documentation.
 *
 * - The user uses a RefinementControl to steer the refinement automation implemented
 *   by this solver. The time stepping is only stopped if refinement is requested.
 *
 * - We perform rollbacks after the mesh refinement as we assume that the
 *   mesh refinement is always triggered to late by the user (a-posteriori refinement).
 *   We thus have to check if we can erase the current and the previous solution.
 *
 * - The refinement automaton uses a RefinementEvent to track the state of
 *   each cell. The automata of different cells interplay when
 *   refining or erasing of cells is performed.
 *
 * - In order to prevent refinement/erasing oscillations, we restrict the
 *   solution and previousSolution to the parent and evaluate the refinement
 *   criterion for both. Only if no refinement is triggered, the
 *   erasing is actually performed.
 *
 * - On the finest level, we evaluate the refinement criterion but we also
 *   spread a status for halo refinement. A cells refinement decision making
 *   is thus depending on two inputs. We combine both information
 *   into a refinement status. We choose always the maximum
 *   value of the inputs.
 *   The hybrid limiter solver further injects the limiter status into
 *   the refinement status.
 *
 * Subsequently, we list refinement "stories" cell descriptions of this solver
 * might "experience".
 *
 * Safe, Oscillation-Free Erasing of Child Cells
 * ---------------------------------------------
 *
 * This procedure can take up to 3 iterations.
 *
 * 1. All children (type: Cell) of a cell (type: Ancestor) want to be erased.
 * They have the refinement status Erase. In the previous iteration, they
 * had this refinement status, too. (adjustSolutionDuringMeshRefinementBody -> markForRefinement)
 * They write their max refinement status to the parent (progressMeshRefinementInEnterCell->updateCoarseGridAncestorRefinementStatus).
 * (If one of them is an Ancestor itself, they set an erasing children veto flag which is reset by
 * the parent in every iteration.)
 *
 * 2. The parent (type: Ancestor) changes its type to Cell and allocates
 * the necessary memory. (progressMeshRefinementInLeaveCell->progressCollectiveRefinementOperationsInLeaveCell)
 * The parent notifies its children by setting the ErasingChildrenRequested RefinementEvent.
 *
 * 3.The children restrict their volume data (solution and previous solution) up to the parent. They
 * do not delete their data yet. (progressMeshRefinementInLeaveCell->restrictVolumeDataIfErasingChildrenRequested)
 *
 * 4. The parent evaluates the refinement criterion for current solution and previous solution.
 * If one of the evaluations indicates that refining is necessary, the children erasing
 * process is stopped and the parent changes its type back to Ancestor. (progressMeshRefinementInLeaveCell->progressCollectiveRefinementOperationsInLeaveCell)
 * During the mesh refinement iterations, it stores the current and previous refinement status in order
 * to prevent the procedure to start again.
 *
 * If no refinement at both time stages is required, the parent keeps its new type (Cell) and signals its children to delete
 * themselves. It does this by changing the RefinementEvent from ErasingChildrenRequested -> ErasingChildren.
 *
 * 5. The parent resets its status to None.
 */
class exahype::solvers::ADERDGSolver : public exahype::solvers::Solver {
  friend class LimitingADERDGSolver;
public:
  /**
   * The maximum helper status.
   * This value is assigned to cell descriptions
   * of type Cell.
   */
  static int CellCommunicationStatus;
  /**
   * The minimum helper status a cell description
   * must have for it allocating boundary data.
   */
  static int MinimumCommunicationStatusForNeighbourCommunication;

  /**
   * The maximum augmentation status.
   * This value is assigned to cell descriptions
   * of type Ancestor.
   */
  static int MaximumAugmentationStatus;
  /**
   * The minimum augmentation status a cell description
   * of type Cell must have for it to refine
   * and add child cells of type Descendant to
   * the grid.
   */
  static int MinimumAugmentationStatusForVirtualRefining;
  /**
   * The minimum augmentation status for refining
   * a cell. Note that there should be at least layer
   * of width 2 between the status for erasing (0)
   * and the one for augmentation && refining (>=3).
   *
   * TODO(Dominic):
   * I have too look a little further into how
   * Peano erases to explain the above experimentally
   * found values better.
   */
  static int MinimumAugmentationStatusForRefining;

  /**
   * Semaphore for fine grid cells restricting
   * volume data to a coarse grid parent.
   */
  static tarch::multicore::BooleanSemaphore RestrictionSemaphore;

  /**
   * Semaphore for fine grid cells accessing a coarse grid parent.
   */
  static tarch::multicore::BooleanSemaphore CoarseGridSemaphore;

  /**
   * Rank-local heap that stores ADERDGCellDescription instances.
   *
   * \note This heap might be shared by multiple ADERDGSolver instances
   * that differ in their solver number and other attributes.
   * @see solvers::Solver::RegisteredSolvers.
   */
  typedef exahype::records::ADERDGCellDescription CellDescription;
  typedef peano::heap::RLEHeap<CellDescription> Heap;

private:

  /**
   * Log device.
   */
  static tarch::logging::Log _log;

  #ifdef Parallel
  DataHeap::HeapEntries _receivedExtrapolatedPredictor;
  DataHeap::HeapEntries _receivedFluctuations;
  DataHeap::HeapEntries _receivedUpdate;
  /**
   * TODO(WORKAROUND): We store these fields in order
   * to use the symmetric boundary exchanger of Peano
   * which does not yet support asymmetric send buffers.
   */
  DataHeap::HeapEntries _invalidExtrapolatedPredictor;
  DataHeap::HeapEntries _invalidFluctuations;
  #endif

  /**
   * Minimum corrector time stamp of all cell descriptions.
   */
  double _previousMinCorrectorTimeStamp;

  /**
   * Minimum corrector time step size of all
   * cell descriptions in the previous iteration.
   *
   * This time step size is necessary for the fused time stepping + limiting
   * to reconstruct the minCorrectorTimeStepSize during a rollback.
   */
  double _previousMinCorrectorTimeStepSize;

  /**
   * Minimum corrector time stamp of all cell descriptions.
   */
  double _minCorrectorTimeStamp;

  /**
   * Minimum corrector time step size of
   * all cell descriptions.
   */
  double _minCorrectorTimeStepSize;

  /**
   * Minimum predictor time stamp of all cell descriptions.
   * Always equal or larger than the minimum corrector time stamp.
   */
  double _minPredictorTimeStamp;

  /**
   * Minimum predictor time step size of
   * all cell descriptions.
   */
  double _minPredictorTimeStepSize;

  /**
   * Minimum next predictor time step size of
   * all cell descriptions.
   */
  double _minNextTimeStepSize;

  /**
   * A flag that is used to track if the
   * CFL condition of a solver was violated.
   */
  bool _stabilityConditionWasViolated;

  /** Special Refinement Status values */
  static constexpr int BoundaryStatus             = -3;
  static constexpr int Pending                    = -2;
  static constexpr int Erase                      = -1; // Erase must be chosen as -1.
  static constexpr int Keep                       =  0;

  int _refineOrKeepOnFineGrid; // can be configured by the user

  /**
   * Number of limiter helper layers in each
   * helper cell subdomain around a troubled cell.
   */
  const int _limiterHelperLayers;
  /**
   * !!! LimitingADERDGSolver functionality !!!
   *
   * The number of observables
   * the discrete maximum principle
   * is applied to.
   */
  const int _DMPObservables;

  /**
   * Cells placed in the separation layers around a troubled cells either project the DG solution
   * onto an FV patch or compute with FV and project onto a DG polynomial depending on
   * if they neighbour a well-behaved or troubled cell, respectively.
   * If a cell in the separation layer becomes troubled, a local recomputation
   * must be performed.
   */
  const int _minRefinementStatusForSeparationCell;

  /**
   * Cells placed in the buffer layers around a troubled cell either project the DG solution
   * onto an FV patch or compute with FV and project onto a DG polynomial depending on
   * if they neighbour a well-behaved or troubled cell, respectively.
   * If a cell in the buffer layer becomes troubled, no local recomputation must be performed
   */
  const int _minRefinementStatusForBufferCell;

  /**
   * Minimum limiter status a troubled cell can have.
   */
  const int _minRefinementStatusForTroubledCell;

  /**
   * Check for NaNs.
   */
  bool _checkForNaNs;

  /**
   * The current mesh update event.
   */
  MeshUpdateEvent _meshUpdateEvent;

  /**
   * The mesh update event which will become active
   * in the next iteration.
   */
  MeshUpdateEvent _nextMeshUpdateEvent;

  /**
   * Different to compress(), this operation is called automatically by
   * mergeNeighbours(). Therefore the routine is private.
   *
   * \note This routine checks if a cell description is
   * compressed. No previous check is necessary.
   */
  void uncompress(CellDescription& cellDescription) const;

  /**
    * Rollback to the previous time step, i.e,
    * overwrite the time step size and time stamp
    * fields of the cell description
    * by previous values.
    */
   void rollbackToPreviousTimeStep(CellDescription& cellDescription) const;

   /**
    * Same as rollbackToPreviousTimeStep
    * but for the fused time stepping scheme.
    *
    * Corrector time stamp and corrector time step size must
    * add up to predictor time stamp after rollback.
    *
    * Corrector time step size is assumed to be used
    * predictor time step size in batch.
    */
   void rollbackToPreviousTimeStepFused(CellDescription& cellDescription) const;

  /**
   * Simply adjust the solution if necessary. Do not modify the time step
   * data or anything else.
   */
  void adjustSolution(CellDescription& cellDescription);

  /**
   * Body of FiniteVolumesSolver::adjustSolutionDuringMeshRefinement(int,int).
   *
   * @note May be called from background task. Do not synchronise time step data here.
   */
  void adjustSolutionDuringMeshRefinementBody(
      CellDescription& cellDescription,
      const bool isInitialMeshRefinement);

  /**
   * Update time step sizes and time stamps for the fused/nonfused time stepping variant.
   *
   * @param cellDescription   a cell description
   * @param fusedTimeStepping if fused time stepping is used
   * @return the new admissible time step size if the cell description is of type Cell. Otherwise, the maximum double value.
   */
  double updateTimeStepSizes(CellDescription& cellDescription,const bool fusedTimeStepping);

  /**
   * Evaluate the refinement criterion and convert the user
   * input to an integer representation, the refinement status.
   *
   * @param cellDescription   a cell description
   * @param solutionHeapIndex solution array heap index. Externalised in order to use the method for current and previous solution.
   * @param timeStamp         time stamp matching the solution array.
   * @return a refinement status (Erase(-1),Keep(0), or _refineOrKeepOnFineGrid (>=1))
   */
  int evaluateRefinementCriterion(
      const CellDescription& cellDescription,
      const double* const solution, const double& timeStamp);

  /**
   * Query the user's refinement criterion and
   * write a refinement request back to the cell description.
   */
  void markForRefinement(CellDescription& cellDescription);

  /**
   * Mark a cell description of Cell for refinement or erasing based
   * on a user supplied physics based refinement criterion.
   *
   * TODO(Dominic): Move docu below to appropriate location.
   *
   * <h2>Erasing</h2> TODO(Dominic): Move docu.
   * Note that we use a not so obvious strategy for performing
   * erasing operations. We first set an erasing request on
   * a parent cell description of type Ancestor or EmptyAncestor,
   * and then let its children of type Cell veto
   * this request if they want to keep their
   * solution or refine even further.
   *
   * No erasing children request can be set on cell descriptions
   * of type Ancestor which have been introduced to the grid during
   * the current mesh update iterations.
   * This prevents races where a refinement criterion has triggered a
   * refinement event on the parent cell but does trigger an erasing
   * event on the children cells.
   *
   * We further veto erasing events if
   * a child of the parent itself is a parent
   * of cell descriptions of type Descendant.
   *
   * <h2>Augmentation</h2>
   * Note that cell descriptions of type Cell are allowed to overwrite an augmentation request
   * by a refinement request if applicable.
   * The refinement event of a cell description of type Cell might be set to
   * an augmentation request in the methods mergeWithNeighbourData(...)
   * as well as in markForAugmentation(...) which is called from within
   * enterCell(...)
   *
   * \note Thread-safe.
   */
  void decideOnRefinement(CellDescription& fineGridCellDescription, const bool stillInRefiningMode);

  /**
   * Performs three operations:
   * 1. Checks if a ErasingVirtualChildrenRequestedTriggered event on the coarse
   * grid parent can be changed to a ErasingVirtualChildrenRequested event.
   * In this case, the triggered request becomes an actual request.
   * The fine grid children can however still veto this request.
   * 2.
   *
   * \note Thread-safe.
   */
  void decideOnVirtualRefinement(CellDescription& fineGridCellDescription);

  /**
   * Change the erasing children request to a change children to descendants
   * request of the coarse grid cell description's parent
   * if the coarse grid cell has children itself (of type Descendant).
   * Rationale: We cannot directly erase a Cell that has children (of type Descendant).
   *
   * Further, reset the erasing virtual children request if a coarse grid
   * Descendant has virtual children itself (of type Descendant). Rationale:
   * We cannot erase a coarse grid cell that has children (of type Descendant)
   * before erasing the children.
   *
   * \note This operation spans over three spacetree levels. Calling
   * it requires that a cell description for
   * the same solver the \p coarseGridCellDescription is associated with
   * is registered on the fine grid cell.
   *
   * \note A more sophisticated procedure has to performed for the refinement event
   * AugmentationRequested. We need to use the taversal's descend event to handle
   * this event. We thus do not rely on fineGridCell.isRefined() in the previous enterCell event
   * to check if we need to reset the deaugmenting request.
   *
   * TODO(Dominic): Make template function as soon as verified.
   *
   * \note Not thread-safe!
   */
  void alterErasingRequestsIfNecessary(
      CellDescription& coarseGridCellDescription,
      const int fineGridCellDescriptionsIndex) const;

  /**
   * Fills the solution and previous solution arrays
   * with zeros.
   */
  void prepareVolumeDataRestriction(
      CellDescription& cellDescription) const;

  /**
   * Start or finish collective operations from a
   * fine cell description point of view.
   *
   * \return true if a new compute cell
   * was allocated as result of an erasing operation.
   */
  void progressCollectiveRefinementOperationsInEnterCell(
      CellDescription& fineGridCellDescription);

  void progressCollectiveRefinementOperationsInLeaveCell(
      CellDescription& fineGridCellDescription,
      const bool stillInRefiningMode);

  /**
   * Checks the current and previous solution on a Cell which was
   * previously an Ancestor if any refinement is/was requested. If so,
   * no erasing must be performed and the Cell must be converted
   * back into an Ancestor.
   *
   * @param  cellDescription the considered cell description.
   * @return if a Cell (which was previously an Ancestor) can
   * erase its children.
   */
  bool markPreviousAncestorForRefinement(CellDescription& cellDescription);

  /**
   * In case, we change the children to a descendant
   * or erase them from the grid, we first restrict
   * volume data up to the parent and further
   * copy the corrector and predictor time stamps.
   *
   * \return true if we erase descendants from
   * the grid. In this case, to call an erase
   * on the grid/Peano cell if no other cell descriptions are
   * registered. Returns false otherwise.
   *
   * TODO(Dominic): More docu.
   *
   * \note This operations is not thread-safe
   */
  void eraseCellDescriptionIfNecessary(
      const int cellDescriptionsIndex,
      const int fineGridCellElement,
      CellDescription& coarseGridCellDescription);

  /**
   * Initialise cell description of type Cell.
   * Initialise the refinement event with None.
   *
   * \note This operations is thread-safe
   */
  void addNewCell(
      exahype::Cell& fineGridCell,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      const int coarseGridCellDescriptionsIndex,
      const int solverNumber);

  /**
   * Initialises helper cell descriptions of type Descendant
   * on the fine level after cell descriptions on the coarse level
   * have been flagged for augmentation and Peano has
   * created the requested new cells.
   *
   * Further sets the refinement event on a coarse grid Descendant to Augmenting
   * if the first new Descendant was initialised on the fine grid.
   *
   * Additionally, copies the information if a face is inside
   * from the parent to the new child cell.
   *
   * Reset an augmentation request if the child cell does hold
   * a Descendant or EmptyDescendant cell description with
   * the same solver number.
   *
   * This scenario occurs if an augmentation request is triggered in
   * enterCell().
   *
   * A similar scenario can never occur for refinement requests
   * since only cell descriptions of type Cell can be refined.
   * Ancestors can never request refinement.
   *
   * \note This operations is not thread-safe
   */
  void addNewDescendantIfVirtualRefiningRequested(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      CellDescription& coarseGridCellDescription,
      const int coarseGridCellDescriptionsIndex);

  /**
   * Initialises compute cell descriptions on the fine level (cell description type is Cell)
   * after coarse grid cell descriptions have been flagged for refinement and Peano has
   * created the requested new cells.
   * Erasing is not performed on cells belonging to the regular initial grid
   * of the solvers (see RegularMesh).
   *
   * Further sets the refinement event on a coarse grid Cell to Refining
   * if the first new Cell was initialised on the fine grid.
   *
   * Additionally, copies the information if a face is inside
   * from the parent to the new child cell.
   *
   * \return True if a cell description of type Cell was allocated
   * on the fineGridCell
   *
   * \note This operations is not thread-safe
   */
  bool addNewCellIfRefinementRequested(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      CellDescription& coarseGridCellDescription,
      const int coarseGridCellDescriptionsIndex);

  /**
   * Change a cell description of type Cell to an Ancestor.
   *
   * @param cellDescription a cell description of type Cell.
   */
  void changeCellToAncestor(CellDescription& cellDescription);

  /**
   * Prolongates Volume data from a parent cell description to
   * \p cellDescription if the fine grid cell associated with
   * \p cellDescription is adjacent to a boundary of the
   * coarse grid cell associated with the parent cell description.
   *
   * Further copy the corrector and predictor time stamp and
   * time step sizes.
   *
   * We prolongate both, the current and the previous
   * solution to the newly created fine grid cell description.
   * This especially important for the LimitingADERDGSolver.
   * Here, the cell descriptions of with LimiterStatus
   * NeighbourOfNeighbourOfTroubledCell need to communicate layers of
   * the previous solution to the neighbour.
   *
   * Furthermore, set the limiterStatus to the value of the
   * coarse grid cell description.
   * Set the value of the mergedLimiterStatus elements to Troubled
   * in case the coarse grid cell descriptions' values are Troubled.
   * Otherwise, set it to Ok.
   *
   * \note No static or const modifiers as kernels are not const.
   */
  void prolongateVolumeData(
      CellDescription& fineGridCellDescription,
      const bool initialGrid);

  /**
   * Restricts Volume data from \p cellDescription to
   * a parent cell description if the fine grid cell associated with
   * \p cellDescription is adjacent to a boundary of the
   * coarse grid cell associated with the parent cell description.
   *
   * \note !!! Currently, we minimise over the time step
   * sizes of the children. Not sure if this makes sense. TODO(Dominic)
   *
   * \note This method makes only sense for real cells.
   */
  void restrictVolumeDataIfErasingRequested(
      const CellDescription& fineGridCellDescription,
      const CellDescription& coarseGridCellDescription);

  /**
   * Ensure that the fine grid cell descriptions's parent index is pointing to the
   * coarse grid cell's cell descriptions index; this is important to re-establish
   * the parent-child relations on a new worker after a fork.
   *
   * Checks if the parent index of a fine grid cell description
   * was set to RemoteAdjacencyIndex during a previous forking event.
   *
   * If so, check if there exists a coarse grid cell description
   * which must have been also received during a previous fork event.
   * If so, update the parent index of the fine grid cell description
   * with the coarse grid cell descriptions index.
   *
   * For cell descriptions of type Descendant, copy offset and
   * level of the top-most parent cell description, which is of type Cell.
   */
  void ensureFineGridCoarseGridConsistency(
      CellDescription& cellDescription,
      const int coarseGridCellDescriptionsIndex);

  /**
   * Update the refinement status of a coarse grid parent of type
   * Ancestor.
   *
   * @param fineGridCellDescription the fine grid cell description
   * @param coarseGridElement       element of the
   */
  void updateCoarseGridAncestorRefinementStatus(
      const CellDescription& fineGridCellDescription,
      CellDescription& coarseGridCellDescription);

  /**
   * Turns checking for NaNs off.
   * If this solver is the main solver
   * of a LimitingADERDGSolver,
   * it does not make sense to check for
   * NaNs as those are cured by the FV limiter.
   */
  void disableCheckForNaNs();

  /**
   * Restrict the obse
   */
  void restrictObservablesMinAndMax(
      const CellDescription& cellDescription,
      const CellDescription& parentCellDescription) const;

  /**
   * Determine if the cell description of type
   * Descendant is on the cell boundary of its parent
   * of type Cell or Descendant with at least one of
   * its faces. If so restrict face data from the parent down
   * to the Descendant for those face(s).
   */
  void prolongateFaceDataToDescendant(
      CellDescription& cellDescription,
      const CellDescription& parentCellDescription,
      const tarch::la::Vector<DIMENSIONS,int>& subcellIndex);

  /**
   * Copies the parent cell descriptions observables'
   * minimum and maximum down to the Descendant.
   */
  void prolongateObservablesMinAndMax(
      const CellDescription& cellDescription,
      const CellDescription& cellDescriptionParent) const;

  /**
   * Solve the Riemann problem at the interface between two cells ("left" and
   * "right"). This method only performs a Riemann solve if at least one of the
   * cell descriptions (per solver) associated with the two cells is of type
   * ::Cell and none of the two cells belongs to the boundary.
   * In case a Riemann problem is solved,
   * the method further sets the ::riemannSolvePerformed
   * flags for the particular faces on both cell descriptions (per solver).
   *
   * This method further synchronises the ADERDGCellDescription
   * with the corresponding solver if this is required by the time stepping
   * scheme.
   * This operation must be performed in mergeWithNeighbour(...) and
   * touchVertexFirstTime(...) since both callbacks touch the
   * ADERDGCellDescriptions before the other callbacks.
   *
   * <h2>Rationale</h2>
   *
   * We did originally split up the boundary condition handling and the Riemann
   * updates into two mappings. This offers a functional decomposition. However,
   * both mappings then need a significiant number of technical administrative
   * code (cmp all the loops in touchVertexFirstTime and the redundant code to
   * manage the semaphores). We thus decided to merge both aspects. This also
   * should make sense from a performance point of view.
   *
   * We could potentially remove the face indices here if we had normals that
   * point outwards. However, we don't evaluate the direction of the normal and
   * thus need these counters as a Riemann problem on a face either could be
   * triggered by the left cell or by the right cell.
   *
   * \note The current implementation might classify cells with vertices that
   * are part of the
   * boundary of the domain or outside to be classified as inside of the domain
   * (volume-ratio based).
   *
   * \note We cannot solely check for indices of value
   * multiscalelinked::HangingVertexBookkepper::DomainBoundaryAdjacencyIndex
   * in vertex.getCellDescriptions() to determine if we are on the boundary of
   * the domain
   * since these values are overwritten by
   * multiscalelinked::HangingVertexBookkepper::RemoteAdjacencyIndex
   * if the domain boundary aligns with an MPI boundary
   * (see
   * multiscalelinkedcell::HangingVertexBookkeeper::updateCellIndicesInMergeWithNeighbour(...)).
   *
   * \note Not thread-safe.
   */
  void solveRiemannProblemAtInterface(
      CellDescription& pLeft,
      CellDescription& pRight,
      Solver::InterfaceInfo& face);

  /**
   * Apply the boundary conditions at the face with index \p faceIndex.
   *
   * This method further synchronises the ADERDGCellDescription
   * with the corresponding solver if this is required by the time stepping
   * scheme.
   * This operation must be performed in mergeWithNeighbour(...) and
   * touchVertexFirstTime(...) since both callbacks touch the
   * ADERDGCellDescriptions before the other callbacks.
   *
   * \note Not thread-safe.
   *
   * \param[in] cellDescription         The cell description
   * \param[in] face                    information about the boundary face
   * \note Not thread-safe.
   */
  void applyBoundaryConditions(CellDescription& p,Solver::BoundaryFaceInfo& face);

#ifdef Parallel
  /**
   * Data messages per neighbour communication.
   * This information is required by the sendEmpty...(...)
   * methods.
   */
  static const int DataMessagesPerNeighbourCommunication;
  /**
   * Data messages per fork/join communication.
   * This information is required by the sendEmpty...(...)
   * methods.
   */
  static const int DataMessagesPerForkOrJoinCommunication;
  /**
   * Data messages per master worker communication.
   * This information is required by the sendEmpty...(...)
   * methods.
   */
  static const int DataMessagesPerMasterWorkerCommunication;

  /**
   * Single-sided version of the other solveRiemannProblemAtInterface(). It
   * works only on one cell and one solver within this cell and in return
   * hands in the F and Q values explicitly through  indexOfQValues and
   * indexOfFValues. The Riemann solver is invoked and the bits are set
   * accordingly no matter of what they did hold before, i.e. different to
   * the standard solveRiemannProblemAtInterface() operation, we do not
   * check whether we shall run a Riemann solver or not.
   *
   * This method further synchronises the ADERDGCellDescription
   * with the corresponding solver if this is required by the time stepping
   * scheme.
   * This operation must be performed in mergeWithNeighbour(...) and
   * touchVertexFirstTime(...) since both callbacks touch the
   * ADERDGCellDescriptions before the other callbacks.
   *
   * \note Not thread-safe.
   */
  void solveRiemannProblemAtInterface(
      CellDescription& cellDescription,
      Solver::BoundaryFaceInfo& face,
      const double* const lQhbnd,
      const double* lFhbnd,
      const int fromRank);

  /**
   * Sets heap indices of an ADER-DG cell description to -1,
   * and the parent index of the cell descriptions to the specified \p
   * parentIndex.
   */
  static void resetIndicesAndFlagsOfReceivedCellDescription(CellDescription& p,const int parentIndex);

  /**
   * Allocate necessary memory and deallocate unnecessary memory.
   */
  static void ensureOnlyNecessaryMemoryIsAllocated(CellDescription& cellDescription);

  /** \copydoc Solver::prepareWorkerCellDescriptionAtMasterWorkerBoundary
   *
   * If the cell description is of type Descendant and
   * is next to a cell description of type Cell
   * or is virtually refined, i.e. has children of type Descendant itself,
   * we set the hasToHoldDataForMasterWorkerCommunication flag
   * on the cell description to true and allocate the required
   * memory.
   */
  static void prepareWorkerCellDescriptionAtMasterWorkerBoundary(
      CellDescription& cellDescription);

  /**
   * As the worker does not know anything about the master's coarse
   * grid cell, we set special child cell based erasing events
   * to notify the worker about the master's coarse grid cell's
   * erasing decision.
   */
  void deduceChildCellErasingEvents(CellDescription& cellDescription) const;

#endif

  /**
   * Determine average of each unknown
   *
   * We run over all sample points (or subcells in a Finite Volume context) and
   * determine the averages of each dof of the PDE. We assume that the arrays
   * first hold all the sample point values of the first PDE unknown. Therefore,
   * our outer loop runs over the PDE unknown and the inner ones run over the
   * sample points.
   *
   * Run over the persistent fields of the ADER-DG cell and determine the
   * average per unknown.' The result is stored within
   *
   * \note The fluctuations and update arrays do not store any material parameters.
   */
  void determineUnknownAverages(CellDescription& cellDescription) const;

  /**
   * Runs over all entries and adds sign times the average value. So if you
   * hand in a -1, you compute the hierarchical transform. If you hand in a +1,
   * you compute the inverse hierarchical transform.
   *
   * \note The fluctuations and update arrays do not store any material parameters.
   */
  void computeHierarchicalTransform(CellDescription& cellDescription, double sign) const;

  /**
   * This routine runs over the unknowns, asks the Heap's compression routines
   * to identify a reasonable compression level, stores this value and
   * afterwrads pipes the dofs into the byte stream. If you don't run with
   * assertions, the code does clear the double heap data afterwards. With
   * assertions, we leave it there and thus allow pullUnknownsFromByteStream()
   * to do quite some validation.
   */
  void putUnknownsIntoByteStream(CellDescription& cellDescription) const;

  /**
   *
   *
   * <h2>Multicore</h2>
   *
   * Unknowns are pulled from the input stream indirectly through
   * touchVertexFirstTime(). It always recycles heap data, so can be triggered
   * while other routines already do something with the cell description. There
   * usually are enough entries to recycle available.
   *
   * However, it may happen that we run out of recycled entries if we run into
   * a large regular subgrid for the first time. We can identify a run out as
   * we get a -1 from the heap. In this case, there are a couple of things to
   * do.
   *
   * - Wait for any background task to finish. Other parts of the grid might
   *   have triggered a compression in the background. So we have to wait for
   *   those guys to finish, as they rely on an invariant heap.
   * - Lock very pessimistically. No two operations (only touchVertexFirstTime
   *   calls should run in parallel, but I'm not 100% sure) should run.
   * - Create additional data.
   */
  void pullUnknownsFromByteStream(CellDescription& cellDescription) const;

  class CompressionJob: public tarch::multicore::jobs::Job {
    private:
      const ADERDGSolver& _solver;
      CellDescription&    _cellDescription;
      const bool          _isSkeletonJob;
    public:
      CompressionJob(
        const ADERDGSolver& solver,
        CellDescription&    cellDescription,
        const bool          isSkeletonJob);

      bool run() override;
  };

  /**
   * A job which performs the prediction and volume integral operations
   * for a cell description.
   */
  class PredictionJob: public tarch::multicore::jobs::Job {
    private:
      ADERDGSolver&    _solver; // not const because of kernels
      CellDescription& _cellDescription;
      const int        _cellDescriptionsIndex; // indices are used for identification of patch
      const int        _element;
      const double     _predictorTimeStamp;
      const double     _predictorTimeStepSize;
      const bool       _uncompressBefore;
      const bool       _isSkeletonJob;
    public:
      PredictionJob(
          ADERDGSolver&     solver,
          CellDescription&  cellDescription,
          const int         cellDescriptionsIndex,
          const int         element,
          const double      predictorTimeStamp,
          const double      predictorTimeStepSize,
          const bool        uncompressBefore,
          const bool        isAtRemoteBoundary);

      bool run() override;

      /**
       * We prefetch the data that is subject to the prediction/updates.
       * As we know that prefetchData on the i+1th tasks is ran before the
       * scheduler does the ith task, we end up with reasonably good
       * prefetching. In theory, we might assume that loading into the L3
       * cache is sufficient, while we found that loading into the L2 cache
       * can destroy the affinity. In practice, it is best to follow Intel's
       * cache recommendations, i.e. we use _MM_HINT_NTA as cache instruction
       * rather than _MM_HINT_T2. This is however empirical evidence and
       * might have to be reevaluated later on.
       */
      void prefetchData() override;
  };

  /**
   * A job which performs prolongation operation
   * for a cell description.
   */
  class ProlongationJob: public tarch::multicore::jobs::Job {
    private:
      ADERDGSolver&                            _solver; // not const because of kernels
      CellDescription&                         _cellDescription;
      const CellDescription&                   _parentCellDescription;
      const tarch::la::Vector<DIMENSIONS,int>  _subcellIndex;
    public:
      ProlongationJob(
          ADERDGSolver&                            solver,
          CellDescription&                         cellDescription,
          const CellDescription&                   parentCellDescription,
          const tarch::la::Vector<DIMENSIONS,int>& subcellIndex);

      bool run() override;
      void prefetchData() override;
  };

  /**
   * A job which performs a fused ADER-DG time step, i.e., it performs the solution update,
   * updates the local time stamp, and finally performs the space-time predictor commputation.
   *
   * TODO(Dominic): Minimise time step sizes and refinement requests per patch
   * (->transpose the typical minimisation order)
   */
  class FusedTimeStepJob: public tarch::multicore::jobs::Job {
    private:
      ADERDGSolver&                                             _solver; // TODO not const because of kernels
      CellDescription&                                          _cellDescription;
      CellInfo                                                  _cellInfo; // copy
      const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char> _neighbourMergePerformed; // copy
      const bool                                                _isFirstTimeStepOfBatch;
      const bool                                                _isLastTimeStepOfBatch;
      const bool                                                _isSkeletonJob;
    public:
      /**
       * Construct a FusedTimeStepJob.
       *
       * @note Job is spawned as high priority job if spawned in the last time step.
       * It further spawns a prediction job in this case in order
       * to overlap work with the reduction of time step size
       * and mesh update events.
       *
       * @note The state of the neighbourMergePerformed flags is used internally by
       * some of the kernels, e.g. in order to determine where to perform a face integral.
       * However, they have to be reset before the next iteration as they indicate on
       * which face a Riemann solve has already been performed or not (their original usage).
       * The flags are thus reset directly after spawning a FusedTimeStepJob.
       * Therefore, we need to copy the neighbourMergePerformed flags when spawning
       * a FusedTimeStep job.
       *
       * @param solver                 the spawning solver
       * @param cellDescription        a cell description
       * @param cellInfo               links to all cell descriptions associated with the cell
       * @param isFirstTimeStepOfBatch if we currently run the first time step of a batch
       * @param isLastTimeStepOfBatch  if we currently run the last time step of a batch
       * @param isSkeletonJob          if the cell is a skeleton cell
       */
      FusedTimeStepJob(
        ADERDGSolver&    solver,
        CellDescription& cellDescription,
        CellInfo&        cellInfo,
        const bool       isFirstTimeStepOfBatch,
        const bool       isLastTimeStepOfBatch,
        const bool       isSkeletonJob);

      bool run() override;
      void prefetchData() override;
  };

  /**
   * A job which performs the solution update and computes a new time step size.
   *
   * \note Spawning these operations as background job makes only sense if you
   * wait in endIteration(...) on the completion of the job.
   * It further important to flag this job as high priority job to
   * ensure completion before the next reduction.
   */
  class UpdateJob: public tarch::multicore::jobs::Job {
    private:
      ADERDGSolver&                                             _solver; // TODO not const because of kernels
      CellDescription&                                          _cellDescription;
      CellInfo                                                  _cellInfo;
      const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char> _neighbourMergePerformed; // copy
      const bool                                                _isAtRemoteBoundary;
    public:
      /**
       * Construct a UpdateJob.
       *
       * @note Job is always spawned as high priority job.
       *
       * @param solver                 the spawning solver
       * @param cellDescription        a cell description
       * @param cellInfo               links to all cell descriptions associated with the cell
       * @param isSkeletonJob          if the cell is a skeleton cell
       */
      UpdateJob(
        ADERDGSolver&    solver,
        CellDescription& cellDescription,
        CellInfo&        cellInfo,
        const bool       isAtRemoteBoundary);

      bool run() override;
      void prefetchData() override;
  };

  /**
   * A job that calls adjustSolutionDuringMeshRefinementBody(...).
   */
  class AdjustSolutionDuringMeshRefinementJob: public tarch::multicore::jobs::Job {
  private:
    ADERDGSolver&    _solver;
    CellDescription& _cellDescription;
    const bool       _isInitialMeshRefinement;
  public:
    AdjustSolutionDuringMeshRefinementJob(
        ADERDGSolver&    solver,
        CellDescription& cellDescription,
        const bool       isInitialMeshRefinement);

    bool run() override;
  };

public:

  /**
   * Compute a load balancing weight for a cell in the mesh.
   */
  static int computeWeight(const int cellDescriptionsIndex);

  /**
   * Push a new cell description to the back
   * of the heap vector at \p cellDescriptionsIndex.
   *
   * !!! Augmentation status
   *
   * For the grid setup, it is important that the
   * previous augmentation status is initialised for new cell
   * descriptions as MaximumAugmentationStatus.
   * This prevents erasing of vertices around newly introduced
   * cell descriptions of type Cell.
   *
   * Note that this is the previous augmentation status.
   * It does not spread.
   */
  static void addNewCellDescription(
      const int                                    solverNumber,
      CellInfo&                                    cellInfo,
      const CellDescription::Type                  cellType,
      const CellDescription::RefinementEvent       refinementEvent,
      const int                                    level,
      const int                                    parentIndex,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize,
      const tarch::la::Vector<DIMENSIONS, double>& cellOffset);

  /**
   * Returns the ADERDGCellDescription heap vector
   * at address \p cellDescriptionsIndex.
   */
  static Heap::HeapEntries& getCellDescriptions(
      const int cellDescriptionsIndex);

  /**
   * Returns the ADERDGCellDescription with index \p element
   * in the heap vector at address \p cellDescriptionsIndex.
   */
  static CellDescription& getCellDescription(
      const int cellDescriptionsIndex,
      const int element);

  /**
   * @return if a cell description (should) hold face data.
   *
   * @param cellDescription a cell description.
   */
  static bool holdsFaceData(const CellDescription& cellDescription);

  /**
   * \return if communication with a (remote) neighbour is necessary.
   *
   * @param cellDescription a cell description.
   * @param faceIndex       an index in the range 0 till DIMENSIONS_TIMES_TWO-1.
   */
  static bool communicateWithNeighbour(const CellDescription& cellDescription,const int faceIndex);

  /**
   * Prefetches Riemann input and output arrays (lQhbnd,LFbhnd) of an ADER-DG cell description.
   *
   * @param cellDescription an ADER-DG cell description
   * @param faceIndex       index in the range of 0 (inclusive) to 2*DIMENSIONS (exclusive)
   *                        numbering the faces of a cell.
   */
  static void prefetchFaceData(CellDescription& cellDescription,const int faceIndex);

  /**
   * Erase all cell descriptions registered for solvers
   * of type Type::ADERDG.
   */
  static void eraseCellDescriptions(const int cellDescriptionsIndex);

  void updateCommunicationStatus(
        exahype::solvers::ADERDGSolver::CellDescription& cellDescription) const;
  /**
   * Determine the communication status of this cell
   * description based on the face wise communication status flags
   * if the cell is of type Descendant.
   *
   * If the cell description is of type Ancestor, return 0.
   * If the cell description of type Cell, return the maximum
   * commmunication status.
   */
  int determineCommunicationStatus(
      exahype::solvers::ADERDGSolver::CellDescription& cellDescription) const;

  /**
   * TODO(Dominic): Add docu.
   */
  void updateAugmentationStatus(
      exahype::solvers::ADERDGSolver::CellDescription& cellDescription) const;

  /**
   * TODO(Dominic): Add docu.
   */
  int determineAugmentationStatus(
      exahype::solvers::ADERDGSolver::CellDescription& cellDescription) const;

  /**
   * Determine a new limiter status for the given direction based on the neighbour's
   * limiter status and the cell's reduced limiter status.
   */
  void mergeWithRefinementStatus(
      CellDescription& cellDescription,
      const int faceIndex,
      const int neighbourLimiterStatus) const;

  /**
   * Determine the refinement status from the face
   * neighbour values.
   *
   * \note It is very important that any troubled cell indicator
   * and any refinement criterion has been evaluated before
   * calling this function.
   */
  void updateRefinementStatus(
      CellDescription&                                           cellDescription,
      const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char>& neighbourMergePerformed) const;

  /**
   * Construct an ADERDGSolver.
   *
   * \param identifier               An identifier for this solver.
   * \param numberOfVariables        the number of variables.
   * \param numberOfParameters       the number of material parameters.
   * \param DOFPerCoordinateAxis     The 1D basis size, i.e. the order + 1.
   * \param maximumMeshSize          The maximum mesh size. From hereon, adaptive mesh refinement is used.
   * \param maximumAdaptiveMeshDepth The maximum depth of the adaptive mesh.
   * \param int DMPObservables       The number of discrete maximum principle observables. Has only
   *                                 a meaning in the context of limiting. Should be set to a value<=0
   *                                 if a pure ADER-DG solver is used.
   * \param timeStepping             the timestepping mode.
   * \param profiler                 a profiler.
   */
  ADERDGSolver(
      const std::string& identifier,
      const int numberOfVariables,
      const int numberOfParameters,
      const int basisSize,
      const double maximumMeshSize,
      const int maximumAdaptiveMeshDepth,
      const int haloCells,
      const int regularisedFineGridLevels,
      const exahype::solvers::Solver::TimeStepping timeStepping,
      const int limiterHelperLayers,
      const int DMPObservables,
      std::unique_ptr<profilers::Profiler> profiler =
          std::unique_ptr<profilers::Profiler>(
              new profilers::simple::NoOpProfiler("")));

  virtual ~ADERDGSolver() {}

  // Disallow copy and assignment
  ADERDGSolver(const ADERDGSolver& other) = delete;
  ADERDGSolver& operator=(const ADERDGSolver& other) = delete;

  /**
   * This operation returns the number of space time
   * unknowns per cell.
   *
   * Note that this operation might only have a meaning for space-time type
   * discretisation methods.
   */
  int getSpaceTimeUnknownsPerCell() const;

  /**
   * This operation returns the number of space time
   * flux unknowns per cell.
   *
   * Note that this operation might only have a meaning for space-time type
   * discretisation methods.
   */
  int getSpaceTimeFluxUnknownsPerCell() const;

  /**
   * This operation returns the number of unknowns per cell located in
   * the interior of a cell.
   */
  int getUnknownsPerCell() const;

  /**
   * This operation returns the number of flux unknowns per cell
   * located in the interior of a cell.
   */
  int getFluxUnknownsPerCell() const;

  /**
   * This operation returns the number of unknowns that are located
   * on or in the vicinity of the boundary of a cell.
   */
  int getUnknownsPerCellBoundary() const;

  /**
   * This operation returns the number of unknowns that are located
   * on or in the vicinity of each face of a cell.
   */
  int getUnknownsPerFace() const;


  /**
   * This operation returns the size of data required
   * to store face area based unknowns and associated parameters.
   *
   * \return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS - 1) * DIMENSIONS_TIMES_TWO;
   */
  int getDataPerCellBoundary() const;

  /**
   * This operation returns the size of data required
   * to store face area based unknowns and associated parameters.
   *
   * \return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS - 1);
   */
  int getDataPerFace() const;
  /**
   * This operation returns the size of data required
   * to store cell volume based unknowns and associated parameters.
   *
   * \return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS + 0);
   */
  int getDataPerCell() const;
  
  /**
   * This operation returns the size of data required
   * to store space-time cell unknowns and associated parameters.
   *
   * \return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS + 1);
   */
  int getSpaceTimeDataPerCell() const;

  /**
   * !!! LimitingADERDGSolver functionality !!!
   *
   * The number of observables
   * the discrete maximum principle
   * is applied to.
   */
  int getDMPObservables() const;

  /**
   * !!! LimitingADERDGSolver functionality !!!
   *
   */
  int getMinRefinementStatusForSeparationCell() const;

  /**
   * !!! LimitingADERDGSolver functionality !!!
   *
   */
  int getMinRefinementStatusForBufferCell() const;

  /**
   * !!! LimitingADERDGSolver functionality !!!
   *
   */
  int getMinRefinementStatusForTroubledCell() const;

  MeshUpdateEvent getNextMeshUpdateEvent() const final override;
  void updateNextMeshUpdateEvent(MeshUpdateEvent meshUpdateEvent) final override;
  void setNextMeshUpdateEvent() final override;
  MeshUpdateEvent getMeshUpdateEvent() const final override;
  void overwriteMeshUpdateEvent(MeshUpdateEvent newMeshUpdateEvent) final override;


  /**
   * Check if the heap array with index \p index could be allocated.
   */
  static void checkDataHeapIndex(const CellDescription& cellDescription, const int arrayIndex, const std::string arrayName);

  /**
   * Checks if no unnecessary memory is allocated for the cell description.
   * If this is not the case, it deallocates the unnecessarily allocated memory.
   *
   * \note This operation is thread safe as we serialise it.
   */
  void ensureNoUnnecessaryMemoryIsAllocated(CellDescription& cellDescription) const;

  /**
   * Checks if all the necessary memory is allocated for the cell description.
   * If this is not the case, it allocates the necessary
   * memory for the cell description.
   *
   * \note This operation is thread safe as we serialise it.
   *
   * \note Heap data creation assumes default policy
   * DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired.
   */
  void ensureNecessaryMemoryIsAllocated(exahype::records::ADERDGCellDescription& cellDescription) const;


  /**
   * Getter for the size of the array allocated that can be overriden
   * to change the allocated size independently of the solver parameters.
   * For example to add padding forthe optimised kernel
   */
  virtual int getTempSpaceTimeUnknownsSize()      const {return getSpaceTimeDataPerCell();} // TODO function should be renamed
  virtual int getTempSpaceTimeFluxUnknowns0Size() const {return getSpaceTimeFluxUnknownsPerCell();}
  virtual int getTempSpaceTimeFluxUnknowns1Size() const {return getSpaceTimeFluxUnknownsPerCell();}
  virtual int getTempUnknownsSize()               const {return getDataPerCell();} // TODO function should be renamed
  virtual int getTempFluxUnknownsSize()           const {return getFluxUnknownsPerCell();}
  virtual int getTempPointForceSourcesSize()      const {return (_nodesPerCoordinateAxis+1)*getUnknownsPerCell();}
  virtual int getBndFaceSize()                    const {return getDataPerFace();} // TODO function should be renamed
  virtual int getBndTotalSize()                   const {return getDataPerCellBoundary();} // TODO function should be renamed
  virtual int getBndFluxSize()                    const {return getUnknownsPerFace();} // TODO function should be renamed
  virtual int getBndFluxTotalSize()               const {return getUnknownsPerCellBoundary();} // TODO function should be renamed
  virtual int getUpdateSize()                     const {return getUnknownsPerCell();}

  virtual bool alignTempArray()                   const {return false;}

  /**
   * False for generic solver, may be true for optimized one
   * Used only for debug assertions
   */
  virtual bool usePaddedData_nVar() const {return false;}
  virtual bool usePaddedData_nDoF() const {return false;}
  

  /**
   * @brief Adds the solution update to the solution.
   *
   * \param[inout] luh    cell-local solution DoF.
   * \param[in]    luhOld the old solution a
   * \param[in]    lduh   cell-local update DoF.
   * \param[dt]    dt     time step size.
   */
  virtual void solutionUpdate(double* luh, const double* const luhOld, const double* const lduh, const double dt) = 0;


  /**
   * @brief Computes a face integral contributions
   * to the cell update.
   *
   * In case of \p levelDelta > 0, the kernel needs to
   * restrict the given boundary-extrapolated flux DoF
   * \p levelDelta levels up before performing the face integral.
   *
   * \param[inout] lduh         Cell-local update DoF.
   * \param[in]    lFhbnd       Cell-local DoF of the boundary extrapolated fluxes for the face
   *                            with the given direction and the given geometry.
   * \param[in]    direction    Coordinate direction the normal vector is aligned with.
   * \param[in]    orientation  Orientation of the normal vector (0: negative sign, 1: positive sign).
   * \param[in[    levelDelta   The difference in levels up to a cell description of type Cell.
   *                            Must be set to zero if we are already performing a face integral for a cell description
   *                            of type Cell. Is greater zero if lFbhnd stems from a Descendant cell description.
   * \param[in]    cellSize     Extent of the cell in each coordinate direction.
   */
  virtual void faceIntegral(
      double* const       lduh,
      const double* const lFhbnd,
      const int direction, const int orientation,
      const tarch::la::Vector<DIMENSIONS-1,int>& subfaceIndex, const int levelDelta,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize) = 0;

  /**
   * @brief Computes the normal fluxes (or fluctuations) at the interface of two
   *cells.
   *
   * \param[inout] FL             Flux DoF belonging to the left cell.
   * \param[inout] FR             Flux DoF belonging the right cell.
   * \param[in]    QL             DoF of the boundary extrapolated predictor
   *                              belonging to the left cell.
   * \param[in]    QR             DoF of the boundary extrapolated predictor
   *                              belonging to the right cell.
   * \param[in]    direction  Index of the nonzero normal vector component,
   *               i.e., 0 for e_x, 1 for e_y, and 2 for e_z.
   */
  virtual void riemannSolver(double* FL, double* FR,
                             const double* const QL,const double* const QR,
                             const double t,
                             const double dt,
                             const int normalNonZero,
                             bool isBoundaryFace,
                             int faceIndex) = 0;

  /**
   * Impose boundary conditions on the fluxes (or fluctuations).
   * The state is only read.
   *
   * \param[inout] fluxIn        boundary-extrapolated (space-time) volume flux.
   *                             Can be overwritten/reused as it is updated anyway after
   *                             the next predictor computation.
   * \param[in]    stateIn       boundary-extraplolated (space-time) predictor
   * \param[in]    luh           the solution values in the patch, for read acces
   * \param[in]    cellCentre    cell centre.
   * \param[in]    cellSize      cell size.
   * \param[in]    t             The time.
   * \param[in]    dt            a time step size.
   * \param[in]    direction     index of the nonzero component of the normal vector
   *                             i.e., 0 for e_x, 1 for e_y, and 2 for e_z.
   * \param[in]    orientation   orientation of the normal vector where 0 means negative and
   *                             1 means positive.
   *
   *
   */
  virtual void boundaryConditions(
      double* const fluxIn,
      const double* const stateIn,
      const double* const luh,
      const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
      const double t,const double dt,
      const int direction,
      const int orientation) = 0;

  /**
   * @brief Computes cell-local space-time predictor, volume, and face DoF
   * and performs volume integral.
   *
   * The space-time predictor computation also includes
   * evaluating the point sources.
   *
   * \param[out] lduh           cell-local update DoF.
   * \param[out] lQhbnd, lFhbnd boundary-extrapolated space-time predictor and volume flux values.
   * \param[int] luh            solution DoF.
   * \param[in]  invDx          inverted extent of the cell per coordinate direction.
   * \param[in]  dt             time step size.
   *
   * \return the number of Picard iterations performed by the
   * space-time predictor computation kernel.
   */
  virtual int fusedSpaceTimePredictorVolumeIntegral(
      double* lduh, double*  lQhbnd, double* lFhbnd,
      const double* const luh,
      const tarch::la::Vector<DIMENSIONS, double>& center,
      const tarch::la::Vector<DIMENSIONS, double>& dx,
      const double t,
      const double dt) = 0;

  /**
   * \brief Returns a stable time step size.
   *
   * \param[in] luh             Cell-local solution DoF.
   * \param[in] cellSize        Extent of the cell in each coordinate direction.
   */
  virtual double stableTimeStepSize(
      const double* const luh,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize) = 0;

  /**
   * This operation allows you to impose time-dependent solution values
   * as well as to add contributions of source terms.
   * Please be aware that this operation is called per time step if
   * the corresponding predicate hasToUpdateSolution() yields true for the
   * region and time interval.
   *
   * \param tNew  The new time stamp after the solution update.
   * \param dtOld The time step size that was used to update the solution.
   *              This time step size was computed based on the old solution.
   *              If we impose initial conditions, i.e, t=0, this value
   *              equals 0.
   */
  virtual void adjustSolution(
      double* luh, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& dx,
      const double t,
      const double dtOld) = 0;

  /**
   * @defgroup AMR Solver routines for adaptive mesh refinement
   */
  ///@{
  /**
   * The refinement criterion that must be defined by the user.
   *
   */
  // @todo: 16/04/06:Dominic Etienne Charrier Consider to correct the level in
  // the invoking code, i.e., level-> level-1
  // since this is was the user expects.
  virtual exahype::solvers::Solver::RefinementControl refinementCriterion(
      const double* luh, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize,
      const double time,
      const int level) = 0;

  /**
   * Project coarse grid face unknowns
   * on level \p coarseGridLevel down to level \p fineGridLevel
   * and writes them to the fine grid unknowns
   *
   * \note For the considered AMR concept, the difference in levels can
   * be larger than one. Let \f$l\f$ be the level difference. The
   * vector \p subfaceIndex does contain values in the range
   * \f$0,1,\ldots,3^l-1\f$.
   */
  virtual void faceUnknownsProlongation(
      double* lQhbndFine, double* lFhbndFine, const double* lQhbndCoarse,
      const double* lFhbndCoarse, const int coarseGridLevel,
      const int fineGridLevel,
      const tarch::la::Vector<DIMENSIONS - 1, int>& subfaceIndex) = 0;

  /**
   * Project coarse grid face unknowns
   * on level \p coarseGridLevel down to level \p fineGridLevel
   * and writes them to the fine grid unknowns
   *
   * \note For the considered AMR concept, the difference in levels is always
   * equal to one. The vector \p subcellIndex does contain values in the range
   * \f$0,1,2\f$.
   */
  virtual void volumeUnknownsProlongation(
      double* luhFine, const double* luhCoarse, const int coarseGridLevel,
      const int fineGridLevel,
      const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) = 0;

  /**
   * Restricts fine grid volume unknowns on level \p fineGridLevel
   * up to level \p coarseGridLevel and adds them to the coarse grid unknowns.
   *
   * \note For the considered AMR concept, the difference in levels is always
   * equal to one. The vector \p subcellIndex does contain values in the range
   * \f$0,1,2\f$.
   */
  virtual void volumeUnknownsRestriction(
      double* luhCoarse, const double* luhFine, const int coarseGridLevel,
      const int fineGridLevel,
      const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) = 0;
  ///@}
  


  /**
   * A criterion determining if the degrees of freedoms of
   * the cell-wise solution are physically admissible.
   *
   * \note We require that the cell-local minimum and maximum
   * of the solution values has been computed
   * a-priori.
   *
   * This operation is required for limiting.
   */
  virtual bool isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,
      const bool wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t, const double dt
      ) const = 0;

  /**
   * Maps the solution values Q to
   * the discrete maximum principle observables.
   *
   * As we can observe all state variables,
   * we interpret an 'observable' here as
   * 'worthy to be observed'.
   *
   *\param[inout] observables The mapped observables.
   *\param[in]                numberOfObservables The number of observables.
   *\param[in]    Q           The state variables.
   */
  virtual void mapDiscreteMaximumPrincipleObservables(
    double* observables,
    const int numberOfObservables,
    const double* const Q) const = 0;

  /**
   * Deduces a time stamp and time step size for the two cells depending
   * on the time stepping strategy and the data on the cell descriptions.
   *
   * @note The arguments are currently unused but they will be used
   * when we interpolate face values between neighbouring cells.
   *
   * @param cellDescription1 one of the cell descriptions
   * @param cellDescription2 one of the cell descriptions
   */
  std::tuple<double,double> getRiemannSolverTimeStepData(
      const CellDescription& cellDescription1,
      const CellDescription& cellDescription2) const;

  /**
   * Copies the time stepping data from the global solver onto the patch's time
   * stepping data.
   */
  void synchroniseTimeStepping(CellDescription& p) const;

  /**
   * \copydoc Solver::startNewTimeStep
   *
   * Update and reset corrector and predictor
   * time stamps and time step sizes according to the chosen
   * time stepping variant.
   *
   * Further reset the minimum and maximum cell sizes
   * to numeric limit values.
   *
   * \note The minimum and maximum cell sizes do
   * not need to be reset to numeric limit values
   * in every time step for uniform mesh refinement
   * static adaptive mesh refinement
   * but we still do it since we want to
   * utilise dynamic adaptive mesh refinement
   * since we want to dynamic adaptive mesh refinement
   * eventually.
   */
  void startNewTimeStep() override;

  void startNewTimeStepFused(
      const bool isFirstTimeStepOfBatch,
      const bool isLastTimeStepOfBatch) final override;

  /**
   * \copydoc Solver::updateTimeStepSizes
   *
   * Does not advance the predictor time stamp in time.
   */
  void updateTimeStepSizes() override;

  /** \copydoc Solver::updateTimeStepSizesFused
   *
   * Does advance the predictor time stamp in time.
   *
   * We further reset the stability condition was violated condition
   * to false for the ADER-DG solver.
   */
  void updateTimeStepSizesFused() override;

  void rollbackToPreviousTimeStep() final override;

  /**
   * \copydoc Solver::rollbackToPreviousTimeStepFused()
   *
   * <h2> Fused ADER-DG Time Stepping </h2>
   *
   * Corrector time stamp and corrector time step size must
   * add up to predictor time stamp after rollback.
   *
   * Corrector time step size is assumed to be used
   * predictor time step size in batch.
   */
  void rollbackToPreviousTimeStepFused() final override;

  /**
   * Update predictor time step size
   *
   * This operation takes the minimum of the current predictor time step size
   * and the argument handed in. The routine is used in
   * TimeStepComputation to determine the subsequent time step size.
   *
   * <h1>Globalfixed timestepping</h1>
   * In case of global fixed timestepping,
   * we rely on the initial condition
   * _predictorTimeStamp==_correctorTimeStamp
   * to detect the first time step size
   * computation.
   *
   * <h1>Thread-safety</h1>
   *
   * This operation is not thread safe.
   *
   */
  void updateMinNextPredictorTimeStepSize(
      const double& nextPredictorTimeStepSize);

  /**
   * Currently, required by TimeStepSizeComputation.
   */
  void setMinPredictorTimeStepSize(const double value);

  double getMinNextPredictorTimeStepSize() const;
  double getMinPredictorTimeStepSize() const;
  double getMinPredictorTimeStamp() const;
  double getMinCorrectorTimeStamp() const;
  double getMinCorrectorTimeStepSize() const;
  double getPreviousMinCorrectorTimeStamp() const;
  double getPreviousMinCorrectorTimeStepSize() const;

  double getMinTimeStamp() const override;
  double getMinTimeStepSize() const override;
  double getMinNextTimeStepSize() const override;

  void updateMinNextTimeStepSize(double value) override;

  /**
   * Set if the CFL condition was violated
   * (by the last fused time step).
   */
  void setStabilityConditionWasViolated(bool state);

  /**
   * \return true if the CFL condition was violated
   * (by the last fused time step).
   */
  bool getStabilityConditionWasViolated() const;

  /**
    * User defined solver initialisation.
    *
    * \param[in] cmdlineargs the command line arguments.
    */
  virtual void init(
      const std::vector<std::string>& cmdlineargs,
      const exahype::parser::ParserView& constants) = 0;

  void initSolver(
      const double timeStamp,
      const tarch::la::Vector<DIMENSIONS,double>& domainOffset,
      const tarch::la::Vector<DIMENSIONS,double>& domainSize,
      const tarch::la::Vector<DIMENSIONS,double>& boundingBoxSize,
      const std::vector<std::string>& cmdlineargs,
      const exahype::parser::ParserView& parserView) override;

  bool isPerformingPrediction(const exahype::State::AlgorithmSection& section) const override;
  bool isMergingMetadata(const exahype::State::AlgorithmSection& section) const override;

  static bool isValidCellDescriptionIndex(const int cellDescriptionsIndex);

  int tryGetElement(
      const int cellDescriptionsIndex,
      const int solverNumber) const override;

  ///////////////////////////////////
  // MODIFY CELL DESCRIPTION
  ///////////////////////////////////

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
      const int solverNumber,
      const bool stillInRefiningMode) override;

  exahype::solvers::Solver::RefinementControl eraseOrRefineAdjacentVertices(
        const int cellDescriptionsIndex,
        const int solverNumber,
        const tarch::la::Vector<DIMENSIONS, double>& cellOffset,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize,
        const bool checkThoroughly) const final override;

  /**\copydoc Solver::attainedStableState
   *
   * Compute flagging gradients in inside cells.
   * If the facewise flags on two opposite sides differ
   * by more than 2, then the flagging has not converged.
   *
   * If this is the case or if the refinement events
   * of a cell are none or the refinement criterion was not
   * evaluated yet, we say the solver has not attained
   * a stable state yet.
   */
  bool attainedStableState(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      const int solverNumber) const override;

  void finaliseStateUpdates(const int solverNumber,CellInfo& cellInfo) override;

  ///////////////////////////////////
  // CELL-LOCAL
  ///////////////////////////////////

  /**
   * Evaluate the refinement criterion after
   * a solution update has been performed and
   * the patch has been advanced in time.
   *
   * We currently only return true if a cell requested refinement.
   * ExaHyPE might then stop the
   * time stepping and update the mesh
   * before continuing. Erasing is here not considered.
   *
   * \note Must be called after startNewTimeStep was called
   *
   * \return True if mesh refinement is requested.
   *
   * \note Has no const modifier since kernels are not const functions yet.
   */
  MeshUpdateEvent evaluateRefinementCriteriaAfterSolutionUpdate(
      CellDescription&                                           cellDescription,
      const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char>& neighbourMergePerformed);

  /*! Perform prediction and volume integral.
   *
   * @note Different to the overloaded method, this method
   * waits for completion of a cell's last operation.
   *
   * @note Requests to uncompress the cell description arrays before computing the predictor
   *
   * @note This function must not be called from a task/thread as it
   * synchronises time step data.
   */
  void performPredictionAndVolumeIntegral(
      const int  solverNumber,
      CellInfo&  cellInfo,
      const bool isAtRemoteBoundary);

  /**
   * Computes the space-time predictor quantities, extrapolates fluxes
   * and (space-time) predictor values to the boundary and
   * computes the volume integral directly afterwards.
   * Furthermore, it restricts face data up to coarser grids
   * and compresses the cell description data again.
   *
   * Can be configured to uncompress the cell description
   * arrays before computing the space-time predictor quantities.
   *
   * <h2>Background Jobs</h2>
   *
   * All FusedTimeStepJob and PredictionJob instances finish here.
   * This is where an ADER-DG time step ends.
   * We thus call cellDescription.setHasCompletedLastStep(true) at
   * the end of this function.
   *
   * \param[in] uncompressBefore             uncompress the cell description arrays before computing
   *                                         the space-time predictor quantities.
   * \param[in] vetoCompressionBackgroundJob veto that the compression is run as a background task
   *
   * @note Might be called by background task. Do not synchronise time step data here.
   */
  void performPredictionAndVolumeIntegralBody(
      CellDescription& cellDescription,
      const double predictorTimeStamp,
      const double predictorTimeStepSize,
      const bool   uncompressBefore,
      const bool   isSkeletonCell );

  /**
   * \note uncompress is not performed in this routine. It must
   * be called before calling this routine if compression is employed.
   *
   * \note Has no const modifier since kernels are not const functions.
   *
   * \param[in] isAtRemoteBoundary indicates that we are at a remote boundary.
   *                               Plays a role in filtering out cells where we cannot
   *                               start backgroudn tasks.
   */
  void performPredictionAndVolumeIntegral(
      const int    solverNumber,
      CellInfo&    cellInfo,
      const double predictorTimeStamp,
      const double predictorTimeStepSize,
      const bool   uncompress,
      const bool   isAtRemoteBoundary);

  /**
   * Validate that the data stored on and for
   * the cell description is valid.
   *
   * \note Must only be called if the compression
   * is currently not in progress, i.e. processed as
   * a background task.
   */
  void validateCellDescriptionData(
      const CellDescription& cellDescription,
      const bool validateTimeStepData,
      const bool afterCompression,
      const bool beforePrediction,
      const std::string& methodTraceOfCaller) const;

  /**
   * Computes the bounds of the inputs and collects
   * them in a message.
   *
   * @param QL left Riemann input state
   * @param QR right Riemann input state
   * @param FL left Riemann input fluxes
   * @param FR right Riemann input fluxes
   */
  std::string riemannDataToString(
      const double* const Q,const double* const F,std::string suffix) const;

  /**
   * Computes a time step size based on the solution
   * values. Does not advance the cell descriptions
   * time stamps forward.
   *
   * \return the newly computed time step size if the cell description is of type Cell or the maximum
   *         double value.
   */
  double computeTimeStepSize(
      CellDescription& cellDescription);


  /**
   * Update the time stamps and time step sizes
   * for a cell description when nonfused time stepping
   * is used.
   *
   * \note No const modifier as kernels are not const yet.
   */
  double startNewTimeStep(CellDescription& cellDescription);

  /**
   * Same as \p startNewTimeStep for the fused time stepping scheme.
   *
   * \param[in] isFirstTimeStepOfBatch indicates that we are in the first iteration
   *                                    of a batch or not. Note that this must be also set to true
   *                                    in case we run a batch of size 1, i.e. no batch at all.
   *
   * \param[in] isLastTimeStepOfBatch indicates that we are in the last iteration
   *                                    of a batch or not. Note that this must be also set to true
   *                                    in case we run a batch of size 1, i.e. no batch at all.
   *
   * \note No const modifier as kernels are not const yet.
   */
  double startNewTimeStepFused(
      CellDescription& cellDescription,
      const bool isFirstTimeStepOfBatch,
      const bool isLastTimeStepOfBatch);

  /** \copydoc Solver::updateTimeStepSizes
   *
   * @param solverNumber      identification number of this solver
   * @param cellInfo          refers to a cell's data
   * @param fusedTimeStepping advances the predictor time stamp in time.
   */
  double updateTimeStepSizes(const int solverNumber,CellInfo& cellInfo,const bool fusedTimeStepping) final override;

  /**
   * Perform a fused time step, i.e. perform the update, update time step data, mark
   * for refinement and then compute the new space-time predictor.
   *
   * Order of operations
   * -------------------
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
   * @param cellDescriptionsIndex cell description heap index
   * @param element               element in the cell description heap array
   *
   * @note Might be called by background task. Do not synchronise time step data here.
   */
  UpdateResult fusedTimeStepBody(
        CellDescription&                                           cellDescription,
        CellInfo&                                                  cellInfo,
        const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char>& neighbourMergePerformed,
        const bool                                                 isFirstTimeStepOfBatch,
        const bool                                                 isLastTimeStepOfBatch,
        const bool                                                 isSkeletonCell,
        const bool                                                 mustBeDoneImmediately);

  /**
   * If the cell description is of type Cell, update the solution, evaluate the refinement criterion,
   * and compute an admissible time step size.
   *
   * @note Not const as kernels are not const.
   *
   * @param cellDescription a cell description
   * @return a struct containing a mesh update event triggered by this cell,
   * and a new time step size.
   *
   * @note Might be called by background task. Do not synchronise time step data here.
   */
  UpdateResult updateBody(
      CellDescription&                                           cellDescription,
      CellInfo&                                                  cellInfo,
      const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char>& neighbourMergePerformed,
      const bool                                                 isAtRemoteBoundary);

  UpdateResult fusedTimeStepOrRestrict(
      const int  solverNumber,
      CellInfo&  cellInfo,
      const bool isFirstTimeStepOfBatch,
      const bool isLastTimeStepOfBatch,
      const bool isAtRemoteBoundary) final override;

  UpdateResult updateOrRestrict(
      const int  solverNumber,
      CellInfo&  cellInfo,
      const bool isAtRemoteBoundary) final override;

  void compress(
      const int solverNumber,
      CellInfo& cellInfo,
      const bool isAtRemoteBoundary) const final override;

  void adjustSolutionDuringMeshRefinement(const int solverNumber,CellInfo& cellInfo) final override;

  /**
   * Computes the surface integral contributions to the
   * cell update and then adds the update degrees
   * to the solution degrees of freedom.
   *
   * Solution adjustments
   * --------------------
   *
   * After the update, the solution is at time
   * cellDescription.getCorrectorTimeStamp() + cellDescription.getCorrectorTimeStepSize().
   * The value cellDescription.getCorrectorTimeStepSize()
   * handed to the solution adjustment function is the one
   * used to update the solution.
   *
   * @todo We will not store the update field anymore
   * but a previous solution. We will thus only perform
   * a solution adjustment and adding of source term contributions here.
   *
   * @param[in] cellDescription         a cell description
   * @param[in] backupPreviousSolution  Set to true if the solution should be backed up before
   *                                    we overwrite it by the updated solution.
   */
  void updateSolution(
      CellDescription&                                           cellDescription,
      const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char>& neighbourMergePerformed,
      const bool                                                 backupPreviousSolution);

  /**
   * TODO(Dominic): Update docu.
   *
   * Old docu:
   *
   * Rolls back the solver's solution on the
   * particular cell description.
   * This method is used by the ADER-DG a-posteriori
   * subcell limiter.
   *
   * Uses the corrector time step size to perform the rollback.
   * Thus make sure to invoke ::rollbackToPreviousTimeStepSize() beforehand
   * if the patch has already advanced to next time step.
   *
   * <h2>Open issues</h2>
   * A rollback is of course not possible if we have adjusted the solution
   * values. Assuming the rollback is invoked by a LimitingADERDGSolver,
   * we should use the adjusted FVM solution as reference solution.
   * A similar issue occurs if we impose initial conditions that
   * include a discontinuity.
   */
  void swapSolutionAndPreviousSolution(CellDescription& cellDescription) const;

  /**
   * Prolongates face data from a parent cell description to
   * the cell description at address (cellDescriptionsIndex,element)
   * in case the fine grid cell associated with the cell description is adjacent to
   * the hull of the coarse grid cell associated with the parent cell description.
   *
   * Further zero out the face data of ancestors.
   *
   * \note This function assumes a top-down traversal of the grid and must thus
   * be called from the enterCell(...) mapping method.
   *
   * \note It is assumed that this operation is applied only to helper cell descriptions
   * of type Descendant and Ancestor. No cell description of type Cell
   * must be touched by this operation. Otherwise, we cannot spawn
   * the prediction and/or the compression as background task.
   *
   * \note Has no const modifier since kernels are not const functions yet.
   */
  void prolongateFaceData(
      const int solverNumber,
      CellInfo& cellInfo,
      const bool isAtRemoteBoundary);

  /**
   * Restrict face data to the top most parent which has allocated face data arrays (Cell)
   * if and only if the fine grid cell (Descendant) has a face which intersects with one of the top most parent
   * cell's faces.
   *
   * \note This function assumes a bottom-up traversal of the grid and must thus
   * be called from the leaveCell(...) or ascend(...) mapping methods.
   *
   * \note Has no const modifier since kernels are not const functions.
   */
  void restrictToTopMostParent(const CellDescription& cellDescription);

  /**
   * Go back to previous time step with
   * time step data, solution, and refinement status.
   *
   * Allocate necessary new limiter patches.
   */
  void rollbackSolutionGlobally(
      const int solverNumber,
      CellInfo& cellInfo,
      const bool fusedTimeStepping) const final override;

  ///////////////////////////////////
  // NEIGHBOUR
  ///////////////////////////////////
  // helper status
  void mergeWithCommunicationStatus(
      CellDescription& cellDescription,
      const int faceIndex,
      const int otherCommunicationStatus) const;

  // augmentation status
  void mergeWithAugmentationStatus(
      CellDescription& cellDescription,
      const int faceIndex,
      const int otherAugmentationStatus) const;

  void mergeNeighboursMetadata(
      const int                                 solverNumber,
      Solver::CellInfo&                         cellInfo1,
      Solver::CellInfo&                         cellInfo2,
      const tarch::la::Vector<DIMENSIONS, int>& pos1,
      const tarch::la::Vector<DIMENSIONS, int>& pos2,
      const tarch::la::Vector<DIMENSIONS,       double>& x,
      const tarch::la::Vector<DIMENSIONS,       double>& h,
      const bool                                checkThoroughly) const;

  void mergeNeighboursData(
      const int                                 solverNumber,
      Solver::CellInfo&                         cellInfo1,
      Solver::CellInfo&                         cellInfo2,
      const tarch::la::Vector<DIMENSIONS, int>& pos1,
      const tarch::la::Vector<DIMENSIONS, int>& pos2);

  void mergeWithBoundaryData(
      const int                                 solverNumber,
      Solver::CellInfo&                         cellInfo,
      const tarch::la::Vector<DIMENSIONS, int>& posCell,
      const tarch::la::Vector<DIMENSIONS, int>& posBoundary);
#ifdef Parallel
  /**
   * Sends all the cell descriptions at address \p
   * cellDescriptionsIndex to the rank \p toRank.
   *
   * <h2>Adaptive mesh refinement</h2>
   * For adaptive meshes, we further fix the type
   * of a descendant to RemoteBoundaryDescendant
   * at both sides of master-worker boundaries.
   *
   * We further fix the type of an Ancestor
   * to RemoteBoundaryAncestor if the parent
   * of the cell description on the master side
   * is also of type RemoteBoundaryAncestor or an
   * Ancestor.
   *
   * \note The data heap indices of the cell descriptions are not
   * valid anymore on rank \p toRank.
   *
   * \param fromWorkerSide Indicates that we sent these cell descriptions from the
   *                       worker side, e.g. during a joining operation.
   */
  static bool sendCellDescriptions(
      const int                                    toRank,
      const int                                    cellDescriptionsIndex,
      const bool                                   fromWorkerSide,
      const peano::heap::MessageType&              messageType,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level);

  /**
   * Sends an empty message to the rank \p toRank.
   */
  static void sendEmptyCellDescriptions(
      const int                                    toRank,
      const peano::heap::MessageType&              messageType,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level);

  /**
   * Receives cell descriptions from rank \p fromRank
   * and resets the data heap indices to -1.
   *
   * If a received cell description has the same
   * solver number as a cell description in the
   * array at address \p cellDescriptionsIndex,
   * we merge the metadata (time stamps, time step size)
   * of both cell descriptions.
   *
   * If no cell description in the array at address
   * \p cellDescriptionsIndex can be found with the
   * same solver number than a received cell description,
   * we push the received cell description to
   * the back of the array at address \p cellDescriptions
   * Index.
   *
   * This operation is intended to be used in combination
   * with the solver method mergeWithWorkerOrMasterDataDueToForkOrJoin(...).
   * Here, we would merge first the cell descriptions sent by the master and worker
   * and then merge the data that is sent out right after.
   */
  static void receiveCellDescriptions(
      const int                                    fromRank,
      exahype::Cell&                               localCell,
      const peano::heap::MessageType&              messageType,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level);

  /**
   * Drop cell descriptions received from \p fromRank.
   */
  static void dropCellDescriptions(
      const int                                    fromRank,
      const peano::heap::MessageType&              messageType,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level);

  ///////////////////////////////////
  // NEIGHBOUR
  ///////////////////////////////////
  /** \copydoc Solver::mergeWithNeighbourMetadata
   *
   * Appends cell type,limiterStatus,augmentationStatus,
   * and communicationStatus to \p metadata.
   */
  void appendNeighbourCommunicationMetadata(
      exahype::MetadataHeap::HeapEntries& metadata,
      const tarch::la::Vector<DIMENSIONS,int>& src,
      const tarch::la::Vector<DIMENSIONS,int>& dest,
      const int cellDescriptionsIndex,
      const int solverNumber) const override;

  /** \copydoc Solver::mergeWithNeighbourMetadata
   *
   * Merges with a metadata message received from
   * a neighbour. The message contains the neighbours
   * cell type,limiterStatus,augmentationStatus,communicationStatus.
   *
   * @note That at metadata merge can happen up to 2 (2D) or 4 (3D) times
   * for the same face as the counter is not used here.
   *
   * <h2>LimitingADERDGSolver</h2>
   * This routine also merges the cell's limiter status
   * with the one of the neighour.
   * We do this here in order to reduce code bloat.
   */
  void mergeWithNeighbourMetadata(
      const int                                 solverNumber,
      CellInfo&                                 cellInfo,
      const MetadataHeap::HeapEntries&          neighbourMetadata,
      const tarch::la::Vector<DIMENSIONS, int>& src,
      const tarch::la::Vector<DIMENSIONS, int>& dest) const;

  /**
   * Sends out two messages, one holding degrees of freedom (DOF)
   * of the boundary-extrapolated space-time predictor and one
   * holding DOF of the boundary-extrapolated space-time flux.
   *
   * <h2>LimitingADERDGSolver's min and max</h2>
   * This method does not send the minimum and
   * maximum values required for the
   * LimitingADERDGSolver's discrete h2>maximum principle.
   * The LimitingADERDGSolver does this in his
   * LimitingADERDGSolver::sendDataToNeighbour method.
   *
   * Min and max have to be merge
   * independent of the limiter status of the cell while
   * a ADER-DG neighbour merge has to be performed
   * only for cells with certain limiter status
   * flags.
   *
   * @param toRank       the adjacent rank we want to send to
   * @param solverNumber identification number for the solver
   * @param cellInfo     links to a cell's data
   * @param src          position of message source relative to vertex
   * @param dest         position of message destination relative to vertex
   * @param x            vertex' position
   * @param level        vertex' level
   */
  void sendDataToNeighbour(
      const int                                     toRank,
      const int                                     solverNumber,
      Solver::CellInfo&                             cellInfo,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level);

  /** \copydoc Solver::mergeWithNeighbourData
   *
   * <h2>LimitingADERDGSolver's min and max</h2>
   * This method does not merge the minimum and
   * maximum values required for the
   * LimitingADERDGSolver's discrete maximum principle.
   * The LimitingADERDGSolver does this in his
   * LimitingADERDGSolver::mergeWithNeighbourData method.
   *
   * Min and max have to be merged
   * independent of the limiter status of the cell while
   * a ADER-DG neighbour merge has to be performed
   * only for cells with certain limiter status
   * flags.
   */
  void mergeWithNeighbourData(
      const int                                    fromRank,
      const int                                    solverNumber,
      Solver::CellInfo&                            cellInfo,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level);

  /** \copydoc Solver::dropNeighbourData
   *
   * <h2>LimitingADERDGSolver's min and max</h2>
   * This method does not drop the minimum and
   * maximum values required for the
   * LimitingADERDGSolver's discrete maximum principle.
   * The LimitingADERDGSolver does this in his
   * LimitingADERDGSolver::dropNeighbourData method.
   */
  void dropNeighbourData(
      const int                                    fromRank,
      const int                                    solverNumber,
      Solver::CellInfo&                            cellInfo,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const;

  ///////////////////////////////////
  // MASTER<=>WORKER
  ///////////////////////////////////
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


  void sendDataToWorkerIfProlongating(
      const int                                     toRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) const final override;

  /**
   * Just receive data depending on the refinement
   * event of a cell description.
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
    * Veto erasing requests of the coarse grid cell if the received
    * cell description has virtual children.
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
       const int                                    level,
       const bool                                   stillInRefiningMode) final override;

  void appendMasterWorkerCommunicationMetadata(
      MetadataHeap::HeapEntries& metadata,
      const int cellDescriptionsIndex,
      const int solverNumber) const override;

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
  /**
   * Compiles a message for the master.
   *
   * Capacity (in byte) of the message vector can be modified
   * in case the calling function wants to push additional
   * entries to the back of the vector.
   *
   * \see LimitingADERDGSolver::sendDataToMaster
   */
  DataHeap::HeapEntries
  compileMessageForMaster(const int capacity=4) const;

  void sendDataToMaster(
      const int masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const override;

  /**
   * Read a message from a worker
   * and adjust solver fields.
   *
   * \see LimitingADERDGSolver::mergeWithWorkerData
   */
  void mergeWithWorkerData(
      const DataHeap::HeapEntries& message);

  void mergeWithWorkerData(
      const int workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  ///////////////////////////////////
  // MASTER->WORKER
  ///////////////////////////////////
  /**
   * Compiles a message for a worker.
   *
   * Capacity of the message vector can be modified
   * in case the calling function wants to push additional
   * entries to the back of the vector.
   *
   * \see LimitingADERDGSolver::sendDataToWorker
   */
  DataHeap::HeapEntries
  compileMessageForWorker(const int capacity=7) const;

  void sendDataToWorker(
      const                                        int workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const override;

  /**
   * Read a message from the master
   * and adjust solver fields.
   *
   * \see LimitingADERDGSolver::mergeWithMasterData
   */
  void mergeWithMasterData(const DataHeap::HeapEntries& message);

  void mergeWithMasterData(
      const                                        int masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;
#endif

  std::string toString() const override;

  void toString (std::ostream& out) const override;

  /**
   * The counterpart of uncompress.
   *
   * <h2> Shared memory parallelisation </h2>
   *
   * Different to the compression, we don't have to take care about any races:
   * the compression is invoked by enterCell or leaveCell respectively, i.e.
   * exactly once per cell. This can happen in parallel for multiple cells
   * which is fine.
   *
   * However, we have to take care about the interplay of compression and
   * uncompression.
   *
   * \param[in] isSkeletonJob decides to which queue we spawn the job if we spawn any
   */
  void compress( CellDescription& cellDescription, const bool isSkeletonCell ) const;
};

#endif
