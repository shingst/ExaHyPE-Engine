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

#ifndef EXAHYPE_MAPPINGS_GlobalRollback_H_
#define EXAHYPE_MAPPINGS_GlobalRollback_H_

#include "tarch/la/Vector.h"
#include "tarch/logging/Log.h"

#include "peano/CommunicationSpecification.h"
#include "peano/MappingSpecification.h"
#include "peano/grid/VertexEnumerator.h"

#include "tarch/multicore/MulticoreDefinitions.h"

#include "exahype/Cell.h"
#include "exahype/State.h"
#include "exahype/Vertex.h"

namespace exahype {
namespace mappings {
class GlobalRollback;
}
}

/**
 * TODO(Dominic): This mapping is unused
 *
 * @author Dominic Charrier
 */
class exahype::mappings::GlobalRollback {
 private:
  /**
   * Logging device for the trace macros.
   */
  static tarch::logging::Log _log;

  /**
   * \return true if we are currently performing
   * a global rollback for the solver.
   */
  static bool performGlobalRollback(exahype::solvers::Solver* solver);

 public:
  /**
   * Mask out data exchange between master and worker.
   * Further let Peano handle heap data exchange internally.
   */
  peano::CommunicationSpecification communicationSpecification() const;

  /**
   * Run through the whole grid. Run concurrently on the fine grid.
   */
  peano::MappingSpecification enterCellSpecification(int level) const;

  /////
  // Below all Specifications are nop
  /////
  /**
   * Nop.
   */
  peano::MappingSpecification touchVertexFirstTimeSpecification(int level) const;
  /**
   * Nop.
   */
  peano::MappingSpecification touchVertexLastTimeSpecification(int level) const;
  /**
   * Nop.
   */
  peano::MappingSpecification leaveCellSpecification(int level) const;
  /**
   * Nop.
   */
  peano::MappingSpecification ascendSpecification(int level) const;
  /**
   * Nop.
   */
  peano::MappingSpecification descendSpecification(int level) const;

  /**
   * Roll the global solver time step data back to the previous
   * state.
   */
  void endIteration(exahype::State& solverState);

  /**
   * For all cells hosting a solver patch of a
   * LimitingADERDGSolver, perform the following operations:
   *
   * 1. Perform a rollback to the previous time step.
   * 2. Reinitialise the solver patch, i.e. perform a rollback in
   *    cells with LimiterStatus other than Ok, and
   *    allocate an additional limiter patch if necessary.
   *    If the solver requested a global recomputation, i.e.
   *    additional mesh refinement,
   *    also perform a rollback in cells with LimiterStatus
   *    Ok.
   * 3. Reset the neighbour merging flags.
   */
  void enterCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  #ifdef Parallel
    /**
     * Drop neighbour data from the previous limiter status spreading.
     */
    void mergeWithNeighbour(exahype::Vertex& vertex,
                            const exahype::Vertex& neighbour, int fromRank,
                            const tarch::la::Vector<DIMENSIONS, double>& x,
                            const tarch::la::Vector<DIMENSIONS, double>& h,
                            int level);

    /**
     * Send the limiter domain change state down to the worker.
     *
     * \note Currently, we perform a full broadcast for all solver, we send up to
     * nine doubles to the worker. This is not necessary, we only
     * need the limiter domain change state.
     */
    bool prepareSendToWorker(
        exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
        int worker);

    /**
     * Receive limiter domain change state from the master.
     *
     * \note Currently, we perform a full broadcast for all solvers, we send up to
     * nine doubles to the worker. This is not necessary, we only
     * need the limiter domain change state.
     */
    void receiveDataFromMaster(
        exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
        const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
        exahype::Vertex* const receivedCoarseGridVertices,
        const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
        exahype::Cell& receivedCoarseGridCell,
        exahype::Vertex* const workersCoarseGridVertices,
        const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
        exahype::Cell& workersCoarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);


    //
    // Below all methods are nop.
    //
    //===================================

    /**
     * Nop.
     */
    void prepareSendToNeighbour(exahype::Vertex& vertex, int toRank,
                                const tarch::la::Vector<DIMENSIONS, double>& x,
                                const tarch::la::Vector<DIMENSIONS, double>& h,
                                int level);


    /**
     * Nop.
     */
    void prepareCopyToRemoteNode(exahype::Vertex& localVertex, int toRank,
                                 const tarch::la::Vector<DIMENSIONS, double>& x,
                                 const tarch::la::Vector<DIMENSIONS, double>& h,
                                 int level);

    /**
     * Nop.
     */
    void prepareCopyToRemoteNode(
        exahype::Cell& localCell, int toRank,
        const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level);

    /**
     * Nop.
     */
    void mergeWithRemoteDataDueToForkOrJoin(
        exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, double>& h, int level);

    /**
     * Nop.
     */
    void mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level);

    /**
     * Nop.
     */
    void mergeWithMaster(
        const exahype::Cell& workerGridCell,
        exahype::Vertex* const workerGridVertices,
        const peano::grid::VertexEnumerator& workerEnumerator,
        exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
        int worker, const exahype::State& workerState,
        exahype::State& masterState);

    /**
     * Nop.
     */
    void prepareSendToMaster(
        exahype::Cell& localCell, exahype::Vertex* vertices,
        const peano::grid::VertexEnumerator& verticesEnumerator,
        const exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        const exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

    /**
     * Nop.
     */
    void mergeWithWorker(exahype::Cell& localCell,
                         const exahype::Cell& receivedMasterCell,
                         const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
                         const tarch::la::Vector<DIMENSIONS, double>& cellSize,
                         int level);

    /**
     * Nop.
     */
    void mergeWithWorker(exahype::Vertex& localVertex,
                         const exahype::Vertex& receivedMasterVertex,
                         const tarch::la::Vector<DIMENSIONS, double>& x,
                         const tarch::la::Vector<DIMENSIONS, double>& h,
                         int level);
  #endif

    /**
     * Nop
     */
    GlobalRollback();

  #if defined(SharedMemoryParallelisation)
    /**
     * Nop.
     */
    GlobalRollback(const GlobalRollback& masterThread);
  #endif

    /**
     * Nop.
     */
    virtual ~GlobalRollback();

  #if defined(SharedMemoryParallelisation)
    /**
     * Nop.
     */
    void mergeWithWorkerThread(const GlobalRollback& workerThread);
  #endif

    /**
     * Nop.
     */
    void beginIteration(exahype::State& solverState);

    /**
     * Nop.
     */
    void touchVertexFirstTime(
        exahype::Vertex& fineGridVertex,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);

    /**
     * Nop.
     */
    void createInnerVertex(
        exahype::Vertex& fineGridVertex,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);

    /**
     * Nop.
     */
    void createBoundaryVertex(
        exahype::Vertex& fineGridVertex,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);

    /**
     * Nop.
     */
    void createHangingVertex(
        exahype::Vertex& fineGridVertex,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);

    /**
     * Nop.
     */
    void destroyHangingVertex(
        const exahype::Vertex& fineGridVertex,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);

    /**
     * Nop.
     */
    void destroyVertex(
        const exahype::Vertex& fineGridVertex,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);

    /**
     * Nop.
     */
    void createCell(
        exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

    /**
     * Nop.
     */
    void destroyCell(
        const exahype::Cell& fineGridCell,
        exahype::Vertex* const fineGridVertices,
        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

    /**
     * Nop.
     */
    void touchVertexLastTime(
        exahype::Vertex& fineGridVertex,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);

    /**
     * Nop.
     */
    void leaveCell(
        exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

    /**
     * Nop.
     */
    void descend(
        exahype::Cell* const fineGridCells,
        exahype::Vertex* const fineGridVertices,
        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell);

    /**
     * Nop.
     */
    void ascend(exahype::Cell* const fineGridCells,
                exahype::Vertex* const fineGridVertices,
                const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
                exahype::Vertex* const coarseGridVertices,
                const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
                exahype::Cell& coarseGridCell);
  };

  #endif
