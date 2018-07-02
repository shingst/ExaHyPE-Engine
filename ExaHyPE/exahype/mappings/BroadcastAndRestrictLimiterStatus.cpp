/**
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
 
#include "exahype/mappings/BroadcastAndRestrictLimiterStatus.h"

#include "tarch/multicore/Loop.h"

#include "peano/utils/Loop.h"
#include "peano/datatraversal/autotuning/Oracle.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/VertexOperations.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "exahype/mappings/LimiterStatusSpreading.h"

bool exahype::mappings::BroadcastAndRestrictLimiterStatus::oneSolverIsLimitingADERDGSolverUsingAMR() {
  static int solversFound = -1;
  if ( solversFound == -1 ) {
    solversFound = 0;
    for (unsigned int solverNumber = 0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      solversFound +=
          solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG &&
          solver->getMaximumAdaptiveMeshDepth()>0;
    }
  }
  return solversFound>0;
}

peano::CommunicationSpecification
exahype::mappings::BroadcastAndRestrictLimiterStatus::communicationSpecification() const {
  peano::CommunicationSpecification::ExchangeWorkerMasterData exchangeWorkerMasterData =
      oneSolverIsLimitingADERDGSolverUsingAMR() ?
          peano::CommunicationSpecification::ExchangeWorkerMasterData::SendDataAndStateAfterLastTouchVertexLastTime :
          peano::CommunicationSpecification::ExchangeWorkerMasterData::MaskOutWorkerMasterDataAndStateExchange;

  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::SendDataAndStateBeforeFirstTouchVertexFirstTime,
      exchangeWorkerMasterData,
      true);
}

peano::MappingSpecification
exahype::mappings::BroadcastAndRestrictLimiterStatus::enterCellSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
}

/* All specifications below are nop. */
peano::MappingSpecification
exahype::mappings::BroadcastAndRestrictLimiterStatus::touchVertexFirstTimeSpecification(int level) const {
  return peano::MappingSpecification(
        peano::MappingSpecification::Nop,
        peano::MappingSpecification::AvoidFineGridRaces,true);
}
peano::MappingSpecification
exahype::mappings::BroadcastAndRestrictLimiterStatus::touchVertexLastTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
}
peano::MappingSpecification
exahype::mappings::BroadcastAndRestrictLimiterStatus::leaveCellSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}
peano::MappingSpecification
exahype::mappings::BroadcastAndRestrictLimiterStatus::ascendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,false);
}
peano::MappingSpecification
exahype::mappings::BroadcastAndRestrictLimiterStatus::descendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,false);
}

tarch::logging::Log exahype::mappings::BroadcastAndRestrictLimiterStatus::_log(
    "exahype::mappings::BroadcastAndRestrictLimiterStatus");

void exahype::mappings::BroadcastAndRestrictLimiterStatus::beginIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);

  if ( exahype::solvers::Solver::SpawnPredictionAsBackgroundJob ) {
    // background threads
    exahype::solvers::Solver::ensureAllJobsHaveTerminated(exahype::solvers::Solver::JobType::SkeletonJob);
    exahype::solvers::Solver::ensureAllJobsHaveTerminated(exahype::solvers::Solver::JobType::EnclaveJob);
  }

  exahype::mappings::LimiterStatusSpreading::IsFirstIteration = true;

  logTraceOutWith1Argument("beginIteration(State)", solverState);
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::restrictLimiterStatus(
    const int fineGridCellDescriptionsIndex,
    const int coarseGridCellDescriptionsIndex
) {
  assertion2(exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(fineGridCellDescriptionsIndex), fineGridCellDescriptionsIndex, coarseGridCellDescriptionsIndex );
  for ( auto& cellDescription : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(fineGridCellDescriptionsIndex) ) {
    cellDescription.setParentIndex(coarseGridCellDescriptionsIndex);
    if (
        exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()]->getType()==exahype::solvers::Solver::Type::LimitingADERDG &&
        cellDescription.getParentIndex() >= 0
    ) {
      auto* limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(
          exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()]);
      const int parentElement =
          limitingADERDGSolver->tryGetElement(cellDescription.getParentIndex(),cellDescription.getSolverNumber());
      if ( parentElement!=exahype::solvers::Solver::NotFound ) {
        limitingADERDGSolver->getSolver().get()->restrictToNextParent(cellDescription,parentElement);
      }
    }
  }
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  if ( fineGridCell.isInitialised() ) {
    if ( oneSolverIsLimitingADERDGSolverUsingAMR() ) {
      restrictLimiterStatus(
          fineGridCell.getCellDescriptionsIndex(),
          coarseGridCell.getCellDescriptionsIndex());
    }

    exahype::Cell::resetNeighbourMergeFlags(
        fineGridCell.getCellDescriptionsIndex(),
        fineGridVertices,fineGridVerticesEnumerator);
  }
}

#ifdef Parallel
///////////////////////////////////////
// NEIGHBOUR
///////////////////////////////////////
void exahype::mappings::BroadcastAndRestrictLimiterStatus::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  if ( exahype::solvers::Solver::FuseADERDGPhases ) {
  vertex.receiveNeighbourData(
      fromRank,false /*no merge*/,true /*no batch*/,
      fineGridX,fineGridH,level);
  }
}

///////////////////////////////////////
// MASTER->WORKER
///////////////////////////////////////
bool exahype::mappings::BroadcastAndRestrictLimiterStatus::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  exahype::Cell::broadcastGlobalDataToWorker(
        worker,
        fineGridVerticesEnumerator.getCellCenter(),
        fineGridVerticesEnumerator.getLevel());

  return true; // see docu
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  exahype::Cell::mergeWithGlobalDataFromMaster(
      tarch::parallel::NodePool::getInstance().getMasterRank(),
      receivedVerticesEnumerator.getCellCenter(),
      receivedVerticesEnumerator.getLevel());
}

//
// Below all methods are nop.
//
// ====================================

void exahype::mappings::BroadcastAndRestrictLimiterStatus::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  if (
      oneSolverIsLimitingADERDGSolverUsingAMR() &&
      localCell.hasToCommunicate( verticesEnumerator.getCellSize())
  ) {
    exahype::solvers::ADERDGSolver::sendCellDescriptions(
        tarch::parallel::NodePool::getInstance().getMasterRank(),
        localCell.getCellDescriptionsIndex(),true/* send data from worker side*/,
        peano::heap::MessageType::MasterWorkerCommunication,
        verticesEnumerator.getCellCenter(),
        verticesEnumerator.getLevel());
  }
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::mergeWithMaster(
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
    exahype::State& masterState) {
  if (
      oneSolverIsLimitingADERDGSolverUsingAMR() &&
      fineGridCell.hasToCommunicate( fineGridVerticesEnumerator.getCellSize())
  ) {
    if ( fineGridCell.isInitialised() ) {
      exahype::solvers::ADERDGSolver::eraseCellDescriptions(fineGridCell.getCellDescriptionsIndex());
    } else {
      fineGridCell.setupMetaData();
    }

    exahype::solvers::ADERDGSolver::receiveCellDescriptions(
        worker,fineGridCell,
        peano::heap::MessageType::MasterWorkerCommunication,
        fineGridVerticesEnumerator.getCellCenter(),
        fineGridVerticesEnumerator.getLevel());

    restrictLimiterStatus( fineGridCell.getCellDescriptionsIndex(), coarseGridCell.getCellDescriptionsIndex() );

    if ( fineGridCell.isEmpty() ) {
      fineGridCell.shutdownMetaDataAndResetCellDescriptionsIndex();
    }
  }
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

#endif


exahype::mappings::BroadcastAndRestrictLimiterStatus::~BroadcastAndRestrictLimiterStatus() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::BroadcastAndRestrictLimiterStatus::BroadcastAndRestrictLimiterStatus(const BroadcastAndRestrictLimiterStatus& masterThread) {
  // do nothing
}
#endif

void exahype::mappings::BroadcastAndRestrictLimiterStatus::endIteration(
    exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
void exahype::mappings::BroadcastAndRestrictLimiterStatus::mergeWithWorkerThread(
    const BroadcastAndRestrictLimiterStatus& workerThread) {
  // do nothing
}
#endif

exahype::mappings::BroadcastAndRestrictLimiterStatus::BroadcastAndRestrictLimiterStatus() {
  // do nothing
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::BroadcastAndRestrictLimiterStatus::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
