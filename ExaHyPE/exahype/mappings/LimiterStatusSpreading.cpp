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
 
#include "exahype/mappings/LimiterStatusSpreading.h"

#include "peano/utils/Globals.h"
#include "peano/utils/Loop.h"

#include "peano/datatraversal/autotuning/Oracle.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "tarch/la/VectorScalarOperations.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/VertexOperations.h"

#ifdef Parallel
bool exahype::mappings::LimiterStatusSpreading::IsFirstIteration = true;
#endif

tarch::logging::Log exahype::mappings::LimiterStatusSpreading::_log("exahype::mappings::LimiterStatusSpreading");

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::LimiterStatusSpreading::communicationSpecification() const {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
          SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
          SendDataAndStateAfterLastTouchVertexLastTime,
      true);
}

// Switched on.
peano::MappingSpecification
exahype::mappings::LimiterStatusSpreading::touchVertexFirstTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}
peano::MappingSpecification
exahype::mappings::LimiterStatusSpreading::enterCellSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}

// Switched off.
peano::MappingSpecification
exahype::mappings::LimiterStatusSpreading::touchVertexLastTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}
peano::MappingSpecification
exahype::mappings::LimiterStatusSpreading::leaveCellSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::Serial,true);
}
peano::MappingSpecification
exahype::mappings::LimiterStatusSpreading::ascendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}
peano::MappingSpecification
exahype::mappings::LimiterStatusSpreading::descendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::LimiterStatusSpreading::LimiterStatusSpreading(
    const LimiterStatusSpreading& masterThread)
  : _localState(masterThread._localState) {
  exahype::solvers::initialiseSolverFlags(_solverFlags);
  exahype::solvers::prepareSolverFlags(_solverFlags);
}
#endif

bool exahype::mappings::LimiterStatusSpreading::spreadLimiterStatus(exahype::solvers::Solver* solver) {
  return
      solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG
      &&
      static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainChange()
      !=exahype::solvers::LimiterDomainChange::Regular;
}

void exahype::mappings::LimiterStatusSpreading::beginIteration(
  exahype::State& solverState
) {
  _localState = solverState;

  exahype::solvers::initialiseSolverFlags(_solverFlags);
  exahype::solvers::prepareSolverFlags(_solverFlags);

  // We memorise the previous request per solver
  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
    if ( spreadLimiterStatus(solver) ) {
      auto* limitingADERDG = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
      limitingADERDG->updateNextMeshUpdateRequest(solver->getMeshUpdateRequest());
      limitingADERDG->updateNextLimiterDomainChange(limitingADERDG->getLimiterDomainChange());
    }
  }

  #ifdef Parallel
  peano::heap::AbstractHeap::allHeapsStartToSendSynchronousData();

  if (! MetadataHeap::getInstance().validateThatIncomingJoinBuffersAreEmpty() ) {
      exit(-1);
  }
  #endif
}

void exahype::mappings::LimiterStatusSpreading::endIteration(exahype::State& solverState) {
  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
    if ( spreadLimiterStatus(solver) ) {
      auto* limitingADERDG = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);

      limitingADERDG->updateNextMeshUpdateRequest(_solverFlags._meshUpdateRequest[solverNumber]);
      limitingADERDG->updateNextAttainedStableState(!limitingADERDG->getNextMeshUpdateRequest());
      limitingADERDG->updateNextLimiterDomainChange(_solverFlags._limiterDomainChange[solverNumber]);

      limitingADERDG->setNextMeshUpdateRequest();
      limitingADERDG->setNextLimiterDomainChange();
      limitingADERDG->setNextAttainedStableState();
    }
  }

  deleteSolverFlags(_solverFlags);

  #ifdef Parallel
  exahype::mappings::LimiterStatusSpreading::IsFirstIteration = false;

  peano::heap::AbstractHeap::allHeapsFinishedToSendSynchronousData();
  #endif
}

void exahype::mappings::LimiterStatusSpreading::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::LimiterStatusSpreading::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  logTraceInWith6Arguments("touchVertexFirstTime(...)", fineGridVertex,
                           fineGridX, fineGridH,
                           coarseGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfVertex);

  fineGridVertex.mergeOnlyNeighboursMetadata(
      exahype::records::State::AlgorithmSection::LimiterStatusSpreading,fineGridX,fineGridH);

  logTraceOutWith1Argument("touchVertexFirstTime(...)", fineGridVertex);
}

void exahype::mappings::LimiterStatusSpreading::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  if (fineGridCell.isInitialised()) {
    for (unsigned int solverNumber=0; solverNumber<exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

      const int element =
          solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
      if (
          spreadLimiterStatus(solver)
          &&
          element!=exahype::solvers::Solver::NotFound
      ) {
        auto* limitingADERDG = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
        limitingADERDG->updateLimiterStatusDuringLimiterStatusSpreading(
            fineGridCell.getCellDescriptionsIndex(),element);

        bool meshUpdateRequest =
            limitingADERDG->
            evaluateLimiterStatusRefinementCriterion(
                fineGridCell.getCellDescriptionsIndex(),element);

        _solverFlags._meshUpdateRequest[solverNumber] |= meshUpdateRequest;
        if (meshUpdateRequest) {
          _solverFlags._limiterDomainChange[solverNumber] =
              exahype::solvers::LimiterDomainChange::IrregularRequiringMeshUpdate;
        }
      }
    }

    exahype::Cell::resetNeighbourMergeFlags(
        fineGridCell.getCellDescriptionsIndex());
    exahype::Cell::resetFaceDataExchangeCounters(
        fineGridCell.getCellDescriptionsIndex(),
        fineGridVertices,fineGridVerticesEnumerator);
  }

  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

#ifdef Parallel
bool exahype::mappings::LimiterStatusSpreading::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  return true;
}

void exahype::mappings::LimiterStatusSpreading::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  logTraceInWith6Arguments("mergeWithNeighbour(...)", vertex, neighbour,
                           fromRank, fineGridX, fineGridH, level);
  if (exahype::mappings::LimiterStatusSpreading::IsFirstIteration) {
    return;
  }

  vertex.mergeOnlyWithNeighbourMetadata(
      fromRank,fineGridX,fineGridH,level,
      _localState.getAlgorithmSection());

  logTraceOut("mergeWithNeighbour(...)");
}

void exahype::mappings::LimiterStatusSpreading::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  logTraceInWith5Arguments("prepareSendToNeighbour(...)", vertex,
                           toRank, x, h, level);

  // the limiter status is part of the metadata
  vertex.sendOnlyMetadataToNeighbour(toRank,x,h,level);

  logTraceOut("prepareSendToNeighbour(...)");
}

void exahype::mappings::LimiterStatusSpreading::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  peano::heap::AbstractHeap::allHeapsStartToSendSynchronousData();

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    solver->sendMeshUpdateFlagsToMaster(
        tarch::parallel::NodePool::getInstance().getMasterRank(),
        verticesEnumerator.getCellCenter(),
        verticesEnumerator.getLevel());
  }

  exahype::sendMasterWorkerCommunicationMetadata(
      tarch::parallel::NodePool::getInstance().getMasterRank(),
      localCell.getCellDescriptionsIndex(),
      verticesEnumerator.getCellCenter(),
      verticesEnumerator.getLevel());

  peano::heap::AbstractHeap::allHeapsFinishedToSendSynchronousData();
}

void exahype::mappings::LimiterStatusSpreading::mergeWithMaster(
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
  // Merge global solver states
  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    solver->mergeWithWorkerMeshUpdateFlags(
        worker,
        fineGridVerticesEnumerator.getCellCenter(),
        fineGridVerticesEnumerator.getLevel());
  }

  // Merge cell states
  if (fineGridCell.isInitialised()) {
    fineGridCell.mergeWithWorkerMetadata(
        worker,
        fineGridVerticesEnumerator.getCellCenter(),
        fineGridVerticesEnumerator.getLevel(),
        _localState.getAlgorithmSection());
  } else {
    exahype::dropMetadata(
        worker,
        peano::heap::MessageType::MasterWorkerCommunication,
        fineGridVerticesEnumerator.getCellCenter(),
        fineGridVerticesEnumerator.getLevel());
  }
}

//
// All methods below are nop,
//
// ==================================

void exahype::mappings::LimiterStatusSpreading::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::LimiterStatusSpreading::mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::LimiterStatusSpreading::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::LimiterStatusSpreading::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::LimiterStatusSpreading::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::LimiterStatusSpreading::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::LimiterStatusSpreading::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

exahype::mappings::LimiterStatusSpreading::LimiterStatusSpreading() {
  // do nothing
}

exahype::mappings::LimiterStatusSpreading::~LimiterStatusSpreading() {
  exahype::solvers::deleteSolverFlags(_solverFlags);
}

#if defined(SharedMemoryParallelisation)
void exahype::mappings::LimiterStatusSpreading::mergeWithWorkerThread(
    const LimiterStatusSpreading& workerThread) {
  for (int i = 0; i < static_cast<int>(exahype::solvers::RegisteredSolvers.size()); i++) {
    _solverFlags._meshUpdateRequest[i]  |= workerThread._solverFlags._meshUpdateRequest[i];
    _solverFlags._limiterDomainChange[i] =
        std::max ( _solverFlags._limiterDomainChange[i],
                   workerThread._solverFlags._limiterDomainChange[i] );
  }
}
#endif



void exahype::mappings::LimiterStatusSpreading::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}


void exahype::mappings::LimiterStatusSpreading::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}


void exahype::mappings::LimiterStatusSpreading::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}


void exahype::mappings::LimiterStatusSpreading::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
 // do nothing
}

void exahype::mappings::LimiterStatusSpreading::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::LimiterStatusSpreading::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::LimiterStatusSpreading::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::LimiterStatusSpreading::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}


void exahype::mappings::LimiterStatusSpreading::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::LimiterStatusSpreading::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
