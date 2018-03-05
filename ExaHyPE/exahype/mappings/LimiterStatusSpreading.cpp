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

tarch::logging::Log exahype::mappings::LimiterStatusSpreading::_log("exahype::mappings::LimiterStatusSpreading");

void exahype::mappings::LimiterStatusSpreading::initialiseLocalVariables(){
  const unsigned int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
  _limiterDomainChanges.resize(numberOfSolvers);
  _meshUpdateRequests.resize(numberOfSolvers);

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    _limiterDomainChanges[solverNumber] = exahype::solvers::LimiterDomainChange::Regular;
    _meshUpdateRequests[solverNumber]   = false;
  }
}

peano::CommunicationSpecification
exahype::mappings::LimiterStatusSpreading::communicationSpecification() const {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
          MaskOutMasterWorkerDataAndStateExchange,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
          SendDataAndStateAfterLastTouchVertexLastTime,
      true);
}

// Switched on.
peano::MappingSpecification
exahype::mappings::LimiterStatusSpreading::touchVertexFirstTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces,true); // TODO(Dominic): false should work in theory
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
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
}
peano::MappingSpecification
exahype::mappings::LimiterStatusSpreading::leaveCellSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::Serial,false);
}
peano::MappingSpecification
exahype::mappings::LimiterStatusSpreading::ascendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,false);
}
peano::MappingSpecification
exahype::mappings::LimiterStatusSpreading::descendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,false);
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::LimiterStatusSpreading::LimiterStatusSpreading(
    const LimiterStatusSpreading& masterThread) {
  initialiseLocalVariables();
}
#endif

exahype::mappings::LimiterStatusSpreading::~LimiterStatusSpreading() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
void exahype::mappings::LimiterStatusSpreading::mergeWithWorkerThread(
    const LimiterStatusSpreading& workerThread) {
  for (int i = 0; i < static_cast<int>(exahype::solvers::RegisteredSolvers.size()); i++) {
    _meshUpdateRequests[i]  =
        _meshUpdateRequests[i] || workerThread._meshUpdateRequests[i];
    _limiterDomainChanges[i] = std::max ( _limiterDomainChanges[i], workerThread._limiterDomainChanges[i] );
  }
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

  if ( exahype::State::isFirstIterationOfBatchOrNoBatch() ) {
    // We memorise the previous request per solver
    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      if ( spreadLimiterStatus(solver) ) {
        auto* limitingADERDG = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
        limitingADERDG->updateNextMeshUpdateRequest(solver->getMeshUpdateRequest());
        limitingADERDG->updateNextLimiterDomainChange(limitingADERDG->getLimiterDomainChange());
      }
    }

    initialiseLocalVariables();
  }

  #ifdef Parallel
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

      limitingADERDG->updateNextMeshUpdateRequest(_meshUpdateRequests[solverNumber]);
      limitingADERDG->updateNextAttainedStableState(!limitingADERDG->getNextMeshUpdateRequest());
      if (limitingADERDG->getNextMeshUpdateRequest()==true) {
        limitingADERDG->updateNextLimiterDomainChange(
            exahype::solvers::LimiterDomainChange::IrregularRequiringMeshUpdate);
      }

      limitingADERDG->setNextMeshUpdateRequest();
      limitingADERDG->setNextLimiterDomainChange();
      limitingADERDG->setNextAttainedStableState();
    }
  }
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
      exahype::State::AlgorithmSection::LimiterStatusSpreading,fineGridX,fineGridH);

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

      const int cellDescriptionsIndex = fineGridCell.getCellDescriptionsIndex();
      const int element = solver->tryGetElement(cellDescriptionsIndex,solverNumber);
      if (
          spreadLimiterStatus(solver) &&
          element!=exahype::solvers::Solver::NotFound
      ) {
        auto* limitingADERDG = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);

        limitingADERDG->updateLimiterStatusDuringLimiterStatusSpreading(
            cellDescriptionsIndex,element);

        auto& solverPatch = exahype::solvers::ADERDGSolver::getCellDescription(
            cellDescriptionsIndex,element);
        _meshUpdateRequests[solverNumber] =
            _meshUpdateRequests[solverNumber] ||
              limitingADERDG->evaluateLimiterStatusRefinementCriterion(solverPatch);
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
  logTraceIn( "prepareSendToWorker(...)" );
  logTraceOutWith1Argument( "prepareSendToWorker(...)", true );
  return true;
}

void exahype::mappings::LimiterStatusSpreading::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  logTraceInWith6Arguments("mergeWithNeighbour(...)", vertex, neighbour,
                           fromRank, fineGridX, fineGridH, level);

  if ( exahype::State::getBatchState()!=exahype::State::BatchState::FirstIterationOfBatch ) {
    vertex.mergeOnlyWithNeighbourMetadata(
        fromRank,fineGridX,fineGridH,level,
        exahype::State::AlgorithmSection::LimiterStatusSpreading);
  }

  logTraceOut("mergeWithNeighbour(...)");
}

void exahype::mappings::LimiterStatusSpreading::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  logTraceInWith5Arguments("prepareSendToNeighbour(...)", vertex,
                           toRank, x, h, level);

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
  logTraceInWith2Arguments( "prepareSendToMaster(...)", localCell, verticesEnumerator.toString() );

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    if ( spreadLimiterStatus(solver) ) {
      solver->sendMeshUpdateFlagsToMaster(
          tarch::parallel::NodePool::getInstance().getMasterRank(),
          verticesEnumerator.getCellCenter(),
          verticesEnumerator.getLevel());
    }
  }

  logTraceOut( "prepareSendToMaster(...)" );
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
  logTraceIn( "mergeWithMaster(...)" );

  // Merge global solver states
  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    if ( spreadLimiterStatus(solver) ) {
      solver->mergeWithWorkerMeshUpdateFlags(
          worker,
          fineGridVerticesEnumerator.getCellCenter(),
          fineGridVerticesEnumerator.getLevel());
    }
  }

  logTraceOut( "mergeWithMaster(...)" );
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
