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
 
#include "exahype/mappings/FusedTimeStep.h"

#include "tarch/multicore/Loop.h"

#include "peano/utils/Globals.h"

#include "peano/datatraversal/autotuning/Oracle.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/VertexOperations.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "exahype/mappings/Prediction.h"

tarch::logging::Log exahype::mappings::FusedTimeStep::_log(
    "exahype::mappings::FusedTimeStep");

void exahype::mappings::FusedTimeStep::prepareLocalTimeStepVariables(){
  const unsigned int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
  _minTimeStepSizes.resize(numberOfSolvers);
  _maxLevels.resize(numberOfSolvers);

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    _minTimeStepSizes[solverNumber] = std::numeric_limits<double>::max();
    _maxLevels    [solverNumber]    = -std::numeric_limits<int>::max(); // "-", min
  }
}

peano::CommunicationSpecification
exahype::mappings::FusedTimeStep::communicationSpecification() const {
  // master->worker
  peano::CommunicationSpecification::ExchangeMasterWorkerData exchangeMasterWorkerData =
      peano::CommunicationSpecification::ExchangeMasterWorkerData::MaskOutMasterWorkerDataAndStateExchange;
  if ( exahype::State::isFirstIterationOfBatchOrNoBatch() ) {
    exchangeMasterWorkerData =
        peano::CommunicationSpecification::ExchangeMasterWorkerData::SendDataAndStateBeforeFirstTouchVertexFirstTime;
  }

  // worker->master
  peano::CommunicationSpecification::ExchangeWorkerMasterData exchangeWorkerMasterData =
      peano::CommunicationSpecification::ExchangeWorkerMasterData::MaskOutWorkerMasterDataAndStateExchange;
  if ( exahype::State::isLastIterationOfBatchOrNoBatch() ) {
    exchangeWorkerMasterData =
        peano::CommunicationSpecification::ExchangeWorkerMasterData::SendDataAndStateAfterLastTouchVertexLastTime;
  }

  return peano::CommunicationSpecification(exchangeMasterWorkerData,exchangeWorkerMasterData,true);
}

peano::MappingSpecification
exahype::mappings::FusedTimeStep::enterCellSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}
peano::MappingSpecification
exahype::mappings::FusedTimeStep::touchVertexFirstTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}
peano::MappingSpecification
exahype::mappings::FusedTimeStep::leaveCellSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}

/**
 * Nop.
 */
peano::MappingSpecification
exahype::mappings::FusedTimeStep::touchVertexLastTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}
peano::MappingSpecification
exahype::mappings::FusedTimeStep::ascendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}
peano::MappingSpecification
exahype::mappings::FusedTimeStep::descendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}



void exahype::mappings::FusedTimeStep::beginIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);

  if ( exahype::State::isFirstIterationOfBatchOrNoBatch() ) {
    _localState = solverState;

    exahype::plotters::startPlottingIfAPlotterIsActive(
        solvers::Solver::getMinTimeStampOfAllSolvers());

    for (auto* solver : exahype::solvers::RegisteredSolvers) {
      solver->setNextMeshUpdateRequest();
      solver->setNextAttainedStableState();

      if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
        static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->setNextLimiterDomainChange();
      }
    }

    // temporary variables
    prepareLocalTimeStepVariables();

    exahype::solvers::initialiseSolverFlags(_solverFlags);
    exahype::solvers::prepareSolverFlags(_solverFlags);
  }

  // background threads
  exahype::solvers::Solver::waitUntilAllBackgroundJobsHaveTerminated();

  logTraceOutWith1Argument("beginIteration(State)", solverState);
}

void exahype::mappings::FusedTimeStep::endIteration(
    exahype::State& state) {
  logTraceInWith1Argument("endIteration(State)", state);

  exahype::plotters::finishedPlotting();

  exahype::solvers::Solver::startNewTimeStepForAllSolvers(
      _solverFlags,_minTimeStepSizes,_maxLevels,
      exahype::State::isFirstIterationOfBatchOrNoBatch(),
      exahype::State::isLastIterationOfBatchOrNoBatch(),
      true);

  // delete temporary variables
  if ( exahype::State::isLastIterationOfBatchOrNoBatch() ) {
    deleteSolverFlags(_solverFlags);
  }

  logTraceOutWith1Argument("endIteration(State)", state);
}

exahype::mappings::FusedTimeStep::~FusedTimeStep() {
  exahype::solvers::deleteSolverFlags(_solverFlags);
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::FusedTimeStep::FusedTimeStep(
    const FusedTimeStep& masterThread)
  : _localState(masterThread._localState) {
  exahype::solvers::initialiseSolverFlags(_solverFlags);
  exahype::solvers::prepareSolverFlags(_solverFlags);

  prepareLocalTimeStepVariables();
}
// Merge over threads
void exahype::mappings::FusedTimeStep::mergeWithWorkerThread(
    const FusedTimeStep& workerThread) {
  for (int i = 0; i < static_cast<int>(exahype::solvers::RegisteredSolvers.size()); i++) {
    // solver flags
    _solverFlags._meshUpdateRequest[i]  |= workerThread._solverFlags._meshUpdateRequest[i];
    _solverFlags._limiterDomainChange[i] =
        std::max ( _solverFlags._limiterDomainChange[i],
            workerThread._solverFlags._limiterDomainChange[i] );
    // time step size and cell sizes
    _minTimeStepSizes[i] =
        std::min(_minTimeStepSizes[i], workerThread._minTimeStepSizes[i]);
    _maxLevels[i] =
        std::max(_maxLevels[i], workerThread._maxLevels[i]);
  }
}
#endif

void exahype::mappings::FusedTimeStep::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,fineGridVerticesEnumerator.toString(),coarseGridCell, fineGridPositionOfCell);

  if ( fineGridCell.isInitialised() ) {
    exahype::Cell::validateThatAllNeighbourMergesHaveBeenPerformed(
        fineGridCell.getCellDescriptionsIndex(),
        fineGridVerticesEnumerator);

    exahype::Cell::resetNeighbourMergeFlags(
        fineGridCell.getCellDescriptionsIndex());
    exahype::Cell::resetFaceDataExchangeCounters(
        fineGridCell.getCellDescriptionsIndex(),
        fineGridVertices,fineGridVerticesEnumerator);

    const int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
    for (int solverNumber=0; solverNumber<numberOfSolvers; solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      const int element = solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
      if (element!=exahype::solvers::Solver::NotFound) {
        exahype::plotters::plotPatchIfAPlotterIsActive(
            solverNumber,fineGridCell.getCellDescriptionsIndex(),element);

        // this operates only on compute cells
        exahype::solvers::Solver::UpdateResult result =
            solver->fusedTimeStep(
                fineGridCell.getCellDescriptionsIndex(),element,
                exahype::State::isFirstIterationOfBatchOrNoBatch(),
                exahype::State::isLastIterationOfBatchOrNoBatch(),
                exahype::Cell::isAtRemoteBoundary(
                    fineGridVertices,fineGridVerticesEnumerator)
            );

        // this operates only on helper cells
        solver->prolongateAndPrepareRestriction(fineGridCell.getCellDescriptionsIndex(),element);

        _solverFlags._meshUpdateRequest  [solverNumber] |= result._refinementRequested;
        _solverFlags._limiterDomainChange[solverNumber]  = std::max( _solverFlags._limiterDomainChange[solverNumber], result._limiterDomainChange );
        assertion(_solverFlags._limiterDomainChange[solverNumber]!=exahype::solvers::LimiterDomainChange::IrregularRequiringMeshUpdate ||
            _solverFlags._meshUpdateRequest[solverNumber]);
        _minTimeStepSizes[solverNumber] = std::min( result._timeStepSize,                 _minTimeStepSizes[solverNumber]);
        _maxLevels       [solverNumber] = std::min( fineGridVerticesEnumerator.getLevel(),_maxLevels       [solverNumber]);
      }
    }
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}


void exahype::mappings::FusedTimeStep::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  fineGridVertex.mergeNeighbours(fineGridX,fineGridH);
}

void exahype::mappings::FusedTimeStep::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("leaveCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  exahype::mappings::Prediction::restriction(
      fineGridCell,exahype::State::AlgorithmSection::TimeStepping);

  logTraceOutWith1Argument("leaveCell(...)", fineGridCell);
}

#ifdef Parallel
void exahype::mappings::FusedTimeStep::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  logTraceInWith6Arguments( "mergeWithNeighbour(...)", vertex, neighbour, fromRank, fineGridX, fineGridH, level );

  vertex.receiveNeighbourData(
        fromRank,true,
        fineGridX,fineGridH,level);

  logTraceOut( "mergeWithNeighbour(...)" );
}

void exahype::mappings::FusedTimeStep::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  logTraceInWith5Arguments( "prepareSendToNeighbour(...)", vertex, toRank, x, h, level );

  vertex.sendToNeighbour(toRank,x,h,level);

  logTraceOut( "prepareSendToNeighbour(...)" );
}

// MASTER->WORKER

bool exahype::mappings::FusedTimeStep::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {

  if ( exahype::State::isFirstIterationOfBatchOrNoBatch() ) {
    exahype::Cell::broadcastGlobalDataToWorker(
        worker,
        fineGridVerticesEnumerator.getCellCenter(),
        fineGridVerticesEnumerator.getLevel());

    fineGridCell.broadcastDataToWorkerPerCell(
        worker,
        fineGridVerticesEnumerator.getCellCenter(),
        fineGridVerticesEnumerator.getCellSize(),
        fineGridVerticesEnumerator.getLevel());
  }

  return exahype::State::isFirstIterationOfBatchOrNoBatch() ||
         exahype::State::isLastIterationOfBatchOrNoBatch();
}

void exahype::mappings::FusedTimeStep::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {

  if ( exahype::State::isFirstIterationOfBatchOrNoBatch() ) {
    exahype::Cell::mergeWithGlobalDataFromMaster(
        tarch::parallel::NodePool::getInstance().getMasterRank(),
        receivedVerticesEnumerator.getCellCenter(),
        receivedVerticesEnumerator.getLevel());

    receivedCell.receiveDataFromMasterPerCell(
        tarch::parallel::NodePool::getInstance().getMasterRank(),
        receivedVerticesEnumerator.getCellCenter(),
        receivedVerticesEnumerator.getCellSize(),
        receivedVerticesEnumerator.getLevel());
  }

}

void exahype::mappings::FusedTimeStep::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  if ( exahype::State::isFirstIterationOfBatchOrNoBatch() ) {
    localCell.mergeWithMasterDataPerCell( cellSize );
  }
}

// WORKER->MASTER

void exahype::mappings::FusedTimeStep::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  if ( exahype::State::isLastIterationOfBatchOrNoBatch() ) {
    exahype::Cell::reduceGlobalDataToMaster(
        tarch::parallel::NodePool::getInstance().getMasterRank(),
        verticesEnumerator.getCellCenter(),
        verticesEnumerator.getLevel());

    localCell.reduceDataToMasterPerCell(
        tarch::parallel::NodePool::getInstance().getMasterRank(),
        verticesEnumerator.getCellCenter(),
        verticesEnumerator.getCellSize(),
        verticesEnumerator.getLevel());
  }
}

void exahype::mappings::FusedTimeStep::mergeWithMaster(
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
  if ( exahype::State::isLastIterationOfBatchOrNoBatch() ) {
    exahype::Cell::mergeWithGlobalDataFromWorker(
        worker,
        fineGridVerticesEnumerator.getCellCenter(),
        fineGridVerticesEnumerator.getLevel());

    fineGridCell.mergeWithDataFromWorkerPerCell(
        worker,
        fineGridVerticesEnumerator.getCellCenter(),
        fineGridVerticesEnumerator.getCellSize(),
        fineGridVerticesEnumerator.getLevel());
  }
}


//
// Below all methods are nop.
//
//=====================================



void exahype::mappings::FusedTimeStep::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::FusedTimeStep::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::FusedTimeStep::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::FusedTimeStep::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::FusedTimeStep::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

exahype::mappings::FusedTimeStep::FusedTimeStep() {
  // do nothing
}

void exahype::mappings::FusedTimeStep::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FusedTimeStep::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FusedTimeStep::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FusedTimeStep::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FusedTimeStep::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FusedTimeStep::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::FusedTimeStep::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::FusedTimeStep::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::FusedTimeStep::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::FusedTimeStep::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
