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
#include "peano/datatraversal/TaskSet.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/VertexOperations.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "exahype/mappings/Prediction.h"

tarch::logging::Log exahype::mappings::FusedTimeStep::_log(
    "exahype::mappings::FusedTimeStep");

void exahype::mappings::FusedTimeStep::updateBatchIterationCounter() {
  if ( exahype::State::isFirstIterationOfBatchOrNoBatch() ) {
    _batchIteration = 0;
  } else {
    _batchIteration++;
  }
  _batchIterationCounterUpdated = true;
}

bool exahype::mappings::FusedTimeStep::issuePredictionJobsInThisIteration() {
  return
      exahype::solvers::Solver::PredictionSweeps==1 ||
      _batchIteration % 2 == 0;
}

bool exahype::mappings::FusedTimeStep::sendOutRiemannDataInThisIteration() {
  return
      exahype::solvers::Solver::PredictionSweeps==1     ||
      exahype::State::isLastIterationOfBatchOrNoBatch() || // covers the NoBatch case
      _batchIteration % 2 != 0;
}


void exahype::mappings::FusedTimeStep::initialiseLocalVariables(){
  const unsigned int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
  _minTimeStepSizes.resize(numberOfSolvers);
  _maxLevels.resize(numberOfSolvers);
  _limiterDomainChanges.resize(numberOfSolvers);
  _meshUpdateRequests.resize(numberOfSolvers);

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    _minTimeStepSizes[solverNumber]     = std::numeric_limits<double>::max();
    _maxLevels[solverNumber]            = -std::numeric_limits<int>::max(); // "-", min
    _limiterDomainChanges[solverNumber] = exahype::solvers::LimiterDomainChange::Regular;
    _meshUpdateRequests[solverNumber]   = false;
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
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
}

/**
 * Nop.
 */
peano::MappingSpecification
exahype::mappings::FusedTimeStep::touchVertexLastTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
}
peano::MappingSpecification
exahype::mappings::FusedTimeStep::ascendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,false);
}
peano::MappingSpecification
exahype::mappings::FusedTimeStep::descendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,false);
}

exahype::mappings::FusedTimeStep::FusedTimeStep() {
}

void exahype::mappings::FusedTimeStep::beginIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);

  if ( exahype::State::isFirstIterationOfBatchOrNoBatch() ) {
    exahype::plotters::startPlottingIfAPlotterIsActive(
        solvers::Solver::getMinTimeStampOfAllSolvers());

    for (auto* solver : exahype::solvers::RegisteredSolvers) {
      solver->setNextMeshUpdateRequest();
      solver->setNextAttainedStableState();

      if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
        static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->setNextLimiterDomainChange();
      }
    }

    initialiseLocalVariables();
  }

  logTraceOutWith1Argument("beginIteration(State)", solverState);
}

void exahype::mappings::FusedTimeStep::endIteration(
    exahype::State& state) {
  logTraceInWith1Argument("endIteration(State)", state);

  if ( sendOutRiemannDataInThisIteration() ) {
    exahype::plotters::finishedPlotting();

    exahype::solvers::Solver::startNewTimeStepForAllSolvers(
        _minTimeStepSizes,_maxLevels,_meshUpdateRequests,_limiterDomainChanges,
        exahype::State::isFirstIterationOfBatchOrNoBatch(),
        exahype::State::isLastIterationOfBatchOrNoBatch(),
        true);
  }

  _batchIterationCounterUpdated = false;

  logTraceOutWith1Argument("endIteration(State)", state);
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::FusedTimeStep::FusedTimeStep(
    const FusedTimeStep& masterThread) {
  _batchIterationCounterUpdated=masterThread._batchIterationCounterUpdated;
  _batchIteration=masterThread._batchIteration;
  initialiseLocalVariables();
}
// Merge over threads
void exahype::mappings::FusedTimeStep::mergeWithWorkerThread(
    const FusedTimeStep& workerThread) {
  for (int i = 0; i < static_cast<int>(exahype::solvers::RegisteredSolvers.size()); i++) {
    _meshUpdateRequests[i]   = _meshUpdateRequests[i] || workerThread._meshUpdateRequests[i];
    _limiterDomainChanges[i] = std::max ( _limiterDomainChanges[i], workerThread._limiterDomainChanges[i] );
    _minTimeStepSizes[i]     = std::min(_minTimeStepSizes[i], workerThread._minTimeStepSizes[i]);
    _maxLevels[i]            = std::max(_maxLevels[i], workerThread._maxLevels[i]);
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
    if ( issuePredictionJobsInThisIteration() ) {
      exahype::Cell::validateThatAllNeighbourMergesHaveBeenPerformed(
          fineGridCell.getCellDescriptionsIndex(),
          fineGridVerticesEnumerator);
    }

    const int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
    for (int solverNumber=0; solverNumber<numberOfSolvers; solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      const int element = solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
      if (element!=exahype::solvers::Solver::NotFound) {
        if ( issuePredictionJobsInThisIteration() ) {
          // this operates only on compute cells
          exahype::plotters::plotPatchIfAPlotterIsActive(
              solverNumber,fineGridCell.getCellDescriptionsIndex(),element); // TODO(Dominic) potential for IO overlap?

          exahype::solvers::Solver::UpdateResult result =
              solver->fusedTimeStep(
                  fineGridCell.getCellDescriptionsIndex(),element,
                  exahype::State::isFirstIterationOfBatchOrNoBatch(),
                  exahype::State::isLastIterationOfBatchOrNoBatch(),
                  exahype::Cell::isAtRemoteBoundary(
                      fineGridVertices,fineGridVerticesEnumerator)
          );

          _meshUpdateRequests    [solverNumber]  =
              _meshUpdateRequests[solverNumber] || result._refinementRequested;
          _limiterDomainChanges  [solverNumber]  = std::max( _limiterDomainChanges[solverNumber], result._limiterDomainChange );
          assertion(_limiterDomainChanges[solverNumber]!=exahype::solvers::LimiterDomainChange::IrregularRequiringMeshUpdate ||
                    _meshUpdateRequests[solverNumber]);
          _minTimeStepSizes[solverNumber] = std::min( result._timeStepSize,                 _minTimeStepSizes[solverNumber]);
          _maxLevels       [solverNumber] = std::min( fineGridVerticesEnumerator.getLevel(),_maxLevels       [solverNumber]);
        }

        if ( sendOutRiemannDataInThisIteration() ) {
          // this operates only on virtual helper cells (pull from below)
          solver->prolongateAndPrepareRestriction(fineGridCell.getCellDescriptionsIndex(),element);
        }
      }
    }

    // Must be performed for all cell descriptions
    exahype::Cell::resetNeighbourMergeFlags(
        fineGridCell.getCellDescriptionsIndex());
    exahype::Cell::resetFaceDataExchangeCounters(
        fineGridCell.getCellDescriptionsIndex(),
        fineGridVertices,fineGridVerticesEnumerator);
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
  logTraceInWith6Arguments("touchVertexFirstTime(...)", fineGridVertex,
                           fineGridX, fineGridH,
                           coarseGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfVertex);

  if ( !_batchIterationCounterUpdated ) {
    updateBatchIterationCounter();
    if ( exahype::solvers::Solver::SpawnPredictionAsBackgroundJob ) {
      if ( issuePredictionJobsInThisIteration() ) {
        exahype::solvers::Solver::ensureAllJobsHaveTerminated(exahype::solvers::Solver::JobType::EnclaveJob);
      }
      if ( sendOutRiemannDataInThisIteration() ) {
        exahype::solvers::Solver::ensureAllJobsHaveTerminated(exahype::solvers::Solver::JobType::SkeletonJob);
        peano::datatraversal::TaskSet::startToProcessBackgroundJobs();
      }
    }
  }

  if ( issuePredictionJobsInThisIteration() ) {
    fineGridVertex.mergeNeighbours(fineGridX,fineGridH);
  }

  logTraceOutWith1Argument("touchVertexFirstTime(...)", fineGridVertex);
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

  if ( sendOutRiemannDataInThisIteration() ) {
    exahype::mappings::Prediction::restriction(
        fineGridCell,exahype::State::AlgorithmSection::TimeStepping);
  }

  logTraceOutWith1Argument("leaveCell(...)", fineGridCell);
}

#ifdef Parallel
void exahype::mappings::FusedTimeStep::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  logTraceInWith6Arguments( "mergeWithNeighbour(...)", vertex, neighbour, fromRank, fineGridX, fineGridH, level );

  if ( !_batchIterationCounterUpdated ) {
    updateBatchIterationCounter();
    if ( exahype::solvers::Solver::SpawnPredictionAsBackgroundJob ) {
      if ( issuePredictionJobsInThisIteration() ) {
        exahype::solvers::Solver::ensureAllJobsHaveTerminated(exahype::solvers::Solver::JobType::EnclaveJob);
      }
      if ( sendOutRiemannDataInThisIteration() ) {
        exahype::solvers::Solver::ensureAllJobsHaveTerminated(exahype::solvers::Solver::JobType::SkeletonJob);
        peano::datatraversal::TaskSet::startToProcessBackgroundJobs();
      }
    }
  }

  if ( issuePredictionJobsInThisIteration() ) {
    vertex.receiveNeighbourData(
        fromRank, true/*merge with data*/,exahype::State::isFirstIterationOfBatchOrNoBatch(),
        fineGridX,fineGridH,level);
  }

  logTraceOut( "mergeWithNeighbour(...)" );
}

void exahype::mappings::FusedTimeStep::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  logTraceInWith5Arguments( "prepareSendToNeighbour(...)", vertex, toRank, x, h, level );

  if ( sendOutRiemannDataInThisIteration() ) {
    vertex.sendToNeighbour(toRank,exahype::State::isLastIterationOfBatchOrNoBatch(),x,h,level);
  }

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
  logTraceIn( "prepareSendToWorker(...)" );

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

  logTraceOutWith1Argument( "prepareSendToWorker(...)", true );

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
  logTraceIn( "receiveDataFromMaster(...)" );

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

  logTraceOut( "receiveDataFromMaster(...)" );
}

void exahype::mappings::FusedTimeStep::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  logTraceInWith2Arguments( "mergeWithWorker(...)", localCell.toString(), receivedMasterCell.toString() );

  if ( exahype::State::isFirstIterationOfBatchOrNoBatch() ) {
    localCell.mergeWithMasterDataPerCell( cellSize );
  }

  logTraceOutWith1Argument( "mergeWithWorker(...)", localCell.toString() );
}

// WORKER->MASTER

void exahype::mappings::FusedTimeStep::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith2Arguments( "prepareSendToMaster(...)", localCell, verticesEnumerator.toString() );

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

  logTraceOut( "prepareSendToMaster(...)" );
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
  logTraceIn( "mergeWithMaster(...)" );

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

  logTraceOut( "mergeWithMaster(...)" );
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


exahype::mappings::FusedTimeStep::~FusedTimeStep() {
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
