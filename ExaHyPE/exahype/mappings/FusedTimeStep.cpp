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

#ifdef DistributedStealing
#include "exahype/stealing/PerformanceMonitor.h"
#include "exahype/stealing/StaticDistributor.h"
#include "exahype/stealing/DiffusiveDistributor.h"
#endif

tarch::logging::Log exahype::mappings::FusedTimeStep::_log(
    "exahype::mappings::FusedTimeStep");

tarch::multicore::BooleanSemaphore exahype::mappings::FusedTimeStep::Semaphore;

bool exahype::mappings::FusedTimeStep::issuePredictionJobsInThisIteration() const {
  return
      exahype::solvers::Solver::PredictionSweeps==1 ||
      _batchIteration % 2 == 0;
}

bool exahype::mappings::FusedTimeStep::sendOutRiemannDataInThisIteration() const {
  return
      exahype::solvers::Solver::PredictionSweeps==1     ||
      _stateCopy.isLastIterationOfBatchOrNoBatch() || // covers the NoBatch case
      _batchIteration % 2 != 0;
}

peano::CommunicationSpecification
exahype::mappings::FusedTimeStep::communicationSpecification() const {
  // master->worker
  peano::CommunicationSpecification::ExchangeMasterWorkerData exchangeMasterWorkerData =
      peano::CommunicationSpecification::ExchangeMasterWorkerData::MaskOutMasterWorkerDataAndStateExchange;
  #ifdef Parallel
  if (
      exahype::solvers::Solver::PredictionSweeps==1 ||
      exahype::State::BroadcastInThisIteration      // must be set in previous iteration
  ) { // must be set in previous iteration
    exchangeMasterWorkerData =
        peano::CommunicationSpecification::ExchangeMasterWorkerData::SendDataAndStateBeforeFirstTouchVertexFirstTime;
  }
  #endif

  // worker->master
  peano::CommunicationSpecification::ExchangeWorkerMasterData exchangeWorkerMasterData =
      peano::CommunicationSpecification::ExchangeWorkerMasterData::MaskOutWorkerMasterDataAndStateExchange;
  #ifdef Parallel
  if (
      exahype::solvers::Solver::PredictionSweeps==1 ||
      exahype::State::ReduceInThisIteration         // must be set in previous iteration
  ) {
    exchangeWorkerMasterData =
        peano::CommunicationSpecification::ExchangeWorkerMasterData::SendDataAndStateAfterLastTouchVertexLastTime;
  }
  #endif

  return peano::CommunicationSpecification(exchangeMasterWorkerData,exchangeWorkerMasterData,true);
}

peano::MappingSpecification
exahype::mappings::FusedTimeStep::enterCellSpecification(int level) {
  #ifdef Parallel
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true); // counter
  #else
  updateBatchIterationCounter(false); // comes after beginIteration in first iteration -> never init counter
  
  const int coarsestSolverLevel = solvers::Solver::getCoarsestMeshLevelOfAllSolvers();
  if ( std::abs(level)>=coarsestSolverLevel && sendOutRiemannDataInThisIteration() ) {
    return peano::MappingSpecification(
          peano::MappingSpecification::WholeTree,
          peano::MappingSpecification::RunConcurrentlyOnFineGrid,true); // counter
  } else {
    return peano::MappingSpecification(
          peano::MappingSpecification::Nop,
          peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
  }
  #endif
}

peano::MappingSpecification
exahype::mappings::FusedTimeStep::leaveCellSpecification(int level) {
  #ifdef Parallel
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true); // counter
  #else
  updateBatchIterationCounter(false); // comes after beginIteration in first iteration -> never init counter
  
  const int coarsestSolverLevel = solvers::Solver::getCoarsestMeshLevelOfAllSolvers();
  if ( std::abs(level)>=coarsestSolverLevel && issuePredictionJobsInThisIteration() ) {
    return peano::MappingSpecification(
          peano::MappingSpecification::WholeTree,
          peano::MappingSpecification::RunConcurrentlyOnFineGrid,true); // performs reductions
  } else {
    return peano::MappingSpecification(
          peano::MappingSpecification::Nop,
          peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
  }
  #endif
}

peano::MappingSpecification
exahype::mappings::FusedTimeStep::touchVertexFirstTimeSpecification(int level) {
  #ifdef Parallel
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true); // counter
  #else
  updateBatchIterationCounter(false); // comes after beginIteration in first iteration -> never init counter

  const int coarsestSolverLevel = solvers::Solver::getCoarsestMeshLevelOfAllSolvers();
  if ( std::abs(level)>=coarsestSolverLevel && issuePredictionJobsInThisIteration() ) {
    return peano::MappingSpecification(
          peano::MappingSpecification::WholeTree,
          peano::MappingSpecification::RunConcurrentlyOnFineGrid,true); // counter
  } else {
    return peano::MappingSpecification(
          peano::MappingSpecification::Nop,
          peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
  }
  #endif
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

void exahype::mappings::FusedTimeStep::updateBatchIterationCounter(bool initialiseBatchIterationCounter) {
  if (!_batchIterationCounterUpdated) {
    _batchIteration = ( initialiseBatchIterationCounter) ? 0 : _batchIteration+1;
    _batchIterationCounterUpdated = true;

    if ( issuePredictionJobsInThisIteration() ) {
      for (unsigned int solverNumber = 0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
        auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
        solver->beginTimeStep(solver->getMinTimeStamp());
      }
    }
    if ( exahype::solvers::Solver::SpawnPredictionAsBackgroundJob && sendOutRiemannDataInThisIteration() ) {
      peano::datatraversal::TaskSet::startToProcessBackgroundJobs();
    }
  }
}

void exahype::mappings::FusedTimeStep::beginIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);
  static bool isFirst = true;

  _stateCopy = solverState;

  if ( _stateCopy.isFirstIterationOfBatchOrNoBatch() ) {
    updateBatchIterationCounter(true);

    exahype::plotters::startPlottingIfAPlotterIsActive(
        solvers::Solver::getMinTimeStampOfAllSolvers());

    for (auto* solver : exahype::solvers::RegisteredSolvers) {
      solver->setNextMeshUpdateEvent();
    }
  }

#ifdef DistributedStealing

#ifdef StealingStrategyDiffusive
  //if(issuePredictionJobsInThisIteration()) {
  //  if(isFirst) {
  //    exahype::stealing::DiffusiveDistributor::getInstance()._iterationTimer.startTimer(); 
  //  }
  //  else {
  //    exahype::stealing::DiffusiveDistributor::getInstance()._iterationTimer.stopTimer();
  //    int elapsed = static_cast<int> (exahype::stealing::DiffusiveDistributor::getInstance()._iterationTimer.getCalendarTime()*1e06);
  //    exahype::stealing::PerformanceMonitor::getInstance().setCurrentLoad(elapsed);
  //    exahype::stealing::DiffusiveDistributor::getInstance().updateLoadDistribution(elapsed);
  //    exahype::stealing::DiffusiveDistributor::getInstance().resetVictimFlag();
  //    exahype::stealing::DiffusiveDistributor::getInstance()._iterationTimer.startTimer();
  //  }
  //}
#endif

  // enable stealing manager job right at the beginning of the very first time step
//  if(isFirst) {
//    for (auto* solver : exahype::solvers::RegisteredSolvers) {
//      // currently only for ADERDG solver supported
//      if (solver->getType()==exahype::solvers::Solver::Type::ADERDG) {
//        static_cast<exahype::solvers::ADERDGSolver*>(solver)->startStealingManager();
//        isFirst=false;
//      }
//      if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
//        static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->startStealingManager();
//        isFirst=false;
//      }
//    }
//  }
    isFirst = false;

#if defined(StealingStrategyStatic) || defined(StealingStrategyStaticHardcoded)
  if(issuePredictionJobsInThisIteration()) {
    exahype::stealing::StaticDistributor::getInstance().resetRemainingTasksToOffload();
  }
#endif

#endif
  logTraceOutWith1Argument("beginIteration(State)", solverState);
}

void exahype::mappings::FusedTimeStep::endIteration(
    exahype::State& state) {
  logTraceInWith1Argument("endIteration(State)", state);

  if ( sendOutRiemannDataInThisIteration() ) {
    exahype::plotters::finishedPlotting();
    
    const int isFirstTimeStep = 
          ( exahype::solvers::Solver::PredictionSweeps==1 ) ? 
          _stateCopy.isFirstIterationOfBatchOrNoBatch() : 
          _stateCopy.isSecondIterationOfBatchOrNoBatch();

    if ( _stateCopy.isLastIterationOfBatchOrNoBatch() ) {
      // background threads
      exahype::solvers::Solver::ensureAllJobsHaveTerminated(exahype::solvers::Solver::JobType::ReductionJob);
    }

    exahype::solvers::Solver::startNewTimeStepForAllSolvers(
        isFirstTimeStep,_stateCopy.isLastIterationOfBatchOrNoBatch(),
        true);
  }

  _batchIterationCounterUpdated = false;

  //logInfo("endIteration(State)", _stateCopy.getBatchIteration() << ", "<<state.getBatchIteration() << ", " << _batchIteration);
  #ifdef Parallel
  // broadcasts - must come after the last commSpec evaluation for beginIteration(...)
  if ( _stateCopy.isFirstIterationOfBatchOrNoBatch() ) { // this is after the broadcast
    assertion1(exahype::State::BroadcastInThisIteration==true,tarch::parallel::Node::getInstance().getRank());
    exahype::State::BroadcastInThisIteration = false;
  }
  if ( _stateCopy.isLastIterationOfBatchOrNoBatch() ) {
    assertion(exahype::State::BroadcastInThisIteration==false);
    exahype::State::BroadcastInThisIteration = true;
  }
  // reduction (must come here; is before first commSpec evaluation for reduction
  if ( _stateCopy.isSecondToLastIterationOfBatchOrNoBatch() ) { // this is after the broadcast
    assertion(exahype::State::ReduceInThisIteration==false);
    exahype::State::ReduceInThisIteration = true;
  }
  if ( _stateCopy.isLastIterationOfBatchOrNoBatch() ) {
    assertion(exahype::State::ReduceInThisIteration==true);
    exahype::State::ReduceInThisIteration = false;
  }
  #endif

  logTraceOutWith1Argument("endIteration(State)", state);
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::FusedTimeStep::FusedTimeStep(
    const FusedTimeStep& masterThread) :
  _stateCopy(masterThread._stateCopy), 
  _batchIterationCounterUpdated(masterThread._batchIterationCounterUpdated),
  _batchIteration(masterThread._batchIteration) {
  // do nothing
}
// Merge over threads
void exahype::mappings::FusedTimeStep::mergeWithWorkerThread(
    const FusedTimeStep& workerThread) {
  // do nothing
}
#endif

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

  updateBatchIterationCounter(false);

  if ( issuePredictionJobsInThisIteration() ) {
    fineGridVertex.mergeNeighbours(fineGridX,fineGridH);
  }

  logTraceOutWith1Argument("touchVertexFirstTime(...)", fineGridVertex);
}

void exahype::mappings::FusedTimeStep::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,fineGridVerticesEnumerator.toString(),coarseGridCell, fineGridPositionOfCell);

  if (
      fineGridCell.isInitialised() &&
      sendOutRiemannDataInThisIteration()
   ) {
     solvers::Solver::CellInfo cellInfo = fineGridCell.createCellInfo();
     const bool isAtRemoteBoundary = exahype::Cell::isAtRemoteBoundary(fineGridVertices,fineGridVerticesEnumerator);

     for (int solverNumber=0; solverNumber<static_cast<int>(solvers::RegisteredSolvers.size()); solverNumber++) {
       auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
       switch ( solver->getType() ) {
         case solvers::Solver::Type::ADERDG:
           static_cast<solvers::ADERDGSolver*>(solver)->prolongateFaceData(solverNumber,cellInfo,isAtRemoteBoundary);
           break;
         case solvers::Solver::Type::LimitingADERDG:
           static_cast<solvers::LimitingADERDGSolver*>(solver)->getSolver()->prolongateFaceData(solverNumber,cellInfo,isAtRemoteBoundary);
           break;
         case solvers::Solver::Type::FiniteVolumes:
           // insert code here
           break;
         default:
           assertionMsg(false,"Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
           logError("mergeWithBoundaryDataIfNotDoneYet(...)","Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
           std::abort();
           break;
       }
     }
  }

  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}



void exahype::mappings::FusedTimeStep::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("leaveCell(...)", fineGridCell,fineGridVerticesEnumerator.toString(),coarseGridCell, fineGridPositionOfCell);

  if (
      issuePredictionJobsInThisIteration() &&
      fineGridCell.isInitialised()
  ) {
    solvers::Solver::CellInfo cellInfo = fineGridCell.createCellInfo();
    const bool isAtRemoteBoundary = exahype::Cell::isAtRemoteBoundary(fineGridVertices,fineGridVerticesEnumerator);
    const int isLastTimeStep =
        ( exahype::solvers::Solver::PredictionSweeps==1 ) ?
            _stateCopy.isLastIterationOfBatchOrNoBatch() :
            _stateCopy.isSecondToLastIterationOfBatchOrNoBatch(); // PredictionSweeps==2

    for (int solverNumber=0; solverNumber<static_cast<int>(solvers::RegisteredSolvers.size()); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

      // this operates only on compute cells
      plotters::plotPatchIfAPlotterIsActive(solverNumber,cellInfo); // TODO(Dominic) potential for IO overlap?

      solvers::Solver::UpdateResult result;
      switch ( solver->getType() ) {
        case solvers::Solver::Type::ADERDG:
          result = static_cast<solvers::ADERDGSolver*>(solver)->fusedTimeStepOrRestrict(
              solverNumber,cellInfo,_stateCopy.isFirstIterationOfBatchOrNoBatch(),isLastTimeStep,isAtRemoteBoundary);
          break;
        case solvers::Solver::Type::LimitingADERDG:
          result = static_cast<solvers::LimitingADERDGSolver*>(solver)->fusedTimeStepOrRestrict(
              solverNumber,cellInfo,_stateCopy.isFirstIterationOfBatchOrNoBatch(),isLastTimeStep,isAtRemoteBoundary);
          break;
        case solvers::Solver::Type::FiniteVolumes:
          result = static_cast<solvers::FiniteVolumesSolver*>(solver)->fusedTimeStepOrRestrict(
              solverNumber,cellInfo,_stateCopy.isFirstIterationOfBatchOrNoBatch(),isLastTimeStep,isAtRemoteBoundary);
          break;
        default:
          assertionMsg(false,"Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
          logError("mergeWithBoundaryDataIfNotDoneYet(...)","Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
          std::abort();
          break;
      }
      // mesh refinement events, cell sizes (for AMR), time
      if ( isLastTimeStep ) {
        tarch::multicore::Lock lock(Semaphore);
        {
          solver->updateNextMeshUpdateEvent(result._meshUpdateEvent);
          solver->updateNextMaxLevel(fineGridVerticesEnumerator.getLevel());
          solver->updateMinNextTimeStepSize(result._timeStepSize);
        }
        lock.free();
      }
    }

    // Must be performed for all cell descriptions
    Cell::resetNeighbourMergePerformedFlags(cellInfo,fineGridVertices,fineGridVerticesEnumerator);
  }

  logTraceOutWith1Argument("leaveCell(...)", fineGridCell);
}

#ifdef Parallel
void exahype::mappings::FusedTimeStep::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  logTraceInWith6Arguments( "mergeWithNeighbour(...)", vertex, neighbour, fromRank, fineGridX, fineGridH, level );

  updateBatchIterationCounter(false);

  if ( issuePredictionJobsInThisIteration() ) {
    vertex.receiveNeighbourData(
        fromRank, true/*merge with data*/,_stateCopy.isFirstIterationOfBatchOrNoBatch(),
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
    vertex.sendToNeighbour(toRank,_stateCopy.isLastIterationOfBatchOrNoBatch(),x,h,level);
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

  if ( _stateCopy.isFirstIterationOfBatchOrNoBatch() ) {
    exahype::Cell::broadcastGlobalDataToWorker(
        worker,
        fineGridVerticesEnumerator.getCellCenter(),
        fineGridVerticesEnumerator.getLevel());
  }

  logTraceInWith1Argument( "prepareSendToWorker(...)", fineGridCell );

  return _stateCopy.isFirstIterationOfBatchOrNoBatch() ||
         _stateCopy.isLastIterationOfBatchOrNoBatch();
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

  if ( _stateCopy.isFirstIterationOfBatchOrNoBatch() ) {
    exahype::Cell::mergeWithGlobalDataFromMaster(
        tarch::parallel::NodePool::getInstance().getMasterRank(),
        receivedVerticesEnumerator.getCellCenter(),
        receivedVerticesEnumerator.getLevel());
  }

  logTraceOut( "receiveDataFromMaster(...)" );
}

void exahype::mappings::FusedTimeStep::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
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

  if ( _stateCopy.isLastIterationOfBatchOrNoBatch() ) {
    exahype::Cell::reduceGlobalDataToMaster(
        tarch::parallel::NodePool::getInstance().getMasterRank(),
        verticesEnumerator.getCellCenter(),
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

  if ( _stateCopy.isLastIterationOfBatchOrNoBatch() ) {
    exahype::Cell::mergeWithGlobalDataFromWorker(
        worker,
        fineGridVerticesEnumerator.getCellCenter(),
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
