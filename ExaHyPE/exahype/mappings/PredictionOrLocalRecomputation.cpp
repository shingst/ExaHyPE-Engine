/**os
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
 
#include "exahype/mappings/PredictionOrLocalRecomputation.h"

#include <algorithm>

#include "tarch/multicore/Loop.h"

#include "peano/utils/Globals.h"
#include "peano/utils/Loop.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/datatraversal/TaskSet.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "exahype/mappings/Prediction.h"

tarch::logging::Log exahype::mappings::PredictionOrLocalRecomputation::_log(
    "exahype::mappings::PredictionOrLocalRecomputation");

bool exahype::mappings::PredictionOrLocalRecomputation::OneSolverRequestedLocalRecomputation = false;

peano::CommunicationSpecification
exahype::mappings::PredictionOrLocalRecomputation::communicationSpecification() const {
  // master->worker
  peano::CommunicationSpecification::ExchangeMasterWorkerData exchangeMasterWorkerData =
      peano::CommunicationSpecification::ExchangeMasterWorkerData::MaskOutMasterWorkerDataAndStateExchange;
  #ifdef Parallel
  if (
      exahype::solvers::Solver::PredictionSweeps==1 ||
      !exahype::solvers::Solver::FuseADERDGPhases   ||
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
      !exahype::solvers::Solver::FuseADERDGPhases   ||
      exahype::State::ReduceInThisIteration         // must be set in previous iteration
  ) {
    exchangeWorkerMasterData =
        peano::CommunicationSpecification::ExchangeWorkerMasterData::SendDataAndStateAfterLastTouchVertexLastTime;
  }
  #endif

  return peano::CommunicationSpecification(exchangeMasterWorkerData,exchangeWorkerMasterData,true);
}

peano::MappingSpecification exahype::mappings::PredictionOrLocalRecomputation::enterCellSpecification(int level) const {
  const int coarsestSolverLevel = solvers::Solver::getCoarsestMeshLevelOfAllSolvers();
  if ( std::abs(level)>=coarsestSolverLevel ) {
    return peano::MappingSpecification(
           peano::MappingSpecification::WholeTree,
           peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);  // performs reduction
  } else {
    return peano::MappingSpecification(
          peano::MappingSpecification::Nop,
          peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
  }
}
peano::MappingSpecification
exahype::mappings::PredictionOrLocalRecomputation::touchVertexFirstTimeSpecification(int level) const {
  return Vertex::getNeighbourMergeSpecification(level);
}

// Below specs are all nop
peano::MappingSpecification
exahype::mappings::PredictionOrLocalRecomputation::leaveCellSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
}

peano::MappingSpecification
exahype::mappings::PredictionOrLocalRecomputation::touchVertexLastTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}
peano::MappingSpecification
exahype::mappings::PredictionOrLocalRecomputation::ascendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}
peano::MappingSpecification
exahype::mappings::PredictionOrLocalRecomputation::descendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}


void exahype::mappings::PredictionOrLocalRecomputation::initialiseLocalVariables(){
  const unsigned int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
  _minTimeStepSizes.resize(numberOfSolvers);
  _maxLevels.resize(numberOfSolvers);

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    _minTimeStepSizes[solverNumber] = std::numeric_limits<double>::max();
    _maxLevels       [solverNumber] = -std::numeric_limits<int>::max(); // "-", min
  }
}

exahype::mappings::PredictionOrLocalRecomputation::PredictionOrLocalRecomputation() {
  // do nothing
}

exahype::mappings::PredictionOrLocalRecomputation::~PredictionOrLocalRecomputation() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::PredictionOrLocalRecomputation::PredictionOrLocalRecomputation(
    const PredictionOrLocalRecomputation& masterThread) :
  _stateCopy(masterThread._stateCopy) {
  initialiseLocalVariables();
}
// Merge over threads
void exahype::mappings::PredictionOrLocalRecomputation::mergeWithWorkerThread(
    const PredictionOrLocalRecomputation& workerThread) {
  for (unsigned int solverNumber = 0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    _minTimeStepSizes[solverNumber] =
        std::min(_minTimeStepSizes[solverNumber], workerThread._minTimeStepSizes[solverNumber]);
    _maxLevels[solverNumber] =
        std::max(_maxLevels[solverNumber], workerThread._maxLevels[solverNumber]);
  }
}
#endif

void exahype::mappings::PredictionOrLocalRecomputation::beginIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);

  _stateCopy = solverState;

  if ( _stateCopy.isFirstIterationOfBatchOrNoBatch() ) {
    OneSolverRequestedLocalRecomputation =
        exahype::solvers::Solver::oneSolverRequestedLocalRecomputation();

    initialiseLocalVariables();
  }

  if (
      exahype::solvers::Solver::SpawnPredictionAsBackgroundJob &&
      _stateCopy.isLastIterationOfBatchOrNoBatch()
  ) {
    peano::datatraversal::TaskSet::startToProcessBackgroundJobs();
  }

  logTraceOutWith1Argument("beginIteration(State)", solverState);
}

bool exahype::mappings::PredictionOrLocalRecomputation::performLocalRecomputation(
    exahype::solvers::Solver* solver) {
  return
      solver->getMeshUpdateEvent()==exahype::solvers::Solver::MeshUpdateEvent::IrregularLimiterDomainChange;
}

bool exahype::mappings::PredictionOrLocalRecomputation::performPrediction(
    exahype::solvers::Solver* solver) {
  return exahype::solvers::Solver::FuseADERDGPhases &&
         solver->hasRequestedMeshRefinement();
}

void exahype::mappings::PredictionOrLocalRecomputation::endIteration(
    exahype::State& state) {
  logTraceInWith1Argument("endIteration(State)", state);

  if (
      _stateCopy.isLastIterationOfBatchOrNoBatch() &&
      OneSolverRequestedLocalRecomputation
  ) {
    for (unsigned int solverNumber = 0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      if ( performLocalRecomputation( solver ) ) {
        assertion1(std::isfinite(_minTimeStepSizes[solverNumber]),_minTimeStepSizes[solverNumber]);
        assertion1(_minTimeStepSizes[solverNumber]>0.0,_minTimeStepSizes[solverNumber]);

        solver->updateNextMaxLevel(_maxLevels[solverNumber]);

        solver->updateMinNextTimeStepSize(_minTimeStepSizes[solverNumber]);

        logDebug("endIteration(state)","[pre] solver="<<solver->toString());
        if (
            exahype::solvers::Solver::FuseADERDGPhases
            #ifdef Parallel
            && tarch::parallel::Node::getInstance().isGlobalMaster()
            #endif
        ) {
          exahype::solvers::Solver::
          reinitialiseTimeStepDataIfLastPredictorTimeStepSizeWasInstable(solver);
        }
        if (exahype::solvers::Solver::FuseADERDGPhases) {
          solver->startNewTimeStepFused(true,true);
        } else {
          solver->startNewTimeStep();
        }

        logDebug("endIteration(state)","[post] updatedTimeStepSize="<<solver->getMinTimeStepSize()<<", solver="<<solver->toString());
      }
    }
  }

  #ifdef Parallel
  // broadcasts
  if ( _stateCopy.isFirstIterationOfBatchOrNoBatch() ) { // this is after the broadcast
    assertion(exahype::State::BroadcastInThisIteration==true);
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

void exahype::mappings::PredictionOrLocalRecomputation::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  if ( fineGridCell.isInitialised() ) {
    solvers::Solver::CellInfo cellInfo = fineGridCell.createCellInfo();
    const bool isAtRemoteBoundary = exahype::Cell::isAtRemoteBoundary(fineGridVertices,fineGridVerticesEnumerator);

    for (unsigned int solverNumber=0; solverNumber<exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      if ( performLocalRecomputation( solver ) && _stateCopy.isFirstIterationOfBatchOrNoBatch() ) {
        auto* limitingADERDG = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
        double admissibleTimeStepSize =
            limitingADERDG->recomputeSolutionLocally(solverNumber,cellInfo,isAtRemoteBoundary,solvers::Solver::FuseADERDGPhases);

        _minTimeStepSizes[solverNumber] = std::min(
            admissibleTimeStepSize, _minTimeStepSizes[solverNumber]);
        _maxLevels[solverNumber] = std::max(
            fineGridVerticesEnumerator.getLevel(),_maxLevels[solverNumber]);

        limitingADERDG->determineMinAndMax(solverNumber,cellInfo); // TODO(Dominic): Optimistation. Do it only in recomputed cells.
      }
      else if ( performPrediction(solver) && _stateCopy.isFirstIterationOfBatchOrNoBatch() ) {
        switch ( solver->getType() ) {
          case solvers::Solver::Type::ADERDG:
            static_cast<solvers::ADERDGSolver*>(solver)->performPredictionAndVolumeIntegral(solverNumber,cellInfo,isAtRemoteBoundary);
            break;
          case solvers::Solver::Type::LimitingADERDG:
            static_cast<solvers::LimitingADERDGSolver*>(solver)->getSolver()->performPredictionAndVolumeIntegral(solverNumber,cellInfo,isAtRemoteBoundary);
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
      else if ( performPrediction(solver) && _stateCopy.isLastIterationOfBatchOrNoBatch() ) {
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
    if ( _stateCopy.isFirstIterationOfBatchOrNoBatch() ) {
      Cell::resetNeighbourMergeFlagsAndCounters(cellInfo,fineGridVertices,fineGridVerticesEnumerator);
    }
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}


void exahype::mappings::PredictionOrLocalRecomputation::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::PredictionOrLocalRecomputation::mergeNeighboursDataDuringLocalRecomputationLoopBody(
    const int                                    pos1Scalar,
    const int                                    pos2Scalar,
    const exahype::Vertex&                       vertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h) {
  tarch::la::Vector<DIMENSIONS,int> pos1 = Vertex::delineariseIndex2(pos1Scalar);
  tarch::la::Vector<DIMENSIONS,int> pos2 = Vertex::delineariseIndex2(pos2Scalar);
  assertion(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1));

  const int  cellDescriptionsIndex1 = vertex.getCellDescriptionsIndex(pos1Scalar);
  const int  cellDescriptionsIndex2 = vertex.getCellDescriptionsIndex(pos2Scalar);
  const bool validIndex1 = cellDescriptionsIndex1 >= 0;
  const bool validIndex2 = cellDescriptionsIndex2 >= 0;

  if ( validIndex1 && validIndex2 ) {
    solvers::Solver::CellInfo cellInfo1 = vertex.createCellInfo(pos1Scalar);
    solvers::Solver::CellInfo cellInfo2 = vertex.createCellInfo(pos2Scalar);
    solvers::Solver::InterfaceInfo face(pos1,pos2);

    const bool mergeWithCell1 = Vertex::hasToMergeAtFace(cellInfo1,face._faceIndex1,false/*prefetchADERDGFace*/);
    const bool mergeWithCell2 = Vertex::hasToMergeAtFace(cellInfo2,face._faceIndex2,false/*prefetchADERDGFace*/);
    assertion2(mergeWithCell1==mergeWithCell2,mergeWithCell1,mergeWithCell2);
    if ( mergeWithCell1 && mergeWithCell2 ) {
      for (int solverNumber=0; solverNumber<static_cast<int>(solvers::RegisteredSolvers.size()); solverNumber++) {
        auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
        if ( performLocalRecomputation(solver) ) {
          static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->
              mergeNeighboursData(solverNumber,cellInfo1,cellInfo2,pos1,pos2,true/* isRecomputation */);
        }
      }
    }
  } else if (
      ((validIndex1 && !validIndex2 &&
        cellDescriptionsIndex2==multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex)
        ||
        (!validIndex1 && validIndex2 &&
        cellDescriptionsIndex1==multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex))
  ) {
    int                               posCellScalar = pos1Scalar;
    tarch::la::Vector<DIMENSIONS,int> posCell       = pos1;
    tarch::la::Vector<DIMENSIONS,int> posBoundary   = pos2;
    if ( cellDescriptionsIndex2 >= 0 ) {
      posCellScalar         = pos2Scalar;
      posCell               = pos2;
      posBoundary           = pos1;
    }
    solvers::Solver::CellInfo         cellInfo = vertex.createCellInfo(posCellScalar);
    solvers::Solver::BoundaryFaceInfo face(posCell,posBoundary);

    if ( !cellInfo.empty() && Vertex::hasToMergeAtFace(cellInfo,face._faceIndex,false/*prefetchADERDGFace*/) ) {
      for (int solverNumber=0; solverNumber<static_cast<int>(solvers::RegisteredSolvers.size()); solverNumber++) {
        auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
        if ( performLocalRecomputation(solver) ) {
          static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->
              mergeWithBoundaryData(solverNumber,cellInfo,posCell,posBoundary,true);
        }
      }
    }
  }
}

void exahype::mappings::PredictionOrLocalRecomputation::touchVertexFirstTime(
  exahype::Vertex& vertex,
  const tarch::la::Vector<DIMENSIONS, double>& x,
  const tarch::la::Vector<DIMENSIONS, double>& h,
  exahype::Vertex* const coarseGridVertices,
  const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
  exahype::Cell& coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>& fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "touchVertexFirstTime(...)", vertex, x, h, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );

  if (
      _stateCopy.isFirstIterationOfBatchOrNoBatch() &&
      OneSolverRequestedLocalRecomputation
  ) {
    for (int i=0; i<2*(DIMENSIONS-1)*(DIMENSIONS); i++) {
      mergeNeighboursDataDuringLocalRecomputationLoopBody(Vertex::pos1Scalar[i],Vertex::pos2Scalar[i],vertex,x,h);
    }
  }

  logTraceOutWith1Argument( "touchVertexFirstTime(...)", vertex );
}


#ifdef Parallel
///////////////////////////////////////
// NEIGHBOUR
///////////////////////////////////////
void exahype::mappings::PredictionOrLocalRecomputation::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  logTraceInWith6Arguments( "mergeWithNeighbour(...)", vertex, neighbour, fromRank, fineGridX, fineGridH, level );
  
  if (
      _stateCopy.isFirstIterationOfBatchOrNoBatch() &&
      OneSolverRequestedLocalRecomputation &&
      vertex.hasToCommunicate(level)
  ) {
    for (unsigned int i = 2*(DIMENSIONS-1)*(DIMENSIONS); i-- > 0;) { // dest and src is swapped & order is swapped
      receiveNeighbourDataLoopBody(fromRank,Vertex::pos1Scalar[i],Vertex::pos2Scalar[i],vertex,fineGridX,level);
      receiveNeighbourDataLoopBody(fromRank,Vertex::pos2Scalar[i],Vertex::pos1Scalar[i],vertex,fineGridX,level);
    }
  }
  logTraceOut( "mergeWithNeighbour(...)" );
}

void exahype::mappings::PredictionOrLocalRecomputation::receiveNeighbourDataLoopBody(
    const int                                    fromRank,
    const int                                    srcScalar,
    const int                                    destScalar,
    const exahype::Vertex&                       vertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  if ( vertex.hasToReceiveMetadata(fromRank,srcScalar,destScalar,vertex.getAdjacentRanks()) ) {
    const int destCellDescriptionsIndex = vertex.getCellDescriptionsIndex(destScalar);
    bool validIndex = destCellDescriptionsIndex >= 0;
    assertion( !validIndex || exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(destCellDescriptionsIndex));

    if ( validIndex ) {
      const tarch::la::Vector<DIMENSIONS,int> src = Vertex::delineariseIndex2(srcScalar);
      const tarch::la::Vector<DIMENSIONS,int> dest = Vertex::delineariseIndex2(destScalar);
      solvers::Solver::CellInfo         cellInfo = vertex.createCellInfo(destScalar);
      solvers::Solver::BoundaryFaceInfo face(dest,src); // dest and src are swapped

      if ( Vertex::hasToReceiveFromNeighbourNow(cellInfo,face,false/*prefetchADERDGFaceData*/) ) {
        for(unsigned int solverNumber = solvers::RegisteredSolvers.size(); solverNumber-- > 0;) {
          auto* solver = solvers::RegisteredSolvers[solverNumber];
          if ( performLocalRecomputation( solver ) ) {
            assertion1( solver->getType()==solvers::Solver::Type::LimitingADERDG, solver->toString() );
            static_cast<solvers::LimitingADERDGSolver*>(solver)->
                mergeWithNeighbourDataBasedOnLimiterStatus(fromRank,solverNumber,cellInfo,src,dest,true/*isRecomputation*/,x,level);
          }
        }
      }
    }
  }
}

void exahype::mappings::PredictionOrLocalRecomputation::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  logTraceInWith3Arguments( "prepareSendToNeighbour(...)", vertex, toRank, level );

  if (
      _stateCopy.isLastIterationOfBatchOrNoBatch() &&
      exahype::solvers::Solver::FuseADERDGPhases
  ) {
   vertex.sendToNeighbour(toRank,true,x,h,level);
  }

  logTraceOut( "prepareSendToNeighbour(...)" );
}

bool exahype::mappings::PredictionOrLocalRecomputation::prepareSendToWorker(
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

  logTraceOutWith1Argument( "prepareSendToWorker(...)", true );
  return true; // this must be sent to worker in first broadcast
}

void exahype::mappings::PredictionOrLocalRecomputation::receiveDataFromMaster(
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

void exahype::mappings::PredictionOrLocalRecomputation::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith2Arguments( "prepareSendToMaster(...)", localCell, verticesEnumerator.toString() );

  if ( _stateCopy.isLastIterationOfBatchOrNoBatch() ) {
    for (auto* solver : exahype::solvers::RegisteredSolvers) {
      if ( performLocalRecomputation(solver) ) {
        solver->sendDataToMaster(
            tarch::parallel::NodePool::getInstance().getMasterRank(),
            verticesEnumerator.getCellCenter(),
            verticesEnumerator.getLevel());
      }
    }
  }

  logTraceOut( "prepareSendToMaster(...)" );
}

void exahype::mappings::PredictionOrLocalRecomputation::mergeWithMaster(
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
    for (auto* solver : exahype::solvers::RegisteredSolvers) {
      if ( performLocalRecomputation(solver) ) {
        solver->mergeWithWorkerData(
            worker,
            fineGridVerticesEnumerator.getCellCenter(),
            fineGridVerticesEnumerator.getLevel());
      }
    }
  }

  logTraceOut( "mergeWithMaster(...)" );
}

//
// Below all methods are nop.
//
//=====================================


void exahype::mappings::PredictionOrLocalRecomputation::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::PredictionOrLocalRecomputation::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::PredictionOrLocalRecomputation::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::PredictionOrLocalRecomputation::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::PredictionOrLocalRecomputation::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::PredictionOrLocalRecomputation::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::PredictionOrLocalRecomputation::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::PredictionOrLocalRecomputation::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::PredictionOrLocalRecomputation::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::PredictionOrLocalRecomputation::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::PredictionOrLocalRecomputation::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::PredictionOrLocalRecomputation::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::PredictionOrLocalRecomputation::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::PredictionOrLocalRecomputation::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::PredictionOrLocalRecomputation::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::PredictionOrLocalRecomputation::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
