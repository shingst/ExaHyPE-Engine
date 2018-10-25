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
 
#include "exahype/mappings/PredictionOrLocalRecomputation.h"

#include <algorithm>

#include "tarch/multicore/Loop.h"

#include "peano/utils/Globals.h"
#include "peano/utils/Loop.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/datatraversal/TaskSet.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/VertexOperations.h"

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

peano::MappingSpecification
exahype::mappings::PredictionOrLocalRecomputation::enterCellSpecification(int level) const {
  if (
      exahype::solvers::Solver::FuseADERDGPhases &&
      exahype::solvers::Solver::oneSolverRequestedMeshRefinement()
  ) {
    return exahype::mappings::Prediction::determineEnterLeaveCellSpecification(level);
  } else {
    return peano::MappingSpecification(
        peano::MappingSpecification::WholeTree,
        peano::MappingSpecification::AvoidFineGridRaces,true); // TODO(Dominic): false should work in theory
  }
}
peano::MappingSpecification
exahype::mappings::PredictionOrLocalRecomputation::touchVertexFirstTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces,true); // TODO(Dominic): false should work in theory
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

exahype::mappings::PredictionOrLocalRecomputation::PredictionOrLocalRecomputation()
  #ifdef Debug
  :
  _interiorFaceMerges(0),
  _boundaryFaceMerges(0)
  #endif
{
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

  #ifdef Debug // TODO(Dominic): And not parallel and not shared memory
  _interiorFaceMerges = 0;
  _boundaryFaceMerges = 0;
  #endif

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

    #if defined(Debug) // TODO(Dominic): Use logDebug if it works with filters
    logInfo("endIteration(...)","interiorFaceSolves: " << _interiorFaceMerges);
    logInfo("endIteration(...)","boundaryFaceSolves: " << _boundaryFaceMerges);
    #endif
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
    const int cellDescriptionsIndex = fineGridCell.getCellDescriptionsIndex();

    const int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
    for (int solverNumber=0; solverNumber<numberOfSolvers; solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      const int element = solver->tryGetElement(cellDescriptionsIndex,solverNumber);
      if ( element!=exahype::solvers::Solver::NotFound ) {
        if (
            performLocalRecomputation( solver ) &&
            _stateCopy.isFirstIterationOfBatchOrNoBatch()
        ) {
          auto* limitingADERDG = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);

          double admissibleTimeStepSize = std::numeric_limits<double>::max();
          if ( exahype::solvers::Solver::FuseADERDGPhases ) {
            admissibleTimeStepSize = limitingADERDG->recomputeSolutionLocallyFused(
                cellDescriptionsIndex,element,
                exahype::Cell::isAtRemoteBoundary(
                    fineGridVertices,fineGridVerticesEnumerator));
          } else {
            admissibleTimeStepSize = limitingADERDG->recomputeSolutionLocally(
                cellDescriptionsIndex,element);
          }

          _minTimeStepSizes[solverNumber] = std::min(
              admissibleTimeStepSize, _minTimeStepSizes[solverNumber]);
          _maxLevels[solverNumber] = std::max(
              fineGridVerticesEnumerator.getLevel(),_maxLevels[solverNumber]);

          limitingADERDG->determineMinAndMax(cellDescriptionsIndex,element);
        }
        else if ( performPrediction(solver) ) {
          if ( _stateCopy.isFirstIterationOfBatchOrNoBatch() ) {
            // this operates only on compute cells
            exahype::solvers::ADERDGSolver::performPredictionAndVolumeIntegral(
                solver,fineGridCell.getCellDescriptionsIndex(),element,
                exahype::Cell::isAtRemoteBoundary(
                    fineGridVertices,fineGridVerticesEnumerator)
            );
          }
          if ( _stateCopy.isLastIterationOfBatchOrNoBatch() ) { // we are sure here that the skeleton STPs have finished
            // this operates only on helper cells
            solver->prolongateFaceData(
                fineGridCell.getCellDescriptionsIndex(),element,
                exahype::Cell::isAtRemoteBoundary(
                    fineGridVertices,fineGridVerticesEnumerator));
          }
        }
      }
    }

    if ( _stateCopy.isFirstIterationOfBatchOrNoBatch() ) {
      exahype::Cell::resetNeighbourMergeFlags(
          fineGridCell.getCellDescriptionsIndex(),
          fineGridVertices,fineGridVerticesEnumerator);
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
    const int pos1Scalar,
    const int pos2Scalar,
    const int cellDescriptionsIndex1,
    const int cellDescriptionsIndex2,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h) {
  tarch::la::Vector<DIMENSIONS,int> pos1 = Vertex::delineariseIndex2(pos1Scalar);
  tarch::la::Vector<DIMENSIONS,int> pos2 = Vertex::delineariseIndex2(pos2Scalar);
  assertion(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1));

  const bool validIndex1 = cellDescriptionsIndex1 >= 0;
  const bool validIndex2 = cellDescriptionsIndex2 >= 0;
  assertion(cellDescriptionsIndex1 < 0 || exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex1));
  assertion(cellDescriptionsIndex2 < 0 || exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex2));

  if ( validIndex1 && validIndex2 ) {
    auto& ADERDGPatches1 = solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex1);
    auto& FVPatches1     = solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex1);

    auto& ADERDGPatches2 = solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex2);
    auto& FVPatches2     = solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex2);

    bool mergeNeighbours = ( !ADERDGPatches1.empty() && !ADERDGPatches2.empty() ) ||
                           ( !FVPatches1.empty() && !FVPatches2.empty() );

    if ( mergeNeighbours ) {
      for (int solverNumber=0; solverNumber<static_cast<int>(solvers::RegisteredSolvers.size()); solverNumber++) {
        auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
        if ( performLocalRecomputation(solver) ) {
          static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->
              mergeNeighboursData(
                  ADERDGPatches1,ADERDGPatches2,FVPatches1,FVPatches2,solverNumber,
                  pos1,pos2,true/* isRecomputation */);
        }
        #ifdef Debug // TODO(Dominic)
        _interiorFaceMerges++;
        #endif
      }
    }
  } else if (
      ((validIndex1 && !validIndex2 &&
        cellDescriptionsIndex2==multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex)
        ||
        (!validIndex1 && validIndex2 &&
        cellDescriptionsIndex1==multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex))
  ) {
    tarch::la::Vector<DIMENSIONS,int> posCell     = pos1;
    tarch::la::Vector<DIMENSIONS,int> posBoundary = pos2;
    int cellDescriptionsIndex                      = cellDescriptionsIndex1;
    if ( cellDescriptionsIndex2 >= 0 ) {
      posCell     = pos2;
      posBoundary = pos1;
      cellDescriptionsIndex = cellDescriptionsIndex2;
    }

    auto& ADERDGPatches = solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex);
    auto& FVPatches     = solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex);

    for (int solverNumber=0; solverNumber<static_cast<int>(solvers::RegisteredSolvers.size()); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      if ( performLocalRecomputation(solver) ) {
        static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->
            mergeWithBoundaryData(ADERDGPatches,FVPatches,solverNumber,posCell,posBoundary,true);
        #ifdef Debug
        _boundaryFaceMerges++;
        #endif
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
    #if DIMENSIONS==2
    mergeNeighboursDataDuringLocalRecomputationLoopBody(0,1,VertexOperations::readCellDescriptionsIndex(vertex,0),VertexOperations::readCellDescriptionsIndex(vertex,1),x,h);
    mergeNeighboursDataDuringLocalRecomputationLoopBody(0,2,VertexOperations::readCellDescriptionsIndex(vertex,0),VertexOperations::readCellDescriptionsIndex(vertex,2),x,h);
    mergeNeighboursDataDuringLocalRecomputationLoopBody(1,3,VertexOperations::readCellDescriptionsIndex(vertex,1),VertexOperations::readCellDescriptionsIndex(vertex,3),x,h);
    mergeNeighboursDataDuringLocalRecomputationLoopBody(2,3,VertexOperations::readCellDescriptionsIndex(vertex,2),VertexOperations::readCellDescriptionsIndex(vertex,3),x,h);
    #elif DIMENSIONS==3
    mergeNeighboursDataDuringLocalRecomputationLoopBody(0,1,VertexOperations::readCellDescriptionsIndex(vertex,0),VertexOperations::readCellDescriptionsIndex(vertex,1),x,h);
    mergeNeighboursDataDuringLocalRecomputationLoopBody(0,2,VertexOperations::readCellDescriptionsIndex(vertex,0),VertexOperations::readCellDescriptionsIndex(vertex,2),x,h);
    mergeNeighboursDataDuringLocalRecomputationLoopBody(0,4,VertexOperations::readCellDescriptionsIndex(vertex,0),VertexOperations::readCellDescriptionsIndex(vertex,4),x,h);
    mergeNeighboursDataDuringLocalRecomputationLoopBody(1,3,VertexOperations::readCellDescriptionsIndex(vertex,1),VertexOperations::readCellDescriptionsIndex(vertex,3),x,h);
    mergeNeighboursDataDuringLocalRecomputationLoopBody(1,5,VertexOperations::readCellDescriptionsIndex(vertex,1),VertexOperations::readCellDescriptionsIndex(vertex,5),x,h);
    mergeNeighboursDataDuringLocalRecomputationLoopBody(2,3,VertexOperations::readCellDescriptionsIndex(vertex,2),VertexOperations::readCellDescriptionsIndex(vertex,3),x,h);
    mergeNeighboursDataDuringLocalRecomputationLoopBody(2,6,VertexOperations::readCellDescriptionsIndex(vertex,2),VertexOperations::readCellDescriptionsIndex(vertex,6),x,h);
    mergeNeighboursDataDuringLocalRecomputationLoopBody(3,7,VertexOperations::readCellDescriptionsIndex(vertex,3),VertexOperations::readCellDescriptionsIndex(vertex,7),x,h);
    mergeNeighboursDataDuringLocalRecomputationLoopBody(4,5,VertexOperations::readCellDescriptionsIndex(vertex,4),VertexOperations::readCellDescriptionsIndex(vertex,5),x,h);
    mergeNeighboursDataDuringLocalRecomputationLoopBody(4,6,VertexOperations::readCellDescriptionsIndex(vertex,4),VertexOperations::readCellDescriptionsIndex(vertex,6),x,h);
    mergeNeighboursDataDuringLocalRecomputationLoopBody(5,7,VertexOperations::readCellDescriptionsIndex(vertex,5),VertexOperations::readCellDescriptionsIndex(vertex,7),x,h);
    mergeNeighboursDataDuringLocalRecomputationLoopBody(6,7,VertexOperations::readCellDescriptionsIndex(vertex,6),VertexOperations::readCellDescriptionsIndex(vertex,7),x,h);
    #endif
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

//  logInfo("mergeWithNeighbour(...)","_stateCopy.isFirstIterationOfBatchOrNoBatch()="<<_stateCopy.isFirstIterationOfBatchOrNoBatch());
//  logInfo("mergeWithNeighbour(...)","OneSolverRequestedLocalRecomputation="<<OneSolverRequestedLocalRecomputation);
  
  if (
      _stateCopy.isFirstIterationOfBatchOrNoBatch() &&
      OneSolverRequestedLocalRecomputation &&
      vertex.hasToCommunicate(level)
  ) {
    dfor2(myDest)
      dfor2(mySrc)
        tarch::la::Vector<DIMENSIONS, int> dest = tarch::la::Vector<DIMENSIONS, int>(1) - myDest;
        tarch::la::Vector<DIMENSIONS, int> src  = tarch::la::Vector<DIMENSIONS, int>(1) - mySrc;

        const int destScalar = TWO_POWER_D - myDestScalar - 1;

        if ( vertex.hasToReceiveMetadata(fromRank,src,dest) ) {
          exahype::MetadataHeap::HeapEntries receivedMetadata;
          receivedMetadata.clear();
          exahype::receiveNeighbourCommunicationMetadata(
              receivedMetadata,fromRank, fineGridX, level);
          assertionEquals(receivedMetadata.size(),
              exahype::NeighbourCommunicationMetadataPerSolver*solvers::RegisteredSolvers.size());

          if( vertex.hasToMergeWithNeighbourData(src,dest) ) { // Only communicate data once per face
            mergeNeighourData(
                fromRank,
                src,dest,
                vertex.getCellDescriptionsIndex()[destScalar],
                fineGridX,level,
                receivedMetadata);
          } else {
            dropNeighbourData(
                fromRank,
                src,dest,
                fineGridX,level,
                receivedMetadata);
          }
        }
      enddforx
    enddforx
  }

  logTraceOut( "mergeWithNeighbour(...)" );
}

void exahype::mappings::PredictionOrLocalRecomputation::dropNeighbourData(
    const int                                    fromRank,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level,
    const exahype::MetadataHeap::HeapEntries& receivedMetadata) {
  for(unsigned int solverNumber = solvers::RegisteredSolvers.size(); solverNumber-- > 0;) {
    auto* solver = solvers::RegisteredSolvers[solverNumber];

    if ( performLocalRecomputation( solver ) ) {
      auto* limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);

      logDebug("dropNeighbourData(...)", "drop data for solver " << solverNumber << " from rank " <<
              fromRank << " at vertex x=" << x << ", level=" << level <<
              ", src=" << src << ", dest=" << dest);

      limitingADERDGSolver->dropNeighbourSolverAndLimiterData(fromRank,src,dest,x,level);
    }
  }
}

void exahype::mappings::PredictionOrLocalRecomputation::mergeNeighourData(
    const int                                    fromRank,
    const tarch::la::Vector<DIMENSIONS,int>&     src,
    const tarch::la::Vector<DIMENSIONS,int>&     dest,
    const int                                    destCellDescriptionIndex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level,
    const exahype::MetadataHeap::HeapEntries&    receivedMetadata) {
  assertion(exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(destCellDescriptionIndex));
  assertion(exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(destCellDescriptionIndex));

  for(unsigned int solverNumber = solvers::RegisteredSolvers.size(); solverNumber-- > 0;) {
    auto* solver = solvers::RegisteredSolvers[solverNumber];

    if ( performLocalRecomputation( solver ) ) {
      const int element = solver->tryGetElement(destCellDescriptionIndex,solverNumber);
      const int offset  = exahype::NeighbourCommunicationMetadataPerSolver*solverNumber;

      if (
          element!=exahype::solvers::Solver::NotFound &&
          receivedMetadata[offset]!=exahype::InvalidMetadataEntry
      ) {
        logDebug("mergeWithNeighbourData(...)", "receive data for solver " << solverNumber << " from rank " <<
                      fromRank << " at vertex x=" << x << ", level=" << level <<
                      ", src=" << src << ", dest=" << dest);

        exahype::MetadataHeap::HeapEntries metadataPortion(
                  receivedMetadata.begin()+offset,
                  receivedMetadata.begin()+offset+exahype::NeighbourCommunicationMetadataPerSolver);

        auto* limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
        limitingADERDGSolver->mergeWithNeighbourDataBasedOnLimiterStatus(
            fromRank,
            destCellDescriptionIndex,element,src,dest,
            true, /* isRecomputation */
            x,level);
      } else {
        logDebug("mergeWithNeighbourData(...)", "drop data for solver " << solverNumber << " from rank " <<
                      fromRank << " at vertex x=" << x << ", level=" << level <<
                      ", src=" << src << ", dest=" << dest);

        auto* limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
        limitingADERDGSolver->dropNeighbourSolverAndLimiterData(fromRank,src,dest,x,level);
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
   vertex.sendToNeighbour(toRank,true,x,level);
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
