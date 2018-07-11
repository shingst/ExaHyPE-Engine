/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
 
#include "exahype/mappings/MeshRefinement.h"

#include "peano/utils/Globals.h"
#include "peano/utils/Loop.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/datatraversal/TaskSet.h"

#include "peano/grid/aspects/VertexStateAnalysis.h"

#include "tarch/la/VectorScalarOperations.h"

#include "tarch/multicore/Lock.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/VertexOperations.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "exahype/mappings/RefinementStatusSpreading.h"

#include <sstream>

bool exahype::mappings::MeshRefinement::IsFirstIteration        = true;
bool exahype::mappings::MeshRefinement::IsInitialMeshRefinement = true;

bool exahype::mappings::MeshRefinement::StillInRefiningMode     = true;

tarch::logging::Log exahype::mappings::MeshRefinement::_log("exahype::mappings::MeshRefinement");

tarch::multicore::BooleanSemaphore exahype::mappings::MeshRefinement::BoundarySemaphore;

void exahype::mappings::MeshRefinement::initialiseLocalVariables(){
  const unsigned int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
  _attainedStableState.resize(numberOfSolvers);
  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    _attainedStableState[solverNumber] = !IsFirstIteration;
  }

  _verticalExchangeOfSolverDataRequired = false;
}

bool exahype::mappings::MeshRefinement::allSolversAttainedStableState() const{
  bool allSolversAttainedStableState =
      #if defined(Parallel)
      !exahype::State::isInvolvedInJoinOrFork() &&
      #endif
      true;

  for (bool isStable : _attainedStableState) {
    allSolversAttainedStableState &= isStable;
  }
  return allSolversAttainedStableState;
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::MeshRefinement::communicationSpecification() const {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
          SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
          SendDataAndStateAfterLastTouchVertexLastTime,
      true);
}

peano::MappingSpecification
exahype::mappings::MeshRefinement::touchVertexFirstTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}

peano::MappingSpecification
exahype::mappings::MeshRefinement::touchVertexLastTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}

peano::MappingSpecification
exahype::mappings::MeshRefinement::enterCellSpecification(int level) const {
  if ( exahype::solvers::Solver::SpawnAMRBackgroundJobs ) {
    return peano::MappingSpecification(
        peano::MappingSpecification::WholeTree,
        peano::MappingSpecification::Serial,true);
  } else {
    return peano::MappingSpecification(
        peano::MappingSpecification::WholeTree,
        peano::MappingSpecification::AvoidFineGridRaces,true);
  }
}

peano::MappingSpecification
exahype::mappings::MeshRefinement::leaveCellSpecification(int level) const {
  if ( exahype::solvers::Solver::SpawnAMRBackgroundJobs ) {
    return peano::MappingSpecification(
        peano::MappingSpecification::WholeTree,
        peano::MappingSpecification::Serial,true);
  } else {
    return peano::MappingSpecification(
        peano::MappingSpecification::WholeTree,
        peano::MappingSpecification::AvoidFineGridRaces,true);
  }
}

// Nop.
peano::MappingSpecification
exahype::mappings::MeshRefinement::ascendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,false);
}
peano::MappingSpecification
exahype::mappings::MeshRefinement::descendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,false);
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::MeshRefinement::MeshRefinement(const MeshRefinement& masterThread):
  _stableIterationsInARow(masterThread._stableIterationsInARow),
  _iterationsSinceLastErasing(masterThread._iterationsSinceLastErasing),
  _localState(masterThread._localState)
{
  initialiseLocalVariables();
}
#endif

#if defined(SharedMemoryParallelisation)
void exahype::mappings::MeshRefinement::mergeWithWorkerThread(
    const MeshRefinement& workerThread) {
  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
    if (solver->hasRequestedMeshRefinement()) {
      _attainedStableState[solverNumber] =
          _attainedStableState[solverNumber] && workerThread._attainedStableState[solverNumber];
    }
  }

  _iterationsSinceLastErasing = std::min( _iterationsSinceLastErasing, workerThread._iterationsSinceLastErasing );
}
#endif

void exahype::mappings::MeshRefinement::beginIteration(
  exahype::State& solverState
) {
  _localState = solverState;

  initialiseLocalVariables();

  if ( exahype::mappings::MeshRefinement::IsFirstIteration ) {
    _stableIterationsInARow     = 0;
    _iterationsSinceLastErasing = (MeshRefinement::IsInitialMeshRefinement) ? 1 : 0;
    StillInRefiningMode = true;
  }
  if (
      // as soon as we have two iterations in a row where
      // nothing has happened we switch to refining mode
      StillInRefiningMode &&
      _stableIterationsInARow>0
  ) {
    StillInRefiningMode = false;
    _iterationsSinceLastErasing = 0;
  }
  if ( StillInRefiningMode ) { // we can use cheap neighbour lookup as long as we are in refining mode
    _iterationsSinceLastErasing = 2;
  }

  _iterationsSinceLastErasing++;

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
    if (solver->hasRequestedMeshRefinement()) {
      solver->zeroTimeStepSizes();
    }
    //assertion2(!solver->getNextMeshUpdateRequest(),solver->toString(),tarch::parallel::Node::getInstance().getRank());
  }

  #ifdef Parallel
  if (! MetadataHeap::getInstance().validateThatIncomingJoinBuffersAreEmpty() ) {
      exit(-1);
  }
  #endif
}

void exahype::mappings::MeshRefinement::endIteration(exahype::State& solverState) {
  logTraceInWith1Argument("endIteration(State)", solverState);

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
    if (solver->hasRequestedMeshRefinement()) {
      if ( allSolversAttainedStableState() ) {
        _stableIterationsInARow++;
      } else {
        _stableIterationsInARow = 0;
      }

      solver->updateNextAttainedStableState( _stableIterationsInARow > 4 ); // Found experimentally
      solver->setNextAttainedStableState();
    }
  }

  solverState.setVerticalExchangeOfSolverDataRequired(_verticalExchangeOfSolverDataRequired);

  exahype::mappings::MeshRefinement::IsFirstIteration=false;

  peano::datatraversal::TaskSet::startToProcessBackgroundJobs();

  // background threads
  exahype::solvers::Solver::ensureAllJobsHaveTerminated(exahype::solvers::Solver::JobType::AMRJob);

  logTraceOutWith1Argument("endIteration(State)", solverState);
}

void exahype::mappings::MeshRefinement::refineSafely(
    exahype::Vertex&                              fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>&  fineGridH,
    int                                           fineGridLevel,
    bool                                          isCalledByCreationalEvent) const {
  if ( fineGridVertex.getRefinementControl()==Vertex::Records::Unrefined ) {
    switch ( _localState.mayRefine(isCalledByCreationalEvent,fineGridLevel) ) {
    case State::RefinementAnswer::DontRefineYet:
      break;
    case State::RefinementAnswer::Refine:
      fineGridVertex.refine();
      break;
    case State::RefinementAnswer::EnforceRefinement:
      fineGridVertex.enforceRefine();
      break;
    }
  }
}

void exahype::mappings::MeshRefinement::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  exahype::solvers::Solver::RefinementControl refinementControl =
      fineGridVertex.evaluateRefinementCriterion(
      fineGridX,
      coarseGridVerticesEnumerator.getLevel()+1,
      fineGridH,
      _iterationsSinceLastErasing < 1);

  if (
      _stableIterationsInARow <= 3 // Found experimentally
      &&
      refinementControl==exahype::solvers::Solver::RefinementControl::Refine
  ) {
    refineSafely(fineGridVertex,fineGridH,coarseGridVerticesEnumerator.getLevel()+1,false);
  } else if (
      refinementControl==exahype::solvers::Solver::RefinementControl::Erase
      &&
      !fineGridVertex.isHangingNode()
      &&
      fineGridVertex.isInside()
      && // otherwise, we compete with ensureRegularityAlongBoundary
      fineGridVertex.getRefinementControl()==
          Vertex::Records::RefinementControl::Refined
      &&
      _stableIterationsInARow > 3 // Found experimentally
  ) {
  // TODO  fineGridVertex.erase(); // TODO(Dominic): vertex erasing is not well understood yet
  }
}


void exahype::mappings::MeshRefinement::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  logTraceInWith6Arguments("createBoundaryVertex(...)", fineGridVertex,
                           fineGridX, fineGridH,
                           coarseGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfVertex);

  if (
      fineGridVertex.evaluateRefinementCriterion(
          fineGridX,
          coarseGridVerticesEnumerator.getLevel()+1,
          fineGridH,
          _iterationsSinceLastErasing < 1)
      == exahype::solvers::Solver::RefinementControl::Refine
  ) {
    refineSafely(fineGridVertex,fineGridH,coarseGridVerticesEnumerator.getLevel()+1,true);
  }

  logTraceOutWith1Argument("createBoundaryVertex(...)", fineGridVertex);
}


void exahype::mappings::MeshRefinement::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  logTraceInWith6Arguments("createInnerVertex(...)", fineGridVertex, fineGridX,
                           fineGridH, coarseGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfVertex);
  if (
      fineGridVertex.evaluateRefinementCriterion(
          fineGridX,
          coarseGridVerticesEnumerator.getLevel()+1,
          fineGridH,
          _iterationsSinceLastErasing < 1)
      == exahype::solvers::Solver::RefinementControl::Refine
  ) {
    refineSafely(fineGridVertex,fineGridH,coarseGridVerticesEnumerator.getLevel()+1,true);
  }

  logTraceOutWith1Argument("createInnerVertex(...)", fineGridVertex);
}


void exahype::mappings::MeshRefinement::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("createCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  // do nothing
  fineGridCell.setCellDescriptionsIndex(
      multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);

  logTraceOutWith1Argument("createCell(...)", fineGridCell);
}

void exahype::mappings::MeshRefinement::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MeshRefinement::touchVertexFirstTime(
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
      exahype::State::AlgorithmSection::MeshRefinement,fineGridX,fineGridH);

  logTraceOutWith1Argument("touchVertexFirstTime(...)", fineGridVertex);
}

void exahype::mappings::MeshRefinement::ensureRegularityAlongBoundary(
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const {
  if (
      _iterationsSinceLastErasing < 1
      &&
      peano::grid::aspects::VertexStateAnalysis::isOneVertexBoundary(
          fineGridVertices,fineGridVerticesEnumerator)
  ) {
    bool oneInnerVertexIsRefined = false;
    bool noInnerVertexIsRefined  = true;
    dfor2(v)
      oneInnerVertexIsRefined |=
          fineGridVertices[fineGridVerticesEnumerator(v)].isInside() &&
          fineGridVertices[fineGridVerticesEnumerator(v)].getRefinementControl()
          ==exahype::Vertex::Records::RefinementControl::Refined;
      noInnerVertexIsRefined &=
          !fineGridVertices[fineGridVerticesEnumerator(v)].isInside() ||
          fineGridVertices[fineGridVerticesEnumerator(v)].getRefinementControl()
          ==exahype::Vertex::Records::RefinementControl::Unrefined;
    enddforx

    if (oneInnerVertexIsRefined) {
      dfor2(v)
        tarch::multicore::Lock lock(BoundarySemaphore);
        if (
            _stableIterationsInARow <= 3 // Found experimentally
            &&
            fineGridVertices[fineGridVerticesEnumerator(v)].isBoundary()
            &&
            fineGridVertices[fineGridVerticesEnumerator(v)].getRefinementControl()==
                exahype::Vertex::Records::RefinementControl::Unrefined
            &&
            fineGridVertices[fineGridVerticesEnumerator(v)].evaluateRefinementCriterion(
                fineGridVerticesEnumerator.getVertexPosition(vScalar),
                fineGridVerticesEnumerator.getLevel(),
                fineGridVerticesEnumerator.getCellSize(),
                false)
            !=exahype::solvers::Solver::RefinementControl::Erase
        ) {
          fineGridVertices[fineGridVerticesEnumerator(v)].refine();
        }
        lock.free();
      enddforx
    } else if (
        noInnerVertexIsRefined
    ) {
      dfor2(v)
        tarch::multicore::Lock lock(BoundarySemaphore);
        if (
            _stableIterationsInARow > 3 // Found experimentally
            &&
            fineGridVertices[fineGridVerticesEnumerator(v)].isBoundary()
            &&
            fineGridVertices[fineGridVerticesEnumerator(v)].getRefinementControl()==
                exahype::Vertex::Records::RefinementControl::Refined
            &&
            fineGridVertices[fineGridVerticesEnumerator(v)].evaluateRefinementCriterion(
              fineGridVerticesEnumerator.getVertexPosition(vScalar),
              fineGridVerticesEnumerator.getLevel(),
              fineGridVerticesEnumerator.getCellSize(),
              false)
            ==exahype::solvers::Solver::RefinementControl::Erase

        ) {
          // TODO
          // fineGridVertices[fineGridVerticesEnumerator(v)].erase(); // TODO(Dominic): vertex erasing is not well understood yet
        }
        lock.free();
      enddforx
    }
  }
}

void exahype::mappings::MeshRefinement::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  assertion(fineGridCell.isInside());

  for (unsigned int solverNumber=0; solverNumber<exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
    if ( solver->hasRequestedMeshRefinement() ) {
      const bool newComputeCell =
          solver->progressMeshRefinementInEnterCell(
              fineGridCell,
              fineGridVertices,
              fineGridVerticesEnumerator,
              coarseGridCell,
              coarseGridVerticesEnumerator,
              solverNumber,
              StillInRefiningMode);

      // Synchronise time stepping, adjust the solution, evaluate refinement criterion if required
      if (
          (
          #ifdef Parallel
          !exahype::State::isNewWorkerDueToForkOfExistingDomain() &&
          #endif
          exahype::mappings::MeshRefinement::IsFirstIteration)     // It has to be the first overall iteration
          ||
          newComputeCell
      ) {
        const int cellDescriptionsIndex = fineGridCell.getCellDescriptionsIndex();
        const int element = solver->tryGetElement(cellDescriptionsIndex,solverNumber);

        if (element!=exahype::solvers::Solver::NotFound) {
          solver->adjustSolutionDuringMeshRefinement(
              cellDescriptionsIndex,element);
        }
      }
    }
  }

  if ( fineGridCell.isInitialised() ) {
    exahype::Cell::resetNeighbourMergeFlags(
        fineGridCell.getCellDescriptionsIndex(),
        fineGridVertices,fineGridVerticesEnumerator);

    // shutdown metadata for empty cells (no cell descriptions)
    if ( fineGridCell.isEmpty() ) {
      fineGridCell.shutdownMetaDataAndResetCellDescriptionsIndex();
      _iterationsSinceLastErasing = 0;
    }
  }

  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::MeshRefinement::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("leaveCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  for (unsigned int solverNumber=0; solverNumber<exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    if (solver->hasRequestedMeshRefinement()) {
      const bool newComputeCell =
          solver->progressMeshRefinementInLeaveCell(
              fineGridCell,
              fineGridVertices,
              fineGridVerticesEnumerator,
              coarseGridCell,
              fineGridPositionOfCell,
              solverNumber);

      const bool isStable =
          solver->attainedStableState(
              fineGridCell,
              fineGridVertices,
              fineGridVerticesEnumerator,
              solverNumber);

      _attainedStableState[solverNumber] =
          _attainedStableState[solverNumber] && isStable;

      // Synchronise time stepping, adjust the solution, evaluate refinement criterion if required
      if ( newComputeCell ) {
        const int cellDescriptionsIndex = fineGridCell.getCellDescriptionsIndex();
        const int element = solver->tryGetElement(cellDescriptionsIndex,solverNumber);

        if (element!=exahype::solvers::Solver::NotFound) {
          solver->adjustSolutionDuringMeshRefinement(
              cellDescriptionsIndex,element);
        }
      }
    }
  }

  ensureRegularityAlongBoundary(fineGridVertices,fineGridVerticesEnumerator);

  logTraceOutWith1Argument("leaveCell(...)", fineGridCell);
}

void exahype::mappings::MeshRefinement::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  if ( fineGridCell.isInitialised() ) {
    fineGridCell.shutdownMetaData();
    _iterationsSinceLastErasing = 0;
  }
}

#ifdef Parallel
void exahype::mappings::MeshRefinement::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  logTraceInWith6Arguments("mergeWithNeighbour(...)", vertex, neighbour,
                           fromRank, fineGridX, fineGridH, level);

  if ( !exahype::mappings::MeshRefinement::IsFirstIteration ) {
    vertex.mergeOnlyWithNeighbourMetadata(
        fromRank,fineGridX,fineGridH,level,
        exahype::State::AlgorithmSection::MeshRefinement);
  }

  logTraceOut("mergeWithNeighbour(...)");
}

void exahype::mappings::MeshRefinement::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  logTraceInWith5Arguments("prepareSendToNeighbour(...)", vertex,
                           toRank, x, h, level);

  vertex.sendOnlyMetadataToNeighbour(toRank,x,h,level);

  logTraceOut("prepareSendToNeighbour(...)");
}

bool exahype::mappings::MeshRefinement::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  logTraceIn( "prepareSendToWorker(...)" );

  if (
      !exahype::State::isForkingRank(worker) &&
      fineGridCell.hasToCommunicate(fineGridVerticesEnumerator.getCellSize())
  ) {
    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      solver->progressMeshRefinementInPrepareSendToWorker(
          worker, fineGridCell, fineGridVertices,fineGridVerticesEnumerator,
          coarseGridCell, coarseGridVerticesEnumerator,
          solverNumber);
    }

    // send out the cell descriptions
    const int cellDescriptionsIndex = fineGridCell.getCellDescriptionsIndex();
    _verticalExchangeOfSolverDataRequired |=
        exahype::solvers::ADERDGSolver::sendCellDescriptions(worker,cellDescriptionsIndex,
            false/* !(send data from worker side) */,
            peano::heap::MessageType::MasterWorkerCommunication,
            fineGridVerticesEnumerator.getCellCenter(),fineGridVerticesEnumerator.getLevel());
    exahype::solvers::FiniteVolumesSolver::sendCellDescriptions(worker,cellDescriptionsIndex,
        peano::heap::MessageType::MasterWorkerCommunication,
        fineGridVerticesEnumerator.getCellCenter(),fineGridVerticesEnumerator.getLevel());

    // possibly send out data
    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      const int element = solver->tryGetElement(cellDescriptionsIndex,solverNumber);
      if ( element!=exahype::solvers::Solver::NotFound ) {
        solver->sendDataToWorkerIfProlongating(
            worker, 
            cellDescriptionsIndex,element,
            fineGridVerticesEnumerator.getCellCenter(),
            fineGridVerticesEnumerator.getLevel());
      }
    }
  }

  logTraceOutWith1Argument( "prepareSendToWorker(...)", true );
  return true;
}

// TODO(Dominic): Add to docu: see documentation in peano/pdt/stdtemplates/MappingHeader.template
// on function receiveDataFromMaster
void exahype::mappings::MeshRefinement::receiveDataFromMaster(
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

  receivedCell.setCellDescriptionsIndex(
      multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
  if (
      !exahype::State::isNewWorkerDueToForkOfExistingDomain() &&
      receivedCell.hasToCommunicate( receivedVerticesEnumerator.getCellSize())
  ) {
    receivedCell.setupMetaData();
    exahype::solvers::ADERDGSolver::receiveCellDescriptions(
        tarch::parallel::NodePool::getInstance().getMasterRank(),
        receivedCell,
        peano::heap::MessageType::MasterWorkerCommunication,
        receivedVerticesEnumerator.getCellCenter(),
        receivedVerticesEnumerator.getLevel());
    exahype::solvers::FiniteVolumesSolver::receiveCellDescriptions(
        tarch::parallel::NodePool::getInstance().getMasterRank(),
        receivedCell, // Two step approach
        peano::heap::MessageType::MasterWorkerCommunication,
        receivedVerticesEnumerator.getCellCenter(),
        receivedVerticesEnumerator.getLevel());

    const int receivedCellDescriptionsIndex = receivedCell.getCellDescriptionsIndex();
    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      const int receivedElement = solver->tryGetElement(receivedCellDescriptionsIndex,solverNumber);
      if ( receivedElement!=exahype::solvers::Solver::NotFound ) {
        solver->receiveDataFromMasterIfProlongating(
            tarch::parallel::NodePool::getInstance().getMasterRank(),
            receivedCellDescriptionsIndex,receivedElement,
            receivedVerticesEnumerator.getCellCenter(),
            receivedVerticesEnumerator.getLevel());
      }
    }
  }

  logTraceOut( "receiveDataFromMaster(...)" );
}

void exahype::mappings::MeshRefinement::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  logTraceInWith2Arguments( "mergeWithWorker(...)", localCell.toString(), receivedMasterCell.toString() );

  if (  receivedMasterCell.isInitialised() ) { // we do not receive anything here
    // Do not merge anything if our worker is on a newly forked part of the mesh
    if ( !exahype::State::isNewWorkerDueToForkOfExistingDomain() ) {
      if ( !localCell.isInitialised() ) { // simply copy the index
        localCell.setupMetaData();
      }

      const int receivedCellDescriptionsIndex = receivedMasterCell.getCellDescriptionsIndex();
      for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
        auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
        const int receivedElement = solver->tryGetElement(receivedCellDescriptionsIndex,solverNumber);
        if ( receivedElement!=exahype::solvers::Solver::NotFound  ) {
          solver->progressMeshRefinementInMergeWithWorker(
              localCell.getCellDescriptionsIndex(),
              receivedCellDescriptionsIndex,receivedElement);
        }
      }
      if ( localCell.isInitialised() && localCell.isEmpty() ) {
        localCell.shutdownMetaDataAndResetCellDescriptionsIndex();
      }
    }
    receivedMasterCell.shutdownMetaData();
  }

  logTraceOutWith1Argument( "mergeWithWorker(...)", localCell.toString() );
}



void exahype::mappings::MeshRefinement::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith2Arguments( "prepareSendToMaster(...)", localCell, verticesEnumerator.toString() );

  // global reductions
  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
    solver->sendMeshUpdateFlagsToMaster(
        tarch::parallel::NodePool::getInstance().getMasterRank(),
        verticesEnumerator.getCellCenter(),
        verticesEnumerator.getLevel());
  }

  if ( localCell.hasToCommunicate( verticesEnumerator.getCellSize()) ) {
    exahype::solvers::ADERDGSolver::sendCellDescriptions(
        tarch::parallel::NodePool::getInstance().getMasterRank(),
        localCell.getCellDescriptionsIndex(),true/* send data from worker side*/,
        peano::heap::MessageType::MasterWorkerCommunication,
        verticesEnumerator.getCellCenter(),
        verticesEnumerator.getLevel());
    exahype::solvers::FiniteVolumesSolver::sendCellDescriptions(
        tarch::parallel::NodePool::getInstance().getMasterRank(),
        localCell.getCellDescriptionsIndex(),
        peano::heap::MessageType::MasterWorkerCommunication,
        verticesEnumerator.getCellCenter(),
        verticesEnumerator.getLevel()); // make collective

    for (unsigned int solverNumber = 0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      const int element = solver->tryGetElement(localCell.getCellDescriptionsIndex(),solverNumber);
      if ( element!=exahype::solvers::Solver::NotFound ) {
        solver->progressMeshRefinementInPrepareSendToMaster(
            tarch::parallel::NodePool::getInstance().getMasterRank(),
            localCell.getCellDescriptionsIndex(),element,
            verticesEnumerator.getCellCenter(),verticesEnumerator.getLevel());
      }
    }
  }
  
  logTraceOut( "prepareSendToMaster(...)" );
}

void exahype::mappings::MeshRefinement::mergeWithMaster(
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
  masterState.merge(workerState);

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
    solver->mergeWithWorkerMeshUpdateFlags(
        worker,
        fineGridVerticesEnumerator.getCellCenter(),
        fineGridVerticesEnumerator.getLevel());
  }

  if ( fineGridCell.hasToCommunicate( fineGridVerticesEnumerator.getCellSize()) ) {
    if ( fineGridCell.isInitialised() ) {
      exahype::solvers::ADERDGSolver::eraseCellDescriptions(fineGridCell.getCellDescriptionsIndex());
      exahype::solvers::FiniteVolumesSolver::eraseCellDescriptions(fineGridCell.getCellDescriptionsIndex());
    } else {
      fineGridCell.setupMetaData();
    }

    exahype::solvers::ADERDGSolver::receiveCellDescriptions(
        worker,fineGridCell,
        peano::heap::MessageType::MasterWorkerCommunication,
        fineGridVerticesEnumerator.getCellCenter(),
        fineGridVerticesEnumerator.getLevel());
    exahype::solvers::FiniteVolumesSolver::receiveCellDescriptions(
        worker,fineGridCell,
        peano::heap::MessageType::MasterWorkerCommunication,
        fineGridVerticesEnumerator.getCellCenter(),
        fineGridVerticesEnumerator.getLevel());

    const int cellDescriptionsIndex = fineGridCell.getCellDescriptionsIndex();
    for (unsigned int solverNumber = 0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      const int element = solver->tryGetElement(cellDescriptionsIndex,solverNumber);
      if ( element!=exahype::solvers::Solver::NotFound ) {
        _verticalExchangeOfSolverDataRequired |=
            solver->progressMeshRefinementInMergeWithMaster(
                worker, cellDescriptionsIndex, element,
                coarseGridCell.getCellDescriptionsIndex(),
                fineGridVerticesEnumerator.getCellCenter(),
                fineGridVerticesEnumerator.getLevel(),
                StillInRefiningMode);
      }
    }

    if ( fineGridCell.isEmpty() ) {
      fineGridCell.shutdownMetaDataAndResetCellDescriptionsIndex();  
    }
  }

  logTraceOut( "mergeWithMaster(...)" );
}


//////////////////
// FORK and JOIN//
//////////////////

void exahype::mappings::MeshRefinement::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  logTraceInWith5Arguments( "prepareCopyToRemoteNode(...)", localCell, toRank, cellCentre, cellSize, level );

  if ( localCell.hasToCommunicate(cellSize) ) {
    const int cellDescriptionsIndex = localCell.getCellDescriptionsIndex();
    exahype::solvers::ADERDGSolver::sendCellDescriptions(toRank,cellDescriptionsIndex,
        exahype::State::isJoiningWithMaster()/* send out data from worker side */,
        peano::heap::MessageType::ForkOrJoinCommunication,cellCentre,level);
    exahype::solvers::FiniteVolumesSolver::sendCellDescriptions(toRank,cellDescriptionsIndex,
        peano::heap::MessageType::ForkOrJoinCommunication,cellCentre,level);

    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      const int element = solver->tryGetElement(cellDescriptionsIndex,solverNumber);
      if( element!=exahype::solvers::Solver::NotFound ) {
        solver->sendDataToWorkerOrMasterDueToForkOrJoin(toRank,cellDescriptionsIndex,element,
            peano::heap::MessageType::ForkOrJoinCommunication,cellCentre,level);
      }
    }

    if ( localCell.isInitialised() && localCell.getRankOfRemoteNode()==toRank ) { // isAsignedToRemoteRank does not work, remeber the halo sends
      localCell.shutdownMetaDataAndResetCellDescriptionsIndex();
      _iterationsSinceLastErasing=0;
    } 
  }

  logTraceOut( "prepareCopyToRemoteNode(...)" );
}

void exahype::mappings::MeshRefinement::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  logTraceInWith6Arguments( "mergeWithRemoteDataDueToForkOrJoin(...)", localVertex, masterOrWorkerVertex, fromRank, x, h, level );

  if ( exahype::State::isNewWorkerDueToForkOfExistingDomain() ) {
    exahype::VertexOperations::writeCellDescriptionsIndex(
        localVertex,multiscalelinkedcell::HangingVertexBookkeeper::getInstance().createVertexLinkMapForNewVertex());
  }

  logTraceOut( "mergeWithRemoteDataDueToForkOrJoin(...)" );
}

void exahype::mappings::MeshRefinement::mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  logTraceInWith3Arguments( "mergeWithRemoteDataDueToForkOrJoin(...)", localCell, masterOrWorkerCell, fromRank );

  if ( localCell.hasToCommunicate(cellSize) ) {
    if ( exahype::State::isNewWorkerDueToForkOfExistingDomain() ) {
      localCell.setCellDescriptionsIndex(
          multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
      localCell.setupMetaData();
    } else if ( exahype::State::isJoiningWithWorker() ) {
      exahype::solvers::ADERDGSolver::eraseCellDescriptions(localCell.getCellDescriptionsIndex());
      exahype::solvers::FiniteVolumesSolver::eraseCellDescriptions(localCell.getCellDescriptionsIndex());
    }

    exahype::solvers::ADERDGSolver::receiveCellDescriptions(
        fromRank,localCell,
        peano::heap::MessageType::ForkOrJoinCommunication,
        cellCentre,level);
    exahype::solvers::FiniteVolumesSolver::receiveCellDescriptions(
        fromRank,localCell,
        peano::heap::MessageType::ForkOrJoinCommunication,
        cellCentre,level);

    // receive accompanying data
    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      const int element = solver->tryGetElement(localCell.getCellDescriptionsIndex(),solverNumber);
      if ( element!=exahype::solvers::Solver::NotFound ) {
        solver->mergeWithWorkerOrMasterDataDueToForkOrJoin(fromRank,
            localCell.getCellDescriptionsIndex(),element,
            peano::heap::MessageType::ForkOrJoinCommunication,cellCentre,level);
      }
    }

    // if no solver was found or the cell does belong to the master,
    // shut down the metadata again
    if (
        localCell.isInitialised() &&
        (
          localCell.isEmpty() ||
          localCell.getRankOfRemoteNode()==fromRank
        )
    ) {
      localCell.shutdownMetaDataAndResetCellDescriptionsIndex();
    }
  }

  logTraceOut( "mergeWithRemoteDataDueToForkOrJoin(...)" );
}

//
// All methods below are nop,
//
// ==================================



void exahype::mappings::MeshRefinement::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::MeshRefinement::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

exahype::mappings::MeshRefinement::MeshRefinement() {
  // do nothing
}

exahype::mappings::MeshRefinement::~MeshRefinement() {
  // do nothing
}

void exahype::mappings::MeshRefinement::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MeshRefinement::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MeshRefinement::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::MeshRefinement::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
