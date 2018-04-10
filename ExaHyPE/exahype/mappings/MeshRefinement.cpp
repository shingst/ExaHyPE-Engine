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
 
#include "exahype/mappings/MeshRefinement.h"

#include "peano/utils/Globals.h"
#include "peano/utils/Loop.h"

#include "peano/datatraversal/autotuning/Oracle.h"

#include "peano/grid/aspects/VertexStateAnalysis.h"

#include "tarch/la/VectorScalarOperations.h"

#include "tarch/multicore/Lock.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/VertexOperations.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "exahype/mappings/LimiterStatusSpreading.h"

bool exahype::mappings::MeshRefinement::IsFirstIteration        = true;
bool exahype::mappings::MeshRefinement::IsInitialMeshRefinement = true;

tarch::logging::Log exahype::mappings::MeshRefinement::_log("exahype::mappings::MeshRefinement");

tarch::multicore::BooleanSemaphore exahype::mappings::MeshRefinement::BoundarySemaphore;

void exahype::mappings::MeshRefinement::initialiseLocalVariables(){
  const unsigned int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
  _attainedStableState.resize(numberOfSolvers);
  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    _attainedStableState[solverNumber] = true;
  }

  _verticalExchangeOfSolverDataRequired = false;
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
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}
peano::MappingSpecification
exahype::mappings::MeshRefinement::leaveCellSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces,true);
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
  _localState(masterThread._localState) {
  initialiseLocalVariables();
}


#if defined(SharedMemoryParallelisation)
void exahype::mappings::MeshRefinement::mergeWithWorkerThread(
    const MeshRefinement& workerThread) {
  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
    if (solver->getMeshUpdateRequest()) {
      _attainedStableState[solverNumber] =
          _attainedStableState[solverNumber] && workerThread._attainedStableState[solverNumber];
    }
  }
}
#endif
#endif

void exahype::mappings::MeshRefinement::beginIteration(
  exahype::State& solverState
) {
  _localState = solverState;

  initialiseLocalVariables();

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    if (solver->getMeshUpdateRequest()) {
      solver->zeroTimeStepSizes();
    }
    //assertion2(!solver->getNextMeshUpdateRequest(),solver->toString(),tarch::parallel::Node::getInstance().getRank());
  }

  // background threads
  exahype::solvers::Solver::ensureAllBackgroundJobsHaveTerminated();

  #ifdef Parallel
  if (! MetadataHeap::getInstance().validateThatIncomingJoinBuffersAreEmpty() ) {
      exit(-1);
  }
  #endif
}

void exahype::mappings::MeshRefinement::endIteration(exahype::State& solverState) {
  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    if (solver->getMeshUpdateRequest()) {
      solver->updateNextAttainedStableState(_attainedStableState[solverNumber]);
      solver->setNextAttainedStableState();
    }
  }

  solverState.setVerticalExchangeOfSolverDataRequired(_verticalExchangeOfSolverDataRequired);

  exahype::mappings::MeshRefinement::IsFirstIteration = false;
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

void exahype::mappings::MeshRefinement::eraseIfInsideAndNotRemote(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>&  fineGridH) const {
  if (
      #ifdef Parallel
      !fineGridVertex.isRemote( _localState, true, false)
      &&
      #endif
      !fineGridVertex.isHangingNode()
      &&
      fineGridVertex.isInside()       && // otherwise, we compete with ensureRegularityAlongBoundary
      fineGridVertex.getRefinementControl()
      == Vertex::Records::RefinementControl::Refined
  ) {
    fineGridVertex.erase();
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
      fineGridVertex.evaluateRefinementCriterion(fineGridH);

  if ( refinementControl==exahype::solvers::Solver::RefinementControl::Refine ) {
    refineSafely(fineGridVertex,fineGridH,coarseGridVerticesEnumerator.getLevel()+1,false);
  } else if ( refinementControl==exahype::solvers::Solver::RefinementControl::Erase ) {
    eraseIfInsideAndNotRemote(fineGridVertex,fineGridH);
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
      fineGridVertex.evaluateRefinementCriterion(fineGridH)==
          exahype::solvers::Solver::RefinementControl::Refine
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
      fineGridVertex.evaluateRefinementCriterion(fineGridH)==
          exahype::solvers::Solver::RefinementControl::Refine
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
        !fineGridVertices[fineGridVerticesEnumerator(v)].isInside() &&
        fineGridVertices[fineGridVerticesEnumerator(v)].getRefinementControl()
        ==exahype::Vertex::Records::RefinementControl::Unrefined;
    enddforx

    if (oneInnerVertexIsRefined) {
      dfor2(v)
        tarch::multicore::Lock lock(BoundarySemaphore);
        if (
            fineGridVertices[fineGridVerticesEnumerator(v)].isBoundary() &&
            fineGridVertices[fineGridVerticesEnumerator(v)].getRefinementControl()
            ==exahype::Vertex::Records::RefinementControl::Unrefined
        ) {
          fineGridVertices[fineGridVerticesEnumerator(v)].refine();
        }
        lock.free();
      enddforx
    }

    if (noInnerVertexIsRefined) {
      dfor2(v)
        tarch::multicore::Lock lock(BoundarySemaphore);
        if (
            fineGridVertices[fineGridVerticesEnumerator(v)].isBoundary() &&
            fineGridVertices[fineGridVerticesEnumerator(v)].getRefinementControl()
            ==exahype::Vertex::Records::RefinementControl::Refined
        ) {
          fineGridVertices[fineGridVerticesEnumerator(v)].erase();
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
    if (solver->getMeshUpdateRequest()) {
      const bool newComputeCell =
          solver->progressMeshRefinementInEnterCell(
              fineGridCell,
              fineGridVertices,
              fineGridVerticesEnumerator,
              coarseGridCell,
              coarseGridVertices,
              coarseGridVerticesEnumerator,
              fineGridPositionOfCell,
              exahype::mappings::MeshRefinement::IsInitialMeshRefinement,
              solverNumber);

      // Synchronise time stepping, adjust the solution, evaluate refinement criterion if required
      if (
          exahype::mappings::MeshRefinement::IsFirstIteration ||
          newComputeCell
      ) {
        const int cellDescriptionsIndex = fineGridCell.getCellDescriptionsIndex();
        const int element = solver->tryGetElement(cellDescriptionsIndex,solverNumber);

        if (element!=exahype::solvers::Solver::NotFound) {
          solver->adjustSolutionDuringMeshRefinement(
              cellDescriptionsIndex,element,IsInitialMeshRefinement);
        }
      }
    }
  }


  if (fineGridCell.isInitialised()) {
    exahype::Cell::resetNeighbourMergeFlags(
        fineGridCell.getCellDescriptionsIndex());
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

    if (solver->getMeshUpdateRequest()) {
      const bool newComputeCell =
          solver->progressMeshRefinementInLeaveCell(
              fineGridCell,
              fineGridVertices,
              fineGridVerticesEnumerator,
              coarseGridCell,
              coarseGridVertices,
              coarseGridVerticesEnumerator,
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
              cellDescriptionsIndex,element,IsInitialMeshRefinement);
        }
      }
    }
  }

  // shutdown metadata for empty cells (no cell descriptions)
  if (fineGridCell.isInitialised() &&
      fineGridCell.isEmpty()) {
    fineGridCell.shutdownMetaData();
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
  // TODO(Dominic): introduce some allSolversDoThis allSolversDoThat...  functions
  if ( fineGridCell.isInitialised() ) {
    exahype::solvers::ADERDGSolver::eraseCellDescriptions(fineGridCell.getCellDescriptionsIndex());
    exahype::solvers::FiniteVolumesSolver::eraseCellDescriptions(fineGridCell.getCellDescriptionsIndex());

    exahype::solvers::ADERDGSolver::Heap::getInstance().
        deleteData(fineGridCell.getCellDescriptionsIndex());
    exahype::solvers::FiniteVolumesSolver::Heap::getInstance().
        deleteData(fineGridCell.getCellDescriptionsIndex());
  }
}

#ifdef Parallel
void exahype::mappings::MeshRefinement::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  logTraceInWith6Arguments("mergeWithNeighbour(...)", vertex, neighbour,
                           fromRank, fineGridX, fineGridH, level);

  if ( exahype::mappings::MeshRefinement::IsFirstIteration==false ) {
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

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
    const int cellDescriptionsIndex = fineGridCell.getCellDescriptionsIndex();
    const int element = solver->tryGetElement(cellDescriptionsIndex,solverNumber);
    if ( element!=exahype::solvers::Solver::NotFound ) {
      _verticalExchangeOfSolverDataRequired |=
          solver->prepareMasterCellDescriptionAtMasterWorkerBoundary(
              cellDescriptionsIndex,element);
    }
  }

  fineGridCell.broadcastMetadataToWorkerPerCell(
      worker,
      fineGridVerticesEnumerator.getCellCenter(),
      fineGridVerticesEnumerator.getCellSize(),
      fineGridVerticesEnumerator.getLevel());

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

  receivedCell.receiveMetadataFromMasterPerCell(
      tarch::parallel::NodePool::getInstance().getMasterRank(),
      receivedVerticesEnumerator.getCellCenter(),
      receivedVerticesEnumerator.getCellSize(),
      receivedVerticesEnumerator.getLevel());

  logTraceOut( "receiveDataFromMaster(...)" );
}

void exahype::mappings::MeshRefinement::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  logTraceInWith2Arguments( "mergeWithWorker(...)", localCell.toString(), receivedMasterCell.toString() );

  localCell.mergeWithMetadataFromMasterPerCell(
      cellSize,
      exahype::State::AlgorithmSection::MeshRefinement);

  logTraceOutWith1Argument( "mergeWithWorker(...)", localCell.toString() );
}

void exahype::mappings::MeshRefinement::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  logTraceInWith5Arguments( "prepareCopyToRemoteNode(...)", localCell, toRank, cellCentre, cellSize, level );

  if (
      localCell.isInside() &&
      localCell.hasToCommunicate(cellSize)
  ) {
    if ( localCell.isInitialised() ) {
      const int cellDescriptionsIndex = localCell.getCellDescriptionsIndex();

      exahype::solvers::ADERDGSolver::sendCellDescriptions(toRank,cellDescriptionsIndex,
          peano::heap::MessageType::ForkOrJoinCommunication,cellCentre,level);
      exahype::solvers::FiniteVolumesSolver::sendCellDescriptions(toRank,cellDescriptionsIndex,
          peano::heap::MessageType::ForkOrJoinCommunication,cellCentre,level);

      for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
        auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
        const int element = solver->tryGetElement(cellDescriptionsIndex,solverNumber);

        if( element!=exahype::solvers::Solver::NotFound ) {
          solver->sendDataToWorkerOrMasterDueToForkOrJoin(
              toRank,cellDescriptionsIndex,element,cellCentre,level);
        } else {
          solver->sendEmptyDataToWorkerOrMasterDueToForkOrJoin(toRank,cellCentre,level);
        }
      }

      if ( exahype::State::isJoiningWithMaster() ) {
        exahype::solvers::ADERDGSolver::eraseCellDescriptions(cellDescriptionsIndex);
        exahype::solvers::FiniteVolumesSolver::eraseCellDescriptions(cellDescriptionsIndex);
        localCell.shutdownMetaData();
      }

    } else {
      exahype::solvers::ADERDGSolver::sendEmptyCellDescriptions(toRank,
          peano::heap::MessageType::ForkOrJoinCommunication,cellCentre,level);
      exahype::solvers::FiniteVolumesSolver::sendEmptyCellDescriptions(toRank,
          peano::heap::MessageType::ForkOrJoinCommunication,cellCentre,level);

      for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
        auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

        solver->sendEmptyDataToWorkerOrMasterDueToForkOrJoin(toRank,cellCentre,level);
      }
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

  if (
      localCell.isInside() &&
      localCell.hasToCommunicate(cellSize)
  ) {
    if ( exahype::State::isNewWorkerDueToForkOfExistingDomain() ) {
      localCell.setCellDescriptionsIndex(multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
    }
    if ( !localCell.isInitialised() ) {
      localCell.setupMetaData();
    }

    const int cellDescriptionsIndex = localCell.getCellDescriptionsIndex();

    exahype::solvers::ADERDGSolver::mergeCellDescriptionsWithRemoteData(
        fromRank,localCell,
        peano::heap::MessageType::ForkOrJoinCommunication,
        cellCentre,level);
    exahype::solvers::FiniteVolumesSolver::mergeCellDescriptionsWithRemoteData(
        fromRank,localCell,
        peano::heap::MessageType::ForkOrJoinCommunication,
        cellCentre,level);

    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      const int element = solver->tryGetElement(cellDescriptionsIndex,solverNumber);

      if ( element!=exahype::solvers::Solver::NotFound ) {
        solver->mergeWithWorkerOrMasterDataDueToForkOrJoin(
            fromRank,cellDescriptionsIndex,element,cellCentre,level);
      } else {
        solver->dropWorkerOrMasterDataDueToForkOrJoin(fromRank,cellCentre,level);
      }
    }
  }

  logTraceOut( "mergeWithRemoteDataDueToForkOrJoin(...)" );
}

void exahype::mappings::MeshRefinement::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith2Arguments( "prepareSendToMaster(...)", localCell, verticesEnumerator.toString() );

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
    solver->sendMeshUpdateFlagsToMaster(
        tarch::parallel::NodePool::getInstance().getMasterRank(),
        verticesEnumerator.getCellCenter(),
        verticesEnumerator.getLevel());

    const int cellDescriptionsIndex = localCell.getCellDescriptionsIndex();
    const int element = solver->tryGetElement(cellDescriptionsIndex,solverNumber);
    if ( element!=exahype::solvers::Solver::NotFound ) {
      solver->prepareWorkerCellDescriptionAtMasterWorkerBoundary(cellDescriptionsIndex,element);
    }
  }

  localCell.reduceMetadataToMasterPerCell(
      tarch::parallel::NodePool::getInstance().getMasterRank(),
      verticesEnumerator.getCellCenter(),
      verticesEnumerator.getCellSize(),
      verticesEnumerator.getLevel());

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

  // Merge cell data
  _verticalExchangeOfSolverDataRequired |=
      fineGridCell.mergeWithMetadataFromWorkerPerCell(
      worker,
      fineGridVerticesEnumerator.getCellCenter(),
      fineGridVerticesEnumerator.getCellSize(),
      fineGridVerticesEnumerator.getLevel(),
      exahype::State::AlgorithmSection::MeshRefinement);

  logTraceOut( "mergeWithMaster(...)" );
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
