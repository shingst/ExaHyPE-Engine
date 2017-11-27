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
 
#include "exahype/mappings/TimeStepSizeComputation.h"

#include "peano/utils/Globals.h"

#include "tarch/multicore/Loop.h"
#include "tarch/parallel/Node.h"

#include "peano/datatraversal/autotuning/Oracle.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include <limits>
#include <iomanip>

tarch::logging::Log exahype::mappings::TimeStepSizeComputation::_log(
    "exahype::mappings::TimeStepSizeComputation");

void exahype::mappings::TimeStepSizeComputation::prepareLocalTimeStepVariables(){
  const unsigned int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
  _minTimeStepSizes.resize(numberOfSolvers);
  _minCellSizes.resize(numberOfSolvers);
  _maxCellSizes.resize(numberOfSolvers);

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    _minTimeStepSizes[solverNumber] = std::numeric_limits<double>::max();
    _minCellSizes    [solverNumber] = std::numeric_limits<double>::max();
    _maxCellSizes    [solverNumber] = -std::numeric_limits<double>::max(); // "-", min
  }
}

peano::CommunicationSpecification
exahype::mappings::TimeStepSizeComputation::communicationSpecification() const {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::SendDataAndStateBeforeDescendIntoLocalSubtree,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::MaskOutWorkerMasterDataAndStateExchange,
      true);
}

peano::MappingSpecification
exahype::mappings::TimeStepSizeComputation::enterCellSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}
/**
 * Nop.
 */
peano::MappingSpecification exahype::mappings::TimeStepSizeComputation::
    touchVertexLastTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}
peano::MappingSpecification exahype::mappings::TimeStepSizeComputation::
    touchVertexFirstTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}
peano::MappingSpecification
exahype::mappings::TimeStepSizeComputation::leaveCellSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}
peano::MappingSpecification
exahype::mappings::TimeStepSizeComputation::ascendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}
peano::MappingSpecification
exahype::mappings::TimeStepSizeComputation::descendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}

exahype::mappings::TimeStepSizeComputation::~TimeStepSizeComputation() {
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::TimeStepSizeComputation::TimeStepSizeComputation(const TimeStepSizeComputation& masterThread)
  : _localState(masterThread._localState) {
  prepareLocalTimeStepVariables();
}

// Merge over threads
void exahype::mappings::TimeStepSizeComputation::mergeWithWorkerThread(
    const TimeStepSizeComputation& workerThread) {
  for (int i = 0; i < static_cast<int>(exahype::solvers::RegisteredSolvers.size()); i++) {
    _minTimeStepSizes[i] =
        std::min(_minTimeStepSizes[i], workerThread._minTimeStepSizes[i]);
    _minCellSizes[i] =
        std::min(_minCellSizes[i], workerThread._minCellSizes[i]);
    _maxCellSizes[i] =
        std::max(_maxCellSizes[i], workerThread._maxCellSizes[i]);
  }
}
#endif

void exahype::mappings::TimeStepSizeComputation::beginIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);

  _localState = solverState;

  prepareLocalTimeStepVariables();

  logTraceOutWith1Argument("beginIteration(State)", solverState);
}

void exahype::mappings::TimeStepSizeComputation::weighMinNextPredictorTimeStepSize(
    exahype::solvers::Solver* solver) {
  exahype::solvers::ADERDGSolver* aderdgSolver = nullptr;

  switch(solver->getType()) {
    case exahype::solvers::Solver::Type::ADERDG:
      aderdgSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);
      break;
    case exahype::solvers::Solver::Type::LimitingADERDG:
      aderdgSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver().get();
      break;
    case exahype::solvers::Solver::Type::FiniteVolumes:
      break;
  }

  if (aderdgSolver!=nullptr) {
    const double stableTimeStepSize = aderdgSolver->getMinNextPredictorTimeStepSize();

    const double timeStepSizeWeight = exahype::State::getTimeStepSizeWeightForPredictionRerun();
    aderdgSolver->updateMinNextPredictorTimeStepSize(
        timeStepSizeWeight * stableTimeStepSize);
    aderdgSolver->setMinPredictorTimeStepSize(
        timeStepSizeWeight * stableTimeStepSize); // This will be propagated to the corrector
  }
}


void exahype::mappings::TimeStepSizeComputation::reinitialiseTimeStepDataIfLastPredictorTimeStepSizeWasInstable(
    exahype::solvers::Solver* solver) {
  exahype::solvers::ADERDGSolver* aderdgSolver = nullptr;

  switch(solver->getType()) {
    case exahype::solvers::Solver::Type::ADERDG:
      aderdgSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);
      break;
    case exahype::solvers::Solver::Type::LimitingADERDG:
      aderdgSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver().get();
      break;
    case exahype::solvers::Solver::Type::FiniteVolumes:
      break;
  }

  if (aderdgSolver!=nullptr) {
    const double stableTimeStepSize = aderdgSolver->getMinNextPredictorTimeStepSize();
    double usedTimeStepSize         = aderdgSolver->getMinPredictorTimeStepSize();

    if (tarch::la::equals(usedTimeStepSize,0.0)) {
      usedTimeStepSize = stableTimeStepSize; // TODO(Dominic): Still necessary?
    }

    bool usedTimeStepSizeWasInstable = usedTimeStepSize > stableTimeStepSize;
    aderdgSolver->setStabilityConditionWasViolated(usedTimeStepSizeWasInstable);

    const double timeStepSizeWeight = exahype::State::getTimeStepSizeWeightForPredictionRerun();
    if (usedTimeStepSizeWasInstable) {
      aderdgSolver->updateMinNextPredictorTimeStepSize(
          timeStepSizeWeight * stableTimeStepSize);
      aderdgSolver->setMinPredictorTimeStepSize(
          timeStepSizeWeight * stableTimeStepSize); // This will be propagated to the corrector
    } else {
      aderdgSolver->updateMinNextPredictorTimeStepSize(
          0.5 * (usedTimeStepSize + timeStepSizeWeight * stableTimeStepSize));
    }
  }
}

void exahype::mappings::TimeStepSizeComputation::endIteration(
    exahype::State& state) {
  logTraceInWith1Argument("endIteration(State)", state);

  for (unsigned int solverNumber = 0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    if (solver->isComputingTimeStepSize(_localState.getAlgorithmSection())) {
      // cell sizes
      solver->updateNextMinCellSize(_minCellSizes[solverNumber]);
      solver->updateNextMaxCellSize(_maxCellSizes[solverNumber]);
      if (tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
        assertion4(solver->getNextMinCellSize()<std::numeric_limits<double>::max(),
            solver->getNextMinCellSize(),_minCellSizes[solverNumber],solver->toString(),
            exahype::records::State::toString(_localState.getAlgorithmSection()));
        assertion4(solver->getNextMaxCellSize()>0,
            solver->getNextMaxCellSize(),_maxCellSizes[solverNumber],solver->toString(),
            exahype::records::State::toString(_localState.getAlgorithmSection()));
      }

      // time
      assertion1(std::isfinite(_minTimeStepSizes[solverNumber]),_minTimeStepSizes[solverNumber]);
      assertion1(_minTimeStepSizes[solverNumber]>0.0,_minTimeStepSizes[solverNumber]);
      solver->updateMinNextTimeStepSize(_minTimeStepSizes[solverNumber]);
      if (exahype::State::fuseADERDGPhases()) {
        #ifdef Parallel
        if (tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
          weighMinNextPredictorTimeStepSize(solver);
        }
        #endif
        solver->updateTimeStepSizesFused();
      } else {
        solver->updateTimeStepSizes();
      }
    }
  }

  logTraceOutWith1Argument("endIteration(State)", state);
}

void exahype::mappings::TimeStepSizeComputation::enterCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  if (fineGridCell.isInitialised()) {
    const int numberOfSolvers = static_cast<int>(exahype::solvers::RegisteredSolvers.size());
    auto grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(numberOfSolvers, peano::datatraversal::autotuning::MethodTrace::UserDefined15);
    pfor(solverNumber, 0, numberOfSolvers, grainSize.getGrainSize())
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      if (solver->isComputingTimeStepSize(_localState.getAlgorithmSection())) {
        const int element = solver->tryGetElement(
            fineGridCell.getCellDescriptionsIndex(),solverNumber);

        double admissibleTimeStepSize = std::numeric_limits<double>::max();

        if (exahype::State::fuseADERDGPhases()) {
          admissibleTimeStepSize = solver->updateTimeStepSizesFused(
              fineGridCell.getCellDescriptionsIndex(),element);
        } else {
          admissibleTimeStepSize = solver->updateTimeStepSizes(
              fineGridCell.getCellDescriptionsIndex(),element);
        }

        _minTimeStepSizes[solverNumber] = std::min(
            admissibleTimeStepSize, _minTimeStepSizes[solverNumber]);
        _minCellSizes[solverNumber] = std::min(
            fineGridVerticesEnumerator.getCellSize()[0],_minCellSizes[solverNumber]);
        _maxCellSizes[solverNumber] = std::max(
            fineGridVerticesEnumerator.getCellSize()[0],_maxCellSizes[solverNumber]);
      }

    endpfor
    grainSize.parallelSectionHasTerminated();
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}


//
// All methods below are nop.
//
//=====================================


#ifdef Parallel

void exahype::mappings::TimeStepSizeComputation::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::mergeWithMaster(
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
  // do nothing
}

bool exahype::mappings::TimeStepSizeComputation::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  // do nothing
  return false;
}


void exahype::mappings::TimeStepSizeComputation::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Vertex& localVertex,
        const exahype::Vertex& masterOrWorkerVertex, int fromRank,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::receiveDataFromMaster(
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

void exahype::mappings::TimeStepSizeComputation::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

exahype::mappings::TimeStepSizeComputation::TimeStepSizeComputation() {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
