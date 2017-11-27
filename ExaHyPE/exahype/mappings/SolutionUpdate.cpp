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
 
#include "exahype/mappings/SolutionUpdate.h"

#include "tarch/multicore/Loop.h"

#include "peano/utils/Globals.h"

#include "peano/datatraversal/autotuning/Oracle.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/VertexOperations.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "exahype/mappings/TimeStepSizeComputation.h"

tarch::logging::Log exahype::mappings::SolutionUpdate::_log(
    "exahype::mappings::SolutionUpdate");

void exahype::mappings::SolutionUpdate::prepareLocalTimeStepVariables(){
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

void exahype::mappings::SolutionUpdate::initialiseTemporaryVariables() {
  if (exahype::State::fuseADERDGPhases()) {
    exahype::solvers::initialiseTemporaryVariables(_predictionTemporaryVariables);
  }
}

void exahype::mappings::SolutionUpdate::deleteTemporaryVariables() {
  if (exahype::State::fuseADERDGPhases()) {
    exahype::solvers::deleteTemporaryVariables(_predictionTemporaryVariables);
  }
}

peano::CommunicationSpecification
exahype::mappings::SolutionUpdate::communicationSpecification() const {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::MaskOutMasterWorkerDataAndStateExchange,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::MaskOutWorkerMasterDataAndStateExchange,
      true);
}



peano::MappingSpecification
exahype::mappings::SolutionUpdate::enterCellSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}
/**
 * Nop.
 */
peano::MappingSpecification
exahype::mappings::SolutionUpdate::touchVertexLastTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}


peano::MappingSpecification
exahype::mappings::SolutionUpdate::touchVertexFirstTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}

peano::MappingSpecification
exahype::mappings::SolutionUpdate::leaveCellSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::SolutionUpdate::ascendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}


peano::MappingSpecification
exahype::mappings::SolutionUpdate::descendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}

exahype::mappings::SolutionUpdate::SolutionUpdate() {
  // do nothing
}

exahype::mappings::SolutionUpdate::~SolutionUpdate() {
  exahype::solvers::deleteSolverFlags(_solverFlags);

  deleteTemporaryVariables();
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::SolutionUpdate::SolutionUpdate(
    const SolutionUpdate& masterThread)
  : _localState(masterThread._localState) {
  exahype::solvers::initialiseSolverFlags(_solverFlags);
  exahype::solvers::prepareSolverFlags(_solverFlags);

  prepareLocalTimeStepVariables();

  initialiseTemporaryVariables();
}
// Merge over threads
void exahype::mappings::SolutionUpdate::mergeWithWorkerThread(
    const SolutionUpdate& workerThread) {
  for (int i = 0; i < static_cast<int>(exahype::solvers::RegisteredSolvers.size()); i++) {
    // solver flags
    _solverFlags._meshUpdateRequest[i]  |= workerThread._solverFlags._meshUpdateRequest[i];
    _solverFlags._limiterDomainChange[i] =
        std::max ( _solverFlags._limiterDomainChange[i],
            workerThread._solverFlags._limiterDomainChange[i] );
    // time step size and cell sizes
    _minTimeStepSizes[i] =
        std::min(_minTimeStepSizes[i], workerThread._minTimeStepSizes[i]);
    _minCellSizes[i] =
        std::min(_minCellSizes[i], workerThread._minCellSizes[i]);
    _maxCellSizes[i] =
        std::max(_maxCellSizes[i], workerThread._maxCellSizes[i]);
  }
}
#endif

void exahype::mappings::SolutionUpdate::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,fineGridVerticesEnumerator.toString(),coarseGridCell, fineGridPositionOfCell);

  if (!exahype::State::EnableNeighbourCommunication) {  // TODO(Dominic): Add to docu.
    return;
  }

  assertion1(_localState.getAlgorithmSection()==exahype::records::State::AlgorithmSection::TimeStepping,_localState.getAlgorithmSection());

  if (fineGridCell.isInitialised()) {
    const int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
    auto grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(numberOfSolvers, peano::datatraversal::autotuning::MethodTrace::UserDefined14);
    pfor(solverNumber, 0, numberOfSolvers, grainSize.getGrainSize())
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      const int element = solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
      if (element!=exahype::solvers::Solver::NotFound) {
        exahype::solvers::Solver::UpdateResult result;

        if (exahype::State::fuseADERDGPhases()) {
          result =
              solver->fusedTimeStep(
                  fineGridCell.getCellDescriptionsIndex(),element,
                  exahype::State::isFirstIterationOfBatchOrNoBatch(),
                  exahype::State::isLastIterationOfBatchOrNoBatch(),
                  _predictionTemporaryVariables._tempSpaceTimeUnknowns    [solverNumber],
                  _predictionTemporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber],
                  _predictionTemporaryVariables._tempUnknowns             [solverNumber],
                  _predictionTemporaryVariables._tempFluxUnknowns         [solverNumber],
                  _predictionTemporaryVariables._tempPointForceSources    [solverNumber]);
        } else {
          solver->updateSolution(
              fineGridCell.getCellDescriptionsIndex(),element,true);
          result._timeStepSize =
              solver->startNewTimeStep(fineGridCell.getCellDescriptionsIndex(),element);

          if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
            auto* limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
            result._limiterDomainChange =  // !!! limiter status must be updated before refinement crit is evaluated
                limitingADERDGSolver->
                updateLimiterStatusAndMinAndMaxAfterSolutionUpdate(
                    fineGridCell.getCellDescriptionsIndex(),element);
            result._refinementRequested |=
                limitingADERDGSolver->evaluateRefinementCriterionAfterSolutionUpdate(
                    fineGridCell.getCellDescriptionsIndex(),element);
          } else {
            result._refinementRequested |=
                solver->evaluateRefinementCriterionAfterSolutionUpdate(
                    fineGridCell.getCellDescriptionsIndex(),element);
          }
        }

        _solverFlags._meshUpdateRequest  [solverNumber] |= result._refinementRequested;
        _solverFlags._limiterDomainChange[solverNumber]  = std::max( _solverFlags._limiterDomainChange[solverNumber], result._limiterDomainChange );
        assertion(_solverFlags._limiterDomainChange[solverNumber]!=exahype::solvers::LimiterDomainChange::IrregularRequiringMeshUpdate ||
            _solverFlags._meshUpdateRequest[solverNumber]);
        _minTimeStepSizes[solverNumber] = std::min( result._timeStepSize,                       _minTimeStepSizes[solverNumber]);
        _minCellSizes    [solverNumber] = std::min( fineGridVerticesEnumerator.getCellSize()[0],_minCellSizes    [solverNumber]);
        _maxCellSizes    [solverNumber] = std::max( fineGridVerticesEnumerator.getCellSize()[0],_maxCellSizes    [solverNumber]);
      }
    endpfor
    grainSize.parallelSectionHasTerminated();

    exahype::Cell::resetNeighbourMergeFlags(
        fineGridCell.getCellDescriptionsIndex());
    exahype::Cell::resetFaceDataExchangeCounters(
        fineGridCell.getCellDescriptionsIndex(),
        fineGridVertices,fineGridVerticesEnumerator);
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::SolutionUpdate::beginIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);

  if ( exahype::State::isFirstIterationOfBatchOrNoBatch() ) {
    _localState = solverState;

    if (_localState.getAlgorithmSection()==exahype::records::State::TimeStepping) {
      for (auto* solver : exahype::solvers::RegisteredSolvers) {
        solver->setNextMeshUpdateRequest();
        solver->setNextAttainedStableState();

        if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
          static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->setNextLimiterDomainChange();
        }
      }
    }

    // temporary variables
    prepareLocalTimeStepVariables();

    initialiseTemporaryVariables();

    exahype::solvers::initialiseSolverFlags(_solverFlags);
    exahype::solvers::prepareSolverFlags(_solverFlags);
  }

  logTraceOutWith1Argument("beginIteration(State)", solverState);
}

void exahype::mappings::SolutionUpdate::endIteration(
    exahype::State& state) {
  logTraceInWith1Argument("endIteration(State)", state);

  assertionEquals(_localState.getAlgorithmSection(),exahype::records::State::AlgorithmSection::TimeStepping);

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    /*
     * Update reduced quantities (over multiple batch iterations)
     */
    // mesh refinement events
    solver->updateNextMeshUpdateRequest(_solverFlags._meshUpdateRequest[solverNumber]);
    solver->updateNextAttainedStableState(!solver->getNextMeshUpdateRequest());
    if (exahype::solvers::RegisteredSolvers[solverNumber]->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
      auto* limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
      limitingADERDGSolver->updateNextLimiterDomainChange(_solverFlags._limiterDomainChange[solverNumber]);
    }
    // cell sizes (for AMR)
    solver->updateNextMinCellSize(_minCellSizes[solverNumber]);
    solver->updateNextMaxCellSize(_maxCellSizes[solverNumber]);

    // time
    assertion1(std::isfinite(_minTimeStepSizes[solverNumber]),_minTimeStepSizes[solverNumber]);
    assertion1(_minTimeStepSizes[solverNumber]>0.0,_minTimeStepSizes[solverNumber]);
    solver->updateMinNextTimeStepSize(_minTimeStepSizes[solverNumber]);

    /*
     * Swap the current values with the next values (in last batch iteration)
     */
    // mesh update events
    if ( exahype::State::isLastIterationOfBatchOrNoBatch() ) {
      solver->setNextMeshUpdateRequest();
      solver->setNextAttainedStableState();
      if (exahype::solvers::RegisteredSolvers[solverNumber]->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
        auto* limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
        limitingADERDGSolver->setNextLimiterDomainChange();
        assertion(
            limitingADERDGSolver->getLimiterDomainChange()
            !=exahype::solvers::LimiterDomainChange::IrregularRequiringMeshUpdate ||
            solver->getMeshUpdateRequest());
      }
    }

    // time
    // TODO(Dominic):
    // only update the time step size in last iteration; just advance with old time step size otherwise
    if (
        exahype::State::fuseADERDGPhases() &&
        exahype::State::isLastIterationOfBatchOrNoBatch() &&
        tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getGlobalMasterRank()
    ) {
      exahype::mappings::TimeStepSizeComputation::
      reinitialiseTimeStepDataIfLastPredictorTimeStepSizeWasInstable(solver);
    }
    if (exahype::State::fuseADERDGPhases()) {
      solver->startNewTimeStepFused(
          exahype::State::isFirstIterationOfBatchOrNoBatch(),
          exahype::State::isLastIterationOfBatchOrNoBatch());
    } else {
      solver->startNewTimeStep();
    }
  }

  // delete temporary variables
  if ( exahype::State::isLastIterationOfBatchOrNoBatch() ) {
    deleteSolverFlags(_solverFlags);
    deleteTemporaryVariables();
  }

  logTraceOutWith1Argument("endIteration(State)", state);
}



//
// Below all methods are nop.
//
//=====================================



void exahype::mappings::SolutionUpdate::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

#ifdef Parallel
void exahype::mappings::SolutionUpdate::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::SolutionUpdate::prepareSendToWorker(
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

void exahype::mappings::SolutionUpdate::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::mergeWithMaster(
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

void exahype::mappings::SolutionUpdate::receiveDataFromMaster(
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

void exahype::mappings::SolutionUpdate::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::SolutionUpdate::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
