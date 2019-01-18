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
 
#include "exahype/mappings/FinaliseMeshRefinement.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "tarch/multicore/Loop.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/datatraversal/TaskSet.h"

#include "peano/parallel/loadbalancing/Oracle.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "exahype/mappings/MeshRefinement.h"
#include "exahype/mappings/RefinementStatusSpreading.h"

#include "exahype/stealing/PerformanceMonitor.h"

#ifdef DistributedStealing
#include "exahype/stealing/AggressiveDistributor.h"
#include "exahype/stealing/AggressiveCCPDistributor.h"
#endif

int exahype::mappings::FinaliseMeshRefinement::NumberOfEnclaveCells = 0;
int exahype::mappings::FinaliseMeshRefinement::NumberOfSkeletonCells = 0;

peano::CommunicationSpecification
exahype::mappings::FinaliseMeshRefinement::communicationSpecification() const {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
          MaskOutMasterWorkerDataAndStateExchange,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
          SendDataAndStateAfterLastTouchVertexLastTime,
      true);
}


peano::MappingSpecification
exahype::mappings::FinaliseMeshRefinement::enterCellSpecification(int level) const {
  const int coarsestSolverLevel = solvers::Solver::getCoarsestMeshLevelOfAllSolvers();
  if ( std::abs(level)>=coarsestSolverLevel ) {
    return peano::MappingSpecification(
          peano::MappingSpecification::WholeTree,
          peano::MappingSpecification::RunConcurrentlyOnFineGrid,true); // performs reduction
  } else {
    return peano::MappingSpecification(
          peano::MappingSpecification::Nop,
          peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
  }
}

// Below all specs are Nop
peano::MappingSpecification exahype::mappings::FinaliseMeshRefinement::
    touchVertexFirstTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
}

peano::MappingSpecification
exahype::mappings::FinaliseMeshRefinement::leaveCellSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
}

peano::MappingSpecification
exahype::mappings::FinaliseMeshRefinement::touchVertexLastTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
}

peano::MappingSpecification
exahype::mappings::FinaliseMeshRefinement::ascendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,false);
}

peano::MappingSpecification
exahype::mappings::FinaliseMeshRefinement::descendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,false);
}

void exahype::mappings::FinaliseMeshRefinement::initialiseLocalVariables(){
  const unsigned int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
  _minTimeStepSizes.resize(numberOfSolvers);
  _maxLevels.resize(numberOfSolvers);

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    _minTimeStepSizes[solverNumber] = std::numeric_limits<double>::max();
    _maxLevels    [solverNumber]    = -std::numeric_limits<int>::max(); // "-", min
  }

  _numberOfEnclaveCells = 0;
  _numberOfSkeletonCells = 0;
}

tarch::logging::Log exahype::mappings::FinaliseMeshRefinement::_log(
    "exahype::mappings::FinaliseMeshRefinement");

bool exahype::mappings::FinaliseMeshRefinement::OneSolverRequestedMeshUpdate = false;

exahype::mappings::FinaliseMeshRefinement::FinaliseMeshRefinement() {}

exahype::mappings::FinaliseMeshRefinement::~FinaliseMeshRefinement() {}

#if defined(SharedMemoryParallelisation)
exahype::mappings::FinaliseMeshRefinement::FinaliseMeshRefinement(const FinaliseMeshRefinement& masterThread) {
  _backgroundJobsHaveTerminated=masterThread._backgroundJobsHaveTerminated;
  initialiseLocalVariables();
}

// Merge over threads
void exahype::mappings::FinaliseMeshRefinement::mergeWithWorkerThread(
    const FinaliseMeshRefinement& workerThread) {
  for (int i = 0; i < static_cast<int>(exahype::solvers::RegisteredSolvers.size()); i++) {
    _minTimeStepSizes[i] =
        std::min(_minTimeStepSizes[i], workerThread._minTimeStepSizes[i]);
    _maxLevels[i] =
        std::max(_maxLevels[i], workerThread._maxLevels[i]);
  }
  _numberOfEnclaveCells += workerThread._numberOfEnclaveCells;
  _numberOfSkeletonCells += workerThread._numberOfSkeletonCells;
}
#endif


void exahype::mappings::FinaliseMeshRefinement::beginIteration(exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);

  tarch::multicore::jobs::Job::setMaxNumberOfRunningBackgroundThreads(
      solvers::Solver::MaxNumberOfRunningBackgroundJobConsumerTasksDuringTraversal); // reset the number of running consumer threads
  peano::parallel::loadbalancing::Oracle::getInstance().activateLoadBalancing(false);

  OneSolverRequestedMeshUpdate =
      exahype::solvers::Solver::oneSolverRequestedMeshRefinement();

  exahype::solvers::Solver::rollbackSolversToPreviousTimeStepIfApplicable();

  exahype::mappings::MeshRefinement::IsFirstIteration = true;

  initialiseLocalVariables();

  logTraceOutWith1Argument("beginIteration(State)", solverState);
}

void exahype::mappings::FinaliseMeshRefinement::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  if (
      OneSolverRequestedMeshUpdate &&
      fineGridCell.isInitialised()
  ) {
    if ( !_backgroundJobsHaveTerminated ) {
      exahype::solvers::Solver::ensureAllJobsHaveTerminated(exahype::solvers::Solver::JobType::AMRJob);
      _backgroundJobsHaveTerminated = true;
    } // TODO(Dominic): Still necessary? Mesh Refinement terminates the background jobs in endIteration now.
    solvers::Solver::CellInfo cellInfo = fineGridCell.createCellInfo();

    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

      if ( solver->hasRequestedMeshRefinement() ) {
        solver->finaliseStateUpdates(solverNumber,cellInfo);

        if ( solver->getMeshUpdateEvent()==exahype::solvers::Solver::MeshUpdateEvent::RefinementRequested ) { // is not the same as the above check
          solver->rollbackSolutionGlobally(solverNumber,cellInfo,exahype::solvers::Solver::FuseADERDGPhases);
        }

        // compute a new time step size
        double admissibleTimeStepSize = solver->updateTimeStepSizes(solverNumber,cellInfo,exahype::solvers::Solver::FuseADERDGPhases);
        _minTimeStepSizes[solverNumber] = std::min(admissibleTimeStepSize, _minTimeStepSizes[solverNumber]);
        _maxLevels[solverNumber]        = std::max(fineGridVerticesEnumerator.getLevel(),_maxLevels[solverNumber]);

        // determine min and max for LimitingADERDGSolver
        if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
          auto* limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
          limitingADERDGSolver->determineMinAndMax(solverNumber,cellInfo);
          assertion(limitingADERDGSolver->getMeshUpdateEvent()!=exahype::solvers::Solver::MeshUpdateEvent::IrregularLimiterDomainChange);
        }
        
        // determine numbers of enclave and skeleton cells
        if (solver->getType()==exahype::solvers::Solver::Type::ADERDG || solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
          //auto* ADERDGSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);
          if( !exahype::Cell::isAtRemoteBoundary(fineGridVertices,fineGridVerticesEnumerator))
            _numberOfEnclaveCells += exahype::solvers::ADERDGSolver::computeWeight(cellInfo._cellDescriptionsIndex);
          else if( exahype::Cell::isAtRemoteBoundary(fineGridVertices,fineGridVerticesEnumerator))
            _numberOfSkeletonCells += exahype::solvers::ADERDGSolver::computeWeight(cellInfo._cellDescriptionsIndex);
        }
        
      }
    }

    Cell::resetNeighbourMergePerformedFlags(cellInfo,fineGridVertices,fineGridVerticesEnumerator);
  }
}

void exahype::mappings::FinaliseMeshRefinement::endIteration(
    exahype::State& solverState) {

  if (OneSolverRequestedMeshUpdate) {
    if (exahype::mappings::MeshRefinement::IsInitialMeshRefinement) {
      peano::performanceanalysis::Analysis::getInstance().enable(true);
    }
    exahype::mappings::MeshRefinement::IsInitialMeshRefinement = false;

    // time stepping
    for (unsigned int solverNumber = 0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

      if (solver->hasRequestedMeshRefinement()) {
        // cell sizes
        solver->updateNextMaxLevel(_maxLevels[solverNumber]);

        // time
        assertion1(std::isfinite(_minTimeStepSizes[solverNumber]),_minTimeStepSizes[solverNumber]);
        assertion1(_minTimeStepSizes[solverNumber]>0.0,_minTimeStepSizes[solverNumber]);
        solver->updateMinNextTimeStepSize(_minTimeStepSizes[solverNumber]);
        if ( exahype::solvers::Solver::FuseADERDGPhases ) {
          #ifdef Parallel
          if (tarch::parallel::Node::getInstance().isGlobalMaster() ) {
            exahype::solvers::Solver::weighMinNextPredictorTimeStepSize(solver);
          }
          #endif
          solver->updateTimeStepSizesFused();
        } else {
          solver->updateTimeStepSizes();
        }
      }
    }
  }

  NumberOfEnclaveCells = _numberOfEnclaveCells;
  NumberOfSkeletonCells = _numberOfSkeletonCells;

#ifdef DistributedStealing
  exahype::stealing::PerformanceMonitor::getInstance().setTasksPerTimestep(_numberOfEnclaveCells + _numberOfSkeletonCells);

#if defined(StealingStrategyAggressive)
  exahype::stealing::AggressiveDistributor::getInstance().enable();
#elif defined(StealingStrategyAggressiveCCP)
  exahype::stealing::AggressiveCCPDistributor::getInstance().enable();
#endif

#endif

  _backgroundJobsHaveTerminated = false;
}

#ifdef Parallel
void exahype::mappings::FinaliseMeshRefinement::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  logTraceInWith6Arguments("mergeWithNeighbour(...)", vertex, neighbour,
                           fromRank, fineGridX, fineGridH, level);

  vertex.dropNeighbourMetadata(fromRank,fineGridX,fineGridH,level);

  logTraceOut("mergeWithNeighbour(...)");
}

bool exahype::mappings::FinaliseMeshRefinement::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  return true;
}

void exahype::mappings::FinaliseMeshRefinement::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {

  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    if ( solver->hasRequestedMeshRefinement() ) {
      solver->sendDataToMaster(
          tarch::parallel::NodePool::getInstance().getMasterRank(),
          verticesEnumerator.getCellCenter(),
          verticesEnumerator.getLevel());
    }
  }
}

void exahype::mappings::FinaliseMeshRefinement::mergeWithMaster(
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
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    if ( solver->hasRequestedMeshRefinement() ) {
      solver->mergeWithWorkerData(
          worker,
          fineGridVerticesEnumerator.getCellCenter(),
          fineGridVerticesEnumerator.getLevel());
    }
  }
}

//
// All methods below are nop.
//
// ====================================



void exahype::mappings::FinaliseMeshRefinement::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {}

void exahype::mappings::FinaliseMeshRefinement::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {}

void exahype::mappings::FinaliseMeshRefinement::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {}

void exahype::mappings::FinaliseMeshRefinement::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Vertex& localVertex,
        const exahype::Vertex& masterOrWorkerVertex, int fromRank,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, double>& h, int level) {}

void exahype::mappings::FinaliseMeshRefinement::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {}

void exahype::mappings::FinaliseMeshRefinement::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {}

void exahype::mappings::FinaliseMeshRefinement::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {}

void exahype::mappings::FinaliseMeshRefinement::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {}
#endif

void exahype::mappings::FinaliseMeshRefinement::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::FinaliseMeshRefinement::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::FinaliseMeshRefinement::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::FinaliseMeshRefinement::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::FinaliseMeshRefinement::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::FinaliseMeshRefinement::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {}

void exahype::mappings::FinaliseMeshRefinement::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {}

void exahype::mappings::FinaliseMeshRefinement::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::FinaliseMeshRefinement::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::FinaliseMeshRefinement::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {}

void exahype::mappings::FinaliseMeshRefinement::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {}

void exahype::mappings::FinaliseMeshRefinement::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {}
