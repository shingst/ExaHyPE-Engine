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
 *
 * \author Dominic E. Charrier, Tobias Weinzierl, Jean-Matthieu Gallard, Fabian Gra, Leonhard Rannabauer, Philipp Samfass
 **/
#include "exahype/solvers/ADERDGSolver.h"

#include "exahype/reactive/ResilienceTools.h"
#include "exahype/reactive/PerformanceMonitor.h"
#include "exahype/reactive/ReactiveContext.h"
#include "exahype/reactive/RequestManager.h"
#include "exahype/reactive/OffloadingProgressService.h"
#include "exahype/reactive/OffloadingProfiler.h"
#include "exahype/reactive/TimeStampAndTriggerTeamHistory.h"
#include "peano/utils/UserInterface.h"
#include "exahype/reactive/ResilienceStatistics.h"

#include "exahype/Cell.h"
#include "exahype/Vertex.h"
#include "exahype/VertexOperations.h"

#include "tarch/la/VectorVectorOperations.h"
#include "tarch/multicore/Lock.h"

#include "exahype/mappings/LevelwiseAdjacencyBookkeeping.h"
#include "exahype/amr/AdaptiveMeshRefinement.h"

#include "peano/heap/CompressedFloatingPointNumbers.h"
#include "peano/datatraversal/TaskSet.h"
#include "peano/grid/aspects/VertexStateAnalysis.h"

#include "exahype/solvers/LimitingADERDGSolver.h"
#include "exahype/solvers/OutcomeDatabase.h"

#include "kernels/KernelUtils.h"

#include "tarch/multicore/Jobs.h"
#include "tarch/multicore/Core.h"
#include "tarch/la/Vector.h"
#include "tarch/timing/Watch.h"

#include <limits>
#include <iomanip>
#include <vector>
#include <chrono>
#include <algorithm> // copy_n

#if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
#include <immintrin.h>
#endif

#if defined(USE_TMPI)
#include "teaMPI.h"
#endif

#if defined(UseSmartMPI)
#include "mpi_offloading.h"
#endif

#ifndef MPI_BLOCKING
#define MPI_BLOCKING false
#endif

namespace {
  constexpr const char* tags[]{"solutionUpdate",
                             "volumeIntegral",
                             "surfaceIntegral",
                             "riemannSolver",
                             "spaceTimePredictor",
                             "stableTimeStepSize",
                             "solutionAdjustment",
                             "faceUnknownsProlongation",
                             "faceUnknownsRestriction",
                             "volumeUnknownsProlongation",
                             "volumeUnknownsRestriction",
                             "boundaryConditions",
                             "deltaDistribution"
                             };
}

#ifdef USE_ITAC
int exahype::solvers::ADERDGSolver::adjustSolutionHandle                  = 0;
int exahype::solvers::ADERDGSolver::fusedTimeStepBodyHandle               = 0;
int exahype::solvers::ADERDGSolver::fusedTimeStepBodyHandleSkeleton       = 0;
int exahype::solvers::ADERDGSolver::predictorBodyHandle                   = 0;
int exahype::solvers::ADERDGSolver::predictorBodyHandleSkeleton           = 0;
int exahype::solvers::ADERDGSolver::updateBodyHandle                      = 0;
int exahype::solvers::ADERDGSolver::mergeNeighboursHandle                 = 0;
int exahype::solvers::ADERDGSolver::prolongateFaceDataToVirtualCellHandle = 0;
int exahype::solvers::ADERDGSolver::restrictToTopMostParentHandle         = 0;

//Todo(Philipp) : remove unused
int exahype::solvers::ADERDGSolver::event_stp = 0;
int exahype::solvers::ADERDGSolver::event_stp_remote = 0;
int exahype::solvers::ADERDGSolver::event_stp_local_replica = 0;
int exahype::solvers::ADERDGSolver::event_offloadingManager = 0;
int exahype::solvers::ADERDGSolver::event_spawn = 0;
int exahype::solvers::ADERDGSolver::event_initial = 0;
int exahype::solvers::ADERDGSolver::event_memory = 0;
int exahype::solvers::ADERDGSolver::event_lock = 0;
int exahype::solvers::ADERDGSolver::event_pack = 0;
int exahype::solvers::ADERDGSolver::event_progress = 0;
int exahype::solvers::ADERDGSolver::event_offload = 0;

#endif

tarch::logging::Log exahype::solvers::ADERDGSolver::_log( "exahype::solvers::ADERDGSolver");

// communication status
int exahype::solvers::ADERDGSolver::LeafCommunicationStatus                             = 2;
int exahype::solvers::ADERDGSolver::MinimumCommunicationStatusForNeighbourCommunication = 1;
// augmentation status
// On-the fly erasing seems to work with those values
int exahype::solvers::ADERDGSolver::MaximumAugmentationStatus                   = 3;
int exahype::solvers::ADERDGSolver::MinimumAugmentationStatusForVirtualRefining = 2;
int exahype::solvers::ADERDGSolver::MinimumAugmentationStatusForRefining        = 1;

/**
 * static constexpr need to declared again when following a
 * C++ standard before C++17.
 */
constexpr int exahype::solvers::ADERDGSolver::EmptyStatus;
constexpr int exahype::solvers::ADERDGSolver::BoundaryStatus;
constexpr int exahype::solvers::ADERDGSolver::Pending;
constexpr int exahype::solvers::ADERDGSolver::Erase; 
constexpr int exahype::solvers::ADERDGSolver::Keep;

tarch::multicore::BooleanSemaphore exahype::solvers::ADERDGSolver::RestrictionSemaphore;
tarch::multicore::BooleanSemaphore exahype::solvers::ADERDGSolver::CoarseGridSemaphore;

#if defined(Parallel)
template class exahype::solvers::OutcomeDatabase<exahype::solvers::ADERDGSolver::MigratablePredictionJobOutcomeKey, exahype::solvers::ADERDGSolver::MigratablePredictionJobData>;
#include "exahype/solvers/OutcomeDatabase.cpph"
#endif

#if defined(SharedTBB) && defined(Parallel)
std::atomic<int> exahype::solvers::ADERDGSolver::MaxIprobesInOffloadingProgress (std::numeric_limits<int>::max());

std::atomic<int> exahype::solvers::ADERDGSolver::MigratablePredictionJob::JobCounter (0);
std::atomic<int> exahype::solvers::ADERDGSolver::NumberOfReceiveJobs (0);
std::atomic<int> exahype::solvers::ADERDGSolver::NumberOfReceiveBackJobs (0);
std::atomic<int> exahype::solvers::ADERDGSolver::LocalStealableSTPCounter (0);
std::atomic<int> exahype::solvers::ADERDGSolver::CompletedSentSTPs(0);
std::atomic<int> exahype::solvers::ADERDGSolver::SentSTPs (0);
std::atomic<int> exahype::solvers::ADERDGSolver::AllocatedSTPs (0);
std::atomic<int> exahype::solvers::ADERDGSolver::AllocatedSTPsSend (0);
std::atomic<int> exahype::solvers::ADERDGSolver::AllocatedSTPsReceive (0);

std::atomic<bool> exahype::solvers::ADERDGSolver::VetoEmergency(false);
const exahype::solvers::ADERDGSolver::CellDescription* exahype::solvers::ADERDGSolver::LastEmergencyCell;
tarch::multicore::BooleanSemaphore exahype::solvers::ADERDGSolver::EmergencySemaphore;

#endif

#ifdef OffloadingUseProgressTask
std::unordered_set<int> exahype::solvers::ADERDGSolver::ActiveSenders;
#endif

int exahype::solvers::ADERDGSolver::computeWeight(const int cellDescriptionsIndex) {
  if ( ADERDGSolver::isValidCellDescriptionIndex(cellDescriptionsIndex) ) {
    int result = 0;
    for ( CellDescription& cellDescription : getCellDescriptions(cellDescriptionsIndex) ) {
      result += ( cellDescription.getType()==CellDescription::Type::Leaf ) ?  1 : 0;
    }
    return result;
  }
  else return 0;
}

/**
 * Returns the ADERDGCellDescription heap vector
 * at address \p cellDescriptionsIndex.
 */
exahype::solvers::ADERDGSolver::Heap::HeapEntries& exahype::solvers::ADERDGSolver::getCellDescriptions(
    const int cellDescriptionsIndex) {
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);

  return Heap::getInstance().getData(cellDescriptionsIndex);
}

/**
 * Returns the ADERDGCellDescription with index \p element
 * in the heap vector at address \p cellDescriptionsIndex.
 */
exahype::solvers::ADERDGSolver::CellDescription& exahype::solvers::ADERDGSolver::getCellDescription(
    const int cellDescriptionsIndex,
    const int element) {
  assertion2(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex,element);
  assertion2(element>=0,cellDescriptionsIndex,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),cellDescriptionsIndex,element);

  return Heap::getInstance().getData(cellDescriptionsIndex)[element];
}

bool exahype::solvers::ADERDGSolver::communicateWithNeighbour(const CellDescription& cellDescription,const int faceIndex) {
  assertion1(cellDescription.getType()!=CellDescription::Type::Leaf ||
            cellDescription.getCommunicationStatus()==LeafCommunicationStatus,cellDescription.toString());
  return
      (cellDescription.getCommunicationStatus()                  == LeafCommunicationStatus &&
      cellDescription.getFacewiseCommunicationStatus(faceIndex)  >= MinimumCommunicationStatusForNeighbourCommunication &&
      cellDescription.getFacewiseAugmentationStatus(faceIndex)   <  MaximumAugmentationStatus)
      ||
      (cellDescription.getFacewiseCommunicationStatus(faceIndex) == LeafCommunicationStatus &&
      cellDescription.getCommunicationStatus()                   >= MinimumCommunicationStatusForNeighbourCommunication &&
      cellDescription.getAugmentationStatus()                    <  MaximumAugmentationStatus);
}

void exahype::solvers::ADERDGSolver::prefetchFaceData(CellDescription& cellDescription,const int faceIndex) {
  #if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
  auto* solver = solvers::RegisteredSolvers[cellDescription.getSolverNumber()];
  int dataPerFace = 0;
  int dofPerFace  = 0;
  switch (solver->getType()) {
    case Solver::Type::ADERDG:
      dataPerFace = static_cast<ADERDGSolver*>(solver)->getBndFaceSize();
      dofPerFace  = static_cast<ADERDGSolver*>(solver)->getBndFluxSize();
      break;
    case Solver::Type::LimitingADERDG:
      dataPerFace = static_cast<LimitingADERDGSolver*>(solver)->getSolver()->getBndFaceSize();
      dofPerFace  = static_cast<LimitingADERDGSolver*>(solver)->getSolver()->getBndFluxSize();
      break;
    default:
      break;
  }

  double* lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor()) + faceIndex * dataPerFace;
  double* lFhbnd = static_cast<double*>(cellDescription.getFluctuation())           + faceIndex * dofPerFace;

  _mm_prefetch(lQhbnd, _MM_HINT_NTA);
  _mm_prefetch(lFhbnd, _MM_HINT_NTA);
  #endif
}

exahype::solvers::ADERDGSolver::ADERDGSolver(
    const std::string& identifier,
    const int numberOfVariables,
    const int numberOfParameters,
    const int numberOfGlobalObservables,
    const int basisSize,
    const double maximumMeshSize,
    const int maximumAdaptiveMeshDepth,
    const int haloCells,
    const int haloBufferCells,
    const int limiterBufferCells,
    const int regularisedFineGridLevels,
    const exahype::solvers::Solver::TimeStepping timeStepping,
    const int DMPObservables,
    std::unique_ptr<profilers::Profiler> profiler)
    : Solver(identifier, Solver::Type::ADERDG, numberOfVariables,
             numberOfParameters, numberOfGlobalObservables, basisSize,
             maximumMeshSize, maximumAdaptiveMeshDepth,
             timeStepping, std::move(profiler)),
     _previousMinTimeStamp( std::numeric_limits<double>::infinity() ),
     _previousMinTimeStepSize( std::numeric_limits<double>::infinity() ),
     _minTimeStamp( std::numeric_limits<double>::infinity() ),
     _minTimeStepSize( std::numeric_limits<double>::infinity() ),
     _estimatedTimeStepSize( std::numeric_limits<double>::infinity() ),
     _admissibleTimeStepSize( std::numeric_limits<double>::infinity() ),
     _stabilityConditionWasViolated( false ),
     _minimumRefinementStatusToRequestMeshRefinementInVirtualCell(1+haloBufferCells),
     _refineOrKeepOnFineGrid(1+haloCells+haloBufferCells),
     _DMPObservables(DMPObservables),
     _minRefinementStatusForTroubledCell(_refineOrKeepOnFineGrid+3),
     _checkForNaNs(true),
     _meshUpdateEvent(MeshUpdateEvent::None)
#if defined(Parallel) // this is not nice at all -> maybe move into constructor body
     ,_lastReceiveTag(tarch::parallel::Node::getInstance().getNumberOfNodes()),
     _lastReceiveBackTag(tarch::parallel::Node::getInstance().getNumberOfNodes()),
     _offloadingManagerJob(nullptr),
     _offloadingManagerJobTerminated(false),
     _offloadingManagerJobStarted(false),
     _offloadingManagerJobTriggerTerminate(false),
     _lastReceiveReplicaTag(tarch::parallel::Node::getInstance().getNumberOfNodes()
                            *exahype::reactive::ReactiveContext::getInstance().getTMPINumTeams()),
     _outcomeDatabase()
#if defined(SharedTBB) //todo(Philipp): this is super ugly, should avoid TBB defines, make design better!
     ,_pendingOutcomesToBeShared(),
     _allocatedOutcomes()
#endif
#endif
{
  // register tags with profiler
  for (const char* tag : tags) {
    _profiler->registerTag(tag);
  }

#ifdef Parallel
#if defined(SharedTBB)
  MigratablePredictionJobMetaData::initDatatype();
  //todo: may need to add support for multiple solvers
  exahype::reactive::OffloadingProgressService::getInstance().setSolver(this);
#endif

#ifdef OffloadingUseProfiler
  exahype::reactive::OffloadingProfiler::getInstance().beginPhase();
#endif

#endif
}

int exahype::solvers::ADERDGSolver::getUnknownsPerFace() const {
  return _numberOfVariables * power(_nodesPerCoordinateAxis, DIMENSIONS - 1);
}

int exahype::solvers::ADERDGSolver::getUnknownsPerCellBoundary() const {
  return DIMENSIONS_TIMES_TWO * getUnknownsPerFace();
}

int exahype::solvers::ADERDGSolver::getUnknownsPerCell() const {
  return _numberOfVariables * power(_nodesPerCoordinateAxis, DIMENSIONS + 0);
}

int exahype::solvers::ADERDGSolver::getFluxUnknownsPerCell() const {
  return (DIMENSIONS + 1) * getUnknownsPerCell(); // +1 for sources
}

int exahype::solvers::ADERDGSolver::getSpaceTimeUnknownsPerCell() const {
  return _numberOfVariables * power(_nodesPerCoordinateAxis, DIMENSIONS + 1);
}

int exahype::solvers::ADERDGSolver::getSpaceTimeFluxUnknownsPerCell() const {
  return (DIMENSIONS + 1) * getSpaceTimeUnknownsPerCell();  // +1 for sources
}

int exahype::solvers::ADERDGSolver::getDataPerFace() const {
  return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS - 1);
}

int exahype::solvers::ADERDGSolver::getDataPerCellBoundary() const {
  return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS - 1) * DIMENSIONS_TIMES_TWO;
}

int exahype::solvers::ADERDGSolver::getDataPerCell() const {
  return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS + 0);
}

int exahype::solvers::ADERDGSolver::getSpaceTimeDataPerCell() const {
  return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS + 1);
}

int exahype::solvers::ADERDGSolver::getDMPObservables() const {
  return _DMPObservables;
}

int exahype::solvers::ADERDGSolver::getMinRefinementStatusForTroubledCell() const {
  return _minRefinementStatusForTroubledCell;
}

void exahype::solvers::ADERDGSolver::resetMeshUpdateEvent() {
  _meshUpdateEvent = MeshUpdateEvent::None;
}

void exahype::solvers::ADERDGSolver::updateMeshUpdateEvent(MeshUpdateEvent meshUpdateEvent) {
  tarch::multicore::Lock lock(_reductionSemaphore);
  _meshUpdateEvent = mergeMeshUpdateEvents(_meshUpdateEvent,meshUpdateEvent);
  lock.free();
}

exahype::solvers::ADERDGSolver::MeshUpdateEvent
exahype::solvers::ADERDGSolver::getMeshUpdateEvent() const {
  return _meshUpdateEvent;
}

std::tuple<double,double> exahype::solvers::ADERDGSolver::getRiemannSolverTimeStepData(
    const CellDescription& cellDescription1,
    const CellDescription& cellDescription2) const {
  auto result  = std::make_tuple<double,double>(0.0,0.0);
  switch (_timeStepping) {
    case TimeStepping::Global:
    case TimeStepping::GlobalFixed:
      std::get<0>(result) = _minTimeStamp;
      std::get<1>(result) = _minTimeStepSize;
      break;
    default:
//      std::get<0>(result) = std::max(cellDescription1.getTimeStamp(),cellDescription2.getTimeStamp());
//      std::get<1>(result) = std::min(cellDescription1.getTimeStepSize(),cellDescription2.getTimeStepSize());
      logError("getRiemannSolverTimeStepData(...)","Unknown time stepping scheme.")
      std::abort();
  }
  return result;
}

std::tuple<double,double> exahype::solvers::ADERDGSolver::getPredictionTimeStepData(
    const CellDescription& cellDescription,const bool duringFusedTimeStep) const {
  auto result  = std::make_tuple<double,double>(0.0,0.0);
  switch (_timeStepping) {
    case TimeStepping::Global:
      if ( duringFusedTimeStep ) {
        std::get<0>(result) = _minTimeStamp+_minTimeStepSize;
        std::get<1>(result) = _estimatedTimeStepSize;
      } else {
        std::get<0>(result) = _minTimeStamp;
        std::get<1>(result) = _minTimeStepSize;
      }
      break;
    case TimeStepping::GlobalFixed:
      if ( duringFusedTimeStep ) {
        std::get<0>(result) = _minTimeStamp+_minTimeStepSize;
        std::get<1>(result) = _minTimeStepSize;
      } else {
        std::get<0>(result) = _minTimeStamp;
        std::get<1>(result) = _minTimeStepSize;
      }
      break;
      break;
    default:
      logError("getRiemannSolverTimeStepData(...)","Unknown time stepping scheme.")
      std::abort();
  }
  return result;
}

void exahype::solvers::ADERDGSolver::synchroniseTimeStepping(
    CellDescription& p) const {
  switch (_timeStepping) {
    case TimeStepping::Global:
    case TimeStepping::GlobalFixed:
      p.setPreviousTimeStamp(_previousMinTimeStamp);
      p.setPreviousTimeStepSize(_previousMinTimeStepSize);
  
      p.setTimeStamp(_minTimeStamp);
      p.setTimeStepSize(_minTimeStepSize);
  }
}

void exahype::solvers::ADERDGSolver::kickOffTimeStep(const bool isFirstTimeStepOfBatchOrNoBatch) {
  if ( isFirstTimeStepOfBatchOrNoBatch ) {
    _meshUpdateEvent               = MeshUpdateEvent::None;
    _admissibleTimeStepSize        = std::numeric_limits<double>::infinity();
    _stabilityConditionWasViolated = false;
    resetGlobalObservables(_nextGlobalObservables.data());
  }

  // call user code
  beginTimeStep(_minTimeStamp,isFirstTimeStepOfBatchOrNoBatch);
}

void exahype::solvers::ADERDGSolver::rollbackTimeStepMetadataToLastConsistentTimeStep() {

  logDebug("rollbackTimeStepMetadataToLastConsistentTimeStep","Limiting solver needs to rollback to team solution of previous consistent time stamp");
  exahype::reactive::TimeStampAndTriggerTeamHistory::getInstance().getLastConsistentTimeStepData(_minTimeStamp, _minTimeStepSize,_estimatedTimeStepSize);
  _admissibleTimeStepSize = _minTimeStepSize;
  logInfo("rollbackTimeStepMetadataToLastConsistentTimeStep","Trying to recover DG solution from other team next timeStamp="
                                               <<_minTimeStamp
                                               <<" time step size = "<<_minTimeStepSize
                                               <<" estimated time step size ="<<_estimatedTimeStepSize);

}

void exahype::solvers::ADERDGSolver::wrapUpTimeStep(const bool isFirstTimeStepOfBatchOrNoBatch,const bool isLastTimeStepOfBatchOrNoBatch) {
  if ( isFirstTimeStepOfBatchOrNoBatch ) {
    _previousMinTimeStepSize  = _minTimeStepSize;
    _previousMinTimeStamp     = _minTimeStamp;
  }

#if defined(Parallel)
  int team = exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber();
#else
  int team = 0;
#endif
  exahype::reactive::TimeStampAndTriggerTeamHistory::getInstance().trackTimeStepAndTriggerActive(team,
                                                                                          _minTimeStamp,
                                                                                          _minTimeStepSize,
                                                                                          exahype::reactive::ResilienceTools::CheckAllMigratableSTPs);


  // with rollback, timestamps should not be adapted, they are set in LimitingADERDGSolver::wrapUpTimeStep and in ADERDGSolver::mergeWithWorkerData
  if(!(_meshUpdateEvent==MeshUpdateEvent::RollbackToTeamSolution)) {
    _minTimeStamp += _minTimeStepSize;

    _stabilityConditionWasViolated = false;
    if (
      tarch::parallel::Node::getInstance().isGlobalMaster() &&
      getTimeStepping() != TimeStepping::GlobalFixed // fix the time step size in intermediate batch iterations
    ) {
      if ( FuseAllADERDGPhases && !isLinear() ) {
        if ( isLastTimeStepOfBatchOrNoBatch ) {
          if ( _estimatedTimeStepSize > _admissibleTimeStepSize ) { // rerun
            _minTimeStepSize       = FusedTimeSteppingRerunFactor * _admissibleTimeStepSize;
            _estimatedTimeStepSize = _minTimeStepSize;
            _stabilityConditionWasViolated = true;
            //logDebug("wrapUpTimeStep","recompute dt_min"<<std::setprecision(30)<<_minTimeStepSize);
          } else {
            _minTimeStepSize       = _estimatedTimeStepSize; // as we have computed the predictor with an estimate, we have to use the estimated time step size to perform the face integral
            _estimatedTimeStepSize = 0.5 * ( FusedTimeSteppingDiffusionFactor * _admissibleTimeStepSize + _estimatedTimeStepSize );
          }
        } else { // use fixed time step size in intermediate batch iterations
          _minTimeStepSize  = _estimatedTimeStepSize;
        }
      } else if ( !isLinear() ) { // non-fused, non-linear
        _minTimeStepSize = _admissibleTimeStepSize;
      } // else if linear do not change the time step size at all
    }
  }
  //rollback to solution of another replica
  else {
    _minTimeStepSize = _admissibleTimeStepSize;
    logError("wrapUpTimeStep","global master determined dt_min="<<_minTimeStepSize<<" after rollback");
  }

  if ( isLastTimeStepOfBatchOrNoBatch ) {
    std::copy(_nextGlobalObservables.begin(),_nextGlobalObservables.end(),_globalObservables.begin());
    wrapUpGlobalObservables(_globalObservables.data());
  }

  if (tarch::parallel::Node::getInstance().isGlobalMaster()
    && _meshUpdateEvent==MeshUpdateEvent::RollbackToTeamSolution) {
    logDebug("wrapUpTimeStep","Rollback activated! Global master uses estimated dt="
      <<_estimatedTimeStepSize
      <<" admissible dt"<<_admissibleTimeStepSize
      <<" min dt"<<_minTimeStepSize
      <<" min stamp"<<_minTimeStamp
      <<" for the next time step.");
  }

  // call user code
  endTimeStep(_minTimeStamp,isLastTimeStepOfBatchOrNoBatch);

  //Todo(Philipp): do this also with local recomp!! OffloadingLocalRecompute
#if defined(SharedTBB) && defined(Parallel)
  if(exahype::reactive::ReactiveContext::getInstance().getResilienceStrategy()!=exahype::reactive::ReactiveContext::ResilienceStrategy::None) {
    exahype::reactive::ResilienceStatistics::getInstance().printStatistics();
    cleanUpStaleTaskOutcomes();
  }

  if(exahype::reactive::ReactiveContext::getInstance().getResilienceStrategy()>=exahype::reactive::ReactiveContext::ResilienceStrategy::TaskSharingResilienceChecks) {
    bool isConsistent = exahype::reactive::TimeStampAndTriggerTeamHistory::getInstance().checkConsistency();
    if(!isConsistent) {
      exahype::reactive::ResilienceTools::getInstance().setCorruptionDetected(true); //todo: reduce to master?
      //MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
#endif
}

void exahype::solvers::ADERDGSolver::updateTimeStepSize() {
  if ( FuseAllADERDGPhases ) {
    _minTimeStepSize        = FusedTimeSteppingRerunFactor * _admissibleTimeStepSize;
    _estimatedTimeStepSize  = _minTimeStepSize;
  } else {
    _minTimeStepSize        = _admissibleTimeStepSize;
  }
  _admissibleTimeStepSize = std::numeric_limits<double>::infinity();
  _stabilityConditionWasViolated = false;
}

void exahype::solvers::ADERDGSolver::updateGlobalObservables() {
  std::copy(_nextGlobalObservables.begin(),_nextGlobalObservables.end(),
            _globalObservables.begin());
}

void exahype::solvers::ADERDGSolver::rollbackToPreviousTimeStep() {
  _minTimeStamp    = _previousMinTimeStamp;
  _minTimeStepSize = _previousMinTimeStepSize;

  _previousMinTimeStamp    = std::numeric_limits<double>::infinity();
  _previousMinTimeStepSize = std::numeric_limits<double>::infinity();

  // TODO(Lukas) Maybe also rollback global observables?
}

double exahype::solvers::ADERDGSolver::getMinTimeStamp() const {
  return _minTimeStamp;
}

double exahype::solvers::ADERDGSolver::getMinTimeStepSize() const {
  return _minTimeStepSize;
}

double exahype::solvers::ADERDGSolver::getEstimatedTimeStepSize() const {
  return _estimatedTimeStepSize;
}

double exahype::solvers::ADERDGSolver::getPreviousMinTimeStepSize() const {
  return _previousMinTimeStepSize;
}

double exahype::solvers::ADERDGSolver::getPreviousMinTimeStamp() const {
  return _previousMinTimeStamp;
}

double exahype::solvers::ADERDGSolver::getAdmissibleTimeStepSize() const {
  return _admissibleTimeStepSize;
}

void exahype::solvers::ADERDGSolver::updateAdmissibleTimeStepSize( double value ) {
  tarch::multicore::Lock lock(_reductionSemaphore);
  _admissibleTimeStepSize = std::min(_admissibleTimeStepSize,value);
  lock.free();
}

void exahype::solvers::ADERDGSolver::resetAdmissibleTimeStepSize() {
  _admissibleTimeStepSize = std::numeric_limits<double>::infinity();
}

void exahype::solvers::ADERDGSolver::initSolver(
    const double                                timeStamp,
    const tarch::la::Vector<DIMENSIONS,double>& domainOffset,
    const tarch::la::Vector<DIMENSIONS,double>& domainSize,
    const double                                boundingBoxSize,
    const double                                boundingBoxMeshSize,
    const std::vector<std::string>&             cmdlineargs,
    const exahype::parser::ParserView&          parserView
) {
  _domainOffset=domainOffset;
  _domainSize=domainSize;
  std::pair<double,int> coarsestMeshInfo =
      exahype::solvers::Solver::computeCoarsestMeshSizeAndLevel(
          std::min(boundingBoxMeshSize,_maximumMeshSize),boundingBoxSize);
  _coarsestMeshSize  = coarsestMeshInfo.first;
  _coarsestMeshLevel = coarsestMeshInfo.second;

  _previousMinTimeStamp = timeStamp;
  _minTimeStamp         = timeStamp;

  _previousMinTimeStepSize = 0.0;
  _minTimeStepSize         = 0.0;

  _estimatedTimeStepSize   = std::numeric_limits<double>::infinity();
  _admissibleTimeStepSize  = std::numeric_limits<double>::infinity();

  _meshUpdateEvent = MeshUpdateEvent::InitialRefinementRequested;

  // global observables
  _globalObservables.resize(_numberOfGlobalObservables);
  _nextGlobalObservables.resize(_numberOfGlobalObservables);
  resetGlobalObservables(_globalObservables.data()    );
  resetGlobalObservables(_nextGlobalObservables.data());

  init(cmdlineargs,parserView); // call user define initalisiation

  
  #ifdef Parallel
  _invalidExtrapolatedPredictor.resize(getBndFaceSize());
  _invalidFluctuations.resize(getBndFluxSize());
  std::fill_n(_invalidExtrapolatedPredictor.data(),_invalidExtrapolatedPredictor.size(),-1);
  std::fill_n(_invalidFluctuations.data(),_invalidFluctuations.size(),-1);

  _receivedExtrapolatedPredictor.resize(getBndFaceSize());
  _receivedFluctuations.resize(getBndFluxSize());

  _receivedUpdate.reserve(getUpdateSize());
  #endif
}

bool exahype::solvers::ADERDGSolver::isPerformingPrediction(
    const exahype::State::AlgorithmSection& section) const {
  bool isPerformingPrediction = false;

  switch (section) {
    case exahype::State::AlgorithmSection::TimeStepping:
      isPerformingPrediction = true;
      break;
    case exahype::State::AlgorithmSection::PredictionRerunAllSend:
      isPerformingPrediction = !hasRequestedAnyMeshRefinement() &&
                               getStabilityConditionWasViolated();
      break;
    case exahype::State::AlgorithmSection::PredictionOrLocalRecomputationAllSend:
      isPerformingPrediction = hasRequestedAnyMeshRefinement();
      break;
    default:
      break;
  }

  return isPerformingPrediction;
}

bool exahype::solvers::ADERDGSolver::isMergingMetadata(
    const exahype::State::AlgorithmSection& section) const {
  return ( exahype::State::AlgorithmSection::MeshRefinement==section ) &&
         hasRequestedAnyMeshRefinement();
}

bool exahype::solvers::ADERDGSolver::getStabilityConditionWasViolated() const {
  return  _stabilityConditionWasViolated;
}

void exahype::solvers::ADERDGSolver::setStabilityConditionWasViolated(const bool state) {
  _stabilityConditionWasViolated = state;
}

bool exahype::solvers::ADERDGSolver::isValidCellDescriptionIndex(
    const int cellDescriptionsIndex) {
  bool result = cellDescriptionsIndex>=0;
  assertion1(!result || Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  return result;
}

int exahype::solvers::ADERDGSolver::tryGetElement(
    const int cellDescriptionsIndex,
    const int solverNumber) const {
  if ( isValidCellDescriptionIndex(cellDescriptionsIndex) ) {
    int element=0;
    for (auto& p : Heap::getInstance().getData(cellDescriptionsIndex)) {
      if (p.getSolverNumber()==solverNumber) {
        return element;
      }
      ++element;
    }
  }
  return NotFound;
}

///////////////////////////////////
// CELL-LOCAL MESH REFINEMENT
///////////////////////////////////


////////////////////////////////////////
// CELL-LOCAL
////////////////////////////////////////
void exahype::solvers::ADERDGSolver::validateCellDescriptionData(
  const CellDescription& cellDescription,
  const bool validateTimeStepData,
  const bool afterCompression,
  const bool beforePrediction,
  const std::string& methodTraceOfCaller) const {
  #ifdef Asserts
  if ( _checkForNaNs && validateTimeStepData ) {
    assertion2(std::isfinite(cellDescription.getTimeStamp()),
        cellDescription.toString(),toString());
    assertion2(cellDescription.getTimeStamp()>=0,
        cellDescription.toString(),toString());
  }
  // TODO(Dominic): Remove assertion1(cellDescription.getRefinementEvent()==CellDescription::None,cellDescription.toString());
  assertion1(getType()==exahype::solvers::Solver::Type::ADERDG,cellDescription.toString());

  if ( _checkForNaNs && afterCompression) {
    // TODO(Dominic)
  } else if ( _checkForNaNs ) {
    assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex()),cellDescription.toString());
    assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getUpdateIndex()),cellDescription.toString());
    assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorIndex()),cellDescription.toString());
    assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getFluctuationIndex()),cellDescription.toString());

    double* luh  = static_cast<double*>(cellDescription.getSolution());
    double* lduh = static_cast<double*>(cellDescription.getUpdate());

    double* lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
    double* lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

    int dataPerCell             = getDataPerCell();
    int updateSize              = getUpdateSize();

    int dataPerCellBoundary     = getBndTotalSize();
    int unknownsPerCellBoundary = getBndFluxTotalSize();

    for (int i=0; i<dataPerCell; i++) {
      assertion4(tarch::la::equals(cellDescription.getTimeStepSize(),0.0) || std::isfinite(luh[i]),
          cellDescription.toString(),toString(),methodTraceOfCaller,i);
    }
    if ( !beforePrediction ) {
      for (int i=0; i<updateSize; i++) {
       assertion4(tarch::la::equals(cellDescription.getTimeStepSize(),0.0) || std::isfinite(lduh[i]),
           cellDescription.toString(),toString(),methodTraceOfCaller,i);
      }

      for (int i=0; i<dataPerCellBoundary; i++) {
        assertion4(tarch::la::equals(cellDescription.getTimeStepSize(),0.0) || std::isfinite(lQhbnd[i]),
            cellDescription.toString(),toString(),methodTraceOfCaller,i);
      }

      for (int i=0; i<unknownsPerCellBoundary; i++) {
        assertion4(tarch::la::equals(cellDescription.getTimeStepSize(),0.0) || std::isfinite(lFhbnd[i]),
            cellDescription.toString(),toString(),methodTraceOfCaller,i);
      }
    }
  }
  #endif
}

exahype::solvers::Solver::MeshUpdateEvent
exahype::solvers::ADERDGSolver::updateRefinementStatusAfterSolutionUpdate(CellDescription& cellDescription) {
  if ( OnlyInitialMeshRefinement ) {
    updateRefinementStatus(cellDescription);
    return MeshUpdateEvent::None;
  }

  cellDescription.setRefinementFlag(false);
  if ( cellDescription.getType()==CellDescription::Type::Leaf ) {
    const double* solution = static_cast<double*>(cellDescription.getSolution());
    RefinementControl refinementControl = refinementCriterion(
                      solution,cellDescription.getOffset()+0.5*cellDescription.getSize(),
                      cellDescription.getSize(),
                      cellDescription.getTimeStamp(), // must be called after advancing in time
                      cellDescription.getLevel());
    if (
        (refinementControl==RefinementControl::Refine ||
        (cellDescription.getLevel()==getMaximumAdaptiveMeshLevel()    &&
        refinementControl==RefinementControl::Keep))
    ) {
      cellDescription.setRefinementStatus(_refineOrKeepOnFineGrid);
      cellDescription.setRefinementFlag(true);
    } else if (
        cellDescription.getRefinementStatus()<_refineOrKeepOnFineGrid &&
        cellDescription.getLevel()<getMaximumAdaptiveMeshLevel()     &&
        refinementControl==RefinementControl::Keep
    ) {
      cellDescription.setRefinementStatus(Keep);
    }
    else if (
        cellDescription.getRefinementStatus()<=Keep &&
        refinementControl==RefinementControl::Erase
    ) {
      cellDescription.setRefinementStatus(Erase);
    }

    // update refinement status after prescribing refinement values
    updateRefinementStatus(cellDescription);

    return
        (cellDescription.getLevel() < getMaximumAdaptiveMeshLevel() &&
         refinementControl==RefinementControl::Refine ) ?
            MeshUpdateEvent::RefinementRequested : MeshUpdateEvent::None;
  } else if ( cellDescription.getType()==CellDescription::Type::Virtual ) {
    // bottom up refinement criterion TODO(Dominic): Add to docu
    // We allow the halo region to diffuse into the virtual subcells
    // up to some point.
    updateRefinementStatus(cellDescription);
    if (
        cellDescription.getLevel()==getMaximumAdaptiveMeshLevel() &&
        cellDescription.getRefinementStatus() >= _minimumRefinementStatusToRequestMeshRefinementInVirtualCell
    ) {
      return MeshUpdateEvent::RefinementRequested;
    } else {
      return MeshUpdateEvent::None;
    }
  } else {
    cellDescription.setRefinementStatus(Erase); // Cannot override the refinement / limiter status in other cells
    return MeshUpdateEvent::None;
  }
}

void exahype::solvers::ADERDGSolver::fusedTimeStepBody(
    CellDescription&                                   cellDescription,
    CellInfo&                                          cellInfo,
    const double                                       predictionTimeStamp,
    const double                                       predictionTimeStepSize,
    const bool                                         isFirstTimeStepOfBatch,
    const bool                                         isLastTimeStepOfBatch,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers,
    const bool                                         isSkeletonCell,
    const bool                                         mustBeDoneImmediately) {
  #ifdef USE_ITAC
  if ( isSkeletonCell ) {
    VT_begin(fusedTimeStepBodyHandleSkeleton);
  } else {
    VT_begin(fusedTimeStepBodyHandle);
  }
  #endif

  correction(cellDescription,boundaryMarkers,isFirstTimeStepOfBatch,isFirstTimeStepOfBatch/*addSurfaceIntegralContributionToUpdate*/);

  UpdateResult result;
  result._timeStepSize    = startNewTimeStep(cellDescription,isFirstTimeStepOfBatch);
  cellDescription.setPreviousRefinementStatus(cellDescription.getRefinementStatus());
  result._meshUpdateEvent = updateRefinementStatusAfterSolutionUpdate(cellDescription);

  reduce(cellDescription,result);

  if (
      SpawnPredictionAsBackgroundJob &&
      !mustBeDoneImmediately 
      && isLastTimeStepOfBatch // only spawned in last iteration if a FusedTimeStepJob was spawned before
  ) {
    const int element = cellInfo.indexOfADERDGCellDescription(cellDescription.getSolverNumber());
    //skeleton cells are not considered for offloading
    if (
       (isSkeletonCell
#if defined(Parallel)
        &&
        exahype::reactive::ReactiveContext::getInstance().getResilienceStrategy()
          < exahype::reactive::ReactiveContext::ResilienceStrategy::TaskSharingResilienceChecks
#endif
       )
#if defined(Parallel)
      || !exahype::reactive::ReactiveContext::getInstance().isEnabled()
#endif
    ) {
      peano::datatraversal::TaskSet( new PredictionJob(
        *this, cellDescription, cellInfo._cellDescriptionsIndex, element,
        predictionTimeStamp,  // corrector time step data is correct; see docu
        predictionTimeStepSize,
        false/*is uncompressed*/, isSkeletonCell, isLastTimeStepOfBatch ));
      exahype::reactive::OffloadingProfiler::getInstance().notifySpawnedTask();
    }
    else {
#if defined(SharedTBB) && defined(Parallel)
      MigratablePredictionJob *migratablePredictionJob = new MigratablePredictionJob(*this,
          cellInfo._cellDescriptionsIndex, element,
          predictionTimeStamp,
          predictionTimeStepSize, false, isSkeletonCell);
      submitOrSendMigratablePredictionJob(migratablePredictionJob);
#else
      peano::datatraversal::TaskSet( new PredictionJob(
        *this, cellDescription, cellInfo._cellDescriptionsIndex, element,
        predictionTimeStamp,  // corrector time step data is correct; see docu
        predictionTimeStepSize,
        false/*is uncompressed*/, isSkeletonCell, isLastTimeStepOfBatch ));
#endif
      exahype::reactive::OffloadingProfiler::getInstance().notifySpawnedTask();
    }
  }
  else {
    predictionAndVolumeIntegralBody(
        cellDescription,
        predictionTimeStamp, predictionTimeStepSize,
        false, isSkeletonCell, isLastTimeStepOfBatch/*addVolumeIntegralResultToUpdate*/);
  }

  #ifdef USE_ITAC
  if ( isSkeletonCell ) {
    VT_end(fusedTimeStepBodyHandleSkeleton);
  } else {
    VT_end(fusedTimeStepBodyHandle);
  }
  #endif
}

void exahype::solvers::ADERDGSolver::fusedTimeStepOrRestrict(
    const int                                          solverNumber,
    CellInfo&                                          cellInfo,
    const bool                                         isFirstTimeStepOfBatch,
    const bool                                         isLastTimeStepOfBatch,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];
    synchroniseTimeStepping(cellDescription);
    cellDescription.setHasCompletedLastStep(false);

    if ( cellDescription.getType()==CellDescription::Type::Leaf ) {
      const bool isAMRSkeletonCell     = belongsToAMRSkeleton(cellDescription);
      const bool isAtRemoteBoundary    = tarch::la::oneEquals(boundaryMarkers,exahype::mappings::LevelwiseAdjacencyBookkeeping::RemoteAdjacencyIndex);
      const bool isSkeletonCell        = isAMRSkeletonCell || isAtRemoteBoundary;
      const bool mustBeDoneImmediately = isSkeletonCell && PredictionSweeps==1;

      if (
          (SpawnUpdateAsBackgroundJob || (SpawnPredictionAsBackgroundJob && !isLastTimeStepOfBatch)) &&
          !mustBeDoneImmediately
      ) {
        const auto predictionTimeStepData = getPredictionTimeStepData(cellDescription,true);
        peano::datatraversal::TaskSet( new FusedTimeStepJob(
            *this, cellDescription, cellInfo,
            std::get<0>(predictionTimeStepData),std::get<1>(predictionTimeStepData),
            isFirstTimeStepOfBatch, isLastTimeStepOfBatch,
            boundaryMarkers, isSkeletonCell) );
      } else {
        const auto predictionTimeStepData = getPredictionTimeStepData(cellDescription,true);
        fusedTimeStepBody(
            cellDescription,cellInfo,
            std::get<0>(predictionTimeStepData),std::get<1>(predictionTimeStepData),
            isFirstTimeStepOfBatch,isLastTimeStepOfBatch,
            boundaryMarkers,isSkeletonCell,mustBeDoneImmediately );
      }
    } else if (
        cellDescription.getType()==CellDescription::Type::Virtual &&
        cellDescription.getCommunicationStatus()>=MinimumCommunicationStatusForNeighbourCommunication
    ) {
      restrictToTopMostParent(cellDescription,isFirstTimeStepOfBatch/*addToCoarseGridUpdate*/);
      updateMeshUpdateEvent( updateRefinementStatusAfterSolutionUpdate(cellDescription) );
      cellDescription.setHasCompletedLastStep(true);
    } else {
      cellDescription.setHasCompletedLastStep(true);
    }
  }
}

void exahype::solvers::ADERDGSolver::reduce(
    const CellDescription& cellDescription,
    const UpdateResult&    result) {
  updateMeshUpdateEvent(result._meshUpdateEvent);
  updateAdmissibleTimeStepSize(result._timeStepSize);

  if ( _numberOfGlobalObservables > 0 ) {
    assertion(cellDescription.getType()==CellDescription::Type::Leaf);
    const double* const luh = static_cast<double*>(cellDescription.getSolution());
    const auto cellCentre   = cellDescription.getOffset() + 0.5 * cellDescription.getSize();
    const auto& cellSize    = cellDescription.getSize();
    const auto t             = cellDescription.getTimeStamp();
    const auto dt            = cellDescription.getTimeStepSize();
    updateGlobalObservables(_nextGlobalObservables.data(),luh,cellCentre,cellSize,t,dt);
  }
}

void exahype::solvers::ADERDGSolver::updateBody(
    CellDescription&                                   cellDescription,
    CellInfo&                                          cellInfo,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers) {
  #ifdef USE_ITAC
  VT_begin(updateBodyHandle);
  #endif

  assertion1(cellDescription.getType()==CellDescription::Type::Leaf,cellDescription.toString());
  uncompress(cellDescription);

  correction(cellDescription,boundaryMarkers,true,false/*effect: add face integral result directly to solution*/);

  UpdateResult result;
  result._timeStepSize    = startNewTimeStep(cellDescription,true);
  cellDescription.setPreviousRefinementStatus(cellDescription.getRefinementStatus());
  result._meshUpdateEvent = updateRefinementStatusAfterSolutionUpdate(cellDescription);

  reduce(cellDescription,result);

  const bool isAtRemoteBoundary = tarch::la::oneEquals(boundaryMarkers,mappings::LevelwiseAdjacencyBookkeeping::RemoteAdjacencyIndex);
  compress(cellDescription,isAtRemoteBoundary);

  cellDescription.setHasCompletedLastStep(true); // required as prediction checks the flag too. Field should be renamed "setHasCompletedLastOperation(...)".

  #ifdef USE_ITAC
  VT_end(updateBodyHandle);
  #endif
}

void exahype::solvers::ADERDGSolver::updateOrRestrict(
    const int                                          solverNumber,
    CellInfo&                                          cellInfo,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];
    synchroniseTimeStepping(cellDescription);
    cellDescription.setHasCompletedLastStep(false);

    if ( cellDescription.getType()==CellDescription::Type::Leaf && SpawnUpdateAsBackgroundJob ) {
      peano::datatraversal::TaskSet ( new UpdateJob(*this,cellDescription,cellInfo,boundaryMarkers) );
    }
    else if ( cellDescription.getType()==CellDescription::Type::Leaf ) {
      updateBody(cellDescription,cellInfo,boundaryMarkers);
    }
    else if (
        cellDescription.getType()==CellDescription::Type::Virtual &&
        cellDescription.getCommunicationStatus()>=MinimumCommunicationStatusForNeighbourCommunication
    ) {
      restrictToTopMostParent(cellDescription,false/*effect: add face integral result directly to solution*/);
      updateMeshUpdateEvent( updateRefinementStatusAfterSolutionUpdate(cellDescription) );
      cellDescription.setHasCompletedLastStep(true);
    }
    else {
      cellDescription.setHasCompletedLastStep(true);
    }
  }
}

void exahype::solvers::ADERDGSolver::compress(
      const int solverNumber,
      CellInfo& cellInfo,
      const bool isAtRemoteBoundary) const {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element!=NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];
    if (cellDescription.getType()==CellDescription::Type::Leaf) {
      const bool isSkeletonCell = belongsToAMRSkeleton(cellDescription);
      compress(cellDescription,isSkeletonCell);
    }
  }
}

int exahype::solvers::ADERDGSolver::predictionAndVolumeIntegralBody(
    CellDescription& cellDescription,
    const double predictorTimeStamp,
    const double predictorTimeStepSize,
    const bool   uncompressBefore,
    const bool   isSkeletonCell,
    const bool   addVolumeIntegralResultToUpdate) {
  #ifdef USE_ITAC
  if ( isSkeletonCell ) {
    VT_begin(predictorBodyHandleSkeleton);
  } else {
    VT_begin(predictorBodyHandle);
  }
  #endif
  if (uncompressBefore) { uncompress(cellDescription); }

  validateCellDescriptionData(cellDescription,true,false,true,"exahype::solvers::ADERDGSolver::performPredictionAndVolumeIntegralBody [pre]");

  double* luh    = static_cast<double*>(cellDescription.getSolution());
  double* lduh   = static_cast<double*>(cellDescription.getUpdate());
  double* lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double* lGradQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient());
  double* lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

  #if !defined(SharedMemoryParallelisation) && !defined(Parallel) && defined(Asserts)
  static int counter = 0;
  static double timeStamp = 0;
  if ( !tarch::la::equals(timeStamp,_minTimeStamp,1e-9) ) {
    logInfo("performPredictionAndVolumeIntegralBody(...)","#predictions="<<counter);
    timeStamp = _minTimeStamp;
    counter=0;
  }
  counter++;
  #endif

  const int numberOfPicardIterations = fusedSpaceTimePredictorVolumeIntegral(
      lduh,lQhbnd,lGradQhbnd,lFhbnd,
      luh,
      cellDescription.getOffset()+0.5*cellDescription.getSize(),
      cellDescription.getSize(),
      predictorTimeStamp,
      predictorTimeStepSize,
      addVolumeIntegralResultToUpdate); // TODO(Dominic): fix 'false' case


  //if(exahype::reactive::ReactiveContext::getInstance().getResilienceStrategy()==exahype::reactive::ReactiveContext::ResilienceStrategy::None)
  exahype::reactive::ResilienceTools::getInstance().corruptDataIfActive(
                                    (cellDescription.getOffset()+0.5*cellDescription.getSize()).data(),
                                    DIMENSIONS, 
                                    predictorTimeStamp,
                                    lduh,
                                    getUpdateSize());

  compress(cellDescription,isSkeletonCell);

  validateCellDescriptionData(cellDescription,true,true,false,"exahype::solvers::ADERDGSolver::performPredictionAndVolumeIntegralBody [post]");

  cellDescription.setHasCompletedLastStep(true);
  
  #ifdef USE_ITAC
  if ( isSkeletonCell ) {
    VT_end(predictorBodyHandleSkeleton);
  } else {
    VT_end(predictorBodyHandle);
  }
  #endif

  return numberOfPicardIterations;
}

void exahype::solvers::ADERDGSolver::predictionAndVolumeIntegral(
    const int    solverNumber,
    CellInfo&    cellInfo,
    const double predictorTimeStamp,
    const double predictorTimeStepSize,
    const bool   uncompressBefore,
    const bool   isAtRemoteBoundary,
    const bool   addVolumeIntegralResultToUpdate) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];

  if ( cellDescription.getType()==CellDescription::Type::Leaf ) {
    cellDescription.setHasCompletedLastStep(false);

    const bool isAMRSkeletonCell     = belongsToAMRSkeleton(cellDescription);
    const bool isSkeletonCell        = isAMRSkeletonCell || isAtRemoteBoundary;
    const bool mustBeDoneImmediately = isSkeletonCell && PredictionSweeps==1;

    if ( SpawnPredictionAsBackgroundJob && !mustBeDoneImmediately ) {
      //skeleton cells are not considered for offloading but for task sharing with resilience checks or correction
      if (
          (isSkeletonCell
#if defined(Parallel)
          && exahype::reactive::ReactiveContext::getInstance().getResilienceStrategy()
             < exahype::reactive::ReactiveContext::ResilienceStrategy::TaskSharingResilienceChecks
#endif
          )
#if defined(Parallel)
       || !exahype::reactive::ReactiveContext::getInstance().isEnabled()
#endif
       ) {
        peano::datatraversal::TaskSet( new PredictionJob(
              *this, cellDescription, cellInfo._cellDescriptionsIndex, element,
              predictorTimeStamp,predictorTimeStepSize,
              uncompressBefore,isSkeletonCell,addVolumeIntegralResultToUpdate) );
      }
      else {
#if defined(SharedTBB) && defined(Parallel)
        MigratablePredictionJob *migratablePredictionJob = new MigratablePredictionJob(*this,
          cellInfo._cellDescriptionsIndex, element,
          predictorTimeStamp,
          predictorTimeStepSize,
          false,
          isSkeletonCell);
        submitOrSendMigratablePredictionJob(migratablePredictionJob);
        exahype::reactive::OffloadingProfiler::getInstance().notifySpawnedTask();
#else
        peano::datatraversal::TaskSet( new PredictionJob(
              *this, cellDescription, cellInfo._cellDescriptionsIndex, element,
              predictorTimeStamp,predictorTimeStepSize,
              uncompressBefore,isSkeletonCell,addVolumeIntegralResultToUpdate) );
#endif
      }
    }
    else {
      predictionAndVolumeIntegralBody(
          cellDescription,
          predictorTimeStamp,predictorTimeStepSize,
          uncompressBefore,isSkeletonCell,addVolumeIntegralResultToUpdate);
    }
  }
}

void exahype::solvers::ADERDGSolver::predictionAndVolumeIntegral(
    const int solverNumber,
    CellInfo& cellInfo,
    const bool isAtRemoteBoundary) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != Solver::NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];
    synchroniseTimeStepping(cellDescription);

    const bool isAMRSkeletonCell = belongsToAMRSkeleton(cellDescription);
    const bool isSkeletonCell    = isAMRSkeletonCell || isAtRemoteBoundary;

    waitUntilCompletedLastStep(cellDescription,isSkeletonCell,false);

    if ( cellDescription.getType()==CellDescription::Type::Leaf ) {
      const auto predictionTimeStepData = getPredictionTimeStepData(cellDescription,false); // this is either the fused scheme or a predictor recomputation

      const bool rollbacksPossible = !OnlyInitialMeshRefinement;
      if ( !FuseAllADERDGPhases && rollbacksPossible ) { // backup previous solution here as prediction already adds a contribution to solution if not all alg. phases are fused.
        std::copy_n(
            static_cast<double*>(cellDescription.getSolution()),getDataPerCell(),
            static_cast<double*>(cellDescription.getPreviousSolution()));
      }

      predictionAndVolumeIntegral(solverNumber,cellInfo,
          std::get<0>(predictionTimeStepData),
          std::get<1>(predictionTimeStepData),
          true,isAtRemoteBoundary,
          FuseAllADERDGPhases/*addVolumeIntegralResultToUpdate*/);
    }
  }
}

double exahype::solvers::ADERDGSolver::computeTimeStepSize(CellDescription& cellDescription) {
  if( cellDescription.getType()==CellDescription::Type::Leaf ) {
    const double* luh = static_cast<double*>(cellDescription.getSolution());

    validateCellDescriptionData(cellDescription,false,false,true,"computeTimeStepSizes(...)");
    double admissibleTimeStepSize = stableTimeStepSize(luh,cellDescription.getSize());

    assertion2(!_checkForNaNs || admissibleTimeStepSize>0,admissibleTimeStepSize,cellDescription.toString());

    assertion3(!_checkForNaNs || admissibleTimeStepSize<std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),admissibleTimeStepSize,cellDescription.toString());
    assertion2(!_checkForNaNs || std::isfinite(admissibleTimeStepSize),admissibleTimeStepSize,cellDescription.toString());

    return admissibleTimeStepSize;
  } else {
    return std::numeric_limits<double>::infinity();
  }
}


double exahype::solvers::ADERDGSolver::startNewTimeStep(
    CellDescription& cellDescription,
    const bool isFirstTimeStepOfBatch) {
  assertion1(cellDescription.getType()==CellDescription::Type::Leaf,cellDescription.toString());
  // n-1
  if ( isFirstTimeStepOfBatch ) {
    cellDescription.setPreviousTimeStamp(cellDescription.getTimeStamp());
    cellDescription.setPreviousTimeStepSize(cellDescription.getTimeStepSize());
  }
  // n
  cellDescription.setTimeStamp(cellDescription.getTimeStamp()+cellDescription.getTimeStepSize());
  if ( getTimeStepping() != TimeStepping::GlobalFixed ) {
    double admissibleTimeStepSize = computeTimeStepSize(cellDescription);
    cellDescription.setTimeStepSize(admissibleTimeStepSize);
    return admissibleTimeStepSize;
  } else {
    return cellDescription.getTimeStepSize();
  }
}

void exahype::solvers::ADERDGSolver::updateTimeStepSize(const int solverNumber,CellInfo& cellInfo) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];
    const double admissibleTimeStepSize = computeTimeStepSize(cellDescription);
    cellDescription.setTimeStepSize( admissibleTimeStepSize );
    cellDescription.setHasCompletedLastStep(true);
    updateAdmissibleTimeStepSize(admissibleTimeStepSize);
  }
}

void exahype::solvers::ADERDGSolver::updateGlobalObservables(const int solverNumber,CellInfo& cellInfo) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];

    if ( _numberOfGlobalObservables > 0 && cellDescription.getType() == CellDescription::Type::Leaf ) {
      const double* const luh  = static_cast<double*>(cellDescription.getSolution());
      const auto cellCentre    = cellDescription.getOffset() + 0.5 * cellDescription.getSize();
      const auto& cellSize     = cellDescription.getSize();
      const auto t             = cellDescription.getTimeStamp();
      const auto dt            = cellDescription.getTimeStepSize();
      updateGlobalObservables(_nextGlobalObservables.data(),luh,cellCentre,cellSize,t,dt);
    }
  }
}

void exahype::solvers::ADERDGSolver::rollbackToPreviousTimeStep(CellDescription& cellDescription) const {
  cellDescription.setTimeStamp   (cellDescription.getPreviousTimeStamp());
  cellDescription.setTimeStepSize(std::numeric_limits<double>::infinity());

  cellDescription.setPreviousTimeStamp(std::numeric_limits<double>::infinity());
  cellDescription.setPreviousTimeStepSize(std::numeric_limits<double>::infinity());
}

void exahype::solvers::ADERDGSolver::adjustSolution(CellDescription& cellDescription) {
  assertion1(cellDescription.getType()==CellDescription::Type::Leaf ||
             cellDescription.getType()==CellDescription::Type::LeafProlongates ||
             cellDescription.getType()==CellDescription::Type::ParentCoarsens
             ,cellDescription.toString());

  double* solution = static_cast<double*>(cellDescription.getSolution());
  adjustSolution(
      solution,
      cellDescription.getOffset()+0.5*cellDescription.getSize(),
      cellDescription.getSize(),
      cellDescription.getTimeStamp(),
      cellDescription.getTimeStepSize());

  double* previousSolution = static_cast<double*>(cellDescription.getPreviousSolution());
  adjustSolution(
      previousSolution,
      cellDescription.getOffset()+0.5*cellDescription.getSize(),
      cellDescription.getSize(),
      cellDescription.getPreviousTimeStamp(),
      cellDescription.getPreviousTimeStepSize());

  #ifdef Asserts
  if ( _checkForNaNs ) {
    for (int i=0; i<getDataPerCell(); i++) {
      assertion3(std::isfinite(solution[i]),cellDescription.toString(),"adjustSolution(...)",i);
    }
  }
  #endif
}


void exahype::solvers::ADERDGSolver::printADERDGSolution2D(const CellDescription& cellDescription)  const {
  #if DIMENSIONS==2
  double* solution = static_cast<double*>(cellDescription.getSolution());
  assertion1(cellDescription.getType()==CellDescription::Type::Leaf,cellDescription.toString());
  assertion1(solution!=nullptr,cellDescription.toString());

  std::cout <<  "solution:" << std::endl;
  const int numberOfData = _numberOfVariables + _numberOfParameters;
  for (int unknown=0; unknown < numberOfData; unknown++) {
    std::cout <<  "unknown=" << unknown << std::endl;
    dfor(i,_nodesPerCoordinateAxis) {
      int iScalar = peano::utils::dLinearisedWithoutLookup(i,_nodesPerCoordinateAxis)*numberOfData+unknown;
      std::cout << std::setprecision(3) << solution[iScalar] << ",";
      if (i(0)==_nodesPerCoordinateAxis-1) {
        std::cout << std::endl;
      }
    }
  }
  std::cout <<  "}" << std::endl;
  #endif
}

void exahype::solvers::ADERDGSolver::printADERDGExtrapolatedPredictor2D(const CellDescription& cellDescription) const {
  #if DIMENSIONS==2
  double* extrapolatedPredictor = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  assertion1(extrapolatedPredictor!=nullptr,cellDescription.toString());

  std::cout <<  "extrapolated predictor:" << std::endl;
  const int numberOfData = _numberOfVariables + _numberOfParameters;
  for (int f=0; f<DIMENSIONS_TIMES_TWO; f++) {
    std::cout <<  "face=" << f << std::endl;
    for (int unknown=0; unknown < numberOfData; unknown++) {
      std::cout <<  "unknown=" << unknown << std::endl;
      for (int i=0; i<_nodesPerCoordinateAxis; i++) {
        int iScalar = f*numberOfData*_nodesPerCoordinateAxis + numberOfData*i + unknown;
        std::cout << std::setprecision(3) << extrapolatedPredictor[iScalar] << ",";
        if (i==_nodesPerCoordinateAxis-1) {
          std::cout << std::endl;
        }
      }
    }
  }
  std::cout <<  "}" << std::endl;
  #endif
}

void exahype::solvers::ADERDGSolver::printADERDGFluctuations2D(const CellDescription& cellDescription) const {
  #if DIMENSIONS==2
  double* fluctuation = static_cast<double*>(cellDescription.getFluctuation());
  assertion1(fluctuation!=nullptr,cellDescription.toString());

  std::cout << "fluctuations:" << std::endl;
  const int numberOfData = _numberOfVariables + _numberOfParameters;
  for (int f=0; f<DIMENSIONS_TIMES_TWO; f++) {
    std::cout <<  "face=" << f << std::endl;
    for (int unknown=0; unknown < numberOfData; unknown++) {
      std::cout <<  "unknown=" << unknown << std::endl;
      for (int i=0; i<_nodesPerCoordinateAxis; i++) {
        int iScalar = f*numberOfData*_nodesPerCoordinateAxis + numberOfData*i + unknown;
        std::cout << std::setprecision(3) << fluctuation[iScalar] << ",";
        if (i==_nodesPerCoordinateAxis-1) {
          std::cout << std::endl;
        }
      }
    }
  }
  std::cout <<  "}" << std::endl;
  #endif
}

void exahype::solvers::ADERDGSolver::surfaceIntegral(
    CellDescription&                                   cellDescription,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers,
    const bool                                         addToUpdate) {
  assertion1(cellDescription.getType()==CellDescription::Type::Leaf,cellDescription.toString());
  double* output= static_cast<double*>(cellDescription.getSolution());
  if ( addToUpdate ) { output = static_cast<double*>(cellDescription.getUpdate()); }

  #ifdef Asserts
  if ( _checkForNaNs ) {
    for (int i=0; i<getUnknownsPerCell(); i++) { // update does not store parameters
      assertion3(tarch::la::equals(cellDescription.getTimeStepSize(),0.0)  || std::isfinite(output[i]),cellDescription.toString(),"surfaceIntegral",i);
    }
  }
  #endif

  // the actual operation
  const int dofsPerFace = getBndFluxSize();
  for (int direction=0; direction<DIMENSIONS; direction++) {
    for (int orientation=0; orientation<2; orientation++) {
      const int faceIndex=2*direction+orientation;
      // impose boundary conditions
      if ( boundaryMarkers[faceIndex]==mappings::LevelwiseAdjacencyBookkeeping::DomainBoundaryAdjacencyIndex ) {
        mergeWithBoundaryData(cellDescription,faceIndex,direction,orientation);
      }
      // perform face integral
      if ( cellDescription.getFacewiseAugmentationStatus(faceIndex)<MaximumAugmentationStatus ) { // ignore Ancestors
        double* const lFhbnd = static_cast<double*>(cellDescription.getFluctuation()) + dofsPerFace * faceIndex;
        faceIntegral(output,lFhbnd,direction,orientation,0/*implicit conversion*/,0,cellDescription.getSize(),
                     cellDescription.getTimeStepSize(),addToUpdate);
      }
    }
  }

  // check that all neighbour merges have been performed; reset the flags
  if ( !tarch::la::equals(cellDescription.getNeighbourMergePerformed(),static_cast<signed char>(true)) && !ProfileUpdate ) {
    logError("surfaceIntegral(...)","Riemann solve was not performed on all faces of cell= "<<cellDescription.toString());
    std::terminate();
  }
  assertion1( tarch::la::equals(cellDescription.getNeighbourMergePerformed(),static_cast<signed char>(true)) || ProfileUpdate,cellDescription.toString());
  cellDescription.setNeighbourMergePerformed(static_cast<signed char>(false));

  #ifdef Asserts
  if ( _checkForNaNs ) {
    for (int i=0; i<getUnknownsPerCell(); i++) { // update does not store parameters
      assertion3(std::isfinite(output[i]),cellDescription.toString(),"updateSolution(...)",i);
    }
  }
  #endif
}

void exahype::solvers::ADERDGSolver::addUpdateToSolution(CellDescription& cellDescription,const bool backupPreviousSolution) {
  assertion1(cellDescription.getType()==CellDescription::Type::Leaf,cellDescription.toString());
  double* update = static_cast<double*>(cellDescription.getUpdate());

  if ( backupPreviousSolution ) {
    // perform the update
    swapSolutionAndPreviousSolution(cellDescription); // solution is overwritten with the current solution plus the update,
    // while current solution is remembered as the previous solution.
    double* solution         = static_cast<double*>(cellDescription.getSolution());
    double* previousSolution = static_cast<double*>(cellDescription.getPreviousSolution());
    addUpdateToSolution(solution,previousSolution,update,cellDescription.getTimeStepSize());
  } else {
    double* solution = static_cast<double*>(cellDescription.getSolution());
    addUpdateToSolution(solution,solution,update,cellDescription.getTimeStepSize());
  }
}

void exahype::solvers::ADERDGSolver::adjustSolutionAfterUpdate(CellDescription& cellDescription) {
  double* solution = static_cast<double*>(cellDescription.getSolution());
  adjustSolution(
      solution,
      cellDescription.getOffset()+0.5*cellDescription.getSize(),
      cellDescription.getSize(),
      cellDescription.getTimeStamp()+cellDescription.getTimeStepSize(),
      cellDescription.getTimeStepSize());

  // only for profiling
  if ( Solver::ProfileUpdate ) { swapSolutionAndPreviousSolution(cellDescription); }

  #ifdef Asserts
  if ( _checkForNaNs ) {
    for (int i=0; i<getUnknownsPerCell(); i++) { // update does not store parameters
      assertion3(std::isfinite(solution[i]),cellDescription.toString(),"updateSolution(...)",i);
    }
  }
  #endif
}

void exahype::solvers::ADERDGSolver::correction(
    CellDescription&                                   cellDescription,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& boundaryMarkers,
    const bool                                         backupPreviousSolution,
    const bool                                         addSurfaceIntegralResultToUpdate) {
  assertion1(cellDescription.getType()==CellDescription::Type::Leaf,cellDescription.toString())
  assertion1(std::isfinite(cellDescription.getTimeStamp()   ),cellDescription.toString());
  assertion1(std::isfinite(cellDescription.getTimeStepSize()),cellDescription.toString());
  #if !defined(SharedMemoryParallelisation) && !defined(Parallel) && defined(Asserts)
  static int counter = 0;
  static double timeStamp = 0;
  if ( !tarch::la::equals(timeStamp,_minTimeStamp,1e-9) ) {
    logInfo("mergeNeighboursData(...)","#updateSolution="<<counter);
    timeStamp = _minTimeStamp;
    counter=0;
  }
  counter++;
  #endif

  surfaceIntegral(cellDescription,boundaryMarkers,addSurfaceIntegralResultToUpdate);
  if ( addSurfaceIntegralResultToUpdate ) {
    addUpdateToSolution(cellDescription,backupPreviousSolution);
  }
  adjustSolutionAfterUpdate(cellDescription);

  // update grid flags
  updateCommunicationStatus(cellDescription);
  updateAugmentationStatus(cellDescription);
}

void exahype::solvers::ADERDGSolver::swapSolutionAndPreviousSolution(CellDescription& cellDescription) const {
  assertion(cellDescription.getType()==CellDescription::Type::Leaf);

  // Simply swap the heap indices
  const int previousSolutionIndex = cellDescription.getPreviousSolutionIndex();
  void* previousSolution          = cellDescription.getPreviousSolution(); // pointer
  cellDescription.setPreviousSolutionIndex(cellDescription.getSolutionIndex());
  cellDescription.setPreviousSolution(cellDescription.getSolution());
  cellDescription.setSolutionIndex(previousSolutionIndex);
  cellDescription.setSolution(previousSolution);
}

void exahype::solvers::ADERDGSolver::prolongateObservablesMinAndMax(
    const CellDescription& cellDescription,
    const CellDescription& parentCellDescription) const {
  const int numberOfObservables = getDMPObservables();
  for (int faceIndex = 0; faceIndex < DIMENSIONS_TIMES_TWO; ++faceIndex) {
    if ( cellDescription.getFacewiseCommunicationStatus(faceIndex)==LeafCommunicationStatus ) { // TODO(Dominic): If the grid changes dynamically during the time steps,
      // fine
      double* minFine = static_cast<double*>(cellDescription.getSolutionMin()) + numberOfObservables * faceIndex;
      double* maxFine = static_cast<double*>(cellDescription.getSolutionMax()) + numberOfObservables * faceIndex;
      // coarse
      const double* minCoarse = static_cast<double*>(parentCellDescription.getSolutionMin()) +  numberOfObservables * faceIndex;
      const double* maxCoarse = static_cast<double*>(parentCellDescription.getSolutionMax()) +  numberOfObservables * faceIndex;

      std::copy_n( minCoarse,numberOfObservables, minFine );
      std::copy_n( maxCoarse,numberOfObservables, maxFine );
    }
  }
}

void exahype::solvers::ADERDGSolver::prolongateFaceDataToVirtualCell(
    CellDescription& cellDescription,
    const CellDescription& parentCellDescription,
    const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) {
  #ifdef USE_ITAC
  VT_begin(prolongateFaceDataToVirtualCellHandle);
  #endif

  assertion(parentCellDescription.getSolverNumber() == cellDescription.getSolverNumber());
  assertion(parentCellDescription.getType() == CellDescription::Type::Leaf ||
            parentCellDescription.getType() == CellDescription::Type::Virtual);

  const int levelFine   = cellDescription.getLevel();
  const int levelCoarse = parentCellDescription.getLevel();
  assertion(levelCoarse < levelFine);

  double* update = static_cast<double*>(cellDescription.getUpdate());
  std::fill_n(update,getUpdateSize(),0.0);

  for (int faceIndex = 0; faceIndex < DIMENSIONS_TIMES_TWO; ++faceIndex) {
    const int direction = faceIndex/2;
    if ( cellDescription.getFacewiseCommunicationStatus(faceIndex)==LeafCommunicationStatus ) { // TODO(Dominic): If the grid changes dynamically during the time steps,
      #ifdef Parallel
      assertion4( exahype::amr::faceIsOnBoundaryOfParent(faceIndex,subcellIndex,levelFine-levelCoarse),
            cellDescription.toString(),parentCellDescription.toString(),tarch::parallel::Node::getInstance().getRank(),
            tarch::parallel::NodePool::getInstance().getMasterRank() ); // necessary but not sufficient
      #else
      assertion2( exahype::amr::faceIsOnBoundaryOfParent(faceIndex,subcellIndex,levelFine-levelCoarse),
                  cellDescription.toString(),parentCellDescription.toString() ); // necessary but not sufficient
      #endif

      logDebug("prolongateFaceDataToDescendant(...)","cell=" << cellDescription.getOffset() <<
               ",level=" << cellDescription.getLevel() <<
               ",face=" << faceIndex <<
               ",subcellIndex" << subcellIndex.toString() <<
               " from " <<
               " parentCell="<<parentCellDescription.getOffset()<<
               " level="<<parentCellDescription.getLevel());

      // extrapolated predictor and flux interpolation
      // extrapolated predictor
      assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorIndex()),cellDescription.toString());
      assertion1(DataHeap::getInstance().isValidIndex(parentCellDescription.getExtrapolatedPredictorIndex()),parentCellDescription.toString());

      const int dataPerFace = getBndFaceSize();
      const int dofsPerFace = getBndFluxSize();

      // fine
      double* lQhbndFine = static_cast<double*>(cellDescription.getExtrapolatedPredictor()) + dataPerFace * faceIndex; // TODO(Dominic ): Pointer should be obtained only once
      double* lFhbndFine = static_cast<double*>(cellDescription.getFluctuation())           + dofsPerFace  * faceIndex ;
      // coarse
      const double* lQhbndCoarse = static_cast<double*>(parentCellDescription.getExtrapolatedPredictor()) + dataPerFace * faceIndex;
      const double* lFhbndCoarse = static_cast<double*>(parentCellDescription.getFluctuation()          ) + dofsPerFace  * faceIndex;

      faceUnknownsProlongation(lQhbndFine,lFhbndFine,lQhbndCoarse,lFhbndCoarse, levelCoarse, levelFine,
                               exahype::amr::getSubfaceIndex(subcellIndex,direction));
    }
  }

  if ( getDMPObservables()>0 ) {
    prolongateObservablesMinAndMax(cellDescription,parentCellDescription);
  }

  cellDescription.setHasCompletedLastStep(true);

  #ifdef USE_ITAC
  VT_end(prolongateFaceDataToVirtualCellHandle);
  #endif
}

void exahype::solvers::ADERDGSolver::prolongateFaceData(
    const int solverNumber,
    CellInfo& cellInfo,
    const bool isAtRemoteBoundary) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);

  if ( element != Solver::NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];

    if (
        cellDescription.getType()==CellDescription::Type::Virtual &&
        cellDescription.getCommunicationStatus()>=MinimumCommunicationStatusForNeighbourCommunication &&
        isValidCellDescriptionIndex(cellDescription.getParentIndex()) // might be at master-worker boundary
    ) {
        Solver::SubcellPosition subcellPosition = amr::computeSubcellPositionOfVirtualCell<CellDescription,Heap>(cellDescription);
        CellDescription& parentCellDescription = getCellDescription(
            subcellPosition.parentCellDescriptionsIndex,subcellPosition.parentElement);
        assertion1(parentCellDescription.getType()==CellDescription::Type::Leaf,parentCellDescription.toString());
        
        waitUntilCompletedLastStep<CellDescription>(parentCellDescription,true,false); // TODO(Dominic): We wait for skeleton jobs here. It might make sense to receiveDanglingMessages here too

        if (
            !SpawnProlongationAsBackgroundJob ||
            isAtRemoteBoundary
        ) {
          prolongateFaceDataToVirtualCell(cellDescription,parentCellDescription,subcellPosition.subcellIndex);
        } else {
          cellDescription.setHasCompletedLastStep(false); // done here in order to skip lookup of cell description in job constructor
          peano::datatraversal::TaskSet spawn( new ProlongationJob( *this, 
              cellDescription, parentCellDescription, subcellPosition.subcellIndex) );
        }
      }
      assertion2(
          cellDescription.getType()!=CellDescription::Type::Virtual ||
          isValidCellDescriptionIndex(cellDescription.getParentIndex()),
          cellDescription.toString(),
          tarch::parallel::Node::getInstance().getRank());
  }
}

void exahype::solvers::ADERDGSolver::restrictObservablesMinAndMax(
    const CellDescription& cellDescription,
    const CellDescription& parentCellDescription) const {
  const int numberOfObservables = getDMPObservables();
  for (int faceIndex=0; faceIndex<DIMENSIONS_TIMES_TWO; faceIndex++) {
    if ( cellDescription.getFacewiseCommunicationStatus(faceIndex)==LeafCommunicationStatus ) {
      // fine
      const double* minFine = static_cast<double*>(cellDescription.getSolutionMin()) + numberOfObservables* faceIndex;
      const double* maxFine = static_cast<double*>(cellDescription.getSolutionMax()) + numberOfObservables* faceIndex;
      // coarse
      double* minCoarse = static_cast<double*>(parentCellDescription.getSolutionMin()) + numberOfObservables * faceIndex;
      double* maxCoarse = static_cast<double*>(parentCellDescription.getSolutionMax()) + numberOfObservables * faceIndex;

      tarch::multicore::Lock lock(RestrictionSemaphore);
      for (int i=0; i<numberOfObservables; i++) {
        *(minCoarse+i) = std::min( *(minFine+i), *(minCoarse+i) );
        *(maxCoarse+i) = std::max( *(maxFine+i), *(maxCoarse+i) );
      }
      lock.free();
    }
  }
}

void exahype::solvers::ADERDGSolver::restrictToTopMostParent(
    CellDescription& cellDescription,
    const bool       addToCoarseGridUpdate) {
  #ifdef USE_ITAC
  VT_begin(restrictToTopMostParentHandle);
  #endif

  // validate and obtain parent
  assertion1( cellDescription.getType()==CellDescription::Type::Virtual &&
              cellDescription.getCommunicationStatus()>=MinimumCommunicationStatusForNeighbourCommunication, cellDescription.toString() );
  assertion1( tryGetElement(cellDescription.getParentIndex(),cellDescription.getSolverNumber()) != NotFound, cellDescription.toString());
  exahype::solvers::Solver::SubcellPosition subcellPosition =
      exahype::amr::computeSubcellPositionOfVirtualCell<CellDescription,Heap>(cellDescription);
  assertion1(subcellPosition.parentElement!=exahype::solvers::Solver::NotFound,cellDescription.toString());

  CellDescription& parentCellDescription =
      getCellDescription(subcellPosition.parentCellDescriptionsIndex,subcellPosition.parentElement);

  assertion(parentCellDescription.getSolverNumber()==cellDescription.getSolverNumber());
  assertion1(cellDescription.getType()==CellDescription::Type::Virtual,cellDescription.toString());
  assertion1(parentCellDescription.getType()==CellDescription::Type::Leaf,parentCellDescription.toString());


  // Perform the face integrals on fine grid cell
  double* updateFine = static_cast<double*>(cellDescription.getUpdate()); // TODO(Dominic): Can be temporary
  std::fill_n(updateFine,getUpdateSize(),0.0);

  const int levelDelta = cellDescription.getLevel() - parentCellDescription.getLevel();
  const tarch::la::Vector<DIMENSIONS,int> subcellIndex =
      exahype::amr::computeSubcellIndex(
          cellDescription.getOffset(),cellDescription.getSize(),
          parentCellDescription.getOffset()); // TODO(Dominic): Maybe, I get can get rid of some variables again

  // gather contributions
  double dt = parentCellDescription.getTimeStepSize();
  switch (getTimeStepping()) {
    case TimeStepping::Global:
    case TimeStepping::GlobalFixed:
      // do nothing
      break;
    default:
      logError("restrictToTopMostParent(...)","Time stepping scheme not supported.");
      break;
  }

  const int dofsPerFace = getBndFluxSize();
  for (int faceIndex=0; faceIndex<DIMENSIONS_TIMES_TWO; faceIndex++) {
    const int direction   = faceIndex / 2;
    const int orientation = faceIndex % 2;
    const tarch::la::Vector<DIMENSIONS-1,int> subfaceIndex =
      exahype::amr::getSubfaceIndex(subcellIndex,direction);
    if ( cellDescription.getFacewiseCommunicationStatus(faceIndex)==LeafCommunicationStatus ) {
      assertion1(exahype::amr::faceIsOnBoundaryOfParent(faceIndex,subcellIndex,levelDelta),cellDescription.toString());
      assertion1(SpawnProlongationAsBackgroundJob || cellDescription.getNeighbourMergePerformed(faceIndex),cellDescription.toString());// necessary but not sufficient

      double* const lFhbnd = static_cast<double*>(cellDescription.getFluctuation()) + dofsPerFace * faceIndex;
      faceIntegral(updateFine,lFhbnd,direction,orientation,subfaceIndex,levelDelta,cellDescription.getSize(),dt,addToCoarseGridUpdate);

      logDebug("restrictToTopMostParent(...)","cell=" << cellDescription.getOffset() <<
             ",level=" << cellDescription.getLevel() <<
             ",face=" << faceIndex <<
             ",subcellIndex" << subcellIndex.toString() <<
             " to " <<
             " parentCell="<<parentCellDescription.getOffset()<<
             " level="<<parentCellDescription.getLevel());
    }
  }
  cellDescription.setNeighbourMergePerformed(static_cast<signed char>(false));

  // Add child contributions to parent
  if ( addToCoarseGridUpdate ) {
    double* updateCoarse = static_cast<double*>(parentCellDescription.getUpdate());
    tarch::multicore::Lock lock(RestrictionSemaphore);
    for (int i = 0; i < getUpdateSize(); ++i) { // TODO(Dominic): There must be a correction somewhere here or to the face integral if LTS is performed
      updateCoarse[i] += updateFine[i];
    }
    lock.free();
  } else {
    double* solutionCoarse = static_cast<double*>(parentCellDescription.getSolution());
    tarch::multicore::Lock lock(RestrictionSemaphore);
    addUpdateToSolution(solutionCoarse,solutionCoarse,updateFine,dt);
    lock.free();
  }

  if ( getDMPObservables()>0 ) {
    restrictObservablesMinAndMax(cellDescription,parentCellDescription);
  }

  #ifdef USE_ITAC
  VT_end(restrictToTopMostParentHandle);
  #endif
}

void exahype::solvers::ADERDGSolver::disableCheckForNaNs() {
  _checkForNaNs = false;
}

// TODO(Dominic): Check that we have rolled back in time as well
void exahype::solvers::ADERDGSolver::rollbackSolutionGlobally(const int solverNumber,CellInfo& cellInfo) const {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];

    // 1. Rollback time step data
    rollbackToPreviousTimeStep(cellDescription);
    // 2. Rollback solution to previous one
    if (cellDescription.getType()==CellDescription::Type::Leaf) {
      swapSolutionAndPreviousSolution(cellDescription);
    }

    // 3. Reset the previous refinement status
    cellDescription.setRefinementStatus(cellDescription.getPreviousRefinementStatus());
  }
}

bool exahype::solvers::ADERDGSolver::updateYieldsPhysicallyAdmissibleSolution(CellDescription& cellDescription) {
  double *tmp_sol = new double[getDataPerCell()];
  double *update = static_cast<double *> (cellDescription.getUpdate());

  double* observablesMin = nullptr;
  double* observablesMax = nullptr;

  const int numberOfObservables = getDMPObservables();
  if (numberOfObservables > 0) {
    observablesMin = static_cast<double*>(cellDescription.getSolutionMin());
    observablesMax = static_cast<double*>(cellDescription.getSolutionMax());
  }

  std::memcpy(tmp_sol, cellDescription.getSolution(), getDataPerCell()*sizeof(double));
  addUpdateToSolution(tmp_sol, tmp_sol, update, cellDescription.getTimeStepSize());

  bool isAdmissible = isPhysicallyAdmissible(
          tmp_sol,
          observablesMin,observablesMax,
          cellDescription.getRefinementStatus()>=_minRefinementStatusForTroubledCell,
          cellDescription.getOffset()+0.5*cellDescription.getSize(),cellDescription.getSize(),
          cellDescription.getTimeStamp());

  delete[] tmp_sol;
  return isAdmissible;
}

void exahype::solvers::ADERDGSolver::mergeWithNeighbourMetadata(
    const int                                    solverNumber,
    Solver::CellInfo&                            cellInfo, // corresponds to dest
    const int                                    neighbourAugmentationStatus,
    const int                                    neighbourCommunicationStatus,
    const int                                    neighbourRefinementStatus,
    //const int                                    neighbourCorruptionStatus,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    const tarch::la::Vector<DIMENSIONS, int>&    posNeighbour,
    const tarch::la::Vector<DIMENSIONS, double>& barycentreFromVertex) {
  assertion(tarch::la::countEqualEntries(posNeighbour,pos)==DIMENSIONS-1); // only consider faces
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element!=Solver::NotFound ) {
    Solver::BoundaryFaceInfo face(pos,posNeighbour); //
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];

    const tarch::la::Vector<DIMENSIONS,double> barycentreFrom1 =
        exahype::Cell::computeFaceBarycentre(cellDescription.getOffset(),cellDescription.getSize(),face._direction,face._orientation);
    if ( Vertex::equalUpToRelativeTolerance(barycentreFrom1,barycentreFromVertex) ) {
      logDebug("mergeWithNeighbourMetadata(...)", "received neighbour metadata="<<neighbourAugmentationStatus<<","<<neighbourCommunicationStatus<<","<<neighbourRefinementStatus);

      mergeWithAugmentationStatus (cellDescription,face._faceIndex,neighbourAugmentationStatus );
      mergeWithCommunicationStatus(cellDescription,face._faceIndex,neighbourCommunicationStatus);
      mergeWithRefinementStatus   (cellDescription,face._faceIndex,neighbourRefinementStatus   );
      //mergeWithCorruptionStatus   (cellDescription,face._faceIndex,neighbourCorruptionStatus   );

      cellDescription.setNeighbourMergePerformed(face._faceIndex,true);
    }
  }
}

void exahype::solvers::ADERDGSolver::mergeNeighboursMetadata(
    const int                                    solverNumber,
    Solver::CellInfo&                            cellInfo1,
    Solver::CellInfo&                            cellInfo2,
    const tarch::la::Vector<DIMENSIONS, int>&    pos1,
    const tarch::la::Vector<DIMENSIONS, int>&    pos2,
    const tarch::la::Vector<DIMENSIONS, double>& barycentreFromVertex,
    const bool                                   checkThoroughly) const {
  const int element1 = cellInfo1.indexOfADERDGCellDescription(solverNumber);
  const int element2 = cellInfo2.indexOfADERDGCellDescription(solverNumber);
  if ( element1 != Solver::NotFound && element2 != Solver::NotFound ) {
    Solver::InterfaceInfo face(pos1,pos2);
    CellDescription& cellDescription1 = cellInfo1._ADERDGCellDescriptions[element1];
    CellDescription& cellDescription2 = cellInfo2._ADERDGCellDescriptions[element2];

    bool mergeMetadata = true;
    if ( checkThoroughly ) {
      const tarch::la::Vector<DIMENSIONS,double> barycentreFrom1 =
          exahype::Cell::computeFaceBarycentre(cellDescription1.getOffset(),cellDescription1.getSize(),face._direction,face._orientation1);
      const tarch::la::Vector<DIMENSIONS,double> barycentreFrom2 =
          exahype::Cell::computeFaceBarycentre(cellDescription2.getOffset(),cellDescription2.getSize(),face._direction,face._orientation2);
      mergeMetadata &= Vertex::equalUpToRelativeTolerance(barycentreFrom1,barycentreFromVertex) &&
                       Vertex::equalUpToRelativeTolerance(barycentreFrom2,barycentreFromVertex);
    }
    if ( mergeMetadata ) {
      mergeWithCommunicationStatus(cellDescription1,face._faceIndex1,cellDescription2.getCommunicationStatus());
      mergeWithAugmentationStatus( cellDescription1,face._faceIndex1,cellDescription2.getAugmentationStatus());
      mergeWithRefinementStatus(   cellDescription1,face._faceIndex1,cellDescription2.getRefinementStatus());
      //mergeWithCorruptionStatus(   cellDescription1,face._faceIndex1,cellDescription2.getCorruptionStatus());

      mergeWithCommunicationStatus(cellDescription2,face._faceIndex2,cellDescription1.getCommunicationStatus());
      mergeWithAugmentationStatus( cellDescription2,face._faceIndex2,cellDescription1.getAugmentationStatus());
      mergeWithRefinementStatus(   cellDescription2,face._faceIndex2,cellDescription1.getRefinementStatus());
      //mergeWithCorruptionStatus(   cellDescription2,face._faceIndex2,cellDescription1.getCorruptionStatus());
    } else {
      mergeWithCommunicationStatus(cellDescription1,face._faceIndex1,EmptyStatus);
      mergeWithAugmentationStatus( cellDescription1,face._faceIndex1,EmptyStatus);
      mergeWithRefinementStatus(   cellDescription1,face._faceIndex1,EmptyStatus);
      //mergeWithCorruptionStatus(   cellDescription1,face._faceIndex1, Uncorrupted); //todo no hardcoding

      mergeWithCommunicationStatus(cellDescription2,face._faceIndex2,EmptyStatus);
      mergeWithAugmentationStatus( cellDescription2,face._faceIndex2,EmptyStatus);
      mergeWithRefinementStatus(   cellDescription2,face._faceIndex2,EmptyStatus);
      //mergeWithCorruptionStatus(   cellDescription2,face._faceIndex1, Uncorrupted); //todo no hardcoding
    }

    cellDescription1.setNeighbourMergePerformed(face._faceIndex1,true); // here we only set, doesn't matter if operation is done twice.
    cellDescription2.setNeighbourMergePerformed(face._faceIndex2,true);
  }
}

// merge compute data
void exahype::solvers::ADERDGSolver::mergeNeighboursData(
    const int                                 solverNumber,
    Solver::CellInfo&                         cellInfo1,
    Solver::CellInfo&                         cellInfo2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2) {
  #ifdef USE_ITAC
  VT_begin(mergeNeighboursHandle);
  #endif

  const int element1 = cellInfo1.indexOfADERDGCellDescription(solverNumber);
  const int element2 = cellInfo2.indexOfADERDGCellDescription(solverNumber);
  if ( element1 != Solver::NotFound && element2 != Solver::NotFound ) {
    Solver::InterfaceInfo face(pos1,pos2);
    CellDescription& cellDescription1 = cellInfo1._ADERDGCellDescriptions[element1];
    CellDescription& cellDescription2 = cellInfo2._ADERDGCellDescriptions[element2];

    if ( ADERDGSolver::communicateWithNeighbour(cellDescription1,face._faceIndex1) ) {
      assertion1( ADERDGSolver::communicateWithNeighbour(cellDescription2,face._faceIndex2),cellDescription2.toString() );

      prefetchFaceData(cellDescription1,face._faceIndex1);
      prefetchFaceData(cellDescription2,face._faceIndex2);

      #if !defined(SharedMemoryParallelisation) && !defined(Parallel) && defined(Asserts)
      static int counter = 0;
      static double timeStamp = 0;
      if ( !tarch::la::equals(timeStamp,_minTimeStamp,1e-9) ) {
        logInfo("mergeNeighboursData(...)","#riemanns="<<counter);
        timeStamp = _minTimeStamp;
        counter=0;
      }
      counter++;
      #endif

      waitUntilCompletedLastStep<CellDescription>(cellDescription1,false,false);  // must be done before any other operation on the patches
      waitUntilCompletedLastStep<CellDescription>(cellDescription2,false,false);

      if ( CompressionAccuracy > 0.0 ) {
        peano::datatraversal::TaskSet uncompression(
            [&] () -> bool {
          uncompress(cellDescription1);
          return false;
        },
        [&] () -> bool {
          uncompress(cellDescription1);
          return false;
        },
        peano::datatraversal::TaskSet::TaskType::Background,
        peano::datatraversal::TaskSet::TaskType::Background,
        true
        );
      }

      //
      // 1. Solve Riemann problem (merge data)
      //
      solveRiemannProblemAtInterface(cellDescription1,cellDescription2,face);
    }

     cellDescription1.setNeighbourMergePerformed(face._faceIndex1,true);
     cellDescription2.setNeighbourMergePerformed(face._faceIndex2,true);
  }

  #ifdef USE_ITAC
  VT_end(mergeNeighboursHandle);
  #endif
}

std::string exahype::solvers::ADERDGSolver::riemannDataToString(
    const double* const Q,const double* const F,std::string suffix) const {
  const int numberOfData = _numberOfVariables + _numberOfParameters;
  const int nodesPerFace = getBndFaceSize()/numberOfData;

  std::ostringstream stream;
  stream << std::endl;
  stream << "riemann states ("<<suffix<<"):" << std::endl;
  for (int i = 0; i<numberOfData; i++) {
    double minQ = std::numeric_limits<double>::infinity();
    double maxQ = -std::numeric_limits<double>::infinity();
    for(int n=0; n<nodesPerFace; ++n) {
      minQ = std::min(Q[n*numberOfData+i],minQ);
      maxQ = std::max(Q[n*numberOfData+i],maxQ);
    }
    stream << "Q"<<suffix<<"["<<i<<"] in ["<<std::setprecision(2)<<minQ<<","<<std::setprecision(2)<<maxQ<<"], ";
    stream << std::endl;
  }
  stream << "riemann fluxes ("<<suffix<<"):" << std::endl;
  for (int i = 0; i<_numberOfVariables; i++) {
    double minF = std::numeric_limits<double>::infinity();
    double maxF = -std::numeric_limits<double>::infinity();
    for(int n=0; n<nodesPerFace; ++n) {
      minF = std::min(F[n*_numberOfVariables+i],minF);
      maxF = std::max(F[n*_numberOfVariables+i],maxF);
    }
    stream << "F"<<suffix<<"["<<i<<"] in  ["<<std::setprecision(2)<<minF<<","<<std::setprecision(2)<<maxF<<"], ";
    stream << std::endl;
  }
  return stream.str();
}

void exahype::solvers::ADERDGSolver::solveRiemannProblemAtInterface(
    CellDescription& cellDescription1,
    CellDescription& cellDescription2,
    Solver::InterfaceInfo& face) {
  CellDescription& pLeft  =
      (face._orientation1==1) ? cellDescription1 : cellDescription2;
  CellDescription& pRight =
      (face._orientation1==1) ? cellDescription2 : cellDescription1;

  assertion2(DataHeap::getInstance().isValidIndex(pLeft.getExtrapolatedPredictorIndex()),pLeft.toString(),pRight.toString());
  assertion2(DataHeap::getInstance().isValidIndex(pLeft.getFluctuationIndex()),pLeft.toString(),pRight.toString());
  assertion2(DataHeap::getInstance().isValidIndex(pRight.getExtrapolatedPredictorIndex()),pLeft.toString(),pRight.toString());
  assertion2(DataHeap::getInstance().isValidIndex(pRight.getFluctuationIndex()),pLeft.toString(),pRight.toString());

  const int dataPerFace = getBndFaceSize();
  const int dofsPerFace  = getBndFluxSize();

  double* QL = static_cast<double*>(pLeft.getExtrapolatedPredictor()) +  dataPerFace * face._faceIndexLeft;
  double* FL = static_cast<double*>(pLeft.getFluctuation()          ) +  dofsPerFace  * face._faceIndexLeft;

  double* QR = static_cast<double*>(pRight.getExtrapolatedPredictor()) + dataPerFace * face._faceIndexRight;
  double* FR = static_cast<double*>(pRight.getFluctuation()          ) + dofsPerFace  * face._faceIndexRight;

  // todo Time step must be interpolated in local time stepping case
  // both time step sizes are the same, so the min has no effect here.
  assertion3(std::isfinite(pLeft.getTimeStepSize()),pLeft.toString(),face._faceIndexLeft,face._direction);
  assertion3(std::isfinite(pRight.getTimeStepSize()),pRight.toString(),face._faceIndexRight,face._direction);
  assertion3(pLeft.getTimeStepSize()>=0.0,pLeft.toString(),face._faceIndexLeft,face._direction);
  assertion3(pRight.getTimeStepSize()>=0.0,pRight.toString(),face._faceIndexRight,face._direction);

  #ifdef Asserts
  std::string inputDataL = riemannDataToString(QL,FL,"L");
  std::string inputDataR = riemannDataToString(QR,FR,"R");
  #endif

  std::tuple<double,double> timeStepData = getRiemannSolverTimeStepData(pLeft,pRight);
  riemannSolver(
      FL,FR,QL,QR,
      std::get<0>(timeStepData),
      std::get<1>(timeStepData),
      pLeft.getSize(),
      face._direction, false, -1); // TODO(Dominic): Merge Riemann solver directly with the face integral and push the result on update
                                   // does not make sense to overwrite the flux when performing local time stepping; coarse grid flux must be constant, or not?

  #ifdef Asserts
  if ( _checkForNaNs ) { // assumes the solver is used as part of the hybrid solver
    std::string outputInformationL = riemannDataToString(QL,FL,"L");
    std::string outputInformationR = riemannDataToString(QR,FR,"R");

    const int nodesPerFace = dofsPerFace/_numberOfVariables;
    for (int i = 0; i<_numberOfVariables; i++) {
      for(int n=0; n<nodesPerFace; ++n) {
        assertion10(tarch::la::equals(pLeft.getTimeStepSize(),0.0) || (std::isfinite(FL[i]) && std::isfinite(FR[i])),
                   pLeft.toString(),pRight.toString(),face._direction,i,FL[i],FR[i],
                   inputDataL,inputDataR,outputInformationL,outputInformationR);
      }
    }
  }
  #endif
}

void exahype::solvers::ADERDGSolver::mergeWithBoundaryData(CellDescription& cellDescription,const int faceIndex,const int direction,const int orientation) {
  assertion1(cellDescription.getType()==CellDescription::Type::Leaf,cellDescription.toString());
  applyBoundaryConditions(cellDescription,faceIndex,direction,orientation);

  mergeWithAugmentationStatus(cellDescription,faceIndex,BoundaryStatus);
  mergeWithCommunicationStatus(cellDescription,faceIndex,BoundaryStatus);
  mergeWithRefinementStatus(cellDescription,faceIndex,BoundaryStatus);

  cellDescription.setNeighbourMergePerformed(faceIndex,true);
}

void exahype::solvers::ADERDGSolver::applyBoundaryConditions(CellDescription& cellDescription,const int faceIndex,const int direction,const int orientation) {
  assertion1(cellDescription.getType()==CellDescription::Type::Leaf,cellDescription.toString());
  assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorIndex()),cellDescription.toString());
  assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getFluctuationIndex()),cellDescription.toString());
  #if !defined(SharedMemoryParallelisation) && !defined(Parallel) && defined(Asserts)
  static int counter = 0;
  static double timeStamp = 0;
  if ( !tarch::la::equals(timeStamp,_minTimeStamp,1e-9) ) {
    logInfo("applyBoundaryConditions(...)","#boundaryConditions="<<counter);
    timeStamp = _minTimeStamp;
    counter=0;
  }
  counter++;
  #endif

  const int dataPerFace = getBndFaceSize();
  const int dofsPerFace = getBndFluxSize();
  const int gradientDataPerFace = getBndGradQSize();
  double* QIn = static_cast<double*>(cellDescription.getExtrapolatedPredictor()) +  dataPerFace * faceIndex;
  double* gradQIn = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient()) +  gradientDataPerFace * faceIndex;
  double* FIn = static_cast<double*>(cellDescription.getFluctuation())           +  dofsPerFace * faceIndex;
  const double* luh = static_cast<double*>(cellDescription.getSolution());

  #ifdef Asserts
  std::string inputData = riemannDataToString(QIn,FIn,"In");
  #endif

  #ifdef Asserts
  assertion4(std::isfinite(cellDescription.getTimeStamp()),cellDescription.toString(),faceIndex,direction,cellDescription.getTimeStamp());
  assertion4(std::isfinite(cellDescription.getTimeStepSize()),cellDescription.toString(),faceIndex,direction,cellDescription.getTimeStepSize());
  assertion4(cellDescription.getTimeStepSize()>=0.0, cellDescription.toString(),faceIndex,direction,cellDescription.getTimeStepSize());
  if ( _checkForNaNs ) {
    for(int i=0; i<dofsPerFace; ++i) {
      assertion6(tarch::la::equals(cellDescription.getTimeStepSize(),0.0) || std::isfinite(FIn[i]),cellDescription.toString(),faceIndex,direction,i,FIn[i],inputData);
    }
  }
  #endif

  // TODO(Dominic): Hand in space-time volume data. Time integrate it afterwards

  std::tuple<double,double> timeStepData = getRiemannSolverTimeStepData(cellDescription,cellDescription);
  boundaryConditions(
      FIn,QIn, gradQIn,
      luh,
      cellDescription.getOffset() + 0.5*cellDescription.getSize(),
      cellDescription.getSize(),
      std::get<0>(timeStepData),
      std::get<1>(timeStepData),
      direction,orientation);

  #ifdef Asserts
  assertion4(std::isfinite(cellDescription.getTimeStamp()),cellDescription.toString(),faceIndex,direction,cellDescription.getTimeStamp());
  assertion4(std::isfinite(cellDescription.getTimeStepSize()),cellDescription.toString(),faceIndex,direction,cellDescription.getTimeStepSize());
  assertion4(cellDescription.getTimeStepSize()>=0.0, cellDescription.toString(),faceIndex,direction,cellDescription.getTimeStepSize());
  if ( _checkForNaNs ) {
    for(int i=0; i<dofsPerFace; ++i) {
      assertion6(tarch::la::equals(cellDescription.getTimeStepSize(),0.0) || std::isfinite(FIn[i]),cellDescription.toString(),faceIndex,direction,i,FIn[i],inputData);
    }
  }
  #endif
}

#ifdef Parallel
const int exahype::solvers::ADERDGSolver::DataMessagesPerNeighbourCommunication    = 2;
const int exahype::solvers::ADERDGSolver::DataMessagesPerMasterWorkerCommunication = 2;

///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////
void
exahype::solvers::ADERDGSolver::appendNeighbourCommunicationMetadata(
    exahype::MetadataHeap::HeapEntries& metadata,
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    const int cellDescriptionsIndex,
    const int solverNumber) const {
  const int element = tryGetElement(cellDescriptionsIndex,solverNumber);

  if (element!=exahype::solvers::Solver::NotFound)  {
    CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

    metadata.push_back(static_cast<int>(cellDescription.getType()));
    metadata.push_back(cellDescription.getAugmentationStatus()); // TODO(Dominic): Add to docu: Might be merged multiple times!
    metadata.push_back(cellDescription.getCommunicationStatus());
    metadata.push_back(cellDescription.getRefinementStatus());
  } else {
    for (int i = 0; i < exahype::NeighbourCommunicationMetadataPerSolver; ++i) {
      metadata.push_back(exahype::InvalidMetadataEntry); // implicit conversion
    }
  }
}

void exahype::solvers::ADERDGSolver::mergeWithNeighbourMetadata(
    const int                                    solverNumber,
    Solver::CellInfo&                            cellInfo,
    const MetadataHeap::HeapEntries&             neighbourMetadata,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    const tarch::la::Vector<DIMENSIONS, int>&    posNeighbour,
    const tarch::la::Vector<DIMENSIONS, double>& barycentre) {
  const int neighbourAugmentationStatus  = neighbourMetadata[exahype::NeighbourCommunicationMetadataAugmentationStatus  ];
  const int neighbourCommunicationStatus = neighbourMetadata[exahype::NeighbourCommunicationMetadataCommunicationStatus ];
  const int neighbourRefinementStatus    = neighbourMetadata[exahype::NeighbourCommunicationMetadataRefinementStatus    ];

  //todo: Corruption status is not yet on metadata heap!

  logDebug("mergeWithNeighbourMetadata(...)", "received neighbour metadata="<<neighbourAugmentationStatus<<","<<neighbourCommunicationStatus<<","<<neighbourRefinementStatus);

  mergeWithNeighbourMetadata(solverNumber,cellInfo,
      neighbourAugmentationStatus,neighbourCommunicationStatus,neighbourRefinementStatus,
      // Uncorrupted, //todo: corruption status needs to be distributed with MPI!
      pos,posNeighbour,barycentre);
}

void exahype::solvers::ADERDGSolver::sendDataToNeighbour(
    const int                                     toRank,
    const int                                     solverNumber,
    Solver::CellInfo&                             cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  barycentre,
    const int                                     level) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != Solver::NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];
    Solver::BoundaryFaceInfo face(src,dest);

    if ( communicateWithNeighbour(cellDescription,face._faceIndex) ) {
      assertion(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorIndex()));
      assertion(DataHeap::getInstance().isValidIndex(cellDescription.getFluctuationIndex()));

      const int dofsPerFace = getBndFluxSize();
      const int dataPerFace = getBndFaceSize();

      const double* lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor()) + dataPerFace * face._faceIndex;
      const double* lFhbnd = static_cast<double*>(cellDescription.getFluctuation())           + dofsPerFace * face._faceIndex;

      waitUntilCompletedLastStep<CellDescription>(cellDescription,true,true);

      // Send order: lQhbnd,lFhbnd,observablesMin,observablesMax
      // Receive order: observablesMax,observablesMin,lFhbnd,lQhbnd
      DataHeap::getInstance().sendData(
          lQhbnd, dataPerFace, toRank, barycentre, level,
          peano::heap::MessageType::NeighbourCommunication);
      DataHeap::getInstance().sendData(
          lFhbnd, dofsPerFace, toRank, barycentre, level,
          peano::heap::MessageType::NeighbourCommunication);
      // TODO(Dominic): If anarchic time stepping send the time step over too.
    }
  }
}

// TODO(Dominic): Add to docu: We only perform a Riemann solve if a Cell is involved.
void exahype::solvers::ADERDGSolver::mergeWithNeighbourData(
    const int                                    fromRank,
    const int                                    solverNumber,
    Solver::CellInfo&                            cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& barycentre,
    const int                                    level) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != NotFound ) {
    Solver::BoundaryFaceInfo face(dest,src);
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];

    if( communicateWithNeighbour(cellDescription,face._faceIndex) ) {
      prefetchFaceData(cellDescription,face._faceIndex);

      // Send order: lQhbnd,lFhbnd
      // Receive order: lFhbnd,lQhbnd
      // TODO(Dominic): If anarchic time stepping, receive the time step too.
      const int dofsPerFace = getBndFluxSize();
      const int dataPerFace = getBndFaceSize();
      DataHeap::getInstance().receiveData(
          const_cast<double*>(_receivedFluctuations.data()),dofsPerFace, // TODO const-correct peano
          fromRank, barycentre, level,peano::heap::MessageType::NeighbourCommunication);
      DataHeap::getInstance().receiveData(                              // TODO const-correct peano
          const_cast<double*>(_receivedExtrapolatedPredictor.data()),dataPerFace,
          fromRank, barycentre, level, peano::heap::MessageType::NeighbourCommunication);
      
      solveRiemannProblemAtInterface(
          cellDescription, face,
          _receivedExtrapolatedPredictor.data(),
          _receivedFluctuations.data(),
          fromRank);
    }

    cellDescription.setNeighbourMergePerformed(face._faceIndex,true);
  }
}

void exahype::solvers::ADERDGSolver::solveRiemannProblemAtInterface(
    CellDescription& cellDescription,
    Solver::BoundaryFaceInfo& face,
    const double* const lQhbnd,
    const double* const lFhbnd,
    const int fromRank) {
  logDebug("solveRiemannProblemAtInterface(...)",
      "cell-description=" << cellDescription.toString());

  assertion(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorIndex()));
  assertion(DataHeap::getInstance().isValidIndex(cellDescription.getFluctuationIndex()));

  const int dataPerFace = getBndFaceSize();
  const int dofsPerFace = getBndFluxSize();
  if ( face._orientation==0 ) {
    const double* const QL = lQhbnd;
    double* FL             = const_cast<double*>(lFhbnd); // TODO const-correct kernels
    const double* const QR = static_cast<double*>( cellDescription.getExtrapolatedPredictor()) + dataPerFace * face._faceIndex;
    double* FR             = static_cast<double*>( cellDescription.getFluctuation())           + dofsPerFace * face._faceIndex;
    // TODO const-correct kernels

    #ifdef Asserts
    std::string inputDataL = riemannDataToString(QL,FL,"L");
    std::string inputDataR = riemannDataToString(QR,FR,"R");
    #endif

    riemannSolver(
        FL, FR, QL, QR,
    cellDescription.getTimeStamp(), cellDescription.getTimeStepSize(),
    cellDescription.getSize(), face._direction,false,face._faceIndex);
    #ifdef Asserts
    if ( _checkForNaNs ) {
      for (int ii = 0; ii<dofsPerFace; ii++) {
        assertion8(std::isfinite(FL[ii]), cellDescription.toString(),
            face._faceIndex, ii, QR[ii], QL[ii], FR[ii], FL[ii],fromRank);
        assertion8(std::isfinite(FR[ii]), cellDescription.toString(),
            face._faceIndex, ii, QR[ii], QL[ii], FR[ii], FL[ii],fromRank);
      }
    }
    #endif
  } else {
    const double* const QR = lQhbnd;
    const double* const QL = static_cast<double*>(cellDescription.getExtrapolatedPredictor()) + dataPerFace * face._faceIndex;
    double* FR = const_cast<double*>(lFhbnd); // TODO const-correct kernels
    double* FL = static_cast<double*>(cellDescription.getFluctuation()) + dofsPerFace * face._faceIndex; // TODO const-correct kernels

    #ifdef Asserts
    std::string inputDataL = riemannDataToString(QL,FL,"L");
    std::string inputDataR = riemannDataToString(QR,FR,"R");
    #endif

    std::tuple<double,double> timeStepData = getRiemannSolverTimeStepData(cellDescription,cellDescription);
    riemannSolver(
        FL, FR, QL, QR,
        std::get<0>(timeStepData),
        std::get<1>(timeStepData),
    cellDescription.getSize(),
        face._direction,false,face._faceIndex);
    
    #ifdef Asserts
    if ( _checkForNaNs ) {
      for (int ii = 0; ii<dofsPerFace; ii++) {
        assertion10(std::isfinite(FL[ii]), cellDescription.toString(),
            face._faceIndex, ii, QR[ii], QL[ii], FR[ii], FL[ii],fromRank,inputDataL,inputDataR);
        assertion10(std::isfinite(FR[ii]), cellDescription.toString(),
            face._faceIndex, ii, QR[ii], QL[ii], FR[ii], FL[ii],fromRank,inputDataL,inputDataR);
      }
    }
    #endif
  }
}

void exahype::solvers::ADERDGSolver::dropNeighbourData(
    const int                                    fromRank,
    const int                                    solverNumber,
    Solver::CellInfo&                            cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != NotFound ) {
    Solver::BoundaryFaceInfo face(dest,src);
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];

    if( communicateWithNeighbour(cellDescription,face._faceIndex) ) {
      for(int receives=0; receives<DataMessagesPerNeighbourCommunication; ++receives)
        DataHeap::getInstance().receiveData(
            fromRank, x, level,
            peano::heap::MessageType::NeighbourCommunication);
    }
  }
}

///////////////////////////////////
// WORKER->MASTER
///////////////////////////////////
void exahype::solvers::ADERDGSolver::sendDataToMaster(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  DataHeap::HeapEntries message = compileMessageForMaster();

  if ( !tarch::parallel::Node::getInstance().isGlobalMaster() ) {
    logDebug("sendDataToMaster(...)","Sending data to master: " <<
        "data[0]=" << message[0] << "," <<
        "data[1]=" << message[1] << "," <<
        "data[2]=" << message[2] << "," <<
        //"data[3]=" << message[3] << "," <<
        "to rank " << masterRank <<
        ", message size="<<message.size()
    );
  }

  DataHeap::getInstance().sendData(
      message.data(), message.size(),
      masterRank,x,level,peano::heap::MessageType::MasterWorkerCommunication);
}

exahype::DataHeap::HeapEntries
exahype::solvers::ADERDGSolver::compileMessageForMaster(const int capacity) const {
  const int messageSize = 3 + _numberOfGlobalObservables;
  DataHeap::HeapEntries message;
  message.reserve(std::max(messageSize,capacity));

  message.push_back(_admissibleTimeStepSize);
  message.push_back(convertToDouble(_meshUpdateEvent));
  message.push_back(_minTimeStamp);
  //message.push_back(_estimatedTimeStepSize);

  for (const auto observable : _globalObservables) {
    logDebug("sendDataToMaster(...)","Sending data to master: " << "entry=" << observable);
    message.push_back(observable);
  }

  assertion1(static_cast<int>(message.size())==messageSize,message.size());
  assertion1(std::isfinite(message[0]),message[0]);
  return message;
}

/**
 * At the time of the merging,
 * the workers and the master have already performed
 * at local update of the next predictor time step size
 * and of the predictor time stamp.
 * We thus need to minimise over both quantities.
 */
void exahype::solvers::ADERDGSolver::mergeWithWorkerData(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  const auto messageSize = 3 + _numberOfGlobalObservables;
  DataHeap::HeapEntries message(messageSize);

  tarch::multicore::RecursiveLock lock( 
        tarch::services::Service::receiveDanglingMessagesSemaphore );
  

  DataHeap::getInstance().receiveData(
      message.data(), message.size(),
      workerRank,x,level,peano::heap::MessageType::MasterWorkerCommunication);

  if ( tarch::parallel::Node::getInstance().isGlobalMaster() ) {
    logDebug("mergeWithWorkerData(...)","Receive data from worker rank: " <<
             "data[0]=" << message[0] << "," <<
             "data[1]=" << message[1] << "," <<
             "data[2]=" << message[2] << "," <<
             //"data[3]=" << message[3] << "," <<
             "from worker " << workerRank << "," <<
             "message size="<<message.size());
   }

  assertion1(static_cast<int>(message.size())==messageSize,message.size());
  mergeWithWorkerData(message);

  if ( tarch::parallel::Node::getInstance().isGlobalMaster() ) {
    logDebug("mergeWithWorkerData(...)","Updated fields: " <<
             "_admissibleTimeStepSize=" << _admissibleTimeStepSize << "," <<
             "_minTimeStamp=" << _minTimeStamp << "," <<
             //"_estimatedTimeStampSize=" << _estimatedTimeStepSize << "," <<
             "_meshUpdateEvent="        << Solver::toString(_meshUpdateEvent) );
  }
}

void exahype::solvers::ADERDGSolver::mergeWithWorkerData(const DataHeap::HeapEntries& message) {
  int index=0; // post update

  if(convertToMeshUpdateEvent(message[index+1])==MeshUpdateEvent::RollbackToTeamSolution) {
    _admissibleTimeStepSize = std::numeric_limits<double>::infinity();
    //_estimatedTimeStepSize = std::numeric_limits<double>::infinity();
    _minTimeStamp = std::numeric_limits<double>::infinity();

    _admissibleTimeStepSize = std::min( _admissibleTimeStepSize, message[index++] );
    _meshUpdateEvent       = mergeMeshUpdateEvents(_meshUpdateEvent,convertToMeshUpdateEvent(message[index++]));
    _minTimeStamp = std::min( _minTimeStamp, message[index++] );
    //_estimatedTimeStepSize = std::min(_estimatedTimeStepSize, message[index++]);
  }
  else {
    _admissibleTimeStepSize = std::min( _admissibleTimeStepSize, message[index++] );
    _meshUpdateEvent       = mergeMeshUpdateEvents(_meshUpdateEvent,convertToMeshUpdateEvent(message[index++]));
    index += 1;
  }
  DataHeap::HeapEntries observablesFromWorker = DataHeap::HeapEntries(_numberOfGlobalObservables);
  for (int i = 0; i < _numberOfGlobalObservables; ++i) {
    observablesFromWorker[i] = message[index++];
  }
  mergeGlobalObservables(_nextGlobalObservables.data(), observablesFromWorker.data()); // !  master hasn't wrapped up time step yet
}

///////////////////////////////////
// MASTER->WORKER
///////////////////////////////////
exahype::DataHeap::HeapEntries
exahype::solvers::ADERDGSolver::compileMessageForWorker(const int capacity) const {
  const auto messageSize = 5 + _numberOfGlobalObservables;
  DataHeap::HeapEntries message;
  message.reserve(std::max(messageSize,capacity));

  message.push_back(_minTimeStamp);
  message.push_back(_minTimeStepSize);
  message.push_back(_estimatedTimeStepSize);
  message.push_back(convertToDouble(_meshUpdateEvent));
  message.push_back(_stabilityConditionWasViolated ? 1.0 : -1.0);

  for (const auto observable : _globalObservables) {
    message.push_back(observable);
  }

  assertion1(static_cast<int>(message.size())==messageSize,message.size());
  return message;
}

void exahype::solvers::ADERDGSolver::sendDataToWorker(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  DataHeap::HeapEntries message = compileMessageForWorker();

  if ( tarch::parallel::Node::getInstance().isGlobalMaster() ) {
    logDebug(
        "sendDataToWorker(...)","Broadcasting time step data: " <<
        "data[0]=" << message[0] << "," <<
        "data[1]=" << message[1] << "," <<
        "data[2]=" << message[2] << "," <<
        "data[3]=" << message[3] << "," <<
        "data[4]=" << message[4]);
  }

  DataHeap::getInstance().sendData(
      message.data(), message.size(),
      workerRank,x,level,peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::ADERDGSolver::mergeWithMasterData(const DataHeap::HeapEntries& message) {
  int index=0;
  _minTimeStamp                  = message[index++];
  _minTimeStepSize               = message[index++];
  _estimatedTimeStepSize         = message[index++];
  _meshUpdateEvent               = convertToMeshUpdateEvent(message[index++]);
  _stabilityConditionWasViolated = (message[index++] > 0.0) ? true : false;

  for (int i = 0; i < _numberOfGlobalObservables; ++i) {
    _globalObservables[i] = message[index++];
  }
}

void exahype::solvers::ADERDGSolver::mergeWithMasterData(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  const auto messageSize = 5 + _numberOfGlobalObservables;
  DataHeap::HeapEntries message(messageSize);

  tarch::multicore::RecursiveLock lock(tarch::services::Service::receiveDanglingMessagesSemaphore, true);


  DataHeap::getInstance().receiveData(
      message.data(), message.size(),
      masterRank,x,level,peano::heap::MessageType::MasterWorkerCommunication);

  assertion1(static_cast<int>(message.size())==messageSize,message.size());
  mergeWithMasterData(message);
  
  lock.free();

  if ( !tarch::parallel::Node::getInstance().isGlobalMaster() ) {
    logDebug(
        "mergeWithMasterData(...)","Received message from master: " <<
        "data[0]=" << message[0] << "," <<
        "data[1]=" << message[1] << "," <<
        "data[2]=" << message[2] << "," <<
        "data[3]=" << message[3] << "," <<
        "data[4]=" << message[4]);
  }
}
#endif

std::string exahype::solvers::ADERDGSolver::toString() const {
  std::ostringstream stringstr;
  toString(stringstr);
  return stringstr.str();
}

void exahype::solvers::ADERDGSolver::toString (std::ostream& out) const {
  out << "(";
  out << "_identifier:"                   << _identifier;
  out << "_type:"                         << Solver::toString(_type) << ",";
  out << "_numberOfVariables:"            << _numberOfVariables << ",";
  out << "_numberOfParameters:"           << _numberOfParameters << ",";
  out << "_nodesPerCoordinateAxis:"       << _nodesPerCoordinateAxis << ",";
  out << "_maximumMeshSize:"              << _maximumMeshSize << ",";
  out << "_timeStepping:"                 << Solver::toString(_timeStepping) << ",";
  out << "_unknownsPerFace:"              << getUnknownsPerFace() << ",";
  out << "_unknownsPerCellBoundary:"      << getUnknownsPerCellBoundary() << ",";
  out << "_unknownsPerCell:"              << getUnknownsPerCell() << ",";
  out << "_fluxUnknownsPerCell:"          << getFluxUnknownsPerCell() << ",";
  out << "_spaceTimeUnknownsPerCell:"     << getSpaceTimeUnknownsPerCell() << ",";
  out << "_spaceTimeFluxUnknownsPerCell:" << getSpaceTimeFluxUnknownsPerCell() << ",";
  out << "_previousMinTimeStamp:"         << _previousMinTimeStamp << ",";
  out << "_previousMinTimeStepSize:"      << _previousMinTimeStepSize << ",";
  out << "_minTimeStamp:"                 << _minTimeStamp << ",";
  out << "_minTimeStepSize:"              << _minTimeStepSize << ",";
  out << "_estimatedTimeStepSize:"        << _estimatedTimeStepSize << ",";
  out <<  ")";
}

///////////////////////////////////
// DISTRIBUTED OFFLOADING
///////////////////////////////////
#if defined(SharedTBB) && defined(Parallel)
//Todo (Philipp) : Is this still needed if we pass this number directly to the progress routine?
void exahype::solvers::ADERDGSolver::setMaxNumberOfIprobesInProgressOffloading(int maxIprobes) {
  MaxIprobesInOffloadingProgress = maxIprobes;
}

///////////////////////////////////
// PROGRESS TASKS
///////////////////////////////////

#ifdef OffloadingUseProgressTask
//TODO: may not be needed but left for now
exahype::solvers::ADERDGSolver::ReceiveJob::ReceiveJob(ADERDGSolver& solver)
  :  tarch::multicore::jobs::Job(tarch::multicore::jobs::JobType::BackgroundTask, 0, tarch::multicore::DefaultPriority*8),
    _solver(solver) {};

bool exahype::solvers::ADERDGSolver::ReceiveJob::run( bool isCalledOnMaster ) {
#if defined(UseSmartMPI)
  MPI_Status_Offload stat;
#else
  MPI_Status stat;
#endif
  MPI_Comm comm = exahype::reactive::ReactiveContext::getInstance().getMPICommunicator();
  //MPI_Comm commStatus = exahype::reactive::OffloadingManager::getInstance().getMPICommunicatorStatus();
  int receivedTask = -1;
  int receivedStatus = -1;
  int lastRecvTag = -1;
  int lastRecvSrc = -1;

  int ierr;

  if(isCalledOnMaster) return true;

  tarch::multicore::Lock lock(OffloadingSemaphore, false);
  bool canRun = lock.tryLock();
  while(!canRun) {
     canRun = lock.tryLock();
  }
  int itcount = 0;

  logDebug("run()","receive job running");

  while(ActiveSenders.size()>0) {
    exahype::reactive::RequestManager::getInstance().progressRequests();
    /* MPI_Iprobe(MPI_ANY_SOURCE, 0, comm, &receivedStatus, &stat2);
    if(receivedStatus) {
      int terminatedSender = stat2.MPI_SOURCE;
      logInfo("run()","active sender "<<terminatedSender<<" has sent termination signal ");
      exahype::reactive::OffloadingManager::getInstance().receiveCompleted(terminatedSender);
      ActiveSenders.erase(terminatedSender);
    }*/

    tarch::multicore::RecursiveLock lock2( tarch::services::Service::receiveDanglingMessagesSemaphore, false );
    if(lock2.tryLock()) {
      tarch::parallel::Node::getInstance().receiveDanglingMessages();
      lock2.free();
    }

#if defined(UseSmartMPI)
    ierr = MPI_Iprobe_offload(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &receivedTask, &stat);
#else
    ierr = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &receivedTask, &stat);
#endif
    assertion(ierr==MPI_SUCCESS);

    // tag 0 is reserved for termination message
    if(receivedTask && stat.MPI_TAG==0) {
      int terminatedSender = stat.MPI_SOURCE;
      logDebug("run()","active sender "<<terminatedSender<<" has sent termination signal ");
      exahype::reactive::ReactiveContext::getInstance().receiveCompleted(terminatedSender); //, stat.rail); //todo: won't work with SmartMPI
      ActiveSenders.erase(terminatedSender);
#if defined(UseSmartMPI)
      ierr = MPI_Iprobe_offload(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &receivedTask, &stat);
#else
      ierr = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &receivedTask, &stat);
#endif
      assertion(ierr==MPI_SUCCESS);
    }

    if(receivedTask) {
      logDebug("run()","adding active sender "<<stat.MPI_SOURCE<< " tag "<<stat.MPI_TAG);
      ActiveSenders.insert(stat.MPI_SOURCE);
      exahype::reactive::ReactiveContext::getInstance().triggerVictimFlag();
      int msgLen = -1;
#if defined(UseSmartMPI)
      MPI_Get_count_offload(&stat, MigratablePredictionJobMetaData::getMPIDatatype(), &msgLen);
#else
      MPI_Get_count(&stat, MigratablePredictionJobMetaData::getMPIDatatype(), &msgLen);
#endif
      // is this message metadata? -> if true, we are about to receive a new STP task
      if(msgLen==MigratablePredictionJobMetaData::getMessageLen() && !(lastRecvTag==stat.MPI_TAG && lastRecvSrc==stat.MPI_SOURCE)) {
        lastRecvTag=stat.MPI_TAG;
        lastRecvSrc=stat.MPI_SOURCE;

        assertion(lastRecvTag!=_solver._lastReceiveTag[lastRecvSrc]);
        _solver._lastReceiveTag[lastRecvSrc] = lastRecvTag;

        MPI_Request receiveRequests[NUM_REQUESTS_MIGRATABLE_COMM+1];
        MigratablePredictionJobData *data = new MigratablePredictionJobData(_solver);
        _solver._mapTagRankToStolenData.insert(std::make_pair(std::make_pair(stat.MPI_SOURCE, stat.MPI_TAG), data));
#if defined(UseSmartMPI)
        _solver.mpiRecvMigratablePredictionJobOffload(
                  data->_luh.data(),
                  stat.MPI_SOURCE,
                  stat.MPI_TAG,
                  exahype::reactive::ReactiveContext::getInstance().getMPICommunicator(),
                  stat.rail,
                  &(data->_metadata[0]));
        MigratablePredictionJob::receiveHandler(solver, stat.MPI_TAG, stats.MPI_SOURCE);
      }
#else
        _solver.mpiIrecvMigratablePredictionJob(
                 data->_luh.data(),
                 stat.MPI_SOURCE,
                 stat.MPI_TAG,
                 exahype::reactive::ReactiveContext::getInstance().getMPICommunicator(),
                 &receiveRequests[0],
                 &(data->_metadata));

        if(tarch::multicore::jobs::getNumberOfWaitingBackgroundJobs()<=1) {
             //logInfo("progressOffloading()","running out of tasks and could not receive stolen task so we just block!");
             double wtime = -MPI_Wtime();
             exahype::reactive::RequestManager::getInstance().submitRequests(
               receiveRequests,
               NUM_REQUESTS_MIGRATABLE_COMM+1,
               stat.MPI_TAG,
               stat.MPI_SOURCE,
               MigratablePredictionJob::receiveHandler,
               exahype::reactive::RequestType::receive,
               &_solver,
               true);
             wtime+= MPI_Wtime();
             if(wtime>0.01)
               logDebug("progressOffloading()","blocking for stolen task took too long:"<<wtime<<"s");
           }
           else {
             exahype::reactive::RequestManager::getInstance().submitRequests(
               receiveRequests,
               NUM_REQUESTS_MIGRATABLE_COMM+1,
               stat.MPI_TAG,
               stat.MPI_SOURCE,
               MigratablePredictionJob::receiveHandler,
               exahype::reactive::RequestType::receive,
               &_solver,
               true);
           }
         }
#endif
      }
#if defined(UseSmartMPI)
      ierr = MPI_Iprobe_offload(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &receivedTask, &stat);
#else
      ierr = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &receivedTask, &stat);
#endif
      itcount++;
  }
  logDebug("run()","terminated receive job after "<<itcount<<" iterations");

  NumberOfReceiveJobs--;
  lock.free();
  return false;
}

#if defined(OffloadingNoEarlyReceiveBacks)
exahype::solvers::ADERDGSolver::ReceiveBackJob::ReceiveBackJob(ADERDGSolver& solver)
  :  tarch::multicore::jobs::Job(tarch::multicore::jobs::JobType::BackgroundTask, 0, tarch::multicore::DefaultPriority*8),
    _solver(solver) {};

bool exahype::solvers::ADERDGSolver::ReceiveBackJob::run( bool isCalledOnMaster ) {
  tarch::multicore::Lock lock(OffloadingSemaphore, false);
  bool canRun = lock.tryLock();
  while(!canRun) {
     canRun = lock.tryLock();
  }

  bool run = true;

  while(run) {
    exahype::reactive::RequestManager::getInstance().progressRequests();

    int tag, srcRank, myRank;
    myRank = tarch::parallel::Node::getInstance().getRank();
    tag = MPI_ANY_TAG;
    srcRank = MPI_ANY_SOURCE;

    int receivedTaskBack = 1;
    MPI_Status statMapped;
    MPI_Comm commMapped = exahype::reactive::ReactiveContext::getInstance().getMPICommunicatorMapped();
    int ierr = MPI_Iprobe(srcRank, tag, commMapped, &receivedTaskBack, &statMapped);
    assertion(ierr==MPI_SUCCESS);
    if(receivedTaskBack) {
      //exahype::reactive::OffloadingManager::getInstance().setRunningAndReceivingBack();
      tbb::concurrent_hash_map<int, CellDescription*>::accessor a_tagToCellDesc;
      bool found = _solver._mapTagToCellDesc.find(a_tagToCellDesc, statMapped.MPI_TAG);
      assertion(found);
      auto cellDescription = a_tagToCellDesc->second;
      a_tagToCellDesc.release();

      double *lduh   = static_cast<double*>(cellDescription->getUpdate());
      double *lQhbnd = static_cast<double*>(cellDescription->getExtrapolatedPredictor());
      double *lFhbnd = static_cast<double*>(cellDescription->getFluctuation());
      double* lGradQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient());

      assertion(statMapped.MPI_TAG!=_solver._lastReceiveBackTag[statMapped.MPI_SOURCE]);
      _solver._lastReceiveBackTag[statMapped.MPI_SOURCE] =  statMapped.MPI_TAG;

      MPI_Request recvRequests[NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME];
      _solver.mpiIrecvMigratablePredictionJobOutcome(
        lduh,
        lQhbnd,
        lFhbnd,
        lGradQhbnd,
        statMapped.MPI_SOURCE,
        statMapped.MPI_TAG,
        commMapped,
        recvRequests);

      exahype::reactive::RequestManager::getInstance().submitRequests(
        recvRequests, NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME, statMapped.MPI_TAG, statMapped.MPI_SOURCE,
      exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveBackHandler,
      exahype::reactive::RequestType::receiveBack, &_solver, false);
    }
    run = receivedTaskBack || exahype::reactive::RequestManager::getInstance().hasOutstandingRequestOfType(exahype::reactive::RequestType::receiveBack);
  }

  NumberOfReceiveBackJobs--;
  lock.free();
  return false;
}
#endif
#endif

exahype::solvers::ADERDGSolver::OffloadingManagerJob::OffloadingManagerJob(ADERDGSolver& solver) :
#ifndef OffloadingUseProgressThread
 tarch::multicore::jobs::Job(tarch::multicore::jobs::JobType::BackgroundTask, 0, tarch::multicore::DefaultPriority-1),
#endif
  _solver(solver),
  _state(State::Running) {}

exahype::solvers::ADERDGSolver::OffloadingManagerJob::~OffloadingManagerJob() {}

/*bool exahype::solvers::ADERDGSolver::OffloadingManagerJob::operator()() {
  return run();
  //return true;
}*/

#if defined(OffloadingUseProgressThread)
tbb::task* exahype::solvers::ADERDGSolver::OffloadingManagerJob::execute() {
   while(run( false )) {};
   return nullptr;
}
#endif

bool exahype::solvers::ADERDGSolver::OffloadingManagerJob::run( bool isCalledOnMaster ) {
// static bool terminated = false;
  bool result=true;
#ifdef USE_ITAC
//  VT_begin(event_offloadingManager);
#endif
  //logInfo("run()", "starting... ");

  switch (_state) {
    case State::Running:
    {
      if( isCalledOnMaster ) {
          return true;
      }

      //if(peano::utils::UserInterface::getMemoryUsageMB()>50000) {
      //    logInfo("run()", "WARNING: memory usage is quite high!");
      //}
      //double time = -MPI_Wtime();
      exahype::solvers::ADERDGSolver::progressOffloading(&_solver, false, std::numeric_limits<int>::max());
      //exahype::solvers::ADERDGSolver::progressOffloading(&_solver, false, 1);
      //time += MPI_Wtime();
      //if(time>0.02)
      //  logInfo("run","took too long"<<time);

      if(_solver._offloadingManagerJobTriggerTerminate) {
          _state = State::Terminate;
      }

      break;
    }
    case State::Resume:
      _state = State::Running;
      break;
    case State::Paused:
      result = false;
      break;
    case State::Terminate:
    {
      exahype::reactive::PerformanceMonitor::getInstance().stop();
      logDebug("offloadingManager", " terminated ");
      _solver._offloadingManagerJobTerminated = true;
      result = false;
      break;
    }
    default:
      result = false;
      break;
  }
#ifdef USE_ITAC
//  VT_end(event_offloadingManager);
#endif
  return result;
}

void exahype::solvers::ADERDGSolver::OffloadingManagerJob::terminate() {
  tarch::multicore::Lock lock(OffloadingSemaphore, true);
  _state = State::Terminate;
  lock.free();
}

void exahype::solvers::ADERDGSolver::OffloadingManagerJob::pause() {
  _state = State::Paused;
}

void exahype::solvers::ADERDGSolver::OffloadingManagerJob::resume() {
  _state = State::Resume;
}

void exahype::solvers::ADERDGSolver::initOffloadingManager() {
  logDebug("startOffloadingManager", " starting ");
  _offloadingManagerJobTerminated = false;
  _offloadingManagerJobTriggerTerminate = false;
#ifdef OffloadingUseProgressThread
  static tbb::task_group_context  backgroundTaskContext(tbb::task_group_context::isolated);
  _offloadingManagerJob = new( backgroundTaskContext ) OffloadingManagerJob(*this);
  //_offloadingManagerJob = new OffloadingManagerJob(*this);
  //assertion(_offloadingManagerJob!=nullptr);
  tbb::task::enqueue(*_offloadingManagerJob);
#else
  _offloadingManagerJob = new OffloadingManagerJob(*this);
  _offloadingManagerJobStarted = true;
#endif
}

#ifndef OffloadingUseProgressThread
void exahype::solvers::ADERDGSolver::pauseOffloadingManager() {
  //tarch::multicore::Lock lock(OffloadingSemaphore, true);
  logDebug("pauseOffloadingManager", "pausing ");
  if(_offloadingManagerJob!=nullptr){
    _offloadingManagerJob->pause();
    _offloadingManagerJob = nullptr;
  }
  //lock.free();
}

void exahype::solvers::ADERDGSolver::resumeOffloadingManager() {
  logDebug("resumeOffloadingManager", "resuming ");
  //old job will be deleted so we create a new one here
  //assertion(_offloadingManagerJob==nullptr);
  if(_offloadingManagerJob==nullptr) {
    _offloadingManagerJob = new OffloadingManagerJob(*this);
    _offloadingManagerJob->resume();
    peano::datatraversal::TaskSet spawnedSet(_offloadingManagerJob);
  }
  //lock.free();
}
#endif

void exahype::solvers::ADERDGSolver::stopOffloadingManager() {
  logDebug("stopOffloadingManager", " stopping ");
  //assertion(_offloadingManagerJob != nullptr);
  _offloadingManagerJobTriggerTerminate = true;

#if defined(OffloadingUseProgressThread)
  while(!_offloadingManagerJobTerminated) {};
  //delete _offloadingManagerJob;
#endif
  //while(!exahype::reactive::PerformanceMonitor::getInstance().isGloballyTerminated()) {tarch::multicore::jobs::finishToProcessBackgroundJobs(); };
  //while(tarch::multicore::jobs::finishToProcessBackgroundJobs()) {};

  //assertion(_offloadingManagerJob != nullptr);
  //delete _offloadingManagerJob;
}

#ifdef OffloadingUseProgressTask
#if defined(OffloadingNoEarlyReceiveBacks)
void exahype::solvers::ADERDGSolver::spawnReceiveBackJob() {
  if(NumberOfReceiveBackJobs==0) {
    NumberOfReceiveBackJobs++;
    peano::datatraversal::TaskSet spawnedSet(new ReceiveBackJob(*this));
  }
}
#endif
#endif


///////////////////////////////////
// TASK OFFLOADING
///////////////////////////////////
exahype::solvers::ADERDGSolver::MigratablePredictionJob* exahype::solvers::ADERDGSolver::createFromData(
  MigratablePredictionJobData *data,
  const int origin,
  const int tag) {
  return new MigratablePredictionJob(*this,
      -1,
      -1,
      data->_metadata._predictorTimeStamp,
      data->_metadata._predictorTimeStepSize,
      data->_luh.data(),
      data->_lduh.data(),
      data->_lQhbnd.data(),
      data->_lFhbnd.data(),
      data->_lGradQhbnd.data(),
      &(data->_metadata._dx[0]),
      &(data->_metadata._center[0]),
      origin,
      tag);
}

int exahype::solvers::ADERDGSolver::getResponsibleRankForCellDescription(const void* cellDescription) {
  int resultRank = -1;

  tbb::concurrent_hash_map<const CellDescription*, std::pair<int, int>>::accessor a_cellDescToTagRank;
  bool found =  _mapCellDescToTagRank.find(a_cellDescToTagRank, static_cast<const CellDescription*>(cellDescription));
  if(found)
    resultRank = a_cellDescToTagRank->second.second;
  a_cellDescToTagRank.release();

  if(!found) return tarch::parallel::Node::getInstance().getRank();

  return resultRank;
}

#if defined(OffloadingLocalRecompute)
tarch::multicore::jobs::Job* exahype::solvers::ADERDGSolver::grabRecomputeJobForCellDescription(const void* cellDescription) {
  tbb::concurrent_hash_map<const CellDescription*, tarch::multicore::jobs::Job* >::accessor a_cellDescToJob;
  bool found = _mapCellDescToRecompJob.find(a_cellDescToJob, static_cast<const CellDescription*>(cellDescription));
  tarch::multicore::jobs::Job *job = nullptr;
  if(found) {
   job = a_cellDescToJob->second;
     _mapCellDescToRecompJob.erase(a_cellDescToJob);
  }
  return job;
}

void exahype::solvers::ADERDGSolver::addRecomputeJobForCellDescription(tarch::multicore::jobs::Job* job, const CellDescription* cellDescription) {
  tbb::concurrent_hash_map<const CellDescription*, tarch::multicore::jobs::Job* >::accessor a_cellDescToJob;
  bool found = _mapCellDescToRecompJob.find(a_cellDescToJob, static_cast<const CellDescription*>(cellDescription));
  assertion(!found);
  _mapCellDescToRecompJob.insert(std::make_pair(cellDescription, job));
}

#endif

void exahype::solvers::ADERDGSolver::getResponsibleRankTagForCellDescription(const void* cellDescription, int& rank, int& tag) {

  tbb::concurrent_hash_map<const CellDescription*, std::pair<int, int>>::accessor a_cellDescToTagRank;
  bool found =  _mapCellDescToTagRank.find(a_cellDescToTagRank, static_cast<const CellDescription*>(cellDescription));
  if(found) {
    rank = a_cellDescToTagRank->second.second;
    tag = a_cellDescToTagRank->second.first;
  }
  else {
    rank = tarch::parallel::Node::getInstance().getRank();
    tag = -1;
  }
  a_cellDescToTagRank.release();

}

int exahype::solvers::ADERDGSolver::getTaskPriorityLocalMigratableJob(int cellDescriptionsIndex, int element, double timeStamp, bool isSkeleton){
  if(exahype::reactive::ReactiveContext::getInstance().getResilienceStrategy()
     !=exahype::reactive::ReactiveContext::ResilienceStrategy::None  && !isSkeleton) {

    int team = exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber();
    int teamSize = exahype::reactive::ReactiveContext::getInstance().getTMPINumTeams();

    CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex, element);

    tarch::la::Vector<DIMENSIONS, double> center;
    center = (cellDescription.getOffset()+0.5*cellDescription.getSize());

    int prio_shuffle = 0;

    // if((exahype::reactive::OffloadingContext::getInstance().getResilienceStrategy()
    //     >=exahype::reactive::OffloadingContext::ResilienceStrategy::TaskSharingResilienceChecks)
    //     && exahype::reactive::ResilienceTools::CheckLimitedCellsOnly) {
    //   prio_shuffle = 0; //rely on shuffling from grid traversal
    // }
    // else {
//     int tasks_per_team = (exahype::mappings::FinaliseMeshRefinement::NumberOfEnclaveCells/teamSize);
//     prio_shuffle = (LocalStealableSTPCounter/tasks_per_team)%teamSize;
//#else
     prio_shuffle = (LocalStealableSTPCounter+team)%teamSize;
     //}
//#endif
     int prio = getTaskPriority(false)+ prio_shuffle;

/*    logDebug("getTaskPriorityLocalStealableJob()", "team = "<<team
                                                <<" center[0] = "<< center[0]
                                                <<" center[1] = "<< center[1]
#if DIMENSIONS==3
                                                <<" center[2] = "<< center[2]
#endif
                                                <<" time stamp = "<<timeStamp
                                                <<" prio = "<<prio);*/

     return prio;
  }
  else {
    return getTaskPriority(isSkeleton);
  }
}

void exahype::solvers::ADERDGSolver::submitOrSendMigratablePredictionJob(MigratablePredictionJob* job) {

   int myRank = tarch::parallel::Node::getInstance().getRank();
   int destRank = myRank;

   bool lastSend = false;
   if(!job->_isSkeleton) {
     exahype::reactive::ReactiveContext::getInstance().selectVictimRank(destRank, lastSend);
     assertion(destRank>=0);
   }

   logDebug("submitOrSendMigratablePredictionJob", "there are "<<NumberOfEnclaveJobs<<" Enclave Jobs and "<<NumberOfRemoteJobs<< " Remote Jobs");

   if(myRank!=destRank) {
    // logInfo("submitOrSendMigratablePredictionJob","element "<<job->_element<<" predictor time stamp"<<job->_predictorTimeStamp<<" predictor time step size "<<job->_predictorTimeStepSize);
     //OffloadEntry entry = {destRank, job->_cellDescriptionsIndex, job->_element, job->_predictorTimeStamp, job->_predictorTimeStepSize};
     //_outstandingOffloads.push( entry );
     auto& cellDescription = getCellDescription(job->_cellDescriptionsIndex, job->_element);

     double *luh    = static_cast<double*>(cellDescription.getSolution());
#if !defined(UseSmartMPI) || defined(SmartMPINB)
     MPI_Request sendRequests[NUM_REQUESTS_MIGRATABLE_COMM+1];
#endif
     int tag = exahype::reactive::ReactiveContext::getInstance().getOffloadingTag(); //cellDescriptionsIndex is not a good idea here, as map entries with key tag may be overwritten if previous sends have not been marked as finished
      //need to create a copy
#if defined(OffloadingLocalRecompute)
     //Todo: we probably don't need this anymore as we don't need a copy
     //costs 70 s!
     //MigratablePredictionJobData *data = new MigratablePredictionJobData(*this);
     //AllocatedSTPsReceive++;

     //std::memcpy(&data->_luh[0], luh, data->_luh.size()*sizeof(double));
     //std::memcpy(&data->_lduh[0], lduh, data->_lduh.size()*sizeof(double));
     //std::memcpy(&data->_lQhbnd[0], lQhbnd, data->_lQhbnd.size()*sizeof(double));
     //std::memcpy(&data->_lFhbnd[0], lFhbnd, data->_lFhbnd.size()*sizeof(double));
     MigratablePredictionJobMetaData *metadata = new MigratablePredictionJobMetaData();
     job->packMetaData(metadata);

     //_mapTagToSTPData.insert(std::make_pair(tag, data));
     _mapTagToCellDesc.insert(std::make_pair(tag, &cellDescription));
     _mapTagToMetaData.insert(std::make_pair(tag, metadata));
     //logInfo("submitOrSendMigratablePredictionJob", "inserting tag"<<tag);

     tbb::concurrent_hash_map<const CellDescription*, std::pair<int,int>>::accessor a_cellDescToTagRank;
     //logInfo("receiveBackHandler", " cleaning up cell desc to tag/rank for "<<cellDescription);
     bool found = _mapCellDescToTagRank.find(a_cellDescToTagRank, &cellDescription);
     assertion(!found);
     a_cellDescToTagRank.release();
     _mapCellDescToTagRank.insert(std::make_pair(&cellDescription, std::make_pair(tag, destRank)));
     //_mapTagToOffloadTime.insert(std::make_pair(tag, -MPI_Wtime()));

     /*mpiIsendMigratablePredictionJobOutcome(
         &data->_luh[0],
         &data->_lduh[0],
         &data->_lQhbnd[0],
         &data->_lFhbnd[0],
         destRank,
         tag,
         exahype::reactive::OffloadingManager::getInstance().getMPICommunicator(),
         sendRequests,
         data->_metadata);*/
#if defined(UseSmartMPI)
     mpiSendMigratablePredictionJobOffload(
              luh,
              destRank,
              tag,
              exahype::reactive::ReactiveContext::getInstance().getMPICommunicator(),
              sendRequests,
              metadata);
     MigratablePredictionJob::sendHandler(this, tag, destRank);
#else
     mpiIsendMigratablePredictionJob(
              luh,
              destRank,
              tag,
              exahype::reactive::ReactiveContext::getInstance().getMPICommunicator(),
              sendRequests,
              metadata);
#endif

#else
     MigratablePredictionJobMetaData *metadata = new MigratablePredictionJobMetaData();
     job->packMetaData(metadata);
#if defined(Asserts)
     tbb::concurrent_hash_map<int, MigratablePredictionJobMetaData*>::accessor a_TagToMetadata;
     bool found = _mapTagToMetaData.find(a_TagToMetadata, tag);
     assert(!found);
#endif

     // we need this info when the task comes back...
     _mapTagToMetaData.insert(std::make_pair(tag, metadata));
     _mapTagToCellDesc.insert(std::make_pair(tag, &cellDescription));
     //logInfo("submitOrSendMigratablePredictionJob", "inserting tag"<<tag);
     _mapCellDescToTagRank.insert(std::make_pair(&cellDescription, std::make_pair(tag, destRank)));
     //_mapTagToOffloadTime.insert(std::make_pair(tag, -MPI_Wtime()));
     logDebug("submitOrSendMigratablePredictionJob()","send away with tag "<<tag<<" to rank "<<destRank<<" job "<<metadata->to_string());
     // send away
#if defined(UseSmartMPI)
#if defined(SmartMPINB)
     mpiIsendMigratablePredictionJobOffload(
         luh,
         destRank,
         tag,
         exahype::reactive::ReactiveContext::getInstance().getMPICommunicator(),
         sendRequests,
         metadata);
#else
     mpiSendMigratablePredictionJobOffload(
         luh,
         destRank,
         tag,
         exahype::reactive::ReactiveContext::getInstance().getMPICommunicator(),
         metadata);
     MigratablePredictionJob::sendHandler(this, tag, destRank);
#endif /*SmartMPINB*/
#else
     mpiIsendMigratablePredictionJob(
         luh,
         destRank,
         tag,
         exahype::reactive::ReactiveContext::getInstance().getMPICommunicator(),
         sendRequests,
         metadata);
#endif

#endif

#ifdef OffloadingLocalRecompute
     //logInfo("submitOrSendMigratablePredictionJob()", "keeping job in local queue");
    /*logInfo("submitOrSendMigratablePredictionJob()", "keeping job in local queue"
                                               <<" center[0] = "<<data->_metadata[0]
                                               <<" center[1] = "<<data->_metadata[1]
                                               <<" center[2] = "<<data->_metadata[2]
                                               <<" time stamp = "<<job->_predictorTimeStamp);*/
     //job->resetPriority(job->getPriority()/2);
     job->_isLocalReplica = true;
     addRecomputeJobForCellDescription(job, &cellDescription);
     //peano::datatraversal::TaskSet spawnedSet( job );
#endif
     //logInfo("submitOrSendMigratablePredictionJob"," there are "<<tarch::multicore::jobs::getNumberOfWaitingBackgroundJobs()<<" background jobs ");

#if !defined(UseSmartMPI) || defined(SmartMPINB)
     exahype::reactive::RequestManager::getInstance().submitRequests(
          sendRequests,
          NUM_REQUESTS_MIGRATABLE_COMM+1,
          tag,
          destRank,
          exahype::solvers::ADERDGSolver::MigratablePredictionJob::sendHandler,
          exahype::reactive::RequestType::send,
          this);
#endif

     // post receive back requests
//     MPI_Request recvRequests[4];
//     mpiIrecvMigratablePredictionJob(
//        luh, lduh, lQhbnd,
//      lFhbnd, destRank, tag, recvRequests);
//
//     exahype::reactive::RequestManager::getInstance().submitRequests(
//         recvRequests, 4, tag, destRank,
//         exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveBackHandler,
//      exahype::reactive::RequestType::receiveBack, this);

     assertion(!job->_isSkeleton); //skeleton jobs should never be sent away!

     NumberOfRemoteJobs++;
#ifndef OffloadingLocalRecompute
     delete job;
#endif

     exahype::reactive::OffloadingProfiler::getInstance().notifyOffloadedTask(destRank);
     exahype::reactive::PerformanceMonitor::getInstance().decCurrentTasks();

#ifdef OffloadingUseProgressTask
     if(lastSend)
        exahype::reactive::ReactiveContext::getInstance().notifyAllVictimsSendCompletedIfNotNotified();
#endif
  }
  else {
    peano::datatraversal::TaskSet spawnedSet( job );
  }
}


///////////////////////////////////
// TASK SHARING
///////////////////////////////////
void exahype::solvers::ADERDGSolver::cleanUpStaleTaskOutcomes(bool isFinal) {
  int unsafe_size = _allocatedOutcomes.unsafe_size();
  assertion(unsafe_size>=0);
  bool gotOne = true;
  int i = 0;

#if defined(OffloadingCheckForSlowOperations)
  double timing = - MPI_Wtime();
#endif

  //Todo (Philipp): refactor and make nice
  logDebug("cleanUpStaleTaskOutcomes()", "before cleanup there are "<<_allocatedOutcomes.unsafe_size()<<" allocated received jobs left, "
                                                                     <<_mapTagToSTPData.size()<<" jobs to send,"
                                                                     <<" allocated jobs send "<<AllocatedSTPsSend
                                                                     <<" allocated jobs receive "<<AllocatedSTPsReceive
                                                                     <<" estimated additional mem consumption "<<(double) getAdditionalCurrentMemoryUsageReplication()/1E9<<"GB"
                                                                     <<" actual mem usage "<<peano::utils::UserInterface::getMemoryUsageMB()
                                                                     <<" memory per stp "<< sizeof(MigratablePredictionJobData) + sizeof(double) * ( getDataPerCell() + getUpdateSize() + getBndTotalSize() + getBndFluxTotalSize() )
                                                                     <<" allocated stps (constructor) "<<AllocatedSTPs
                                                                     <<" entries in hash map "<<_outcomeDatabase.size()
                                                                     <<" sent STPs "<<SentSTPs
                                                                     <<" completed sends "<<CompletedSentSTPs
                                                                     <<" outstanding requests (for outcomes) "<<exahype::reactive::RequestManager::getInstance().getNumberOfOutstandingRequests(exahype::reactive::RequestType::sendOutcome)
                                                                                            +exahype::reactive::RequestManager::getInstance().getNumberOfOutstandingRequests(exahype::reactive::RequestType::receiveOutcome));

  double lastconsistentTimeStamp, lastconsistentTimeStepSize, lastconsistentEstimatedSize;
  exahype::reactive::TimeStampAndTriggerTeamHistory::getInstance().getLastConsistentTimeStepData(lastconsistentTimeStamp, lastconsistentTimeStepSize, lastconsistentEstimatedSize);

  double minTimeStampToKeep = _previousMinTimeStamp;

  if(exahype::reactive::ReactiveContext::getInstance().getResilienceStrategy()
    >= exahype::reactive::ReactiveContext::ResilienceStrategy::TaskSharingResilienceCorrection)
    minTimeStampToKeep = std::min(minTimeStampToKeep, lastconsistentTimeStepSize); //keep even older outcomes to be able to rollback

  while( (i< unsafe_size || isFinal) && gotOne) {
    MigratablePredictionJobOutcomeKey key;
    gotOne = _allocatedOutcomes.try_pop_front(&key);

    if(!gotOne) break;

    logDebug("cleanUpStaleTaskOutcomes()", " time stamp of key ="<<key._timestamp);
    i++;

    assertion(key._center!=nullptr);
    //logInfo("cleanUpStaleReplicatedSTPs()", " trying to find key - "
    //                                        <<" center[0] = "<<key._center[0]
    //                                       <<" center[1] = "<<key._center[1]
    //                                       <<" center[2] = "<<key._center[2]
    //                                       <<" time stamp = "<<key._timestamp);

    if(key._timestamp>=minTimeStampToKeep) {
      _allocatedOutcomes.push_front(key);
      logDebug("cleanUpStaleTaskOutcomes()", " breaking out of loop: time stamp of key ="<<key._timestamp)
      break;
    }

    MigratablePredictionJobData *data = nullptr;
    DeliveryStatus status;
    bool found = _outcomeDatabase.tryFindAndExtractOutcome(key, &data, status);

    if(found) {
      //logInfo("cleanUpStaleReplicatedSTPs()", " time stamp "<<a_jobToData->first.timestamp<< " _minTimeStamp "<<_minTimeStamp);
      logDebug("cleanUpStaleTaskOutcomes()",   data->_metadata.to_string());
      assertion(data!=nullptr);
      delete data;
      AllocatedSTPsReceive--;
      exahype::reactive::ResilienceStatistics::getInstance().notifyLateTask();
    }
  }

  if(isFinal) {
      for(auto & elem: _mapTagToSTPData) {
          delete elem.second;
      }
  }

  logDebug("cleanUpStaleTaskOutcomes()", " there are "<<_allocatedOutcomes.unsafe_size()<<" allocated received jobs left, "<<_mapTagToSTPData.size()<<" jobs to send,"
                                          <<" allocated jobs send "<<AllocatedSTPsSend<<" allocated jobs receive "<<AllocatedSTPsReceive<< " iterated through "<<i<<" keys");
#if defined(OffloadingCheckForSlowOperations)
  timing += MPI_Wtime();
  if(timing > OFFLOADING_SLOW_OPERATION_THRESHOLD)
    logError("cleanUpStaleTaskOutcomes()", " took "<<timing<<"s");
#endif

}

size_t exahype::solvers::ADERDGSolver::getAdditionalCurrentMemoryUsageReplication() {
  size_t sizePerSTP = sizeof(MigratablePredictionJobData)
                    + sizeof(double) * ( getDataPerCell() + getUpdateSize() + getBndTotalSize() + getBndFluxTotalSize() );
  return (AllocatedSTPsSend + AllocatedSTPsReceive) * sizePerSTP;
}

bool exahype::solvers::ADERDGSolver::tryToFindAndExtractOutcome(
    CellDescription& cellDescription,
    double predictionTimeStamp,
    double predictorTimeStepSize,
    DeliveryStatus &status,
    MigratablePredictionJobData **outcome) {


  //Caution: calls progress and may delay calling thread significantly
#if !defined(OffloadingUseProgressThread)
  exahype::solvers::ADERDGSolver::progressOffloading(this, false, MAX_PROGRESS_ITS);
#endif

  tarch::la::Vector<DIMENSIONS, double> center;
  center = cellDescription.getOffset() + 0.5 * cellDescription.getSize();

  //logInfo("tryToFindAndExtractOutcome()", "looking for center[0] = "<<center[0]
  //                                       <<" center[1] = "<<center[1]
  //                                       <<" timestamp = "<<predictionTimeStamp
  //                                       <<" time step = "<<predictorTimeStepSize);

  MigratablePredictionJobOutcomeKey key(center.data(), predictionTimeStamp, predictorTimeStepSize, 0); //todo: verify that element is always 0
  bool found = _outcomeDatabase.tryFindAndExtractOutcome(key, outcome, status);

  if(found && status==DeliveryStatus::Transit) {
    _outcomeDatabase.insertOutcome(key, *outcome, DeliveryStatus::Transit);
    outcome = nullptr;
    found = false;
  }
  else if(found){
    logDebug("tryToFindAndExtractOutcome()",
        "team "<<exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber()
        <<" found STP in received jobs:"
        <<(*outcome)->_metadata.to_string());
  }
  return found;
}

bool exahype::solvers::ADERDGSolver::tryToFindAndExtractOutcome(
    int cellDescriptionsIndex,
    int element,
    double predictionTimeStamp,
    double predictorTimeStepSize,
    DeliveryStatus &status,
    MigratablePredictionJobData **outcome) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,
      element);

  return tryToFindAndExtractOutcome(cellDescription, predictionTimeStamp, predictorTimeStepSize, status, outcome);
}

void exahype::solvers::ADERDGSolver::storePendingOutcomeToBeShared(MigratablePredictionJob *job) {
  //create copy
  MigratablePredictionJobData *data = new MigratablePredictionJobData(*this);
  AllocatedSTPsSend++;

  auto& cellDescription = getCellDescription(job->_cellDescriptionsIndex, job->_element);
  double *luh   = static_cast<double*>(cellDescription.getSolution());
  double *lduh   = static_cast<double*>(cellDescription.getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());
#if defined(OffloadingGradQhbnd)
  double *lGradQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient());
#endif

  logDebug("sendTaskOutcomeToOtherTeams","allocated STPs send "<<AllocatedSTPsSend );
  //logInfo("sendFullReplicatedSTPToOtherTeams", "allocated "<<sizeof(MigratablePredictionJobData)
  //                                                         +sizeof(double)*(data->_luh.size()+data->_lduh.size()+data->_lQhbnd.size()+data->_lFhbnd.size())<<" bytes ");
  std::memcpy(&data->_luh[0], luh, data->_luh.size()*sizeof(double));
  std::memcpy(&data->_lduh[0], lduh, data->_lduh.size()*sizeof(double));
  std::memcpy(&data->_lQhbnd[0], lQhbnd, data->_lQhbnd.size()*sizeof(double));
  std::memcpy(&data->_lFhbnd[0], lFhbnd, data->_lFhbnd.size()*sizeof(double));
#if OffloadingGradQhbnd
  std::memcpy(&data->_lGradQhbnd[0], lGradQhbnd, data->_lGradQhbnd.size()*sizeof(double));
#endif
  job->packMetaData(&data->_metadata);

  //outcome may already be available due to predictor re-run -> replace old data
  tbb::concurrent_hash_map<std::pair<int,int>, MigratablePredictionJobData*>::accessor accessor;
  bool found = _pendingOutcomesToBeShared.find(accessor, std::make_pair(job->_cellDescriptionsIndex,job->_element));
  //assert(!found);
  if(found) {
      MigratablePredictionJobData *data2;
      data2 = accessor->second;
      delete data2;
      assertion(data2!=nullptr);
      _pendingOutcomesToBeShared.erase(accessor);
      accessor.release();
      logInfo("storePendingOutcomeToBeShared", "replacing cellDesc ="<<job->_cellDescriptionsIndex
                                                <<" time step size "<<std::setprecision(30)<<job->_predictorTimeStepSize);
  }

  _pendingOutcomesToBeShared.insert(std::make_pair(std::make_pair(job->_cellDescriptionsIndex,job->_element), data));
}

void exahype::solvers::ADERDGSolver::storePendingOutcomeToBeShared(int cellDescriptionsIndex, int element, double timestamp, double timestep) {
  MigratablePredictionJobData *data = new MigratablePredictionJobData(*this);

  AllocatedSTPsSend++;
  auto& cellDescription = getCellDescription(cellDescriptionsIndex, element);
  double *luh   = static_cast<double*>(cellDescription.getSolution());
  double *lduh   = static_cast<double*>(cellDescription.getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

  std::memcpy(&data->_luh[0], luh, data->_luh.size()*sizeof(double));
  std::memcpy(&data->_lduh[0], lduh, data->_lduh.size()*sizeof(double));
  std::memcpy(&data->_lQhbnd[0], lQhbnd, data->_lQhbnd.size()*sizeof(double));
  std::memcpy(&data->_lFhbnd[0], lFhbnd, data->_lFhbnd.size()*sizeof(double));

  tarch::la::Vector<DIMENSIONS, double> center;
  center = cellDescription.getOffset() + 0.5 * cellDescription.getSize();
  tarch::la::Vector<DIMENSIONS, double> dx;
  dx =  cellDescription.getSize();

  for(int i = 0; i<DIMENSIONS; i++) {
    data->_metadata._dx[i] = dx[i];
    data->_metadata._center[i] = center[i];
  }

  data->_metadata._predictorTimeStamp = timestamp;
  data->_metadata._predictorTimeStepSize = timestep;
  data->_metadata._isCorrupted = false;
  data->_metadata._isPotSoftErrorTriggered = true; //will be re-set once limiter status is known

  //outcome may already be available due to predictor re-run -> replace old data
  tbb::concurrent_hash_map<std::pair<int,int>, MigratablePredictionJobData*>::accessor accessor;
  bool found = _pendingOutcomesToBeShared.find(accessor, std::make_pair(cellDescriptionsIndex, element));
  //assert(!found);
  if(found) {
    MigratablePredictionJobData *data2;
    data2 = accessor->second;
    delete data2;
    assertion(data2!=nullptr);
    _pendingOutcomesToBeShared.erase(accessor);
    accessor.release();
    logInfo("storePendingOutcomeToBeShared", "replacing cellDesc ="<<cellDescriptionsIndex
                                               <<" time step size "<<std::setprecision(30)<<timestep);
  }

  _pendingOutcomesToBeShared.insert(std::make_pair(std::make_pair( cellDescriptionsIndex, element), data));

}

void exahype::solvers::ADERDGSolver::releaseDummyOutcomeAndShare(int cellDescriptionsIndex, int element, double timestamp, double timestep, bool isTroubled) {
  MigratablePredictionJobData *data = new MigratablePredictionJobData(*this);

  AllocatedSTPsSend++;
  auto& cellDescription = getCellDescription(cellDescriptionsIndex, element);
  double *luh   = static_cast<double*>(cellDescription.getSolution());
  double *lduh   = static_cast<double*>(cellDescription.getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

  std::memcpy(&data->_luh[0], luh, data->_luh.size()*sizeof(double));
  std::memcpy(&data->_lduh[0], lduh, data->_lduh.size()*sizeof(double));
  std::memcpy(&data->_lQhbnd[0], lQhbnd, data->_lQhbnd.size()*sizeof(double));
  std::memcpy(&data->_lFhbnd[0], lFhbnd, data->_lFhbnd.size()*sizeof(double));

  logDebug("releaseDummyOutcomeAndShare","releasing initial dummy outcome");

  //set trigger -> need to reset limiter trigger, might wanna pull this out
  logDebug("releaseDummyOutcomeAndShare", " celldesc ="<<cellDescriptionsIndex<<" isPotCorrupted "<<isTroubled<<" time stamp "<<timestamp<<std::setprecision(30)<<" time step "<<timestep);

  tarch::la::Vector<DIMENSIONS, double> center;
  center = cellDescription.getOffset() + 0.5 * cellDescription.getSize();
  tarch::la::Vector<DIMENSIONS, double> dx;
  dx =  cellDescription.getSize();

  for(int i = 0; i<DIMENSIONS; i++) {
    data->_metadata._dx[i] = dx[i];
    data->_metadata._center[i] = center[i];
  }

  data->_metadata._predictorTimeStamp = timestamp;
  data->_metadata._predictorTimeStepSize = timestep;
  data->_metadata._isCorrupted = false;
  if(exahype::reactive::ResilienceTools::CheckLimitedCellsOnly)
    data->_metadata._isPotSoftErrorTriggered = isTroubled;

  //Share now
  int teams = exahype::reactive::ReactiveContext::getInstance().getTMPINumTeams();
  int interCommRank = exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber();
  MPI_Comm teamInterComm = exahype::reactive::ReactiveContext::getInstance().getTMPIInterTeamCommunicatorData();
  MPI_Request *sendRequests = new MPI_Request[(NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME+1)*(teams-1)];
  //int tag = job->_cellDescriptionsIndex; // exahype::reactive::OffloadingContext::getInstance().getOffloadingTag();
  int tag = exahype::reactive::ReactiveContext::getInstance().getOffloadingTag();
  _mapTagToSTPData.insert(std::make_pair(tag, data));

  if(data->_metadata._isCorrupted) {
    logWarning("releaseDummyOutcomeAndShare", "Caution: a corrupted outcome is shared. SDC should be detected...");
    if(!data->_metadata._isPotSoftErrorTriggered)
      logError("releaseDummyOutcomeAndShare","has not been detected by SDC mechanism, softErrorTriggered="<<data->_metadata._isPotSoftErrorTriggered);
  }

  int j = 0;
  for(int i=0; i<teams; i++) {
    if(i!=interCommRank) {
      logDebug("releaseDummyOutcomeAndShare"," team "<< interCommRank
                                               <<" send replica job: "
                                               << data->_metadata.to_string()
                                               <<" to team "<<i);
      mpiIsendMigratablePredictionJobOutcomeSolution(
                                    &(data->_luh[0]),
                                    &(data->_lduh[0]),
                                    &(data->_lQhbnd[0]),
                                    &(data->_lFhbnd[0]),
                                    &(data->_lGradQhbnd[0]),
                                    i,
                                    tag,
                                    teamInterComm,
                                    &sendRequests[(NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME+1)*j],
                                    &(data->_metadata));
      j++;
    }
  }
  SentSTPs++;
  exahype::reactive::RequestManager::getInstance().submitRequests(sendRequests,
                                                                       (teams-1)*(NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME+1),
                                                                       tag,
                                                                       exahype::reactive::RequestManager::MULTIPLE_SOURCES,
                                                                       MigratablePredictionJob::sendHandlerTaskSharing,
                                                                       exahype::reactive::RequestType::sendOutcome,
                                                                       this, MPI_BLOCKING);
  delete[] sendRequests;
}

void exahype::solvers::ADERDGSolver::releasePendingOutcomeAndShare(int cellDescriptionsIndex, int element, double timeStamp, double timeStepSize, bool isTroubled) {
  MigratablePredictionJobData *data = nullptr;
  tbb::concurrent_hash_map<std::pair<int,int>, MigratablePredictionJobData*>::accessor accessor;
  bool found = _pendingOutcomesToBeShared.find(accessor, std::make_pair(cellDescriptionsIndex, element));
  if(found) {
    data = accessor->second;
    assertion(data!=nullptr);
    _pendingOutcomesToBeShared.erase(accessor);
    accessor.release();

    auto& cellDescription = getCellDescription(cellDescriptionsIndex, element);
    //set trigger -> need to reset limiter trigger, might wanna pull this out
    logDebug("releasePendingOutcomeAndShare", "team = "<<exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber()
        <<" releasing cellDescriptionsIndex = "<<cellDescriptionsIndex<<" "<<cellDescription.toString())

    //need to reset this: with predictor re-runs, the solution belongs to a different time step size
    if(data->_metadata._predictorTimeStamp != timeStamp)
      logError("releasePendingOutcomeAndShare", "time stamp="<<timeStamp<<" predictorTimeStamp="<<data->_metadata._predictorTimeStamp);

    assert(data->_metadata._predictorTimeStamp == timeStamp);
    assert(data->_metadata._predictorTimeStepSize == timeStepSize);
    //data->_metadata._predictorTimeStamp = timeStamp;
    //data->_metadata._predictorTimeStepSize = timeStepSize;
    logDebug("releasePendingOutcomeAndShare", " celldesc ="<<cellDescriptionsIndex<<" isPotentiallyCorrupted "<<isTroubled
                                              <<" timestepsize "<<std::setprecision(30)<<data->_metadata._predictorTimeStepSize);
    if(exahype::reactive::ResilienceTools::CheckLimitedCellsOnly)
      data->_metadata._isPotSoftErrorTriggered = isTroubled;

    //update solution
    double *luh   = static_cast<double*>(cellDescription.getSolution());
    std::memcpy(&data->_luh[0], luh, data->_luh.size()*sizeof(double));

    //Share now
    int teams = exahype::reactive::ReactiveContext::getInstance().getTMPINumTeams();
    int interCommRank = exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber();
    MPI_Comm teamInterComm = exahype::reactive::ReactiveContext::getInstance().getTMPIInterTeamCommunicatorData();
    MPI_Request *sendRequests = new MPI_Request[(NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME+1)*(teams-1)];
    //int tag = job->_cellDescriptionsIndex; // exahype::reactive::OffloadingContext::getInstance().getOffloadingTag();
    int tag = exahype::reactive::ReactiveContext::getInstance().getOffloadingTag();
    _mapTagToSTPData.insert(std::make_pair(tag, data));

    if(data->_metadata._isCorrupted) {
      logWarning("releasePendingOutcomeAndShare", "Caution: a corrupted outcome is shared. SDC should be detected...");
      if(!data->_metadata._isPotSoftErrorTriggered)
        logError("releasePendingOutcomeAndShare","has not been detected by SDC mechanism, softErrorTriggered="<<data->_metadata._isPotSoftErrorTriggered);
    }

    int j = 0;
    for(int i=0; i<teams; i++) {
      if(i!=interCommRank) {
        logDebug("releasePendingOutcomeAndShare"," team "<< interCommRank
                                                 <<" send replica job: "
                                                 << data->_metadata.to_string()
                                                 <<" to team "<<i
                                                 <<" timestepsize "<<std::setprecision(30)<<data->_metadata._predictorTimeStepSize);
        mpiIsendMigratablePredictionJobOutcomeSolution(
                                    &(data->_luh[0]),
                                    &(data->_lduh[0]),
                                    &(data->_lQhbnd[0]),
                                    &(data->_lFhbnd[0]),
                                    &(data->_lGradQhbnd[0]),
                                    i,
                                    tag,
                                    teamInterComm,
                                    &sendRequests[(NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME+1)*j],
                                    &(data->_metadata));

        j++;
      }
    }
    SentSTPs++;
    exahype::reactive::RequestManager::getInstance().submitRequests(sendRequests,
                                                                    (teams-1)*(NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME+1),
                                                                    tag,
                                                                    exahype::reactive::RequestManager::MULTIPLE_SOURCES,
                                                                    MigratablePredictionJob::sendHandlerTaskSharing,
                                                                    exahype::reactive::RequestType::sendOutcome,
                                                                    this, MPI_BLOCKING);
    delete[] sendRequests;
  }
}

void exahype::solvers::ADERDGSolver::correctCellDescriptionWithOutcome(CellDescription& cellDescription, MigratablePredictionJobData *outcome) {
  double *lduh = static_cast<double*>(cellDescription.getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double *lGradQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient());
  double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

  //correct here
  std::memcpy(lduh, &outcome->_lduh[0], outcome->_lduh.size() * sizeof(double));
  std::memcpy(lQhbnd, &outcome->_lQhbnd[0], outcome->_lQhbnd.size() * sizeof(double));
  std::memcpy(lFhbnd, &outcome->_lFhbnd[0], outcome->_lFhbnd.size() * sizeof(double));
#if OffloadingGradQhbnd
  std::memcpy(lGradQhbnd, &outcome->_lGradQhbnd[0], outcome->_lGradQhbnd.size() * sizeof(double));
#endif

  logError("correctWithOutcome()","Corrected an error in STP.");

  exahype::reactive::ResilienceStatistics::getInstance().notifyHealedTask();
}

exahype::solvers::ADERDGSolver::SDCCheckResult exahype::solvers::ADERDGSolver::checkCellDescriptionAgainstOutcome(CellDescription& cellDescription, MigratablePredictionJobData *data){
  double *lduh = static_cast<double*>(cellDescription.getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
#if defined(OffloadingGradQhbnd)
  double *lGradQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient());
#endif
  double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

  tarch::la::Vector<DIMENSIONS, double> center = cellDescription.getOffset()+0.5*cellDescription.getSize();

  bool equal = true;
  bool tmp;

  tmp = data->_metadata._predictorTimeStamp == cellDescription.getTimeStamp(); equal &= tmp;
  tmp = data->_metadata._predictorTimeStepSize == cellDescription.getTimeStepSize(); equal &= tmp;

  tmp = exahype::reactive::ResilienceTools::getInstance().isAdmissibleNumericalError(data->_lQhbnd.data(), lQhbnd, data->_lQhbnd.size()); equal&=tmp;
  if(!tmp) {
    logError("checkAgainstOutcome", "lQhbnd is not (numerically) equal for cell "<<"center[0]="<<center[0]<<" center[1]="<<center[1]<<" timestamp "<<cellDescription.getTimeStamp());
  }

#if defined(OffloadingGradQhbnd)
  tmp = exahype::reactive::ResilienceTools::getInstance().isAdmissibleNumericalError(data->_lGradQhbnd.data(), lGradQhbnd, data->_lGradQhbnd.size()); equal&=tmp;
  if(!tmp) {
    logError("checkAgainstOutcome", "lGradQhbnd is not (numerically) equal for cell "<<"center[0]="<<center[0]<<" center[1]="<<center[1]<<" timestamp "<<cellDescription.getTimeStamp());
  }
#endif
  tmp = exahype::reactive::ResilienceTools::getInstance().isAdmissibleNumericalError(data->_lFhbnd.data(), lFhbnd, data->_lFhbnd.size()); equal&=tmp;
  if(!tmp) {
    logError("checkAgainstOutcome", "lFhbnd is not  (numerically) equal for cell "<<"center[0]="<<center[0]<<" center[1]="<<center[1]<<" timestamp "<<cellDescription.getTimeStamp());
  }
  tmp = exahype::reactive::ResilienceTools::getInstance().isAdmissibleNumericalError(data->_lduh.data(), lduh, data->_lduh.size()); equal&=tmp;
  if(!tmp) {
    logError("checkAgainstOutcome", "lduh is not  (numerically) equal for cell "<<"center[0]="<<center[0]<<" center[1]="<<center[1]<<" timestamp "<<cellDescription.getTimeStamp());
  }

  if(!equal) {
    logError("checkAgainstOutcome", "soft error detected: "<<data->_metadata.to_string());
    exahype::reactive::ResilienceStatistics::getInstance().notifyDetectedError();
  }

  exahype::reactive::ResilienceStatistics::getInstance().notifyDoubleCheckedTask();

  if(!equal) {
    if(data->_metadata._isPotSoftErrorTriggered)
      return SDCCheckResult::UncorrectableSoftError; //don't know which outcome is sane
    else {
      return SDCCheckResult::OutcomeSaneAsTriggerNotActive; //we assume that the trigger is good enough to tell us that the outcome is ok
    }
  }
  else
    return SDCCheckResult::NoCorruption;
}


///////////////////////////////////
// COMMUNICATION_ROUTINES
///////////////////////////////////
void exahype::solvers::ADERDGSolver::pollForOutstandingCommunicationRequests(exahype::solvers::ADERDGSolver *solver, bool calledOnMaster, int maxIts) {

#if defined(OffloadingCheckForSlowOperations)
  double timing = - MPI_Wtime();
#endif

#if defined(UseSmartMPI)
  MPI_Status_Offload stat, statMapped;
#else
  MPI_Status stat;
#if defined(OffloadingNoEarlyReceiveBacks)
  MPI_Status statMapped;
#endif
#endif
  int receivedTask = 0;
  int receivedTaskBack = 0;
  int msgLen = -1;
  //int lastTag = -1;
  //int lastSrc = -1;
  int lastRecvTag = -1;
  int lastRecvSrc = -1;
#if defined(OffloadingNoEarlyReceiveBacks)
  int lastRecvBackTag = -1;
  int lastRecvBackSrc = -1;
  MPI_Comm commMapped = exahype::reactive::ReactiveContext::getInstance().getMPICommunicatorMapped();
#endif
  MPI_Comm comm = exahype::reactive::ReactiveContext::getInstance().getMPICommunicator();
  int iprobesCounter = 0;
  int ierr;

  MPI_Comm interTeamComm = exahype::reactive::ReactiveContext::getInstance().getTMPIInterTeamCommunicatorData();
  int receivedReplicaTask = 0;
  MPI_Status statRepData;

#if defined(UseSmartMPI)
  MPI_CHECK("pollForOutstandingCommunicationRequests", MPI_Iprobe_offload(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &receivedTask, &stat));
  assertion(ierr==MPI_SUCCESS);
#else
  MPI_CHECK("pollForOutstandingCommunicationRequests", MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &receivedTask, &stat));
  assertion(ierr==MPI_SUCCESS);
#endif

#if defined(OffloadingNoEarlyReceiveBacks)
#if defined(UseSmartMPI)
  MPI_CHECK("pollForOutstandingCommunicationRequests", MPI_Iprobe_offload(MPI_ANY_SOURCE, MPI_ANY_TAG, commMapped, &receivedTaskBack, &statMapped));
  assertion(ierr==MPI_SUCCESS);
#else
  MPI_CHECK("pollForOutstandingCommunicationRequests", MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, commMapped, &receivedTaskBack, &statMapped));
  assertion(ierr==MPI_SUCCESS);
#endif
#endif

#if defined(UseSmartMPI)
  MPI_Status_Offload statRepDataOffload;
  MPI_CHECK("pollForOutstandingCommunicationRequests", MPI_Iprobe_offload(MPI_ANY_SOURCE, MPI_ANY_TAG, interTeamComm, &receivedReplicaTask, &statRepDataOffload));
  logDebug("progressOffloading", "Iprobe for replica task "<<receivedReplicaTask<<" statRepDataOffload.MPI_TAG="<<statRepDataOffload.MPI_TAG<<" statRepDataOffload.size = "<<statRepDataOffload.size);
#else
  MPI_CHECK("pollForOutstandingCommunicationRequests", MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, interTeamComm, &receivedReplicaTask, &statRepData));
#endif

  bool terminateImmediately = false;

  bool warningPrinted = false;
  double timing = -MPI_Wtime();

  while(
      (receivedTask || receivedTaskBack || receivedReplicaTask )
      && (iprobesCounter<MaxIprobesInOffloadingProgress || receivedReplicaTask)
      && !terminateImmediately
      && iprobesCounter<maxIts)
  {
    if((timing+MPI_Wtime()) >2 && !warningPrinted) {
      logError("pollForOutstanding", " warning: polling very long iprobes counter "<<iprobesCounter);
      warningPrinted = true;
    }

    iprobesCounter++;
    // RECEIVE TASK BACK
#if defined(OffloadingNoEarlyReceiveBacks)
    if(receivedTaskBack && (statMapped.MPI_TAG!=lastRecvBackTag || statMapped.MPI_SOURCE!=lastRecvBackSrc)) {
      lastRecvBackSrc = statMapped.MPI_SOURCE;
      lastRecvBackTag = statMapped.MPI_TAG;

      assertion(lastRecvBackTag!=solver->_lastReceiveBackTag[lastRecvBackSrc]);
      solver->_lastReceiveBackTag[lastRecvBackSrc]=lastRecvBackTag; 

#if defined(UseSmartMPI)
      receiveBackMigratableJob(statMapped.MPI_TAG, statMapped.MPI_SOURCE, solver, statMapped.rail);
#else
      receiveBackMigratableJob(statMapped.MPI_TAG, statMapped.MPI_SOURCE, solver);
#endif
    }
#if defined(UseSmartMPI)
    MPI_CHECK("pollForOutstandingCommunicationRequests", MPI_Iprobe_offload(MPI_ANY_SOURCE, MPI_ANY_TAG, commMapped, &receivedTaskBack, &statMapped));
#else
    MPI_CHECK("pollForOutstandingCommunicationRequests", MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, commMapped, &receivedTaskBack, &statMapped));
#endif
    assertion(ierr==MPI_SUCCESS);
#endif

#ifdef OffloadingUseProgressTask
    // special termination signal
    if(receivedTask && stat.MPI_TAG==0) {
       int terminatedSender = stat.MPI_SOURCE;
       //logInfo("progressOffloading()","active sender "<<terminatedSender<<" has sent termination signal ");
       exahype::reactive::ReactiveContext::getInstance().receiveCompleted(terminatedSender);
       ActiveSenders.erase(terminatedSender);
#if defined(UseSmartMPI)
       MPI_CHECK("pollForOutstandingCommunicationRequests", MPI_Iprobe_offload(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &receivedTask, &stat));
       assertion(ierr==MPI_SUCCESS);
#else
       MPI_CHECK("pollForOutstandingCommunicationRequests", MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &receivedTask, &stat));
       assertion(ierr==MPI_SUCCESS);
#endif
    }
#endif

    // RECEIVE TASK
    if(receivedTask) {
#ifdef OffloadingUseProgressTask
       //logInfo("progressOffloading()","inserting active sender "<<stat.MPI_SOURCE);
       ActiveSenders.insert(stat.MPI_SOURCE);
       if(NumberOfReceiveJobs==0) {
          NumberOfReceiveJobs++;
          assertion(NumberOfReceiveJobs<=1);
          //logInfo("progressOffloading()","spawning receive job, receive jobs "<<NumberOfReceiveJobs);
          ReceiveJob *receiveJob = new ReceiveJob(*solver);
          peano::datatraversal::TaskSet spawnedSet(receiveJob);     
          terminateImmediately = true; // we'll receive this task but then terminate to give the receive job the opportunity to run
       }
#endif

      exahype::reactive::ReactiveContext::getInstance().triggerVictimFlag();
      msgLen = -1;
#if defined(UseSmartMPI)
      MPI_Get_count_offload(&stat, MigratablePredictionJobMetaData::getMPIDatatype(), &msgLen);
#else
      MPI_Get_count(&stat, MigratablePredictionJobMetaData::getMPIDatatype(), &msgLen);
#endif
      // is this message metadata? -> if true, we are about to receive a new STP task
      if((size_t) msgLen==MigratablePredictionJobMetaData::getMessageLen() && !(lastRecvTag==stat.MPI_TAG && lastRecvSrc==stat.MPI_SOURCE)) {
        lastRecvTag = stat.MPI_TAG;
        lastRecvSrc = stat.MPI_SOURCE;
      
        assertion(solver->_lastReceiveTag[lastRecvSrc]!=lastRecvTag); //Todo: still necessary?
        solver->_lastReceiveTag[lastRecvSrc] = lastRecvTag;

#if defined(UseSmartMPI)
        receiveMigratableJob(stat.MPI_TAG, stat.MPI_SOURCE, solver, stat.rail);
#else
        receiveMigratableJob(stat.MPI_TAG, stat.MPI_SOURCE, solver);
#endif
      }
    }
#if defined(UseSmartMPI)
    MPI_CHECK("pollForOutstandingCommunicationRequests", MPI_Iprobe_offload(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &receivedTask, &stat));
#else
    MPI_CHECK("pollForOutstandingCommunicationRequests", MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &receivedTask, &stat));
#endif
    assertion(ierr==MPI_SUCCESS);
//#if defined(TaskSharing)
#if !defined(TaskSharingUseHandshake)
    if(receivedReplicaTask) {
      int msgLenDouble = -1;
#if defined UseSmartMPI
      MPI_Get_count_offload(&statRepDataOffload, MigratablePredictionJobMetaData::getMPIDatatype(), &msgLenDouble);
      if(msgLenDouble==MigratablePredictionJobMetaData::getMessageLen()) {
        logDebug("progressOffloading","received replica task from "<<statRepDataOffload.MPI_SOURCE<<" , tag "<<statRepDataOffload.MPI_TAG);
        assert(solver->_lastReceiveReplicaTag[statRepDataOffload.MPI_SOURCE]!=statRepDataOffload.MPI_TAG);
        solver->_lastReceiveReplicaTag[statRepDataOffload.MPI_SOURCE] = statRepDataOffload.MPI_TAG;
        receiveTaskOutcome(statRepDataOffload.MPI_TAG, statRepDataOffload.MPI_SOURCE, solver, statRepDataOffload.rail);
      }
#else
      MPI_Get_count(&statRepData, MigratablePredictionJobMetaData::getMPIDatatype(), &msgLenDouble);
      // is this message metadata? -> if true, we are about to receive a new STP task
      if((size_t) msgLenDouble==MigratablePredictionJobMetaData::getMessageLen()) {
        assertion(solver->_lastReceiveReplicaTag[statRepData.MPI_SOURCE]!=statRepData.MPI_TAG);
        solver->_lastReceiveReplicaTag[statRepData.MPI_SOURCE] = statRepData.MPI_TAG;
        logDebug("progressOffloading","team "<<exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber()<<" received replica task from "<<statRepData.MPI_SOURCE<<" , tag "<<statRepData.MPI_TAG);
        receiveTaskOutcome(statRepData.MPI_TAG, statRepData.MPI_SOURCE, solver);
      }
#endif /*UseSmartMPI */
    }
#endif /*UseSmartMPI */
#if defined(UseSmartMPI)
    MPI_CHECK("pollForOutstandingCommunicationRequests", MPI_Iprobe_offload(MPI_ANY_SOURCE, MPI_ANY_TAG, interTeamComm, &receivedReplicaTask, &statRepDataOffload));
    logDebug("progressOffloading", "Iprobe for replica task "<<receivedReplicaTask<<" statRepDataOffload.MPI_TAG="<<statRepDataOffload.MPI_TAG<<" statRepDataOffload.size = "<<statRepDataOffload.size);
#else
    MPI_CHECK("pollForOutstandingCommunicationRequests", MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, interTeamComm, &receivedReplicaTask, &statRepData));
#endif
    assertion(ierr==MPI_SUCCESS);
  }

#if defined(OffloadingCheckForSlowOperations)
  timing += MPI_Wtime();
  if(timing > OFFLOADING_SLOW_OPERATION_THRESHOLD)
    logError("pollForOutstandingCommunicationRequests()", " took "<<timing<<"s with iprobes = "<<iprobesCounter);
#endif
}

void exahype::solvers::ADERDGSolver::progressOffloading(exahype::solvers::ADERDGSolver* solver, bool runOnMaster, int maxIts) {
#if defined(OffloadingCheckForSlowOperations)
  double timing = -MPI_Wtime();
#endif

#if !defined(UseMPIThreadSplit) //skip spin lock when MPI thread split model is used
  bool canRun;
  tarch::multicore::Lock lock(OffloadingSemaphore, false);
  
#if defined(OffloadingUseProgressThread)
  if(runOnMaster)
    canRun = false;
  else
    canRun = lock.tryLock();
#else
  //assert(!runOnMaster);
  //if(tarch::multicore::Core::getInstance().getThreadNum()==0) {
  //  logError("progressOffloading", "Progress about to run on master thread!");
 // }
  //assertion(!runOnMaster);
  // First, we ensure here that only one thread at a time progresses offloading
  // this avoids multithreaded MPI problems
  canRun = lock.tryLock();
#endif

  if(!canRun) {
    return;
  }
#endif

#ifdef USE_ITAC
  //VT_begin(event_progress);
#endif

#ifdef OffloadingUseProfiler
  double timing = -MPI_Wtime();
  exahype::reactive::OffloadingProfiler::getInstance().beginProgress();
#endif

  /*tarch::timing::Watch watch("ADERDGSolver","progress", false, false);
  if(tarch::multicore::Core::getInstance().getThreadNum()==0) {
     watch.startTimer(); 
     std::cout<<"Error!!!!!!"<<std::endl;
     lock.free();
     return;
  }*/

  // 2. make progress on any outstanding MPI communication
  //if(!runOnMaster)
#ifdef OffloadingUseProfiler
  double timing_requests = -MPI_Wtime();
  exahype::reactive::OffloadingProfiler::getInstance().beginProgressRequests();
#endif
  exahype::reactive::RequestManager::getInstance().progressRequests();
#ifdef OffloadingUseProfiler
  timing_requests += MPI_Wtime();
  exahype::reactive::OffloadingProfiler::getInstance().endProgressRequests(timing_requests);
#endif

  // 3. progress on performance monitor
  exahype::reactive::PerformanceMonitor::getInstance().run();

  // 4. detect whether local rank should receive anything
  //if(!runOnMaster)
#ifdef OffloadingUseProfiler
  double timing_poll = -MPI_Wtime();
  exahype::reactive::OffloadingProfiler::getInstance().beginPolling();
#endif
  pollForOutstandingCommunicationRequests(solver, runOnMaster, maxIts);
#ifdef OffloadingUseProfiler
  timing_poll += MPI_Wtime();
  exahype::reactive::OffloadingProfiler::getInstance().endPolling(timing_poll);
#endif

#ifdef OffloadingUseProfiler
  timing += MPI_Wtime();
  exahype::reactive::OffloadingProfiler::getInstance().endProgress(timing);
#endif

#if !defined(UseMPIThreadSplit)
  lock.free();
#endif
  
  /*if(tarch::multicore::Core::getInstance().getThreadNum()==0) {
    watch.stopTimer();
    if(watch.getCalendarTime()>1)
      std::cout<<"Master blocked for "<<watch.getCalendarTime()<< " s "<<std::endl;
  }*/
 
#ifdef USE_ITAC
  //VT_end(event_progress);
#endif

#if defined(OffloadingCheckForSlowOperations)
  timing += MPI_Wtime();
  if(timing > OFFLOADING_SLOW_OPERATION_THRESHOLD)
    logError("progressOffloading()", " took "<<timing<<"s");
#endif
}

bool exahype::solvers::ADERDGSolver::tryToReceiveTaskBack(exahype::solvers::ADERDGSolver* solver, const void* cellDescription) {

  // First, we ensure here that only one thread at a time progresses offloading
  // this attempts to avoid multithreaded MPI problems
  tarch::multicore::Lock lock(OffloadingSemaphore, false);
  bool canRun = lock.tryLock();
  if(!canRun) {
#if defined(PerformanceAnalysisOffloadingDetailed)
    watch.stopTimer();
    if(watch.getCalendarTime() >= 0.0) {
      logDebug(
          "progressOffloading() ",
          "couldn't run "<<
          "time=" << std::fixed <<
          watch.getCalendarTime() <<
          ", cpu time=" <<
          watch.getCPUTime()
      );
    }
#endif
    return true;
  }
  //Todo (Philipp): fix when no early receive backs are active
  return exahype::reactive::RequestManager::getInstance().progressReceiveBackRequests();

//todo(Philipp): this code won't be executed anymore
#if defined(OffloadingNoEarlyReceiveBacks)
  int tag, srcRank, myRank;
  myRank = tarch::parallel::Node::getInstance().getRank();
  if(cellDescription==nullptr) {
    tag = MPI_ANY_TAG;
    srcRank = MPI_ANY_SOURCE; 
  }
  else {
    solver->getResponsibleRankTagForCellDescription(cellDescription, srcRank, tag);
    if(srcRank==myRank) {
      lock.free();
      return false;
    }
    //logInfo("tryToReceiveTaskBack()","probing for tag "<<tag<<" from rank "<<srcRank);
  }
  
  MPI_Comm commMapped = exahype::reactive::ReactiveContext::getInstance().getMPICommunicatorMapped();
  int receivedTaskBack = 0;
#if defined(UseSmartMPI)
  //todo: implement
#else
  MPI_Status statMapped;
  int ierr = MPI_Iprobe(srcRank, tag, commMapped, &receivedTaskBack, &statMapped);
  assertion(ierr==MPI_SUCCESS);
  if(receivedTaskBack) {
      //exahype::reactive::OffloadingManager::getInstance().setRunningAndReceivingBack();
      tbb::concurrent_hash_map<int, CellDescription*>::accessor a_tagToCellDesc;
      bool found = solver->_mapTagToCellDesc.find(a_tagToCellDesc, statMapped.MPI_TAG);
      assertion(found);
      auto cellDescription = a_tagToCellDesc->second;
      a_tagToCellDesc.release();
      double *lduh   = static_cast<double*>(cellDescription->getUpdate());
      double *lQhbnd = static_cast<double*>(cellDescription->getExtrapolatedPredictor());
      double *lFhbnd = static_cast<double*>(cellDescription->getFluctuation());
      double *lGradQhbnd = static_cast<double*>(cellDescription->getExtrapolatedPredictorGradient());

      assertion(statMapped.MPI_TAG!=solver->_lastReceiveBackTag[statMapped.MPI_SOURCE]);
      solver->_lastReceiveBackTag[statMapped.MPI_SOURCE] = statMapped.MPI_TAG;

      MPI_Request recvRequests[NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME];
      solver->mpiIrecvMigratablePredictionJobOutcome(
        lduh, lQhbnd,
        lFhbnd, lGradQhbnd, statMapped.MPI_SOURCE, statMapped.MPI_TAG, commMapped, recvRequests);

      exahype::reactive::RequestManager::getInstance().submitRequests(
        recvRequests,
        NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME,
        statMapped.MPI_TAG,
        statMapped.MPI_SOURCE,
        exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveBackHandler,
        exahype::reactive::RequestType::receiveBack, solver, false);
      lock.free();
      return true;
  }
#endif
  // now, a different thread can progress the offloading
  lock.free();
  return false;
#else
  lock.free();
  return true;
#endif
  //return exahype::reactive::OffloadingManager::getInstance().hasOutstandingRequestOfType(exahype::reactive::RequestType::receiveBack);
}

void exahype::solvers::ADERDGSolver::sendTaskOutcomeToOtherTeams(MigratablePredictionJob *job) {

    int teams = exahype::reactive::ReactiveContext::getInstance().getTMPINumTeams();
    int interCommRank = exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber();
    MPI_Comm teamInterComm = exahype::reactive::ReactiveContext::getInstance().getTMPIInterTeamCommunicatorData();

    auto& cellDescription = getCellDescription(job->_cellDescriptionsIndex, job->_element);
    double *luh   = static_cast<double*>(cellDescription.getSolution());
    double *lduh   = static_cast<double*>(cellDescription.getUpdate());
    double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
    double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());
#if defined(OffloadingGradQhbnd)
    double *lGradQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient());
#endif

#if defined(UseSmartMPI)
    logDebug("sendTaskOutcomeToOtherTeams","allocated STPs send "<<AllocatedSTPsSend );

    MigratablePredictionJobMetaData *metadata = new MigratablePredictionJobMetaData();
    job->packMetaData(metadata);

    int tag = exahype::reactive::ReactiveContext::getInstance().getOffloadingTag();
    //int tag = job->_cellDescriptionsIndex; //exahype::reactive::OffloadingContext::getInstance().getOffloadingTag();
    //_mapTagToReplicationSendData.insert(std::make_pair(tag, data));

    bool hasSent = false; //indicates whether at least one send was successful
    int j = 0;
    for(int i=0; i<teams; i++) {
      if(i!=interCommRank) {
        logDebug("sendTaskOutcomeToOtherTeams"," team "<< interCommRank
                                                         <<" send replica job: "
                                                         <<metadata->to_string()
                                                         <<" time stamp = "<<job->_predictorTimeStamp
                                                         <<" to team "<<i);

        hasSent |= mpiSendMigratablePredictionJobOutcomeOffload(&lduh[0],
                                           &lQhbnd[0],
                                           &lFhbnd[0],
                                           &lGradQhbnd[0],
                                           i,
                                           tag,
                                           teamInterComm,
                                           metadata);
        j++;
      }
    }
    if(hasSent) {
      SentSTPs++;
      CompletedSentSTPs++;
      exahype::reactive::ResilienceStatistics::getInstance().notifySentTask();
    }

    delete metadata;
#else
    //create copy
    MigratablePredictionJobData *data = new MigratablePredictionJobData(*this);
    AllocatedSTPsSend++;

    logDebug("sendTaskOutcomeToOtherTeams","allocated STPs send "<<AllocatedSTPsSend );
    //logInfo("sendFullReplicatedSTPToOtherTeams", "allocated "<<sizeof(MigratablePredictionJobData)
    //                                                         +sizeof(double)*(data->_luh.size()+data->_lduh.size()+data->_lQhbnd.size()+data->_lFhbnd.size())<<" bytes ");
    std::memcpy(&data->_luh[0], luh, data->_luh.size()*sizeof(double));
    std::memcpy(&data->_lduh[0], lduh, data->_lduh.size()*sizeof(double));
    std::memcpy(&data->_lQhbnd[0], lQhbnd, data->_lQhbnd.size()*sizeof(double));
    std::memcpy(&data->_lFhbnd[0], lFhbnd, data->_lFhbnd.size()*sizeof(double));
#if OffloadingGradQhbnd
    std::memcpy(&data->_lGradQhbnd[0], lGradQhbnd, data->_lGradQhbnd.size()*sizeof(double));
#endif
    //double *metadata = new double[2*DIMENSIONS+2];
    job->packMetaData(&data->_metadata);

    MPI_Request *sendRequests = new MPI_Request[(NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME+1)*(teams-1)];

    //int tag = job->_cellDescriptionsIndex; // exahype::reactive::OffloadingContext::getInstance().getOffloadingTag();
    int tag = exahype::reactive::ReactiveContext::getInstance().getOffloadingTag();

    _mapTagToSTPData.insert(std::make_pair(tag, data));

    double time = -MPI_Wtime();

    int j = 0;
    for(int i=0; i<teams; i++) {
      if(i!=interCommRank) {
          logDebug("sendTaskOutcomeToOtherTeams"," team "<< interCommRank
                                                   <<" send replica job: "
                                                   << data->_metadata.to_string()
                                                   <<" time stamp = "<<job->_predictorTimeStamp
                                                   <<" to team "<<i);
          mpiIsendMigratablePredictionJobOutcomeSolution(
                                      &(data->_luh[0]),
                                      &(data->_lduh[0]),
                                      &(data->_lQhbnd[0]),
                                      &(data->_lFhbnd[0]),
                                      &(data->_lGradQhbnd[0]),
                                      i,
                                      tag,
                                      teamInterComm,
                                      &sendRequests[(NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME+1)*j],
                                      &(data->_metadata));
                                      j++;
      }
    }

    time += MPI_Wtime();
    if(time>0.02)
      logError("sendTaskOutcome","took too long "<<time<<" AllocatedSTPsSend "<<AllocatedSTPsSend);
    SentSTPs++;
    exahype::reactive::RequestManager::getInstance().submitRequests(sendRequests,
                                                                         (teams-1)*(NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME+1),
                                                                         tag,
                                                                         exahype::reactive::RequestManager::MULTIPLE_SOURCES,
                                                                         MigratablePredictionJob::sendHandlerTaskSharing,
                                                                         exahype::reactive::RequestType::sendOutcome,
                                                                         this, MPI_BLOCKING);
   delete[] sendRequests;
#endif
}

void exahype::solvers::ADERDGSolver::receiveMigratableJob(int tag, int src, exahype::solvers::ADERDGSolver *solver, int rail) {
#if !defined(UseSmartMPI) || defined(SmartMPINB)
  MPI_Request receiveRequests[NUM_REQUESTS_MIGRATABLE_COMM+1];
#endif
  int ierr;
  MigratablePredictionJobData *data = new MigratablePredictionJobData(*solver);
  solver->_mapTagRankToStolenData.insert(std::make_pair(std::make_pair(src, tag), data));
#if defined(UseSmartMPI)
#if !defined(SmartMPINB)
  solver->mpiRecvMigratablePredictionJobOffload(
       data->_luh.data(),
       src,
       tag,
       exahype::reactive::ReactiveContext::getInstance().getMPICommunicator(),
       rail,
       &(data->_metadata));
  MigratablePredictionJob::receiveHandler(solver, tag, src);
#else
  solver->mpiIrecvMigratablePredictionJobOffload(
       data->_luh.data(),
       src,
       tag,
       exahype::reactive::ReactiveContext::getInstance().getMPICommunicator(),
       rail,
       &receiveRequests[0],
       &(data->_metadata));

  exahype::reactive::RequestManager::getInstance().submitRequests(
       receiveRequests,
       NUM_REQUESTS_MIGRATABLE_COMM+1,
       tag,
       src,
       MigratablePredictionJob::receiveHandler,
       exahype::reactive::RequestType::receive,
       solver,
       false);
#endif
#else
  solver->mpiIrecvMigratablePredictionJob(
       data->_luh.data(),
       src,
       tag,
       exahype::reactive::ReactiveContext::getInstance().getMPICommunicator(),
       &receiveRequests[0],
       &(data->_metadata));
   //double wtime = -MPI_Wtime();
  int canComplete = 0;
  MPI_CHECK("receiveMigratableJob",MPI_Testall(NUM_REQUESTS_MIGRATABLE_COMM+1, &receiveRequests[0], &canComplete, MPI_STATUSES_IGNORE));
  assertion(ierr==MPI_SUCCESS);
  if(canComplete)
    MigratablePredictionJob::receiveHandler(solver, tag, src);
  else {
    if(tarch::multicore::jobs::getNumberOfWaitingBackgroundJobs()<=1) {
       //logInfo("progressOffloading()","running out of tasks and could not receive stolen task so we just block!");
       double wtime = -MPI_Wtime();
       exahype::reactive::RequestManager::getInstance().submitRequests(
           receiveRequests,
           NUM_REQUESTS_MIGRATABLE_COMM+1,
           tag,
           src,
           MigratablePredictionJob::receiveHandler,
           exahype::reactive::RequestType::receive,
           solver,
           true);
       wtime+= MPI_Wtime();
       if(wtime>0.01)
         logDebug("progressOffloading()","blocking for stolen task took too long:"<<wtime<<"s");
    }
    else {
       exahype::reactive::RequestManager::getInstance().submitRequests(
           receiveRequests,
           NUM_REQUESTS_MIGRATABLE_COMM+1,
           tag,
           src,
           MigratablePredictionJob::receiveHandler,
           exahype::reactive::RequestType::receive,
           solver,
           false);
    }
  }
#endif
}

void exahype::solvers::ADERDGSolver::receiveBackMigratableJob(int tag, int src, exahype::solvers::ADERDGSolver *solver, int rail) {

  MPI_Comm commMapped = exahype::reactive::ReactiveContext::getInstance().getMPICommunicatorMapped();

#if defined (OffloadingLocalRecompute)
  //Todo (Philipp): we actually may not need to transfer back metadata as it may be available locally
  tbb::concurrent_hash_map<int, MigratablePredictionJobData*>::accessor a_tagToData;
  //bool found = solver->_mapTagToSTPData.find(a_tagToData, tag);
  assertion(found);
  //MigratablePredictionJobData *data = a_tagToData->second;
  a_tagToData.release();

  MigratablePredictionJobData *data = new MigratablePredictionJobData(*solver);
  AllocatedSTPsReceive++;
  solver->_mapTagToSTPData.insert(std::make_pair(tag, data));

  tbb::concurrent_hash_map<int, CellDescription*>::accessor a_tagToCellDesc;
  bool found = solver->_mapTagToCellDesc.find(a_tagToCellDesc, tag);
  assertion(found);
  auto cellDescription = a_tagToCellDesc->second;
  a_tagToCellDesc.release();
  double *lduh   = static_cast<double*>(cellDescription->getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription->getExtrapolatedPredictor());
  double *lFhbnd = static_cast<double*>(cellDescription->getFluctuation());
  double *lGradQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient());

#if defined(UseSmartMPI)
#if defined(SmartMPINB)
  MPI_Request recvRequests[NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME];

  solver->mpiIrecvMigratablePredictionJobOutcomeOffload(
      &(data->_lduh[0]),
      &(data->_lQhbnd[0]),
      &(data->_lFhbnd[0]),
      &(data->_lGradQhbnd[0]),
      src,
      tag,
      commMapped,
      rail,
      recvRequests);

  exahype::reactive::RequestManager::getInstance().submitRequests(
      recvRequests,
      NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME, //5,
      tag,
      src,
      exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveBackHandler,
      exahype::reactive::RequestType::receiveBack,
      solver,
      false);
#else
  solver->mpiRecvMigratablePredictionJobOutcomeOffload(
      &(data->_lduh[0]),
      &(data->_lQhbnd[0]),
      &(data->_lFhbnd[0]),
      &(data->_lGradQhbnd[0]),
      src,
      tag,
      commMapped,
      rail);
  MigratablePredictionJob::receiveBackHandler(solver, tag, src);
#endif /*SmartMPINB*/
#else
  MPI_Request recvRequests[NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME];
  solver->mpiIrecvMigratablePredictionJobOutcome(
      &(data->_lduh[0]),
      &(data->_lQhbnd[0]),
      &(data->_lFhbnd[0]),
      &(data->_lGradQhbnd[0]),
      src,
      tag,
      commMapped,
      recvRequests);
   //   &(data->_metadata[0]));

  exahype::reactive::RequestManager::getInstance().submitRequests(
      recvRequests,
      NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME, //5,
      tag,
      src,
      exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveBackHandler,
      exahype::reactive::RequestType::receiveBack,
      solver,
      false);
#endif

  /*tbb::concurrent_hash_map<int, CellDescription*>::accessor a_tagToCellDesc;
  found = solver->_mapTagToCellDesc.find(a_tagToCellDesc, tag);
  assertion(found);
  auto cellDescription = a_tagToCellDesc->second;
  a_tagToCellDesc.release();
  double *luh    = static_cast<double*>(cellDescription->getSolution());
  double *lduh   = static_cast<double*>(cellDescription->getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription->getExtrapolatedPredictor());
  double *lFhbnd = static_cast<double*>(cellDescription->getFluctuation());

  MPI_Request recvRequests[5];
   solver->mpiIrecvMigratablePredictionJob(
       luh,
       lduh,
       lQhbnd,
       lFhbnd,
       src,
       tag,
       commMapped,
       recvRequests,
       &(data->_metadata[0]));

   exahype::reactive::RequestManager::getInstance().submitRequests(
       recvRequests,
       5,
       tag,
       src,
       exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveBackHandler,
       exahype::reactive::RequestType::receiveBack,
       solver,
       false);*/


#else
  tbb::concurrent_hash_map<int, CellDescription*>::accessor a_tagToCellDesc;
  //logInfo("receiveBackMigratableJob", "looking for tag"<<tag);
  bool found = solver->_mapTagToCellDesc.find(a_tagToCellDesc, tag);
  if(!found)
    logError("receiveBackMigratableJob","Inconsistent maps, couldn't find cell descriptions\n");
  assertion(found);
  auto cellDescription = a_tagToCellDesc->second;
  a_tagToCellDesc.release();
  double *lduh   = static_cast<double*>(cellDescription->getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription->getExtrapolatedPredictor());
  double *lFhbnd = static_cast<double*>(cellDescription->getFluctuation());
  double *lGradQhbnd = static_cast<double*>(cellDescription->getExtrapolatedPredictorGradient());

#if defined(UseSmartMPI)
#if defined(SmartMPINB)
  MPI_Request recvRequests[NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME];

  solver->mpiIrecvMigratablePredictionJobOutcomeOffload(
      lduh,
      lQhbnd,
      lFhbnd,
      lGradQhbnd,
      src,
      tag,
      commMapped,
      rail,
      recvRequests);

  exahype::reactive::RequestManager::getInstance().submitRequests(
      recvRequests,
      NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME, //5,
      tag,
      src,
      exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveBackHandler,
      exahype::reactive::RequestType::receiveBack,
      solver,
      false);
#else
  solver->mpiRecvMigratablePredictionJobOutcomeOffload(
    lduh,
    lQhbnd,
    lFhbnd,
    lGradQhbnd,
    src,
    tag,
    commMapped,
    rail);
  MigratablePredictionJob::receiveBackHandler(solver, tag, src);
#endif
#else
  MPI_Request recvRequests[NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME];
  solver->mpiIrecvMigratablePredictionJobOutcome(
    lduh,
    lQhbnd,
          lFhbnd,
    lGradQhbnd,
    src,
    tag,
    commMapped,
    recvRequests);
  exahype::reactive::RequestManager::getInstance().submitRequests(
      recvRequests,
      NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME,
      tag,
      src,
      exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveBackHandler,
      exahype::reactive::RequestType::receiveBack, solver, false);
#endif

#endif
}

void exahype::solvers::ADERDGSolver::receiveTaskOutcome(int tag, int src, exahype::solvers::ADERDGSolver *solver, int rail) {
  MPI_Comm interTeamComm = exahype::reactive::ReactiveContext::getInstance().getTMPIInterTeamCommunicatorData();

  MigratablePredictionJobData *data = new MigratablePredictionJobData(*solver);
  AllocatedSTPsReceive++;

#ifdef UseSmartMPI

#if defined(UseSmartMPINB) || defined(ResilienceHealing)
#error "Not implemented yet!"
#endif

  solver->mpiRecvMigratablePredictionJobOutcomeOffload(
         data->_lduh.data(),
         data->_lQhbnd.data(),
         data->_lFhbnd.data(),
         data->_lGradQhbnd.data(),
         src,
         tag,
         interTeamComm,
         rail,
         &(data->_metadata));

  data->_metadata.unpackContiguousBuffer();

  MigratablePredictionJobOutcomeKey key(data->_metadata.getCenter(), data->_metadata.getPredictorTimeStamp(),
                                         data->_metadata.getPredictorTimeStepSize(), data->_metadata.getElement());
  if(key._timestamp<solver->getMinTimeStamp()) {
    exahype::reactive::ResilienceStatistics::getInstance().notifyLateTask();
    delete data;
    AllocatedSTPsReceive--;
  }
  else {
    MigratablePredictionJobData *data2 = nullptr;
    DeliveryStatus status;
    bool found = static_cast<exahype::solvers::ADERDGSolver*>(solver)->_outcomeDatabase.tryFindAndExtractOutcome(key, &data2, status);

    if(found && status==DeliveryStatus::Transit) {
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_outcomeDatabase.insertOutcome(key, data2, DeliveryStatus::Received);
    }
    else {
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_outcomeDatabase.insertOutcome(key, data, DeliveryStatus::Received);
    }
    static_cast<exahype::solvers::ADERDGSolver*>(solver)->_allocatedOutcomes.push_back(key);
  }
  exahype::reactive::ResilienceStatistics::getInstance().notifyReceivedTask();
#else
   //logInfo("progressOffloading", "allocated stps receive"<<AllocatedSTPsReceive);
  MPI_Request receiveReplicaRequests[NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME+1];

  solver->mpiIrecvMigratablePredictionJobOutcomeSolution(
         data->_luh.data(),
         data->_lduh.data(),
         data->_lQhbnd.data(),
         data->_lFhbnd.data(),
         data->_lGradQhbnd.data(),
         src,
         tag,
         interTeamComm,
         &receiveReplicaRequests[0],
         &(data->_metadata));

  solver->_mapTagRankToReplicaData.insert(std::make_pair(std::make_pair(src, tag), data));
  exahype::reactive::RequestManager::getInstance().submitRequests(
         receiveReplicaRequests,
         NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME+1,
         tag,
         src,
         MigratablePredictionJob::receiveHandlerTaskSharing,
         exahype::reactive::RequestType::receiveOutcome,
         solver,
         MPI_BLOCKING);
#endif
}

void exahype::solvers::ADERDGSolver::mpiIsendMigratablePredictionJob(
  double *luh,
  int dest,
  int tag,
  MPI_Comm comm,
  MPI_Request *requests,
  MigratablePredictionJobMetaData *metadata) {

#if defined(OffloadingCheckForSlowOperations)
  double timing = - MPI_Wtime();
#endif

  int i = 0;
  int ierr;
  //MPI_Comm comm = exahype::reactive::OffloadingManager::getInstance().getMPICommunicator();

  if(metadata != nullptr) {
    MPI_CHECK("mpiIsendMigratablePredictionJob", MPI_Isend(metadata->getMPIBuffer(), MigratablePredictionJobMetaData::getMessageLen(), MigratablePredictionJobMetaData::getMPIDatatype(), dest, tag, comm, &requests[i++]));
    assertion(ierr==MPI_SUCCESS);
    assertion(requests[i-1]!=MPI_REQUEST_NULL);
  }

  assertion(luh!=NULL);
  MPI_CHECK("mpiIsendMigratablePredictionJob", MPI_Isend(luh, getDataPerCell(), MPI_DOUBLE, dest, tag, comm, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);


#if defined(OffloadingCheckForSlowOperations)
  timing += MPI_Wtime();
  if(timing > OFFLOADING_SLOW_OPERATION_THRESHOLD)
    logError("mpiIsendMigratablePredictionJobOutcome()", " took "<<timing<<"s");
#endif

}

void exahype::solvers::ADERDGSolver::mpiIsendMigratablePredictionJobOutcome(
  double *lduh,
  double *lQhbnd,
  double *lFhbnd,
  double *lGradQhbnd,
  int dest,
  int tag,
  MPI_Comm comm,
  MPI_Request *requests,
  MigratablePredictionJobMetaData *metadata) {

#if defined(OffloadingCheckForSlowOperations)
  double timing = - MPI_Wtime();
#endif

  int i = 0;
  int ierr;
  //MPI_Comm comm = exahype::reactive::OffloadingManager::getInstance().getMPICommunicator();

  if(metadata != nullptr) {
    logDebug("mpiIsendDebug", metadata->to_string()<<" , len = "<<MigratablePredictionJobMetaData::getMessageLen()<<" dest "<<dest<<" tag "<<tag);
    MPI_CHECK("mpiIsendMigratablePredictionJobOutcome", MPI_Isend(metadata->getMPIBuffer(), MigratablePredictionJobMetaData::getMessageLen(), MigratablePredictionJobMetaData::getMPIDatatype(), dest, tag, comm, &requests[i++]));
    
    assertion(ierr==MPI_SUCCESS);
    assertion(requests[i-1]!=MPI_REQUEST_NULL);
  }

  assertion(lduh!=NULL);
  MPI_CHECK("mpiIsendMigratablePredictionJobOutcome", MPI_Isend(lduh, getUpdateSize(), MPI_DOUBLE, dest, tag, comm, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

  assertion(lQhbnd!=NULL);
  MPI_CHECK("mpiIsendMigratablePredictionJobOutcome", MPI_Isend(lQhbnd, getBndTotalSize(), MPI_DOUBLE, dest, tag, comm, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

  assertion(lFhbnd!=NULL);
  MPI_CHECK("mpiIsendMigratablePredictionJobOutcome", MPI_Isend(lFhbnd, getBndFluxTotalSize(), MPI_DOUBLE, dest, tag, comm, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

#if defined(OffloadingGradQhbnd)
  assertion(lGradQhbnd!=NULL);
  MPI_CHECK("mpiIsendMigratablePredictionJobOutcome", MPI_Isend(lGradQhbnd, getBndGradQTotalSize(), MPI_DOUBLE, dest, tag, comm, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);
#endif

#if defined(OffloadingCheckForSlowOperations)
  timing += MPI_Wtime();
  if(timing > OFFLOADING_SLOW_OPERATION_THRESHOLD)
    logError("mpiIsendMigratablePredictionJobOutcome()", " took "<<timing<<"s");
#endif
}

void exahype::solvers::ADERDGSolver::mpiIsendMigratablePredictionJobOutcomeSolution(
  double *luh,
  double *lduh,
  double *lQhbnd,
  double *lFhbnd,
  double *lGradQhbnd,
  int dest,
  int tag,
  MPI_Comm comm,
  MPI_Request *requests,
  MigratablePredictionJobMetaData *metadata) {

#if defined(OffloadingCheckForSlowOperations)
  double timing = - MPI_Wtime();
#endif

  int i = 0;
  int ierr;
  //MPI_Comm comm = exahype::reactive::OffloadingManager::getInstance().getMPICommunicator();

  if(metadata != nullptr) {
    logDebug("mpiIsendMigratablePredictionJobOutcomeSolution", metadata->to_string()<<" , len = "<<MigratablePredictionJobMetaData::getMessageLen()<<" dest "<<dest<<" tag "<<tag);
    MPI_CHECK("mpiIsendMigratablePredictionJobOutcomeSolution", MPI_Isend(metadata->getMPIBuffer(), MigratablePredictionJobMetaData::getMessageLen(), MigratablePredictionJobMetaData::getMPIDatatype(), dest, tag, comm, &requests[i++]));

    assertion(ierr==MPI_SUCCESS);
    assertion(requests[i-1]!=MPI_REQUEST_NULL);
  }

  assertion(luh!=NULL);
  MPI_CHECK("mpiIsendMigratablePredictionJobOutcomeSolution", MPI_Isend(luh, getDataPerCell(), MPI_DOUBLE, dest, tag, comm, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

  assertion(lduh!=NULL);
  MPI_CHECK("mpiIsendMigratablePredictionJobOutcomeSolution", MPI_Isend(lduh, getUpdateSize(), MPI_DOUBLE, dest, tag, comm, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

  assertion(lQhbnd!=NULL);
  MPI_CHECK("mpiIsendMigratablePredictionJobOutcomeSolution", MPI_Isend(lQhbnd, getBndTotalSize(), MPI_DOUBLE, dest, tag, comm, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

  assertion(lFhbnd!=NULL);
  MPI_CHECK("mpiIsendMigratablePredictionJobOutcomeSolution", MPI_Isend(lFhbnd, getBndFluxTotalSize(), MPI_DOUBLE, dest, tag, comm, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

#if defined(OffloadingGradQhbnd)
  assertion(lGradQhbnd!=NULL);
  MPI_CHECK("mpiIsendMigratablePredictionJobOutcomeSolution", MPI_Isend(lGradQhbnd, getBndGradQTotalSize(), MPI_DOUBLE, dest, tag, comm, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);
#endif

#if defined(OffloadingCheckForSlowOperations)
  timing += MPI_Wtime();
  if(timing > OFFLOADING_SLOW_OPERATION_THRESHOLD)
    logError("mpiIsendMigratablePredictionJobOutcomeSolution()", " took "<<timing<<"s");
#endif
}

void exahype::solvers::ADERDGSolver::mpiIrecvMigratablePredictionJob(
    double *luh,
    int srcRank,
    int tag,
    MPI_Comm comm,
    MPI_Request *requests,
    MigratablePredictionJobMetaData *metadata ) {
  int ierr;
  //MPI_Comm comm = exahype::reactive::OffloadingManager::getInstance().getMPICommunicator();
  int i = 0;

#if defined(OffloadingCheckForSlowOperations)
  double timing = - MPI_Wtime();
#endif
  //logInfo("mpiIrecvMigratablePredictionJob", "receiving job "<<tag<<" from srcRank "<<srcRank);

  if(metadata != nullptr) {
    MPI_CHECK("mpiIrecvMigratablePredictionJob", MPI_Irecv(metadata->getMPIBuffer(), MigratablePredictionJobMetaData::getMessageLen(), MigratablePredictionJobMetaData::getMPIDatatype(), srcRank, tag, comm, &requests[i++]));
    assertion(ierr==MPI_SUCCESS);
    assertion(requests[i-1]!=MPI_REQUEST_NULL);
  }

  assertion(luh!=NULL);
  MPI_CHECK("mpiIrecvMigratablePredictionJob", MPI_Irecv(luh, getDataPerCell(), MPI_DOUBLE, srcRank, tag, comm, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

#if defined(OffloadingCheckForSlowOperations)
  timing += MPI_Wtime();
  if(timing > OFFLOADING_SLOW_OPERATION_THRESHOLD)
    logError("mpiIrecvMigratablePredictionJob()", " took "<<timing<<"s");
#endif
}

void exahype::solvers::ADERDGSolver::mpiIrecvMigratablePredictionJobOutcomeSolution(
    double *luh,
    double *lduh,
    double *lQhbnd,
    double *lFhbnd,
    double *lGradQhbnd,
    int srcRank,
    int tag,
    MPI_Comm comm,
    MPI_Request *requests,
    MigratablePredictionJobMetaData *metadata ) {
  int ierr;
  //MPI_Comm comm = exahype::reactive::OffloadingManager::getInstance().getMPICommunicator();
  int i = 0;

#if defined(OffloadingCheckForSlowOperations)
  double timing = - MPI_Wtime();
#endif
  //logInfo("mpiIrecvMigratablePredictionJob", "receiving job "<<tag<<" from srcRank "<<srcRank);

  if(metadata != nullptr) {
    MPI_CHECK("mpiIrecvMigratablePredictionJobOutcomeSolution",  MPI_Irecv(metadata->getMPIBuffer(), MigratablePredictionJobMetaData::getMessageLen(), MigratablePredictionJobMetaData::getMPIDatatype(), srcRank, tag, comm, &requests[i++]));
    assertion(ierr==MPI_SUCCESS);
    assertion(requests[i-1]!=MPI_REQUEST_NULL);
  }

  assertion(luh!=NULL);
  MPI_CHECK("mpiIrecvMigratablePredictionJobOutcomeSolution", MPI_Irecv(luh, getDataPerCell(), MPI_DOUBLE, srcRank, tag, comm, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

  assertion(lduh!=NULL);
  MPI_CHECK("mpiIrecvMigratablePredictionJobOutcomeSolution", MPI_Irecv(lduh, getUpdateSize(), MPI_DOUBLE, srcRank, tag, comm, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

  assertion(lQhbnd!=NULL);
  MPI_CHECK("mpiIrecvMigratablePredictionJobOutcomeSolution", MPI_Irecv(lQhbnd, getBndTotalSize(), MPI_DOUBLE, srcRank, tag, comm, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

  assertion(lFhbnd!=NULL);
  MPI_CHECK("mpiIrecvMigratablePredictionJobOutcomeSolution", MPI_Irecv(lFhbnd, getBndFluxTotalSize(), MPI_DOUBLE, srcRank, tag, comm, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

#if defined(OffloadingGradQhbnd)
  assertion(lGradQhbnd!=NULL);
  MPI_CHECK("mpiIrecvMigratablePredictionJobOutcomeSolution", MPI_Irecv(lGradQhbnd, getBndGradQTotalSize(), MPI_DOUBLE, srcRank, tag, comm, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);
#endif

#if defined(OffloadingCheckForSlowOperations)
  timing += MPI_Wtime();
  if(timing > OFFLOADING_SLOW_OPERATION_THRESHOLD)
    logError("mpiIrecvMigratablePredictionJobOutcomeSolution()", " took "<<timing<<"s");
#endif
}

void exahype::solvers::ADERDGSolver::mpiIrecvMigratablePredictionJobOutcome(
    double *lduh,
    double *lQhbnd,
    double *lFhbnd,
	  double *lGradQhbnd,
    int srcRank,
    int tag,
    MPI_Comm comm,
    MPI_Request *requests,
    MigratablePredictionJobMetaData *metadata ) {
  int ierr;
  //MPI_Comm comm = exahype::reactive::OffloadingManager::getInstance().getMPICommunicator();
  int i = 0;

#if defined(OffloadingCheckForSlowOperations)
  double timing = - MPI_Wtime();
#endif
  //logInfo("mpiIrecvMigratablePredictionJob", "receiving job "<<tag<<" from srcRank "<<srcRank);

  if(metadata != nullptr) {
    MPI_CHECK("mpiIrecvMigratablePredictionJobOutcome",  MPI_Irecv(metadata->getMPIBuffer(), MigratablePredictionJobMetaData::getMessageLen(), MigratablePredictionJobMetaData::getMPIDatatype(), srcRank, tag, comm, &requests[i++]));
    assertion(ierr==MPI_SUCCESS);
    assertion(requests[i-1]!=MPI_REQUEST_NULL);
  }

  assertion(lduh!=NULL);
  MPI_CHECK("mpiIrecvMigratablePredictionJobOutcome", MPI_Irecv(lduh, getUpdateSize(), MPI_DOUBLE, srcRank, tag, comm, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

  assertion(lQhbnd!=NULL);
  MPI_CHECK("mpiIrecvMigratablePredictionJobOutcome", MPI_Irecv(lQhbnd, getBndTotalSize(), MPI_DOUBLE, srcRank, tag, comm, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

  assertion(lFhbnd!=NULL);
  MPI_CHECK("mpiIrecvMigratablePredictionJobOutcome", MPI_Irecv(lFhbnd, getBndFluxTotalSize(), MPI_DOUBLE, srcRank, tag, comm, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

#if defined(OffloadingGradQhbnd)
  assertion(lGradQhbnd!=NULL);
  MPI_CHECK("mpiIrecvMigratablePredictionJobOutcome", MPI_Irecv(lGradQhbnd, getBndGradQTotalSize(), MPI_DOUBLE, srcRank, tag, comm, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);
#endif

#if defined(OffloadingCheckForSlowOperations)
  timing += MPI_Wtime();
  if(timing > OFFLOADING_SLOW_OPERATION_THRESHOLD)
    logError("mpiIrecvMigratablePredictionJobOutcome()", " took "<<timing<<"s");
#endif
}

void exahype::solvers::ADERDGSolver::mpiRecvMigratablePredictionJobOutcome(
  double *lduh,
  double *lQhbnd,
  double *lFhbnd,
  double *lGradQhbnd,
  int srcRank,
  int tag,
  MPI_Comm comm,
  MigratablePredictionJobMetaData *metadata ) {
  int ierr;
  //MPI_Comm comm = exahype::reactive::OffloadingManager::getInstance().getMPICommunicator();
#if defined(OffloadingCheckForSlowOperations)
  double timing = - MPI_Wtime();
#endif

  if(metadata != nullptr) {
    MPI_CHECK("mpiRecvMigratablePredictionJobOutcome", MPI_Recv(metadata->getMPIBuffer(),  MigratablePredictionJobMetaData::getMessageLen(), MigratablePredictionJobMetaData::getMPIDatatype(), srcRank, tag, comm, MPI_STATUS_IGNORE));
    assertion(ierr==MPI_SUCCESS);
  }

  assertion(lduh!=NULL);
  MPI_CHECK("mpiRecvMigratablePredictionJobOutcome", MPI_Recv(lduh, getUpdateSize(), MPI_DOUBLE, srcRank, tag, comm, MPI_STATUS_IGNORE));
  assertion(ierr==MPI_SUCCESS);

  assertion(lQhbnd!=NULL);
  MPI_CHECK("mpiRecvMigratablePredictionJobOutcome", MPI_Recv(lQhbnd, getBndTotalSize(), MPI_DOUBLE, srcRank, tag, comm, MPI_STATUS_IGNORE));
  assertion(ierr==MPI_SUCCESS);

  assertion(lFhbnd!=NULL);
  MPI_CHECK("mpiRecvMigratablePredictionJobOutcome", MPI_Recv(lFhbnd, getBndFluxTotalSize(), MPI_DOUBLE, srcRank, tag, comm, MPI_STATUS_IGNORE));
  assertion(ierr==MPI_SUCCESS);

#if defined(OffloadingGradQhbnd)
  assertion(lGradQhbnd!=NULL);
  MPI_CHECK("mpiRecvMigratablePredictionJobOutcome", MPI_Recv(lGradQhbnd, getBndGradQTotalSize(), MPI_DOUBLE, srcRank, tag, comm, MPI_STATUS_IGNORE));
  assertion(ierr==MPI_SUCCESS);
#endif

#if defined(OffloadingCheckForSlowOperations)
  timing += MPI_Wtime();
  if(timing > OFFLOADING_SLOW_OPERATION_THRESHOLD)
    logError("mpiRecvMigratablePredictionJobOutcome()", " took "<<timing<<"s");
#endif

}

#if defined(UseSmartMPI)
void exahype::solvers::ADERDGSolver::mpiIrecvMigratablePredictionJobOffload(
    double *luh,
    int srcRank,
    int tag,
    MPI_Comm comm,
    int rail,
    MPI_Request *requests,
    MigratablePredictionJobMetaData *metadata ) {
  int ierr;
  //MPI_Comm comm = exahype::reactive::OffloadingManager::getInstance().getMPICommunicator();
  int i = 0;
  
#if defined(OffloadingCheckForSlowOperations)
  double timing = - MPI_Wtime();
#endif

  //todo: won't work as SmartMPI doesn't support derived datatypes
  if(metadata != nullptr) {
    MPI_CHECK("mpiIrecvMigratablePredictionJobOffload", MPI_Irecv_offload(metadata->getContiguousBuffer(), MigratablePredictionJobMetaData::getMessageLen(), MigratablePredictionJobMetaData::getMPIDatatype(), srcRank, tag, comm, rail, &requests[i++]));
    assertion(ierr==MPI_SUCCESS);
  }

  assertion(luh!=NULL);
  MPI_CHECK("mpiIrecvMigratablePredictionJobOffload", MPI_Irecv_offload(luh, getDataPerCell(), MPI_DOUBLE, srcRank, tag, comm, rail, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);

#if defined(OffloadingCheckForSlowOperations)
  timing += MPI_Wtime();
  if(timing > OFFLOADING_SLOW_OPERATION_THRESHOLD)
    logError("mpiIrecvMigratablePredictionJobOutcomeOffload()", " took "<<timing<<"s");
#endif

}

void exahype::solvers::ADERDGSolver::mpiRecvMigratablePredictionJobOffload(
    double *luh,
    int srcRank,
    int tag,
    MPI_Comm comm,
    int rail,
    MigratablePredictionJobMetaData *metadata ) {
  int ierr;
  MPI_Status_Offload stat;

#if defined(OffloadingCheckForSlowOperations)
  double timing = - MPI_Wtime();
#endif

  //todo: won't work as SmartMPI doesn't support derived datatypes
  if(metadata != nullptr) {
    MPI_CHECK("mpiRecvMigratablePredictionJobOffload", MPI_Recv_offload(metadata->getContiguousBuffer(), MigratablePredictionJobMetaData::getMessageLen(), MigratablePredictionJobMetaData::getMPIDatatype(), srcRank, tag, comm, &stat, rail));
    assertion(ierr==MPI_SUCCESS);
  }

  assertion(luh!=NULL);
  MPI_CHECK("mpiRecvMigratablePredictionJobOffload", MPI_Recv_offload(luh, getDataPerCell(), MPI_DOUBLE, srcRank, tag, comm, &stat, rail));
  assertion(ierr==MPI_SUCCESS);

#if defined(OffloadingCheckForSlowOperations)
  timing += MPI_Wtime();
  if(timing > OFFLOADING_SLOW_OPERATION_THRESHOLD)
    logError("mpiRecvMigratablePredictionJobOutcomeOffload()", " took "<<timing<<"s");
#endif
}

void exahype::solvers::ADERDGSolver::mpiSendMigratablePredictionJobOffload(
  double *luh,
  int dest,
  int tag,
  MPI_Comm comm,
  MigratablePredictionJobMetaData *metadata) {

  int ierr;
  //MPI_Comm comm = exahype::reactive::OffloadingManager::getInstance().getMPICommunicator();
  //int tid = tarch::multicore::Core::getInstance().getThreadNum();

#if defined(OffloadingUseProfiler)
  double time = -MPI_Wtime();
  exahype::reactive::OffloadingProfiler::getInstance().beginCommunication();
#endif
#if defined(OffloadingCheckForSlowOperations)
  double timing = - MPI_Wtime();
#endif

  int rail = get_next_rail();

  if(metadata != nullptr) {
    //ierr = MPI_Send_offload(metadata, 2*DIMENSIONS+3, MPI_DOUBLE, dest, tag, comm, tid);
    MPI_CHECK("mpiSendMigratablePredictionJobOffload", MPI_Send_offload(metadata->getContiguousBuffer(), MigratablePredictionJobMetaData::getMessageLen(), MigratablePredictionJobMetaData::getMPIDatatype(), dest, tag, comm, rail));
    assertion(ierr==MPI_SUCCESS);
  }

  assertion(luh!=NULL);
  //ierr = MPI_Send_offload(lduh, getUpdateSize(), MPI_DOUBLE, dest, tag, comm, tid);
  MPI_CHECK("mpiSendMigratablePredictionJobOffload", MPI_Send_offload(luh, getDataPerCell(), MPI_DOUBLE, dest, tag, comm, rail));
  assertion(ierr==MPI_SUCCESS);

#if defined(OffloadingUseProfiler)
  time += MPI_Wtime();
  exahype::reactive::OffloadingProfiler::getInstance().endCommunication(true, time);
#endif
#if defined(OffloadingCheckForSlowOperations)
  timing += MPI_Wtime();
  if(timing > OFFLOADING_SLOW_OPERATION_THRESHOLD)
    logError("mpiSendMigratablePredictionJobOutcomeOffload()", " took "<<timing<<"s");
#endif
}

void exahype::solvers::ADERDGSolver::mpiIsendMigratablePredictionJobOffload(
  double *luh,
  int dest,
  int tag,
  MPI_Comm comm,
  MPI_Request *requests,
  MigratablePredictionJobMetaData *metadata) {

  int i = 0;
  int ierr;
  //MPI_Comm comm = exahype::reactive::OffloadingManager::getInstance().getMPICommunicator();
  //int tid = tarch::multicore::Core::getInstance().getThreadNum();

#if defined(OffloadingUseProfiler)
  double time = -MPI_Wtime();
  exahype::reactive::OffloadingProfiler::getInstance().beginCommunication();
#endif
#if defined(OffloadingCheckForSlowOperations)
  double timing = - MPI_Wtime();
#endif

  int rail = get_next_rail();

  if(metadata != nullptr) {
    //ierr = MPI_Send_offload(metadata, 2*DIMENSIONS+3, MPI_DOUBLE, dest, tag, comm, tid);
    MPI_CHECK("mpiIsendMigratablePredictionJobOffload", MPI_Isend_offload(metadata->getContiguousBuffer(), MigratablePredictionJobMetaData::getMessageLen(), MigratablePredictionJobMetaData::getMPIDatatype(), dest, tag, comm, rail, &requests[i++]));
    assertion(ierr==MPI_SUCCESS);
  }

  assertion(luh!=NULL);
  //ierr = MPI_Send_offload(lduh, getUpdateSize(), MPI_DOUBLE, dest, tag, comm, tid);
  MPI_CHECK("mpiIsendMigratablePredictionJobOffload", MPI_Isend_offload(luh, getDataPerCell(), MPI_DOUBLE, dest, tag, comm, rail, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);

#if defined(OffloadingUseProfiler)
  time += MPI_Wtime();
  exahype::reactive::OffloadingProfiler::getInstance().endCommunication(true, time);
#endif
#if defined(OffloadingCheckForSlowOperations)
  timing += MPI_Wtime();
  if(timing > OFFLOADING_SLOW_OPERATION_THRESHOLD)
    logError("mpiIsendMigratablePredictionJobOffload()", " took "<<timing<<"s");
#endif
}

void exahype::solvers::ADERDGSolver::mpiRecvMigratablePredictionJobOutcomeOffload(
    double *lduh,
    double *lQhbnd,
    double *lFhbnd,
    double *lGradQhbnd,
    int srcRank,
    int tag,
    MPI_Comm comm,
    int rail,
    MigratablePredictionJobMetaData *metadata ) {
  int ierr;
  //MPI_Comm comm = exahype::reactive::OffloadingManager::getInstance().getMPICommunicator();
  MPI_Status_Offload stat;

#if defined(OffloadingCheckForSlowOperations)
  double timing = - MPI_Wtime();
#endif

#if defined(OffloadingUseProfiler)
  double time = -MPI_Wtime();
  exahype::reactive::OffloadingProfiler::getInstance().beginCommunication();
#endif

  if(metadata != nullptr) {
    MPI_CHECK("mpiRecvMigratablePredictionJobOutcomeOffload", MPI_Recv_offload(metadata->getContiguousBuffer(),  MigratablePredictionJobMetaData::getMessageLen(), MigratablePredictionJobMetaData::getMPIDatatype(), srcRank, tag, comm, &stat, rail));
    assertion(ierr==MPI_SUCCESS);
  }

  assertion(lduh!=NULL);
  MPI_CHECK("mpiRecvMigratablePredictionJobOutcomeOffload", MPI_Recv_offload(lduh, getUpdateSize(), MPI_DOUBLE, srcRank, tag, comm, &stat, rail));
  assertion(ierr==MPI_SUCCESS);

  assertion(lQhbnd!=NULL);
  MPI_CHECK("mpiRecvMigratablePredictionJobOutcomeOffload", MPI_Recv_offload(lQhbnd, getBndTotalSize(), MPI_DOUBLE, srcRank, tag, comm, &stat, rail));
  assertion(ierr==MPI_SUCCESS);

  assertion(lFhbnd!=NULL);
  MPI_CHECK("mpiRecvMigratablePredictionJobOutcomeOffload", MPI_Recv_offload(lFhbnd, getBndFluxTotalSize(), MPI_DOUBLE, srcRank, tag, comm, &stat, rail));
  assertion(ierr==MPI_SUCCESS);

#if defined(OffloadingGradQhbnd)
  assertion(lGradQhbnd!=NULL);
  MPI_CHECK("mpiRecvMigratablePredictionJobOutcomeOffload", MPI_Recv_offload(lGradQhbnd, getBndGradQTotalSize(), MPI_DOUBLE, srcRank, tag, comm, &stat, rail));
  assertion(ierr==MPI_SUCCESS);
#endif

#if defined(OffloadingUseProfiler)
  time += MPI_Wtime();
  exahype::reactive::OffloadingProfiler::getInstance().endCommunication(true, time);
#endif
  
#if defined(OffloadingCheckForSlowOperations)
  timing += MPI_Wtime();
  if(timing > OFFLOADING_SLOW_OPERATION_THRESHOLD)
    logError("mpiRecvMigratablePredictionJobOutcomeOffload()", " took "<<timing<<"s");
#endif
}

void exahype::solvers::ADERDGSolver::mpiIrecvMigratablePredictionJobOutcomeOffload(
    double *lduh,
    double *lQhbnd,
    double *lFhbnd,
    double *lGradQhbnd,
    int srcRank,
    int tag,
    MPI_Comm comm,
    int rail,
    MPI_Request *requests,
    MigratablePredictionJobMetaData *metadata ) {
  int ierr;
  //MPI_Comm comm = exahype::reactive::OffloadingManager::getInstance().getMPICommunicator();
  int i = 0;

#if defined(OffloadingCheckForSlowOperations)
  double timing = - MPI_Wtime();
#endif
  //logInfo("mpiIrecvMigratablePredictionJob", "receiving job "<<tag<<" from srcRank "<<srcRank);

  if(metadata != nullptr) {
    MPI_CHECK("mpiIrecvMigratablePredictionJobOutcomeOffload",  MPI_Irecv_offload(metadata->getMPIBuffer(), MigratablePredictionJobMetaData::getMessageLen(), MigratablePredictionJobMetaData::getMPIDatatype(), srcRank, tag, comm, rail, &requests[i++]));
    assertion(ierr==MPI_SUCCESS);
    assertion(requests[i-1]!=MPI_REQUEST_NULL);
  }

  assertion(lduh!=NULL);
  MPI_CHECK("mpiIrecvMigratablePredictionJobOutcomeOffload", MPI_Irecv_offload(lduh, getUpdateSize(), MPI_DOUBLE, srcRank, tag, comm, rail, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

  assertion(lQhbnd!=NULL);
  MPI_CHECK("mpiIrecvMigratablePredictionJobOutcomeOffload", MPI_Irecv_offload(lQhbnd, getBndTotalSize(), MPI_DOUBLE, srcRank, tag, comm, rail, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

  assertion(lFhbnd!=NULL);
  MPI_CHECK("mpiIrecvMigratablePredictionJobOutcomeOffload", MPI_Irecv_offload(lFhbnd, getBndFluxTotalSize(), MPI_DOUBLE, srcRank, tag, comm, rail, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

#if defined(OffloadingGradQhbnd)
  assertion(lGradQhbnd!=NULL);
  MPI_CHECK("mpiIrecvMigratablePredictionJobOutcomeOffload", MPI_Irecv_offload(lGradQhbnd, getBndGradQTotalSize(), MPI_DOUBLE, srcRank, tag, comm, rail, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);
#endif

#if defined(OffloadingCheckForSlowOperations)
  timing += MPI_Wtime();
  if(timing > OFFLOADING_SLOW_OPERATION_THRESHOLD)
    logError("mpiIrecvMigratablePredictionJobOutcomeOffload()", " took "<<timing<<"s");
#endif
}

void exahype::solvers::ADERDGSolver::mpiIsendMigratablePredictionJobOutcomeOffload(
  double *lduh,
  double *lQhbnd,
  double *lFhbnd,
  double *lGradQhbnd,
  int dest,
  int tag,
  MPI_Comm comm,
  MPI_Request *requests,
  MigratablePredictionJobMetaData *metadata) {

#if defined(OffloadingCheckForSlowOperations)
  double timing = - MPI_Wtime();
#endif

  int rail = get_next_rail();
  int i = 0;
  int ierr;
  //MPI_Comm comm = exahype::reactive::OffloadingManager::getInstance().getMPICommunicator();

  if(metadata != nullptr) {
    logDebug("mpiIsendDebug", metadata->to_string()<<" , len = "<<MigratablePredictionJobMetaData::getMessageLen()<<" dest "<<dest<<" tag "<<tag);
    MPI_CHECK("mpiIsendMigratablePredictionJobOutcomeOffload", MPI_Isend_offload(metadata->getMPIBuffer(), MigratablePredictionJobMetaData::getMessageLen(), MigratablePredictionJobMetaData::getMPIDatatype(), dest, tag, comm, rail, &requests[i++]));
    
    assertion(ierr==MPI_SUCCESS);
    assertion(requests[i-1]!=MPI_REQUEST_NULL);
  }

  assertion(lduh!=NULL);
  MPI_CHECK("mpiIsendMigratablePredictionJobOutcomeOffload", MPI_Isend_offload(lduh, getUpdateSize(), MPI_DOUBLE, dest, tag, comm, rail, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

  assertion(lQhbnd!=NULL);
  MPI_CHECK("mpiIsendMigratablePredictionJobOutcomeOffload", MPI_Isend_offload(lQhbnd, getBndTotalSize(), MPI_DOUBLE, dest, tag, comm, rail, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

  assertion(lFhbnd!=NULL);
  MPI_CHECK("mpiIsendMigratablePredictionJobOutcomeOffload", MPI_Isend_offload(lFhbnd, getBndFluxTotalSize(), MPI_DOUBLE, dest, tag, comm, rail, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);

#if defined(OffloadingGradQhbnd)
  assertion(lGradQhbnd!=NULL);
  MPI_CHECK("mpiIsendMigratablePredictionJobOutcomeOffload", MPI_Isend_offload(lGradQhbnd, getBndGradQTotalSize(), MPI_DOUBLE, dest, tag, comm, rail, &requests[i++]));
  assertion(ierr==MPI_SUCCESS);
  assertion(requests[i-1]!=MPI_REQUEST_NULL);
#endif

#if defined(OffloadingCheckForSlowOperations)
  timing += MPI_Wtime();
  if(timing > OFFLOADING_SLOW_OPERATION_THRESHOLD)
    logError("mpiIsendMigratablePredictionJobOutcomeOffload()", " took "<<timing<<"s");
#endif
}

bool exahype::solvers::ADERDGSolver::mpiSendMigratablePredictionJobOutcomeOffload(
  double *lduh,
  double *lQhbnd,
  double *lFhbnd,
  double *lGradQhbnd,
  int dest,
  int tag,
  MPI_Comm comm,
  MigratablePredictionJobMetaData *metadata) {

  int ierr;
  //int tid = tarch::multicore::Core::getInstance().getThreadNum();

#if defined(OffloadingCheckForSlowOperations)
  double timing = - MPI_Wtime();
#endif

#if defined(OffloadingUseProfiler)
  double time = -MPI_Wtime();
  exahype::reactive::OffloadingProfiler::getInstance().beginCommunication();
#endif

  int rail = get_next_rail();
  //int rail = tid;
  //int rail = 0;

  if(!is_allowed_to_send_to_destination(dest, rail)) {
    logInfo("mpiSendMigratablePredictionJobOutcomeOffload", " Decided to not send a task outcome as BlueField advised against it!");
    return false;
  }

  if(metadata != nullptr) {
    //ierr = MPI_Send_offload(metadata, 2*DIMENSIONS+3, MPI_DOUBLE, dest, tag, comm, tid);
    ierr = MPI_Send_offload(metadata->getContiguousBuffer(), MigratablePredictionJobMetaData::getMessageLen(), MigratablePredictionJobMetaData::getMPIDatatype(), dest, tag, comm, rail);
    assertion(ierr==MPI_SUCCESS);
  }

  assertion(lduh!=NULL);
  //ierr = MPI_Send_offload(lduh, getUpdateSize(), MPI_DOUBLE, dest, tag, comm, tid);
  MPI_CHECK("mpiSendMigratablePredictionJobOutcomeOffload", MPI_Send_offload(lduh, getUpdateSize(), MPI_DOUBLE, dest, tag, comm, rail));
  assertion(ierr==MPI_SUCCESS);

  assertion(lQhbnd!=NULL);
  //ierr = MPI_Send_offload(lQhbnd, getBndTotalSize(), MPI_DOUBLE, dest, tag, comm, tid);
  MPI_CHECK("mpiSendMigratablePredictionJobOutcomeOffload", MPI_Send_offload(lQhbnd, getBndTotalSize(), MPI_DOUBLE, dest, tag, comm, rail));
  assertion(ierr==MPI_SUCCESS);

  assertion(lFhbnd!=NULL);
  //ierr = MPI_Send_offload(lFhbnd, getBndFluxTotalSize(), MPI_DOUBLE, dest, tag, comm, tid);
  MPI_CHECK("mpiSendMigratablePredictionJobOutcomeOffload", MPI_Send_offload(lFhbnd, getBndFluxTotalSize(), MPI_DOUBLE, dest, tag, comm, rail));
  assertion(ierr==MPI_SUCCESS);

#if defined(OffloadingGradQhbnd)
  assertion(lGradQhbnd!=NULL);
  MPI_CHECK("mpiSendMigratablePredictionJobOutcomeOffload", MPI_Send_offload(lGradQhbnd, getBndGradQTotalSize(), MPI_DOUBLE, dest, tag, comm, rail));
  assertion(ierr==MPI_SUCCESS);
#endif

#if defined(OffloadingUseProfiler)
  time += MPI_Wtime();
  exahype::reactive::OffloadingProfiler::getInstance().endCommunication(true, time);
#endif
  
#if defined(OffloadingCheckForSlowOperations)
  timing += MPI_Wtime();
  if(timing > OFFLOADING_SLOW_OPERATION_THRESHOLD)
    logError("mpiSendMigratablePredictionJobOutcomeOffload()", " took "<<timing<<"s");
#endif

  return true;
}

#endif

void exahype::solvers::ADERDGSolver::finishOutstandingInterTeamCommunication () {
  MPI_Comm interTeamComm = exahype::reactive::ReactiveContext::getInstance().getTMPIInterTeamCommunicatorData();

  while(exahype::reactive::RequestManager::getInstance().hasOutstandingRequestOfType(exahype::reactive::RequestType::sendOutcome)
    || exahype::reactive::RequestManager::getInstance().hasOutstandingRequestOfType(exahype::reactive::RequestType::receiveOutcome) ) {
    progressOffloading(this, false, std::numeric_limits<int>::max());
  }
  MPI_Request request;

  MPI_Ibarrier(interTeamComm, &request);
  int finished = 0;
  while(!finished) {
    progressOffloading(this, false, std::numeric_limits<int>::max());
    MPI_Test(&request, &finished, MPI_STATUS_IGNORE);
  }
}
#endif

exahype::solvers::ADERDGSolver::CompressionJob::CompressionJob(
  const ADERDGSolver& solver,
  CellDescription&    cellDescription,
  const bool          isSkeletonJob):
  tarch::multicore::jobs::Job(tarch::multicore::jobs::JobType::BackgroundTask,0,getTaskPriority(isSkeletonJob)),
  _solver(solver),
  _cellDescription(cellDescription),
  _isSkeletonJob(isSkeletonJob) {
  if (_isSkeletonJob) {
    NumberOfSkeletonJobs.fetch_add(1);
  } else {
    NumberOfEnclaveJobs.fetch_add(1);
  }
}

bool exahype::solvers::ADERDGSolver::CompressionJob::run(bool runOnMasterThread) {
  _solver.determineUnknownAverages(_cellDescription);
  _solver.computeHierarchicalTransform(_cellDescription,-1.0);
  _solver.putUnknownsIntoByteStream(_cellDescription);

  if (_isSkeletonJob) {
    NumberOfSkeletonJobs.fetch_sub(1);
    assertion( NumberOfSkeletonJobs.load()>=0 );
  } else {
    NumberOfEnclaveJobs.fetch_sub(1);
    assertion( NumberOfEnclaveJobs.load()>=0 );
  }
  return false;
}

void exahype::solvers::ADERDGSolver::compress( CellDescription& cellDescription, const bool isSkeletonCell ) const {
  assertion1( cellDescription.getCompressionState() ==  CellDescription::Uncompressed, cellDescription.toString() );
  if (CompressionAccuracy>0.0) {
    if ( SpawnCompressionAsBackgroundJob ) {
      cellDescription.setCompressionState(CellDescription::CurrentlyProcessed);
      CompressionJob compressionJob( *this, cellDescription, isSkeletonCell );
      assertionMsg( false, "this call is invalid" );
/*
      if ( isSkeletonCell ) {
        peano::datatraversal::TaskSet spawnedSet( compressionJob, peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible  );
      } else {
        peano::datatraversal::TaskSet spawnedSet( compressionJob, peano::datatraversal::TaskSet::TaskType::Background  );
      }
*/
    }
    else {
      determineUnknownAverages(cellDescription);
      computeHierarchicalTransform(cellDescription,-1.0);
      putUnknownsIntoByteStream(cellDescription);
      cellDescription.setCompressionState(CellDescription::Compressed);
    }
  }
}

void exahype::solvers::ADERDGSolver::uncompress(CellDescription& cellDescription) const {
  #ifdef SharedMemoryParallelisation
  bool madeDecision = CompressionAccuracy<=0.0;
  bool uncompress   = false;

  while (!madeDecision) {
    peano::datatraversal::TaskSet::finishToProcessBackgroundJobs();

    tarch::multicore::Lock lock(exahype::HeapSemaphore);
    madeDecision = cellDescription.getCompressionState() != CellDescription::CurrentlyProcessed;
    uncompress   = cellDescription.getCompressionState() == CellDescription::Compressed;
    if (uncompress) {
      cellDescription.setCompressionState( CellDescription::CurrentlyProcessed );
    }
    lock.free();
  }
  #else
  bool uncompress = CompressionAccuracy>0.0
      && cellDescription.getCompressionState() == CellDescription::Compressed;
  #endif

  if (uncompress) {
    pullUnknownsFromByteStream(cellDescription);
    computeHierarchicalTransform(cellDescription,1.0);

    tarch::multicore::Lock lock(exahype::HeapSemaphore);
      cellDescription.setCompressionState(CellDescription::Uncompressed);
    lock.free();
  }
}

void exahype::solvers::ADERDGSolver::determineUnknownAverages(
  CellDescription& cellDescription) const {
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getUpdateIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getFluctuationIndex()), cellDescription.toString() );

  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getSolutionAveragesIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionAveragesIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getUpdateAveragesIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorAveragesIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getFluctuationAveragesIndex()), cellDescription.toString() );

  const int dataPerNode  = getNumberOfParameters()+getNumberOfVariables();
  const int nodesPerCell = getDataPerCell()/ dataPerNode;
  const int nodesPerFace = getDataPerFace() / dataPerNode;

  double* solutionAverages              = static_cast<double*>(cellDescription.getSolutionAverages());
  double* previousSolutionAverage       = static_cast<double*>(cellDescription.getPreviousSolutionAverages());
  double* updateAverages                = static_cast<double*>(cellDescription.getUpdateAverages());
  double* extrapolatedPredictorAverages = static_cast<double*>(cellDescription.getExtrapolatedPredictorAverages());
  double* fluctuationAverages           = static_cast<double*>(cellDescription.getFluctuationAverages());

  double* solution              = static_cast<double*>(cellDescription.getSolution());
  double* previousSolution      = static_cast<double*>(cellDescription.getPreviousSolution());
  double* update                = static_cast<double*>(cellDescription.getUpdate());
  double* extrapolatedPredictor = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double* fluctuation           = static_cast<double*>(cellDescription.getFluctuation());

  // patch data
  kernels::idx2 idx_cellData    (nodesPerCell,dataPerNode);
  kernels::idx2 idx_cellUnknowns(nodesPerCell,getNumberOfVariables());
  for (int i=0; i<nodesPerCell; i++) {
    for (int variableNumber=0; variableNumber<dataPerNode; variableNumber++) { // variables+parameters
      solutionAverages[variableNumber]        += solution        [idx_cellData(i,variableNumber)];
      previousSolutionAverage[variableNumber] += previousSolution[idx_cellData(i,variableNumber)];
    }
    for (int variableNumber=0; variableNumber<getNumberOfVariables(); variableNumber++) { // variables
      updateAverages[variableNumber]          += update[idx_cellUnknowns(i,variableNumber)];
    }
  }
  for (int variableNumber=0; variableNumber<dataPerNode; variableNumber++) { // variables+parameters
    solutionAverages[variableNumber]        = solutionAverages[variableNumber]        / (double) nodesPerCell;
    previousSolutionAverage[variableNumber] = previousSolutionAverage[variableNumber] / (double) nodesPerCell;
  }
  for (int variableNumber=0; variableNumber<getNumberOfVariables(); variableNumber++) { // variables
    updateAverages[variableNumber]          = updateAverages[variableNumber]          / (double) nodesPerCell;
  }

  // face data
  kernels::idx2 idx_faceDataAvg    (DIMENSIONS_TIMES_TWO,dataPerNode);
  kernels::idx2 idx_faceUnknownsAvg(DIMENSIONS_TIMES_TWO,getNumberOfVariables());
  kernels::idx3 idx_faceData       (DIMENSIONS_TIMES_TWO,nodesPerFace,dataPerNode);
  kernels::idx3 idx_faceUnknowns   (DIMENSIONS_TIMES_TWO,nodesPerFace,getNumberOfVariables());
  for (int face=0; face<2*DIMENSIONS; face++) {
    for (int i=0; i<nodesPerFace; i++) {
      for (int variableNumber=0; variableNumber<dataPerNode; variableNumber++) { // variables+parameters
        extrapolatedPredictorAverages[idx_faceDataAvg(face,variableNumber)] +=
            extrapolatedPredictor[idx_faceData(face,i,variableNumber)];
      }
      for (int variableNumber=0; variableNumber<getNumberOfVariables(); variableNumber++) { // variables
        fluctuationAverages[idx_faceUnknownsAvg(face,variableNumber)] +=
            fluctuation[idx_faceUnknowns(face,i,variableNumber)];
      }
    }
    for (int variableNumber=0; variableNumber<dataPerNode; variableNumber++) { // variables+parameters
      extrapolatedPredictorAverages[idx_faceDataAvg(face,variableNumber)] =
          extrapolatedPredictorAverages[idx_faceDataAvg(face,variableNumber)] / (double) nodesPerFace;
    }
    for (int variableNumber=0; variableNumber<getNumberOfVariables(); variableNumber++) { // variables
      fluctuationAverages[idx_faceUnknownsAvg(face,variableNumber)] =
          fluctuationAverages[idx_faceUnknownsAvg(face,variableNumber)]       / (double) nodesPerFace;
    }
  }
}

void exahype::solvers::ADERDGSolver::computeHierarchicalTransform(
    CellDescription& cellDescription, double sign) const {
  const int dataPerNode  = getNumberOfParameters()+getNumberOfVariables();
  const int nodesPerCell = getDataPerCell()/ dataPerNode;
  const int nodesPerFace = getDataPerFace() / dataPerNode;

  double* solutionAverages              = static_cast<double*>(cellDescription.getSolutionAverages());
  double* previousSolutionAverage       = static_cast<double*>(cellDescription.getPreviousSolutionAverages());
  double* updateAverages                = static_cast<double*>(cellDescription.getUpdateAverages());
  double* extrapolatedPredictorAverages = static_cast<double*>(cellDescription.getExtrapolatedPredictorAverages());
  double* fluctuationAverages           = static_cast<double*>(cellDescription.getFluctuationAverages());

  double* solution              = static_cast<double*>(cellDescription.getSolution());
  double* previousSolution      = static_cast<double*>(cellDescription.getPreviousSolution());
  double* update                = static_cast<double*>(cellDescription.getUpdate());
  double* extrapolatedPredictor = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double* fluctuation           = static_cast<double*>(cellDescription.getFluctuation());

  // patch data
  kernels::idx2 idx_cellData    (nodesPerCell,dataPerNode);
  kernels::idx2 idx_cellUnknowns(nodesPerCell,getNumberOfVariables());
  for (int i=0; i<nodesPerCell; i++) {
    for (int variableNumber=0; variableNumber<dataPerNode; variableNumber++) { // variables+parameters
      solution        [idx_cellData(i,variableNumber)] += sign * solutionAverages[variableNumber];
      previousSolution[idx_cellData(i,variableNumber)] += sign * previousSolutionAverage[variableNumber];
    }
    for (int variableNumber=0; variableNumber<getNumberOfVariables(); variableNumber++) { // variables
      update[idx_cellUnknowns(i,variableNumber)] += sign * updateAverages[variableNumber];
    }
  }

  // face data
  kernels::idx2 idx_faceDataAvg    (DIMENSIONS_TIMES_TWO,dataPerNode);
  kernels::idx2 idx_faceUnknownsAvg(DIMENSIONS_TIMES_TWO,getNumberOfVariables());
  kernels::idx3 idx_faceData       (DIMENSIONS_TIMES_TWO,nodesPerFace,dataPerNode);
  kernels::idx3 idx_faceUnknowns   (DIMENSIONS_TIMES_TWO,nodesPerFace,getNumberOfVariables());
  for (int face=0; face<DIMENSIONS_TIMES_TWO; face++) {
    for (int i=0; i<nodesPerFace; i++) {
      for (int variableNumber=0; variableNumber<dataPerNode; variableNumber++) {  // variables+parameters
        extrapolatedPredictor[idx_faceData(face,i,variableNumber)] +=
            sign * extrapolatedPredictorAverages[idx_faceDataAvg(face,variableNumber)];
      }
      for (int variableNumber=0; variableNumber<getNumberOfVariables(); variableNumber++) {  // variables
        fluctuation[idx_faceUnknowns(face,i,variableNumber)] +=
            sign * fluctuationAverages[idx_faceUnknownsAvg(face,variableNumber)];
      }
    }
  }
}

void exahype::solvers::ADERDGSolver::putUnknownsIntoByteStream(
    CellDescription& cellDescription) const {
  assertion(CompressionAccuracy>0.0);

  assertion( cellDescription.getPreviousSolutionCompressedIndex()==-1 );
  assertion( cellDescription.getSolutionCompressedIndex()==-1 );
  assertion( cellDescription.getUpdateCompressedIndex()==-1 );
  assertion( cellDescription.getExtrapolatedPredictorCompressedIndex()==-1 );
  assertion( cellDescription.getFluctuationCompressedIndex()==-1 );

  int compressionOfPreviousSolution;
  int compressionOfSolution;
  int compressionOfUpdate;
  int compressionOfExtrapolatedPredictor;
  int compressionOfFluctuation;

  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getSolutionIndex() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionIndex() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getUpdateIndex() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedPredictorIndex() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getFluctuationIndex() ));

  peano::datatraversal::TaskSet compressionFactorIdentification(
    [&]() -> bool { compressionOfPreviousSolution = peano::heap::findMostAgressiveCompression(
      static_cast<double*>(cellDescription.getPreviousSolution()),
      getDataPerCell(),
      CompressionAccuracy,true
      );
      return false;
      },
    [&] () -> bool  { compressionOfSolution = peano::heap::findMostAgressiveCompression(
      static_cast<double*>(cellDescription.getSolution()),
      getDataPerCell(),
      CompressionAccuracy,true
      );
      return false;
      },
    [&]() -> bool  { compressionOfUpdate = peano::heap::findMostAgressiveCompression(
      static_cast<double*>(cellDescription.getUpdate()),
      getUnknownsPerCell(),
      CompressionAccuracy,true
      );
      return false;
      },
    [&]() -> bool  { compressionOfExtrapolatedPredictor = peano::heap::findMostAgressiveCompression(
      static_cast<double*>(cellDescription.getExtrapolatedPredictor()),
      getDataPerCellBoundary(),
      CompressionAccuracy,true
      );
      return false;
      },
    [&]() -> bool  { compressionOfFluctuation = peano::heap::findMostAgressiveCompression(
      static_cast<double*>(cellDescription.getFluctuation()),
      getUnknownsPerCellBoundary(),
      CompressionAccuracy,true
      );
      return false;
      },
    peano::datatraversal::TaskSet::TaskType::Background,
    peano::datatraversal::TaskSet::TaskType::Background,
    peano::datatraversal::TaskSet::TaskType::Background,
    peano::datatraversal::TaskSet::TaskType::Background,
    peano::datatraversal::TaskSet::TaskType::Background,
    true
  );

  assertion(1<=compressionOfPreviousSolution);
  assertion(1<=compressionOfSolution);
  assertion(1<=compressionOfUpdate);
  assertion(1<=compressionOfExtrapolatedPredictor);
  assertion(1<=compressionOfFluctuation);

  assertion(compressionOfPreviousSolution<=7);
  assertion(compressionOfSolution<=7);
  assertion(compressionOfUpdate<=7);
  assertion(compressionOfExtrapolatedPredictor<=7);
  assertion(compressionOfFluctuation<=7);

  peano::datatraversal::TaskSet runParallelTasks(
    [&]() -> bool {
      cellDescription.setBytesPerDoFInPreviousSolution(compressionOfPreviousSolution);
      if (compressionOfPreviousSolution<7) {
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          cellDescription.setPreviousSolutionCompressedIndex( CompressedDataHeap::getInstance().createData(0,0) );
          assertion( cellDescription.getPreviousSolutionCompressedIndex()>=0 );
          cellDescription.setPreviousSolutionCompressed( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getPreviousSolutionCompressedIndex()).data() ) );
        lock.free();

        const int numberOfEntries = getDataPerCell();
        tearApart(numberOfEntries, cellDescription.getPreviousSolutionIndex(), cellDescription.getPreviousSolutionCompressedIndex(), compressionOfPreviousSolution);

        #if defined(TrackGridStatistics)
        lock.lock();
          PipedUncompressedBytes += getDataHeapEntries(cellDescription.getPreviousSolutionIndex()).size() * 8.0;
          PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getPreviousSolutionCompressedIndex() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
          DataHeap::getInstance().deleteData( cellDescription.getPreviousSolutionIndex(), true );
          cellDescription.setPreviousSolutionIndex(-1);
          cellDescription.setPreviousSolution(nullptr);
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          PipedUncompressedBytes += getDataHeapEntries(cellDescription.getPreviousSolutionIndex()).size() * 8.0;
          PipedCompressedBytes   += getDataHeapEntries(cellDescription.getPreviousSolutionIndex()).size() * 8.0;
        lock.free();
        #endif
      }
      return false;
    },
    [&]() -> bool {
      cellDescription.setBytesPerDoFInSolution(compressionOfSolution);
      if (compressionOfSolution<7) {
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          cellDescription.setSolutionCompressedIndex(CompressedDataHeap::getInstance().createData(0,0));
          assertion1( cellDescription.getSolutionCompressedIndex()>=0, cellDescription.getSolutionCompressedIndex() );
          cellDescription.setSolutionCompressed( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getSolutionCompressedIndex()).data() ) );
        lock.free();

        const int numberOfEntries = getDataPerCell();

        tearApart(numberOfEntries, cellDescription.getSolutionIndex(), cellDescription.getSolutionCompressedIndex(), compressionOfSolution);

        #if defined(TrackGridStatistics)
        lock.lock();
          PipedUncompressedBytes += getDataHeapEntries(cellDescription.getSolutionIndex()).size() * 8.0;
          PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getSolutionCompressedIndex() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
          DataHeap::getInstance().deleteData( cellDescription.getSolutionIndex(), true );
          cellDescription.setSolutionIndex(-1);
          cellDescription.setSolution(nullptr);
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          PipedUncompressedBytes += getDataHeapEntries(cellDescription.getSolutionIndex()).size() * 8.0;
          PipedCompressedBytes   += getDataHeapEntries(cellDescription.getSolutionIndex()).size() * 8.0;
        lock.free();
        #endif
      }
      return false;
    },
    [&]() -> bool {
      cellDescription.setBytesPerDoFInUpdate(compressionOfUpdate);
      if (compressionOfUpdate<7) {
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          cellDescription.setUpdateCompressedIndex( CompressedDataHeap::getInstance().createData(0,0) );
          assertion( cellDescription.getUpdateCompressedIndex()>=0 );
          cellDescription.setUpdateCompressed( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getUpdateCompressedIndex()).data() ) );
        lock.free();

        const int numberOfEntries = getUnknownsPerCell();
        tearApart(numberOfEntries, cellDescription.getUpdateIndex(), cellDescription.getUpdateCompressedIndex(), compressionOfUpdate);

        #if defined(TrackGridStatistics)
        lock.lock();
          PipedUncompressedBytes += getDataHeapEntries(cellDescription.getUpdateIndex()).size() * 8.0;
          PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getUpdateCompressedIndex() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
          DataHeap::getInstance().deleteData( cellDescription.getUpdateIndex(), true );
          cellDescription.setUpdateIndex(-1);
          cellDescription.setUpdate(nullptr);
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          PipedUncompressedBytes += getDataHeapEntries(cellDescription.getUpdateIndex()).size() * 8.0;
          PipedCompressedBytes   += getDataHeapEntries(cellDescription.getUpdateIndex()).size() * 8.0;
        lock.free();
        #endif
      }
      return false;
    },
    [&]() -> bool {
      cellDescription.setBytesPerDoFInExtrapolatedPredictor(compressionOfExtrapolatedPredictor);
      if (compressionOfExtrapolatedPredictor<7) {
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          cellDescription.setExtrapolatedPredictorCompressedIndex( CompressedDataHeap::getInstance().createData(0,0) );
          assertion( cellDescription.getExtrapolatedPredictorCompressedIndex()>=0 );
          cellDescription.setExtrapolatedPredictorCompressed( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictorCompressedIndex()).data() ) );
        lock.free();

        const int numberOfEntries = getDataPerCellBoundary();
        tearApart(numberOfEntries, cellDescription.getExtrapolatedPredictorIndex(), cellDescription.getExtrapolatedPredictorCompressedIndex(), compressionOfExtrapolatedPredictor);

        #if defined(TrackGridStatistics)
        lock.lock();
          PipedUncompressedBytes += getDataHeapEntries(cellDescription.getExtrapolatedPredictorIndex()).size() * 8.0;
          PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getExtrapolatedPredictorCompressedIndex() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
          DataHeap::getInstance().deleteData( cellDescription.getExtrapolatedPredictorIndex(), true );
          cellDescription.setExtrapolatedPredictorIndex(-1);
          cellDescription.setExtrapolatedPredictor(nullptr);
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          PipedUncompressedBytes += getDataHeapEntries(cellDescription.getExtrapolatedPredictorIndex()).size() * 8.0;
          PipedCompressedBytes   += getDataHeapEntries(cellDescription.getExtrapolatedPredictorIndex()).size() * 8.0;
        lock.free();
        #endif
      }
      return false;
    },
    [&]() -> bool {
      cellDescription.setBytesPerDoFInFluctuation(compressionOfFluctuation);
      if (compressionOfFluctuation<7) {
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          cellDescription.setFluctuationCompressedIndex( CompressedDataHeap::getInstance().createData(0,0) );
          assertion( cellDescription.getFluctuationCompressedIndex()>=0 );
          cellDescription.setFluctuationCompressed( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getFluctuationCompressedIndex()).data() ) );
        lock.free();

        const int numberOfEntries = getUnknownsPerCellBoundary();
        tearApart(numberOfEntries, cellDescription.getFluctuationIndex(), cellDescription.getFluctuationCompressedIndex(), compressionOfFluctuation);

        #if defined(TrackGridStatistics)
        lock.lock();
          PipedUncompressedBytes += getDataHeapEntries(cellDescription.getFluctuationIndex()).size() * 8.0;
          PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getFluctuationCompressedIndex() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
          DataHeap::getInstance().deleteData( cellDescription.getFluctuationIndex(), true );
          cellDescription.setFluctuationIndex(-1);
          cellDescription.setFluctuation(nullptr);
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          PipedUncompressedBytes += getDataHeapEntries(cellDescription.getFluctuationIndex()).size() * 8.0;
          PipedCompressedBytes   += getDataHeapEntries(cellDescription.getFluctuationIndex()).size() * 8.0;
        lock.free();
        #endif
      }
      return false;
    },
    peano::datatraversal::TaskSet::TaskType::Background,
    peano::datatraversal::TaskSet::TaskType::Background,
    peano::datatraversal::TaskSet::TaskType::Background,
    peano::datatraversal::TaskSet::TaskType::Background,
    peano::datatraversal::TaskSet::TaskType::Background,
    true
  );
}

void exahype::solvers::ADERDGSolver::pullUnknownsFromByteStream(
    CellDescription& cellDescription) const {
  assertion(CompressionAccuracy>0.0);

  #if !defined(ValidateCompressedVsUncompressedData)
  const int dataPointsPerCell       = getDataPerCell();
  const int unknownsPerCellBoundary = getUnknownsPerCellBoundary();

  {
    tarch::multicore::Lock lock(exahype::HeapSemaphore);
      cellDescription.setPreviousSolutionIndex( DataHeap::getInstance().createData( dataPointsPerCell, dataPointsPerCell ) );
      cellDescription.setSolutionIndex( DataHeap::getInstance().createData(         dataPointsPerCell, dataPointsPerCell ) );
      cellDescription.setUpdateIndex( DataHeap::getInstance().createData(           getUpdateSize(),   getUpdateSize() ) );

      cellDescription.setExtrapolatedPredictorIndex( DataHeap::getInstance().createData( unknownsPerCellBoundary, unknownsPerCellBoundary ) );
      cellDescription.setFluctuationIndex( DataHeap::getInstance().createData(           unknownsPerCellBoundary, unknownsPerCellBoundary ) );
    lock.free();

    if (cellDescription.getPreviousSolutionIndex()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
        cellDescription.setPreviousSolutionIndex( DataHeap::getInstance().createData( dataPointsPerCell, dataPointsPerCell ) );
      lock.free();
    }
    if (cellDescription.getSolutionIndex()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
        cellDescription.setSolutionIndex( DataHeap::getInstance().createData( dataPointsPerCell, dataPointsPerCell ) );
      lock.free();
    }
    if (cellDescription.getUpdateIndex()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
        cellDescription.setUpdateIndex( DataHeap::getInstance().createData( getUpdateSize(), getUpdateSize() ) );
      lock.free();
    }
    if (cellDescription.getExtrapolatedPredictorIndex()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
        cellDescription.setExtrapolatedPredictorIndex( DataHeap::getInstance().createData(unknownsPerCellBoundary ) );
      lock.free();
    }
    if (cellDescription.getFluctuationIndex()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
        cellDescription.setFluctuationIndex( DataHeap::getInstance().createData( unknownsPerCellBoundary, unknownsPerCellBoundary ) );
      lock.free();
    }

    cellDescription.setPreviousSolution     ( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getPreviousSolutionIndex()).data() ) );
    cellDescription.setSolution             ( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getSolutionIndex()).data() ) );
    cellDescription.setExtrapolatedPredictor( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictorIndex()).data() ) );

  }
  #else
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolution() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getSolution() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getUpdate() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedPredictor() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getFluctuation() ));
  #endif

  assertion1(
      CompressedDataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionCompressedIndex() ),
      cellDescription.getPreviousSolutionCompressedIndex()
    );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getSolutionCompressedIndex() ),
    cellDescription.getSolutionCompressedIndex()
  );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getUpdateCompressedIndex() ),
    cellDescription.getUpdateCompressedIndex()
  );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedPredictorCompressedIndex() ),
    cellDescription.getExtrapolatedPredictorCompressedIndex()
  );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getFluctuationCompressedIndex() ),
    cellDescription.getFluctuationCompressedIndex()
  );

  peano::datatraversal::TaskSet glueTasks(
    [&]() -> bool {
      if (cellDescription.getBytesPerDoFInPreviousSolution()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionIndex() ), cellDescription.getPreviousSolutionIndex());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionCompressedIndex() ));
        const int numberOfEntries = getDataPerCell();
        glueTogether(numberOfEntries, cellDescription.getPreviousSolutionIndex(), cellDescription.getPreviousSolutionCompressedIndex(), cellDescription.getBytesPerDoFInPreviousSolution());
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          CompressedDataHeap::getInstance().deleteData( cellDescription.getPreviousSolutionCompressedIndex(), true );
          cellDescription.setPreviousSolutionCompressedIndex(-1);
          cellDescription.setPreviousSolutionCompressed(nullptr);
        lock.free();
      }
      return false;
    },
    [&]() -> bool {
      if (cellDescription.getBytesPerDoFInSolution()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getSolutionIndex() ), cellDescription.getSolutionIndex() );
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getSolutionCompressedIndex() ));
        const int numberOfEntries = getDataPerCell();
        glueTogether(numberOfEntries, cellDescription.getSolutionIndex(), cellDescription.getSolutionCompressedIndex(), cellDescription.getBytesPerDoFInSolution());
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          CompressedDataHeap::getInstance().deleteData( cellDescription.getSolutionCompressedIndex(), true );
          cellDescription.setSolutionCompressedIndex(-1);
          cellDescription.setSolutionCompressed(nullptr);
        lock.free();
      }
      return false;
    },
    [&]() -> bool {
      if (cellDescription.getBytesPerDoFInUpdate()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getUpdateIndex() ), cellDescription.getUpdateIndex());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getUpdateCompressedIndex() ));
        const int numberOfEntries = getUnknownsPerCell();
        glueTogether(numberOfEntries, cellDescription.getUpdateIndex(), cellDescription.getUpdateCompressedIndex(), cellDescription.getBytesPerDoFInUpdate());
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          CompressedDataHeap::getInstance().deleteData( cellDescription.getUpdateCompressedIndex(), true );
          cellDescription.setUpdateCompressedIndex(-1);
          cellDescription.setUpdateCompressed(nullptr);
        lock.free();
      }
      return false;
    },
    [&]() -> bool {
      if (cellDescription.getBytesPerDoFInExtrapolatedPredictor()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedPredictorIndex() ), cellDescription.getExtrapolatedPredictorIndex());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedPredictorCompressedIndex() ));
        const int numberOfEntries = getDataPerCellBoundary();
        glueTogether(numberOfEntries, cellDescription.getExtrapolatedPredictorIndex(), cellDescription.getExtrapolatedPredictorCompressedIndex(), cellDescription.getBytesPerDoFInExtrapolatedPredictor());
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          CompressedDataHeap::getInstance().deleteData( cellDescription.getExtrapolatedPredictorCompressedIndex(), true );
          cellDescription.setExtrapolatedPredictorCompressedIndex(-1);
          cellDescription.setExtrapolatedPredictorCompressed(nullptr);
        lock.free();
      }
      return false;
    },
    [&]() -> bool {
      if (cellDescription.getBytesPerDoFInFluctuation()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getFluctuationIndex() ), cellDescription.getFluctuationIndex());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getFluctuationCompressedIndex() ));
        const int numberOfEntries = getUnknownsPerCellBoundary();
        glueTogether(numberOfEntries, cellDescription.getFluctuationIndex(), cellDescription.getFluctuationCompressedIndex(), cellDescription.getBytesPerDoFInFluctuation());
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          CompressedDataHeap::getInstance().deleteData( cellDescription.getFluctuationCompressedIndex(), true );
          cellDescription.setFluctuationCompressedIndex(-1);
          cellDescription.setFluctuationCompressed(nullptr);
        lock.free();
      }
      return false;
    },
    peano::datatraversal::TaskSet::TaskType::Background,
    peano::datatraversal::TaskSet::TaskType::Background,
    peano::datatraversal::TaskSet::TaskType::Background,
    peano::datatraversal::TaskSet::TaskType::Background,
    peano::datatraversal::TaskSet::TaskType::Background,
    true
  );
}

///////////////////////
// PROFILING
///////////////////////

exahype::solvers::Solver::CellProcessingTimes exahype::solvers::ADERDGSolver::measureCellProcessingTimes(const int numberOfRuns) {
  // Setup
  const int cellDescriptionsIndex = ADERDGSolver::Heap::getInstance().createData(0,1);
  FiniteVolumesSolver::Heap::getInstance().createDataForIndex(cellDescriptionsIndex,0,1); // needs to be done

  Solver::CellInfo cellInfo(cellDescriptionsIndex);
  addNewCellDescription(
      0,cellInfo,CellDescription::Type::Leaf,
      getMaximumAdaptiveMeshLevel(), /* needs to be on the fine grid for the limiter cells */-1,
      getCoarsestMeshSize(),
      _domainOffset);

  CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[0];

  ensureNecessaryMemoryIsAllocated(cellDescription);

  adjustSolutionDuringMeshRefinementBody(cellDescription,true);
  updateTimeStepSize(0,cellInfo);
  const double dt = cellDescription.getTimeStepSize();

  // ADER-DG specific setup ( all Riemanns have been performed, cell is surrounded by other Cell type cells )
  cellDescription.setAugmentationStatus(0);
  cellDescription.setFacewiseAugmentationStatus(0);
  cellDescription.setCommunicationStatus(ADERDGSolver::LeafCommunicationStatus);
  cellDescription.setFacewiseCommunicationStatus(ADERDGSolver::LeafCommunicationStatus);

  // MEASUREMENTS
  CellProcessingTimes result;

  // measure ADERDG STP
  {
    const std::chrono::high_resolution_clock::time_point timeStart = std::chrono::high_resolution_clock::now();
    int numberOfPicardIterations = std::numeric_limits<int>::max();
    for (int it=0; it<numberOfRuns; it++) {
      numberOfPicardIterations = predictionAndVolumeIntegralBody(cellDescription,cellDescription.getTimeStamp(),cellDescription.getTimeStepSize(),false,true,true);
    }
    const double time_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-timeStart).count() * 1e-9;
    result._minTimePredictor = time_sec / numberOfRuns / std::abs(numberOfPicardIterations);
    if ( isLinear() ) {
      result._maxTimePredictor = result._minTimePredictor;
    } else {
      result._maxTimePredictor = result._minTimePredictor * getNodesPerCoordinateAxis(); // * (order+1)
    }
  }

  // measure ADER-DG cell update
  cellDescription.setRefinementStatus(_refineOrKeepOnFineGrid);
  cellDescription.setFacewiseRefinementStatus(_refineOrKeepOnFineGrid);
  {
    const std::chrono::high_resolution_clock::time_point timeStart = std::chrono::high_resolution_clock::now();
    for (int it=0; it<numberOfRuns; it++) {
      cellDescription.setNeighbourMergePerformed(static_cast<unsigned char>(true));
      const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int> boundaryMarkers(0); // >= 0 indicates no remote/domain boundary
      updateBody(cellDescription,cellInfo,boundaryMarkers);

      swapSolutionAndPreviousSolution(cellDescription); // assumed to be very cheap
      rollbackToPreviousTimeStep(cellDescription);
      cellDescription.setTimeStepSize(dt);
    }
    const double time_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-timeStart).count() * 1e-9;
    result._timeADERDGUpdate = time_sec / numberOfRuns;
  }

  // measure Riemann solve
  {
    const tarch::la::Vector<DIMENSIONS,int> pos1(0);
    tarch::la::Vector<DIMENSIONS,int> pos2(0); pos2[0]=1;
    Solver::InterfaceInfo face(pos1,pos2);
    const std::chrono::high_resolution_clock::time_point timeStart = std::chrono::high_resolution_clock::now();
    for (int it=0; it<numberOfRuns; it++) {
      solveRiemannProblemAtInterface(cellDescription,cellDescription,face);
    }
    const double time_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-timeStart).count() * 1e-9;
    result._timeADERDGRiemann = time_sec / numberOfRuns;
  }

  // Clean up
  cellDescription.setType(CellDescription::Type::Erased);
  ensureNoUnnecessaryMemoryIsAllocated(cellDescription);

  DataHeap::getInstance().deleteAllData();
  ADERDGSolver::Heap::getInstance().deleteAllData();
  FiniteVolumesSolver::Heap::getInstance().deleteAllData();

  return result;
}

//#undef assertion
//#define assertion(expr) 
