#if defined(Parallel) && defined(SharedTBB)
#include "ADERDGSolver.h"

#if defined(ScoreP)
#include "scorep/SCOREP_User.h"
#endif

#if defined(FileTrace)
#include "exahype/reactive/STPStatsTracer.h"
#endif

#include "exahype/reactive/PerformanceMonitor.h"
#include "exahype/reactive/OffloadingAnalyser.h"
#include "exahype/reactive/OffloadingProfiler.h"
#include "exahype/reactive/ResilienceStatistics.h"
#include "exahype/reactive/MemoryMonitor.h"
#include "exahype/reactive/NoiseGenerator.h"
#include "exahype/reactive/ResilienceTools.h"
#include "exahype/reactive/TimeStampAndLimiterTeamHistory.h"
#include "exahype/reactive/RequestManager.h"

#include "tarch/multicore/Core.h"

#if defined(USE_TMPI)
#include "teaMPI.h"
#endif

//#undef assertion
//#define assertion assert

MPI_Datatype exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::_datatype;

exahype::solvers::ADERDGSolver::MigratablePredictionJob::MigratablePredictionJob(
    ADERDGSolver& solver,
    const int cellDescriptionsIndex,
    const int element,
    const double predictorTimeStamp,
    const double predictorTimeStepSize,
    const bool isPotSoftErrorTriggered,
    const bool isSkeleton) :
      tarch::multicore::jobs::Job(
          tarch::multicore::jobs::JobType::BackgroundTask, 0,
          getTaskPriorityLocalMigratableJob(cellDescriptionsIndex, element,
              predictorTimeStamp, isSkeleton)),  //this is a locally executed job
      _solver(solver),
      _cellDescriptionsIndex(cellDescriptionsIndex),
      _element(element),
      _predictorTimeStamp(predictorTimeStamp),
      _predictorTimeStepSize(predictorTimeStepSize),
      _originRank(tarch::parallel::Node::getInstance().getRank()),
      _tag(-1),
      _isSkeleton(isSkeleton),
      _luh(nullptr),
      _lduh(nullptr),
      _lQhbnd(nullptr),
      _lFhbnd(nullptr),
      _lGradQhbnd(nullptr),
      _isLocalReplica(false),
      _isPotSoftErrorTriggered(isPotSoftErrorTriggered),
      _isCorrupted(false),
      _currentState(State::INITIAL)
{
  if (_isSkeleton) {
    NumberOfSkeletonJobs.fetch_add(1);
  } else {
    LocalStealableSTPCounter++;
    NumberOfEnclaveJobs.fetch_add(1);
  }
  exahype::reactive::ResilienceStatistics::getInstance().notifySpawnedTask();
  exahype::reactive::PerformanceMonitor::getInstance().incCurrentTasks();

  auto& cellDescription = getCellDescription(cellDescriptionsIndex, element);
  logDebug("MigratablePredictionJob","team "<<exahype::reactive::OffloadingContext::getInstance().getTMPITeamNumber()<<" spawning STP for "
        <<cellDescription.toString());

  tarch::la::Vector<DIMENSIONS, double> center = cellDescription.getOffset()+0.5*cellDescription.getSize();
  tarch::la::Vector<DIMENSIONS, double> dx = cellDescription.getSize();

  for (int i = 0; i < DIMENSIONS; i++) {
    _center[i] = center[i];
    _dx[i] = dx[i];
  }

}

exahype::solvers::ADERDGSolver::MigratablePredictionJob::MigratablePredictionJob(
    ADERDGSolver& solver,
    const int cellDescriptionsIndex,
    const int element,
    const double predictorTimeStamp,
    const double predictorTimeStepSize,
    double *luh, double *lduh,
    double *lQhbnd, double *lFhbnd,
    double *lGradQhbnd, double *dx,
    double *center, const int originRank,
    const int tag) :
      tarch::multicore::jobs::Job(
        tarch::multicore::jobs::JobType::BackgroundTask, 0,
        tarch::multicore::DefaultPriority * 4), //this is a remotely executed job -> high priority
      _solver(solver),
      _cellDescriptionsIndex(cellDescriptionsIndex),
      _element(element),
      _predictorTimeStamp(predictorTimeStamp),
      _predictorTimeStepSize(predictorTimeStepSize),
      _originRank(originRank),
      _tag(tag),
      _isSkeleton(false),  //this constructor should never be invoked for skeletons!
      _luh(luh),
      _lduh(lduh),
      _lQhbnd(lQhbnd),
      _lFhbnd(lFhbnd),
      _lGradQhbnd(lGradQhbnd),
      _isLocalReplica(false),
      _isPotSoftErrorTriggered(false),
      _isCorrupted(false),
      _currentState(State::INITIAL)
{
  assertion(!_isSkeleton);
  assertion(element==0); //todo: Is element still used somewhere? Dominic's code seems to assume it to be zero...

  for (int i = 0; i < DIMENSIONS; i++) {
    _center[i] = center[i];
    _dx[i] = dx[i];
  }

  if (isRemoteJob()) {
    NumberOfStolenJobs++;
  }
  else
    NumberOfEnclaveJobs++;
  exahype::reactive::PerformanceMonitor::getInstance().incCurrentTasks();
}

exahype::solvers::ADERDGSolver::MigratablePredictionJob::~MigratablePredictionJob() {
}

//Caution: Compression is not supported yet!
bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::run(
    bool isCalledOnMaster) {

  bool reschedule = false;

  switch(_currentState)
  {
    case State::INITIAL:
      reschedule = runExecution(isCalledOnMaster);
      break;
    case State::CHECK_REQUIRED:
      reschedule = runCheck();
      break;
    /*case State::CHECK_PREVIOUS:
      reschedule = tryFindPreviousOutcomeAndCheck();
      if(!reschedule)
        NumberOfEnclaveJobs--;
      break;*/
    /*case State::HEALING_REQUIRED:
      reschedule = tryFindOutcomeAndHeal();
      if(!reschedule)
        NumberOfEnclaveJobs--;
      break;*/
    default:
      logError("run","inconsistent state");
  }
  return reschedule;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::runCheck() {
  //caution, this may delay setting the completed flag and cause slowdowns for the master!
#if !defined(OffloadingUseProgressThread)
  exahype::solvers::ADERDGSolver::progressOffloading(&_solver, false, MAX_PROGRESS_ITS);
#endif

  DeliveryStatus status;
  MigratablePredictionJobData *outcome;
  bool found = tryToFindAndExtractEquivalentSharedOutcome(false, status, &outcome);
  bool reschedule;

  if(found && status==DeliveryStatus::Received) {
    if(matches(outcome)) {
      CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,_element);
      cellDescription.setHasCompletedLastStep(true);
      reschedule = false;
     }
     else {
       //soft error detected
       reschedule = false;
       MPI_Abort(MPI_COMM_WORLD, -1);
     }
  }
  else {
    reschedule = true;
  }

  return reschedule;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::runExecution(bool isCalledOnMaster) {
#ifdef USE_ITAC
  VT_begin(event_stp);
#endif

#if defined ScoreP
  SCOREP_USER_REGION( "exahype::solvers::ADERDGSolver::MigratablePredictionJob::run", SCOREP_USER_REGION_TYPE_FUNCTION )
#endif

  bool result = false;
  bool hasComputed = false;
  int curr = std::atomic_fetch_add(&JobCounter, 1);

  if (curr % 1000 == 0) {
    tarch::timing::Watch watch("exahype::MigratablePredictionJob::", "-", false, false);
    watch.startTimer();
    result = handleExecution(isCalledOnMaster, hasComputed);
    watch.stopTimer();

    logDebug("run()","measured time per STP "<<watch.getCalendarTime());
    if (hasComputed) {
      exahype::reactive::OffloadingAnalyser::getInstance().setTimePerSTP(
          watch.getCalendarTime());
    }
  }
  else
    result = handleExecution(isCalledOnMaster, hasComputed);

#if defined(GenerateNoise)
    exahype::reactive::NoiseGenerator::getInstance().generateNoiseSTP();
#endif

#ifdef USE_ITAC
  VT_end(event_stp);
#endif
  return result;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::handleExecution(
    bool isRunOnMaster, bool& hasComputed) {

  bool result = false;
#if defined(OffloadingUseProfiler)
  exahype::reactive::OffloadingProfiler::getInstance().beginComputation();
  double time = -MPI_Wtime();
#endif
  //local treatment if this job belongs to the local rank
  if (!isRemoteJob()) {
    result = handleLocalExecution(isRunOnMaster, hasComputed);
    if (!_isLocalReplica) {
      if (_isSkeleton) {
        NumberOfSkeletonJobs.fetch_sub(1);
      } else {
        NumberOfEnclaveJobs.fetch_sub(1);
      }
    }
    assertion( NumberOfEnclaveJobs>=0 );
    assertion( NumberOfSkeletonJobs>=0 );
#if !defined(OffloadingUseProgressThread)
    //double time = -MPI_Wtime();
    if (!isRunOnMaster)
      exahype::solvers::ADERDGSolver::progressOffloading(&_solver,
          isRunOnMaster, std::numeric_limits<int>::max());
      //exahype::solvers::ADERDGSolver::progressOffloading(&_solver,
      //    isRunOnMaster, MAX_PROGRESS_ITS);
   //time += MPI_Wtime();
   //if(time>0.02) {
   //  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
   //                                          _element);
   //  logError("handleExecution","progress took too long "<<time<<" cellDesc "<<cellDescription.toString());

   //}
#endif
  }
  //remote task, need to execute and send it back
  else {
    result = handleRemoteExecution(hasComputed);
    sendBackOutcomeToOrigin();
  }
#if defined(OffloadingUseProfiler)
  time += MPI_Wtime();
  exahype::reactive::OffloadingProfiler::getInstance().endComputation(time);
#endif
  return result;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::handleLocalExecution(bool isRunOnMaster, bool& hasComputed) {
  //fast path
  if(exahype::reactive::OffloadingContext::getInstance().getResilienceStrategy()
      ==exahype::reactive::OffloadingContext::ResilienceStrategy::None) {
    executeLocally();
    hasComputed = true;
    corruptIfActive();
    setFinished();
    return false;
  }

  //got some task sharing here
  bool isEqual = false;
  bool hasFlipped = false;
  //check for outcome and get it + corruption status
  MigratablePredictionJobData *outcome = nullptr;
  DeliveryStatus status;
  if(tryToFindAndExtractEquivalentSharedOutcome(false, status, &outcome)) {
    executeOrCopySTPOutcome(outcome, hasComputed, hasFlipped); //corruption can also happen here if executed
    setSTPPotCorrupted(hasFlipped);
    if(needToCheckThisSTP(hasComputed)) {
      isEqual = matches(outcome);
      if(!isEqual) {
        logError("handleLocalExecution", "Soft error detected in job execution: "<<to_string());
        MPI_Abort(MPI_COMM_WORLD, -1);
      }
    }
    if(needToShare(true /*hasOutcome*/, outcome->_metadata._isPotSoftErrorTriggered)) {
      shareSTPImmediatelyOrLater();
    }
    if(needToPutBackOutcome()) {
      MigratablePredictionJobOutcomeKey key(outcome->_metadata.getCenter(), outcome->_metadata.getPredictorTimeStamp(),
                                             outcome->_metadata.getPredictorTimeStepSize(), outcome->_metadata.getElement());
      _solver._outcomeDatabase.insertOutcome(key, outcome, DeliveryStatus::Received);
    }
    setFinished(); //have already checked
    return false;
  }
  //don't have an outcome
  else {
    executeLocally();
    hasComputed = true;
    hasFlipped = corruptIfActive();
    setSTPPotCorrupted(hasFlipped);

    if(needToShare(false, false)) {
      shareSTPImmediatelyOrLater();
    }

    if(needToCheckThisSTP(hasComputed)) {
      logInfo("handleLocalExecution","going into check mode "<<to_string());
      setState(State::CHECK_REQUIRED);
      return true; //re-enqueue
    }

    setFinished();
    return false;
  }
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::shareSTPImmediatelyOrLater() {
  if(exahype::reactive::ResilienceTools::CheckLimitedCellsOnly
     && (exahype::reactive::OffloadingContext::getInstance().getResilienceStrategy()
        >=exahype::reactive::OffloadingContext::ResilienceStrategy::TaskSharingResilienceChecks)) {
    logInfo("handleLocalExecution", "Delaying outcome, as we need to compute corrector first!");
    _solver.storePendingOutcomeToBeShared(this); //delay sharing until we can be sure that trigger has been set (correction must happen first)
  }
  else {
    if(_isCorrupted && !_isPotSoftErrorTriggered)
      logWarning("handleLocalExecution","Sharing corrupted outcome but soft error has not been triggered!");
    _solver.sendTaskOutcomeToOtherTeams(this);
  }
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::needToPutBackOutcome() {
  return (exahype::reactive::OffloadingContext::getInstance().getResilienceStrategy()
        >=exahype::reactive::OffloadingContext::ResilienceStrategy::TaskSharingResilienceChecks)
      &&  exahype::reactive::ResilienceTools::CheckLimitedCellsOnly;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::needToCheckThisSTP(bool hasComputed) {
  return (exahype::reactive::OffloadingContext::getInstance().getResilienceStrategy()
        >=exahype::reactive::OffloadingContext::ResilienceStrategy::TaskSharingResilienceChecks)
      &&  hasComputed
      && _isPotSoftErrorTriggered;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::needToShare(bool hasOutcome, bool isOutcomePotCorrupt) {
  if(exahype::reactive::OffloadingContext::getInstance().getResilienceStrategy()
      ==exahype::reactive::OffloadingContext::ResilienceStrategy::None) {
    return false;
  }
  else if(exahype::reactive::OffloadingContext::getInstance().getResilienceStrategy()
          ==exahype::reactive::OffloadingContext::ResilienceStrategy::TaskSharing){
    return !hasOutcome
        && (AllocatedSTPsSend
            <= exahype::reactive::PerformanceMonitor::getInstance().getTasksPerTimestep()/((exahype::reactive::OffloadingContext::getInstance().getTMPINumTeams()-1)*tarch::multicore::Core::getInstance().getNumberOfThreads()));
  }
  else {
    return !hasOutcome
         || exahype::reactive::ResilienceTools::CheckLimitedCellsOnly
         || isOutcomePotCorrupt;
  }
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::corruptIfActive() {
  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
   _element);

  double *lduh = static_cast<double*>(cellDescription.getUpdate());
  bool hasFlipped = exahype::reactive::ResilienceTools::getInstance().corruptDataIfActive(lduh, _solver.getUpdateSize());

  if(hasFlipped) {
    logInfo("corruptIfActive","Has corrupted STP job "<<to_string());
  }

  return hasFlipped;
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::copyOutcome(MigratablePredictionJobData *outcome) {
  assert(outcome!=nullptr);
  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
   _element);

  double *lduh = static_cast<double*>(cellDescription.getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
#if OffloadingGradQhbnd
  double *lGradQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient());
#endif
  double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

  if(outcome->_metadata._isCorrupted)
    logError("copyOutcome","Error: we're using a corrupted outcome!"); assertion(false);

  assertion(outcome!=nullptr);
  std::memcpy(lduh, &outcome->_lduh[0], outcome->_lduh.size() * sizeof(double));
  std::memcpy(lQhbnd, &outcome->_lQhbnd[0], outcome->_lQhbnd.size() * sizeof(double));
  std::memcpy(lFhbnd, &outcome->_lFhbnd[0], outcome->_lFhbnd.size() * sizeof(double));
#if OffloadingGradQhbnd
  std::memcpy(lGradQhbnd, &outcome->_lGradQhbnd[0], outcome->_lGradQhbnd.size() * sizeof(double));
#endif
  logInfo("copyOutcome","team "<<exahype::reactive::OffloadingContext::getInstance().getTMPITeamNumber()<<" copied STP outcome for "
        <<to_string());

  exahype::reactive::ResilienceStatistics::getInstance().notifySavedTask();
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::executeOrCopySTPOutcome(MigratablePredictionJobData *outcome, bool& hasComputed, bool& hasFlipped) {
  if(exahype::reactive::OffloadingContext::getInstance().getResilienceStrategy()==exahype::reactive::OffloadingContext::ResilienceStrategy::TaskSharing) {
    copyOutcome(outcome);
    hasComputed = false;
  }
  else {
    if(exahype::reactive::OffloadingContext::getInstance().getResilienceStrategy()>=exahype::reactive::OffloadingContext::ResilienceStrategy::TaskSharingResilienceChecks) {
      //if(outcome->_metadata._isPotSoftErrorTriggered) {
        executeLocally();
        hasFlipped = corruptIfActive();
        hasComputed = true;
      //}
      //else {
      //  copyOutcome(outcome);
      //  hasComputed = false;
       // }
    }
  }
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::executeLocally() {
  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
   _element);

  double *luh = static_cast<double*>(cellDescription.getSolution());
  double *lduh = static_cast<double*>(cellDescription.getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double *lGradQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient());
  double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

  logInfo("handleLocalExecution","team "<<exahype::reactive::OffloadingContext::getInstance().getTMPITeamNumber()<<" computing STP for "
        <<to_string()
        <<std::setprecision(30)<<" time step size "<<_predictorTimeStepSize);

#if defined FileTrace
    auto start = std::chrono::high_resolution_clock::now();
#endif

  int iterations = _solver.fusedSpaceTimePredictorVolumeIntegral(lduh,
          lQhbnd,
          lGradQhbnd,
          lFhbnd,
          luh,
          cellDescription.getOffset() + 0.5 * cellDescription.getSize(),
          cellDescription.getSize(),
          _predictorTimeStamp,
          _predictorTimeStepSize,
          true);

#if defined(FileTrace)
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    exahype::reactive::STPStatsTracer::getInstance().writeTracingEventRunIterations(duration.count(), iterations, exahype::reactive::STPTraceKey::ADERDGOwnMigratable);
#endif

  exahype::reactive::ResilienceStatistics::getInstance().notifyExecutedTask();

}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::setFinished() {
  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
   _element);

  cellDescription.setHasCompletedLastStep(true);
  exahype::reactive::PerformanceMonitor::getInstance().decRemainingTasks();
  exahype::reactive::PerformanceMonitor::getInstance().decCurrentTasks();
}

/*bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::handleLocalExecution(
    bool isRunOnMaster, bool& hasComputed) {

  bool hasFlipped = false;
  bool reschedule = false;
  bool needToCompute = true;
//#if defined(TaskSharing)
  bool needToCheck = false;
  bool copyResult = false;
  bool hasResult = false;
  bool needToShare = false;
//#endif

  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
	   _element);

  double *luh = static_cast<double*>(cellDescription.getSolution());
  double *lduh = static_cast<double*>(cellDescription.getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double *lGradQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient());
  double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

  //JobOutcomeStatus status;
  //MigratablePredictionJobData *outcomes[2];
  //int required = exahype::reactive::OffloadingContext::getInstance().getTMPINumTeams()-1;
  DeliveryStatus status;
  MigratablePredictionJobData *outcome = nullptr;
  bool found = false;

//#if defined(TaskSharing)
  if(exahype::reactive::OffloadingContext::getInstance().getResilienceStrategy()
      != exahype::reactive::OffloadingContext::ResilienceStrategy::None) {
    needToShare = (AllocatedSTPsSend
                   <= exahype::reactive::PerformanceMonitor::getInstance().getTasksPerTimestep()/tarch::multicore::Core::getInstance().getNumberOfThreads())
             || (exahype::reactive::OffloadingContext::getInstance().getResilienceStrategy()>=exahype::reactive::OffloadingContext::ResilienceStrategy::TaskSharingResilienceChecks
                  && exahype::reactive::ResilienceTools::CheckLimitedCellsOnly);

    if(!exahype::reactive::ResilienceTools::CheckLimitedCellsOnly)
      found = tryToFindAndExtractEquivalentSharedOutcome(false, status, outcome);
    else
      found = false;

    if(found) {
      switch (status) {
        case DeliveryStatus::Received:
      	  hasResult = true;
#if defined(ResilienceChecks)
          needToCompute = outcomes[0]->_metadata._isPotSoftErrorTriggered;
          copyResult = !outcomes[0]->_metadata._isPotSoftErrorTriggered;
          //needToCheck = outcomes[0]->_metadata._isPotSoftErrorTriggered;
          needToCheck = false;
          needToShare = outcomes[0]->_metadata._isPotSoftErrorTriggered;
#else
      	  needToShare = false;
      	  needToCompute = false;
          copyResult = true;
#endif
          break;
      case DeliveryStatus::Transit:
#ifdef TaskSharingRescheduleIfInTransit
        reschedule = true;
        needToCompute = false;
#endif
        break;
      }
    }

    if(copyResult) {
      if(outcome->_metadata._isCorrupted)
        logError("handleLocalExecution","Error: we're using a corrupted outcome!");

      assertion(outcome!=nullptr);
      std::memcpy(lduh, &outcome->_lduh[0], outcome->_lduh.size() * sizeof(double));
      std::memcpy(lQhbnd, &outcome->_lQhbnd[0], outcome->_lQhbnd.size() * sizeof(double));
      std::memcpy(lFhbnd, &outcome->_lFhbnd[0], outcome->_lFhbnd.size() * sizeof(double));
#if OffloadingGradQhbnd
      std::memcpy(lGradQhbnd, &outcome->_lGradQhbnd[0], outcome->_lGradQhbnd.size() * sizeof(double));
#endif
      exahype::reactive::JobTableStatistics::getInstance().notifySavedTask();
      AllocatedSTPsReceive--;
      delete outcome;
    }
  }

  if(needToCompute) {
#if defined(USE_ITAC)
    if(_isLocalReplica)
      VT_begin(event_stp_local_replica);
#endif
#if defined FileTrace
    auto start = std::chrono::high_resolution_clock::now();
#endif

    tarch::la::Vector<DIMENSIONS, double> center;
    center = cellDescription.getOffset() + 0.5 * cellDescription.getSize();
    logInfo("handleLocalExecution","team "<<exahype::reactive::OffloadingContext::getInstance().getTMPIInterTeamRank()<<" computing STP for "
          <<to_string());

    int iterations = _solver.fusedSpaceTimePredictorVolumeIntegral(lduh,
       	lQhbnd,
        lGradQhbnd,
        lFhbnd,
        luh,
        cellDescription.getOffset() + 0.5 * cellDescription.getSize(),
        cellDescription.getSize(),
        _predictorTimeStamp,
        _predictorTimeStepSize,
        true);

    hasComputed = true;

    exahype::reactive::JobTableStatistics::getInstance().notifyExecutedTask();

    hasFlipped = exahype::reactive::ResilienceTools::getInstance().corruptDataIfActive(lduh, _solver.getUpdateSize());
    setTrigger(hasFlipped);

    if(hasFlipped)  {
      logError("handleLocalExecution","team  "<<exahype::reactive::OffloadingContext::getInstance().getTMPIInterTeamRank()<<" has corrupted STP for "
                                      <<to_string());

      _isCorrupted = true; //indicates that this outcome has corrupt data!
    }
  
    logDebug("handleLocalExecution", "celldesc ="<<_cellDescriptionsIndex<<" _isPotSoftErrorTriggered ="<<(int) _isPotSoftErrorTriggered);

#if defined(ResilienceChecks)
    needToCheck = needToCheck || _isPotSoftErrorTriggered;
    needToShare = needToShare || _isPotSoftErrorTriggered;
#endif

#if defined(FileTrace)
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    exahype::reactive::STPStatsTracer::getInstance().writeTracingEventRunIterations(duration.count(), iterations, exahype::reactive::STPTraceKey::ADERDGOwnMigratable);
#endif

#if defined(USE_ITAC)
      if(_isLocalReplica)
        VT_end(event_stp_local_replica);
#endif
  }

#if defined(OffloadingLocalRecompute)
  if (_isLocalReplica) {
    tbb::concurrent_hash_map<const CellDescription*, std::pair<int, int>>::accessor a_cellDescToTagRank;
    //logInfo("handleExecution()", "cleaning up cell desc "<<&cellDescription);
    bool found = _solver._mapCellDescToTagRank.find(a_cellDescToTagRank,
        &cellDescription);
    assert(found);
    _solver._mapCellDescToTagRank.erase(a_cellDescToTagRank);
    a_cellDescToTagRank.release();
  }
#endif

  if((exahype::reactive::OffloadingContext::getInstance().getResilienceStrategy()
    != exahype::reactive::OffloadingContext::ResilienceStrategy::None)) {

    if(!hasResult && !exahype::reactive::ResilienceTools::CheckLimitedCellsOnly) {
      hasResult = tryToFindAndExtractEquivalentSharedOutcomes(false, required, status, outcomes);
    }

    if(needToCheck) {
      if(hasResult) {
#if defined(ResilienceChecks)
        int saneIdx = -1;
        bool hasNoError = checkAgainstOutcomesAndFindSaneOne(outcomes, required, saneIdx);
        if(!hasNoError) {
#if defined(ResilienceHealing)
          //logInfo("handleLocalExecution", "could switch on healing mode now!");
          recoverWithOutcome(outcomes[saneIdx]);
          reschedule = false;
          //_solver.switchToHealingMode();
          //MPI_Abort(MPI_COMM_WORLD, -1);
#else
          //MPI_Abort(MPI_COMM_WORLD, -1); //for now, abort immediately
          exahype::reactive::ResilienceTools::getInstance().setCorruptionDetected(true);
          logError("handleLocalExecution", "soft error detected but we ignore it...");
          reschedule = false;
#endif
        }
#endif
      }
      else {
        _currentState = State::CHECK_REQUIRED; //check later
         reschedule = true;
      }
    }
    else {
     if(hasResult)
      needToShare = false;
    }

    if(needToShare) {
      if(exahype::reactive::ResilienceTools::CheckLimitedCellsOnly) {
        logInfo("handleLocalExecution", "Delaying outcome, as we need to compute corrector first!");
        _solver.storePendingOutcomeToBeShared(this); //delay sharing until we can be sure that trigger has been set (correction must happen first)
      }
      else {
        if(_isCorrupted)
          logError("handleLocalExecution","Sharing corrupted outcome");
        _solver.sendTaskOutcomeToOtherTeams(this);
      }
    }
  }

  if(hasFlipped) {
    tarch::la::Vector<DIMENSIONS, double> center;
    center = cellDescription.getOffset() + 0.5 * cellDescription.getSize();
    logInfo("handleLocalExecution","Corrupted data in "
          << to_string()
          << " hasResult "<<hasResult
          << " needToCheck "<<needToCheck
          << " needToShare "<<needToShare
          << " reschedule  "<<reschedule);
  }

  if (!reschedule) {
    cellDescription.setHasCompletedLastStep(true);
    exahype::reactive::PerformanceMonitor::getInstance().decRemainingTasks();
    exahype::reactive::PerformanceMonitor::getInstance().decCurrentTasks();
  }

  return reschedule;
}*/

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::handleRemoteExecution(bool& hasComputed) {
#if defined(USE_ITAC)
  VT_begin(event_stp_remote);
#endif
#if defined FileTrace
  auto start = std::chrono::high_resolution_clock::now();
#endif
  bool result = false;

  assertion(_lduh!=nullptr);
  assertion(_lQhbnd!=nullptr);
  assertion(_lGradQhbnd!=nullptr);
  assertion(_lFhbnd!=nullptr);
  assertion(_luh!=nullptr);

  logDebug("handleRemoteExecution",
        " processJob: "<<to_string());

  int iterations=_solver.fusedSpaceTimePredictorVolumeIntegral(
     _lduh,
     _lQhbnd,
     _lGradQhbnd,
     _lFhbnd,
     _luh,
     _center,
     _dx,
     _predictorTimeStamp,
     _predictorTimeStepSize,
     true);

  logDebug("handleRemoteExecution",
        " finished job: "<<to_string());
   
  hasComputed = true;
#if defined(USE_ITAC)
  VT_end(event_stp_remote);
#endif
#if defined(FileTrace)
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

  // exahype::reactive::STPStatsTracer::getInstance().writeTracingEventIteration(iterations, exahype::reactive::STPTraceKey::ADERDGRemoteMigratable);
  exahype::reactive::STPStatsTracer::getInstance().writeTracingEventRunIterations(duration.count(), iterations, exahype::reactive::STPTraceKey::ADERDGRemoteMigratable);
#endif
  return result;
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::sendBackOutcomeToOrigin() {
  logDebug("sendBackOutcomeToOrigin",
        " send job outcome: "<<to_string());

  assertion(!_isSkeleton); //skeleton jobs should not be offloaded to a different rank
    //logInfo("handleLocalExecution()", "postSendBack");
#if defined(UseSmartMPI)
#if defined(SmartMPINB)
  MPI_Request sendBackRequests[NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME];
  _solver.mpiIsendMigratablePredictionJobOutcomeOffload(
       _lduh,
       _lQhbnd,
       _lFhbnd,
       _lGradQhbnd,
       _originRank,
       _tag,
       exahype::reactive::OffloadingContext::getInstance().getMPICommunicatorMapped(),
       sendBackRequests);
  exahype::reactive::OffloadingContext::getInstance().submitRequests(
       sendBackRequests,
       NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME,
       _tag,
       _originRank,
       sendBackHandler,
       exahype::reactive::RequestType::sendBack,
       &_solver);
#else
  _solver.mpiSendMigratablePredictionJobOutcomeOffload(
      _lduh,
      _lQhbnd,
      _lFhbnd,
      _lGradQhbnd,
      _originRank,
      _tag,
      exahype::reactive::OffloadingContext::getInstance().getMPICommunicatorMapped()
  );
  MigratablePredictionJob::sendBackHandler(&_solver, _tag, _originRank);
#endif
#else
  MPI_Request sendBackRequests[NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME];
  _solver.mpiIsendMigratablePredictionJobOutcome(
       _lduh,
       _lQhbnd,
       _lFhbnd,
       _lGradQhbnd,
       _originRank,
       _tag,
       exahype::reactive::OffloadingContext::getInstance().getMPICommunicatorMapped(),
       sendBackRequests);
  exahype::reactive::RequestManager::getInstance().submitRequests(
       sendBackRequests,
       NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME,
       _tag,
       _originRank,
       sendBackHandler,
       exahype::reactive::RequestType::sendBack,
       &_solver);
#endif
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::setSTPPotCorrupted(bool flipped) {

  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
                                                        _element);

  //trigger is set later if we are using limiter as a trigger
  _isPotSoftErrorTriggered =  (exahype::reactive::ResilienceTools::CheckFlipped && flipped)
                            ||  exahype::reactive::ResilienceTools::CheckAllMigratableSTPs;

  logInfo("setSTPPotCorrupted", " celldesc ="<<_cellDescriptionsIndex<<" isPotCorrupted "<<_isPotSoftErrorTriggered);
}

//Caution: Compression is not supported yet!
void exahype::solvers::ADERDGSolver::MigratablePredictionJob::setState(State newState) {
  _currentState = newState;
}
/*
void exahype::solvers::ADERDGSolver::MigratablePredictionJob::recoverWithOutcome(MigratablePredictionJobData *outcome) {
  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
        _element);

  double *luh = static_cast<double*>(cellDescription.getSolution());
  double *lduh = static_cast<double*>(cellDescription.getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double *lGradQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient());
  double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

  logError("recoverWithOutcome", "recovering sane solution and correct STP");
  std::memcpy(luh, outcome->_luh.data(), outcome->_luh.size() * sizeof(double));
  std::memcpy(lduh, outcome->_lduh.data(), outcome->_lduh.size() * sizeof(double));
  std::memcpy(lQhbnd, outcome->_lQhbnd.data(), outcome->_lQhbnd.size() * sizeof(double));
  std::memcpy(lFhbnd, outcome->_lFhbnd.data(), outcome->_lFhbnd.size() * sizeof(double));
#if OffloadingGradQhbnd
  std::memcpy(_lGradQhbnd, &outcome->_lGradQhbnd.data(), outcome->_lGradQhbnd.size() * sizeof(double));
#endif
  exahype::reactive::JobTableStatistics::getInstance().notifyHealedTask();
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::tryFindOutcomeAndHeal() {
#if !defined(OffloadingUseProgressThread)
  exahype::solvers::ADERDGSolver::progressOffloading(&_solver, false, MAX_PROGRESS_ITS);
#endif
  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
        _element);

  //std::cout<<"try to heal..."<<std::endl;

  JobOutcomeStatus status;
  MigratablePredictionJobData *outcomes[2];
  int requiredOutcomes = exahype::reactive::OffloadingContext::getInstance().getTMPINumTeams()-1;

  bool found = tryToFindAndExtractEquivalentSharedOutcomes(false, requiredOutcomes, status, outcomes);

  if(found) {
    int saneIdx = -1;
    if(!checkAgainstOutcomesAndFindSaneOne(outcomes, requiredOutcomes, saneIdx)) {
      logInfo("tryFindOutcomeAndHeal", "Healing a cell!");
      assert(saneIdx>=0 && saneIdx<=2); //otherwise, we have a problem
      recoverWithOutcome(outcomes[saneIdx]);
    }

    for(int i = 0; i<requiredOutcomes; i++)
      delete outcomes[i];
    cellDescription.setHasCompletedLastStep(true);
    logError("tryFindOutcomeAndHeal","has healed!");
    return false;
  }
  else {
    //std::cout<<"still waiting!"<<std::endl;
    return true;
  }
}*/

/*bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::tryFindPreviousOutcomeAndCheck() {
  //caution, this may delay setting the completed flag and cause slowdowns for the master!
#if !defined(OffloadingUseProgressThread)
  exahype::solvers::ADERDGSolver::progressOffloading(&_solver, false, MAX_PROGRESS_ITS);
#endif
  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
            _element);

  logInfo("tryFindPreviousOutcomeAndCheck", "Trying to check previously troubled solution..."
          <<_cellDescriptionsIndex
          <<" center[0] = "
          << _center[0]
          <<" center[1] = "
          << _center[1]
          <<" stamp "<<_predictorTimeStamp
          <<" step "<<_predictorTimeStepSize
          <<" previous stamp "<<cellDescription.getPreviousTimeStamp()
          <<" previous step "<<cellDescription.getPreviousTimeStepSize());

  DeliveryStatus status;
  MigratablePredictionJobData *outcomes[2]; //todo: avoid hardcoding
  int required = exahype::reactive::OffloadingContext::getInstance().getTMPINumTeams()-1;
  bool found = tryToFindAndExtractEquivalentSharedOutcomes(true, required, status, outcomes);
  //bool reschedule;

  if(found) {
    //todo: check other data, too?
    double *luh = static_cast<double*>(cellDescription.getSolution());
    double *lduh = static_cast<double*>(cellDescription.getUpdate());
    double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
    double *lGradQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient());
    double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

    if(_isPotSoftErrorTriggered!=outcomes[0]->_metadata._isPotSoftErrorTriggered) {
      logError("tryFindPreviousOutcomeAndCheck","Limiter was not active previously in other team, we must have a soft error!"
          <<_cellDescriptionsIndex
          <<" center[0] = "
          << _center[0]
          <<" center[1] = "
          << _center[1]
          <<" stamp "<<_predictorTimeStamp
          <<" step "<<_predictorTimeStepSize
          <<" previous stamp "<<cellDescription.getPreviousTimeStamp()
          <<" previous step "<<cellDescription.getPreviousTimeStepSize());

      //correct here
      std::memcpy(luh, &outcomes[0]->_luh[0], outcomes[0]->_luh.size() * sizeof(double));
      std::memcpy(lduh, &outcomes[0]->_lduh[0], outcomes[0]->_lduh.size() * sizeof(double));
      std::memcpy(lQhbnd, &outcomes[0]->_lQhbnd[0], outcomes[0]->_lQhbnd.size() * sizeof(double));
      std::memcpy(lFhbnd, &outcomes[0]->_lFhbnd[0], outcomes[0]->_lFhbnd.size() * sizeof(double));
#if OffloadingGradQhbnd
      std::memcpy(lGradQhbnd, &outcomes[0]->_lGradQhbnd[0], outcomes[0]->_lGradQhbnd.size() * sizeof(double));
#endif

      cellDescription.setRefinementStatus(Keep);
      cellDescription.setCorruptionStatus(CorruptedAndCorrected);
      //logError("tryFindPreviousOutcomeAndCheck", " soft error corrected and detected, but we compute limiter next, as error has propagated!");

      cellDescription.setHasCompletedLastStep(true); //let's use the limiter for now to correct on faulty team, the other one will continue as usual
      //MPI_Abort(MPI_COMM_WORLD, -1);

      return false; //run limiter next, todo: for the faulty outcome we could try to use sane solution from other team
    }
    else {
      bool equalSolution = exahype::reactive::ResilienceTools::getInstance().isAdmissibleNumericalError(outcomes[0]->_luh.data(), luh, outcomes[0]->_luh.size());
      if(!equalSolution)
        logError("tryFindPreviousOutcomeAndCheck"," solutions to not match, we must have a soft error but we can't correct!");
      cellDescription.setHasCompletedLastStep(true);

      //return false; //run limiter next for both
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    //limiter was activated for both cells correctly -> run limiter next
    cellDescription.setHasCompletedLastStep(true);
    return false;
  }
  else
    return true;

}*/

//bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::tryFindOutcomeAndCheck() {
//
//  //caution, this may delay setting the completed flag and cause slowdowns for the master!
//#if !defined(OffloadingUseProgressThread)
//  exahype::solvers::ADERDGSolver::progressOffloading(&_solver, false, MAX_PROGRESS_ITS);
//#endif
//
//  DeliveryStatus status;
//  MigratablePredictionJobData *outcome;
//  bool found = tryToFindAndExtractEquivalentSharedOutcome(false, status, &outcome);
//  bool reschedule;
//
//  if(found && status==DeliveryStatus::Received) {
//    if(matches(outcome)) {
//      CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,_element);
//      cellDescription.setHasCompletedLastStep(true);
//      reschedule = false;
//     }
//     else {
//       //soft error detected
//       reschedule = false;
//       MPI_Abort(MPI_COMM_WORLD, -1);
//     }
//  }
//  else {
//    reschedule = true;
//  }
//
//  return reschedule;
//
///*  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
//        _element);
//  tarch::la::Vector<DIMENSIONS, double> center;
//  center = cellDescription.getOffset() + 0.5 * cellDescription.getSize();
//
//  int saneIdx = -1;
//
//  //todo: can we delete here already?
//  if(found && status==JobOutcomeStatus::received) {
//    if(checkAgainstOutcomesAndFindSaneOne(outcomes, required, saneIdx)) {
//	    cellDescription.setHasCompletedLastStep(true);
//	    reschedule = false;
//    }
//    else {
//      //soft error detected
//#if defined(ResilienceHealing)
//      assert(saneIdx>=0 && saneIdx<=2);
//      recoverWithOutcome(outcomes[saneIdx]);
//      cellDescription.setHasCompletedLastStep(true);
//      reschedule = false;
//      logError("tryFindOutcomeAndCheck", "switch on healing mode now!");
//      //_solver.switchToHealingMode();
//      logInfo("tryFindOutcomeAndCheck","tried to recover cell "
//             <<_cellDescriptionsIndex
//             <<" center[0] = "
//             << center[0]
//             <<" center[1] = "
//             << center[1]
//             <<" center[2] = "
//             << center[2]);
//#else
//      reschedule = false;
//      logError("tryFindOutcomeAndCheck", "soft error detected but we are ignoring it...");
//      exahype::reactive::ResilienceTools::getInstance().setCorruptionDetected(true);
//      cellDescription.setHasCompletedLastStep(true);
//      //MPI_Abort(MPI_COMM_WORLD, -1);
//#endif
//    }
//    for(int i=0; i<required; i++)
//      delete outcomes[i];
//  }
//  else {
//    logInfo("tryFindOutcomeAndCheck","Looking for cell "
//             <<_cellDescriptionsIndex
//             <<" center[0] = "
//             << center[0]
//             <<" center[1] = "
//             << center[1]
//#if DIMENSIONS==3 //fixme: need a to string method for a migratable job!
//             <<" center[2] = "
//             << center[2]
//#endif
//             );
//    reschedule = true;
//  }
//  return reschedule;
//*/
//
//}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::tryToFindAndExtractEquivalentSharedOutcome(bool previous, DeliveryStatus &status, MigratablePredictionJobData **outcome) {
  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
      _element);

  logInfo("tryToFindAndExtractEquivalentSharedOutcome", "team = "<<exahype::reactive::OffloadingContext::getInstance().getTMPITeamNumber()
      <<" looking for "<<cellDescription.toString())

  return _solver.tryToFindAndExtractOutcome(_cellDescriptionsIndex, _element, _predictorTimeStamp, _predictorTimeStepSize, status, outcome);
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::matches(MigratablePredictionJobData *data) {
  if(_isCorrupted)
    logError("checkAgainstOutcome", "checking against outcome");
#if defined(USE_TMPI)
  logInfo("matches", "team "<<TMPI_GetTeamNumber()<<" comparing center[0]="<<_center[0]<<" center[1]="<<_center[1]<<" timestamp "<<_predictorTimeStamp<<" with received task outcome "<<data->_metadata.to_string());
#endif

  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
      _element);

  //todo: check full metadata, too?

  double *lduh = static_cast<double*>(cellDescription.getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
#if defined(OffloadingGradQhbnd)
  double *lGradQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient());
#endif
  double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

  bool equal = true;
  bool tmp;

  tmp = data->_metadata._predictorTimeStamp == _predictorTimeStamp; equal &= tmp;
  tmp = data->_metadata._predictorTimeStepSize == _predictorTimeStepSize; equal &= tmp;

  tmp = exahype::reactive::ResilienceTools::getInstance().isAdmissibleNumericalError(data->_lQhbnd.data(), lQhbnd, data->_lQhbnd.size()); equal&=tmp;
  if(!tmp) {
    logError("matchesOtherOutcome", "lQhbnd is not (numerically) equal for cell "<<"center[0]="<<_center[0]<<" center[1]="<<_center[1]<<" timestamp "<<_predictorTimeStamp);
  }

#if defined(OffloadingGradQhbnd)
  tmp = exahype::reactive::ResilienceTools::getInstance().isAdmissibleNumericalError(data->_lGradQhbnd.data(), lGradQhbnd, data->_lGradQhbnd.size()); equal&=tmp;
  if(!tmp) {
    logError("matchesOtherOutcome", "lGradQhbnd is not (numerically) equal for cell "<<"center[0]="<<_center[0]<<" center[1]="<<_center[1]<<" timestamp "<<_predictorTimeStamp);
  }
#endif
  tmp = exahype::reactive::ResilienceTools::getInstance().isAdmissibleNumericalError(data->_lFhbnd.data(), lFhbnd, data->_lFhbnd.size()); equal&=tmp;
  if(!tmp) {
    logError("matchesOtherOutcome", "lFhbnd is not  (numerically) equal for cell "<<"center[0]="<<_center[0]<<" center[1]="<<_center[1]<<" timestamp "<<_predictorTimeStamp);
  }
  tmp = exahype::reactive::ResilienceTools::getInstance().isAdmissibleNumericalError(data->_lduh.data(), lduh, data->_lduh.size()); equal&=tmp;
  if(!tmp) {
    logError("matchesOtherOutcome", "lduh is not  (numerically) equal for cell "<<"center[0]="<<_center[0]<<" center[1]="<<_center[1]<<" timestamp "<<_predictorTimeStamp);
  }

  if(!equal) {
    logError("matchesOtherOutcome", "soft error detected: "<<data->_metadata.to_string());
    exahype::reactive::ResilienceStatistics::getInstance().notifyDetectedError();
  }

  exahype::reactive::ResilienceStatistics::getInstance().notifyDoubleCheckedTask();
  return equal;
}

/*
bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::checkAgainstOutcomesAndFindSaneOne(MigratablePredictionJobData **data, int outcomes, int &idxSane) {

  bool equal = true;
  assertion(outcomes>0);

  equal = matches(data[0]);
  if(equal || !data[0]->_metadata._isPotSoftErrorTriggered) {
    idxSane = 0;
  }

  if(outcomes==1 || equal) {
    logDebug("checkAgainstOutcomesAndFindSaneOne", "checked duplicate executions for soft errors, result = "<<equal);
    exahype::reactive::JobTableStatistics::getInstance().notifyDoubleCheckedTask();
    return equal; //already found two matching outcomes or we only have a single one to check against
  }

  assertion(!equal);
  //A!=B  -> there must be some error somewhere
  //majority vote to find sane one
  if(matches(data[1])) {
    idxSane = 1;
    logDebug("checkAgainstOutcomesAndFindSaneOne", "checked duplicate executions for soft errors, result = "<<equal);
    exahype::reactive::JobTableStatistics::getInstance().notifyDoubleCheckedTask();
    return true; // A != B && A==C -> C is correct
  }
  else { // A!=B && A!=C
    if(data[0]->matches(*data[1])) {
      //  A!=B && A!=C && B==C
      idxSane = 1;
      return false; //A is wrong
    }
    else {
      logError("checkAgainstOutcomesAndFindSaneOne", "We have three different outcomes, there must be something wrong here!");
      assertion(false);
      MPI_Abort(MPI_COMM_WORLD, -1);
      return false;
    }
  }

  return equal;
}*/

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::sendHandler(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {
#if defined(PerformanceAnalysisOffloadingDetailed)
  static std::atomic<int> cnt=0;
  cnt++;
#endif
  logDebug("sendHandler","successful send request tag = "<<tag<<" remoteRank = "<<remoteRank);
#if !defined(OffloadingNoEarlyReceiveBacks) || defined(OffloadingLocalRecompute)
  ADERDGSolver::receiveBackMigratableJob(tag, remoteRank,
      static_cast<exahype::solvers::ADERDGSolver*>(solver));
#endif
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveHandler(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {
  //logInfo("receiveHandler","successful receive request");

  tbb::concurrent_hash_map<std::pair<int, int>, MigratablePredictionJobData*>::accessor a_tagRankToData;
  MigratablePredictionJobData *data;
  bool found = static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToStolenData.find(
      a_tagRankToData, std::make_pair(remoteRank, tag));
  if(!found)
    logError("receiveHandler","Didn't find migrated data, maps are inconsistent!");
  assertion(found);
  data = a_tagRankToData->second;
#if defined(UseSmartMPI) || defined(OffloadingMetadataPacked)
  // hack: when sending back a job outcome, remoteRank is different (it is then the actual client rank and not a SmartMPI server rank)
  // therefore, we adapt the rank here
  data->_metadata.unpackContiguousBuffer();
#endif
#if defined(UseSmartMPI)
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToStolenData.erase(a_tagRankToData);
  a_tagRankToData.release();
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToStolenData.insert(std::make_pair(
                                                                                       std::make_pair(data->_metadata.getOrigin(), tag),
                                                                                       data));
#else
  a_tagRankToData.release();
#endif

  exahype::reactive::OffloadingAnalyser::getInstance().notifyReceivedSTPJob();
  MigratablePredictionJob *job =
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->createFromData(data,
          data->_metadata.getOrigin(), tag);
  peano::datatraversal::TaskSet spawnedSet(job);

  logDebug("receiveHandler",
      " received task : "<< data->_metadata.to_string()<<" from "<<data->_metadata.getOrigin()<<" tag = "<<tag);

  exahype::reactive::OffloadingProfiler::getInstance().notifyReceivedTask(data->_metadata.getOrigin());
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::sendBackHandler(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {
  //logInfo("sendBackHandler","successful sendBack request");
  tbb::concurrent_hash_map<std::pair<int, int>, MigratablePredictionJobData*>::accessor a_tagRankToData;

  logDebug("sendBackHandler", "looking for tag = "<<tag<<" remoteRank = "<<remoteRank);

  MigratablePredictionJobData *data;
  bool found =
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToStolenData.find(
          a_tagRankToData, std::make_pair(remoteRank, tag));
  if(!found)
    logError("sendBackHandler","Couldn't find sent back task data. Inconsistent maps! Tag = "<<tag<<" rank = "<<remoteRank);
  assertion(found);
  data = a_tagRankToData->second;
  a_tagRankToData.release();
  delete data;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToStolenData.erase(
      std::make_pair(remoteRank, tag));

  NumberOfStolenJobs--;
  assertion( NumberOfStolenJobs>=0 );
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveBackHandler(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {

  // logInfo("receiveBackHandler","successful receiveBack request");
  logDebug("receiveBackHandler", "received back STP job tag="<<tag<<" rank="<<remoteRank);
  tbb::concurrent_hash_map<int, CellDescription*>::accessor a_tagToCellDesc;
  bool found =
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToCellDesc.find(
          a_tagToCellDesc, tag);
  if(!found)
    logError("receiveBackHandler", "Received back task but couldn't find cell, inconsistent maps! Tag = "<<tag<<" remote rank = "<<remoteRank);
  assertion(found);
  auto cellDescription = a_tagToCellDesc->second;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToCellDesc.erase(
      a_tagToCellDesc);
  a_tagToCellDesc.release();

  tbb::concurrent_hash_map<int, MigratablePredictionJobMetaData*>::accessor a_tagToMetaData;
  found =
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToMetaData.find(
          a_tagToMetaData, tag);
  assertion(found);
  MigratablePredictionJobMetaData *metadata = a_tagToMetaData->second;
  delete metadata;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToMetaData.erase(
      a_tagToMetaData);
  a_tagToMetaData.release();

#ifndef OffloadingLocalRecompute
  tbb::concurrent_hash_map<const CellDescription*, std::pair<int,int>>::accessor a_cellDescToTagRank;
  found = static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapCellDescToTagRank.find(a_cellDescToTagRank, cellDescription);
  assertion(found);
  // do not erase for local recompute as we need this information later on
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapCellDescToTagRank.erase(a_cellDescToTagRank);
  a_cellDescToTagRank.release();
#endif

  /*tbb::concurrent_hash_map<int, double>::accessor a_tagToOffloadTime;
  found =
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToOffloadTime.find(
          a_tagToOffloadTime, tag);
  double elapsed = MPI_Wtime() + a_tagToOffloadTime->second;
  assertion(found);
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToOffloadTime.erase(
      a_tagToOffloadTime);
  a_tagToOffloadTime.release();*/

  NumberOfRemoteJobs--;
  NumberOfEnclaveJobs--;
#ifndef OffloadingLocalRecompute
  cellDescription->setHasCompletedLastStep(true);
#else
  logDebug("receiveBackHandler", "received back STP job tag="<<tag<<" rank="<<remoteRank);
  MigratablePredictionJobData *data = nullptr;
  tbb::concurrent_hash_map<int, MigratablePredictionJobData*>::accessor a_tagToData;
  found = static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToSTPData.find(
          a_tagToData, tag);
  assertion(found);
  data = a_tagToData->second;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToSTPData.erase(
      a_tagToData);
  a_tagToData.release();

  tbb::concurrent_hash_map<int, double*>::accessor a_tagToMetaData;
  found =
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToMetaData.find(
          a_tagToMetaData, tag);
  assertion(found);
  double *metadata = a_tagToMetaData->second;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToMetaData.erase(
      a_tagToMetaData);
  a_tagToMetaData.release();

  logDebug("receiveBackHandler",
      " received task outcome: center[0] = "<<metadata[0]
      <<" center[1] = "<<metadata[1]
#if DIMENSIONS==3
      <<" center[2] = "<<metadata[2]
#endif
      <<" time stamp = "<<metadata[2*DIMENSIONS] <<" element = "<<(int) metadata[2*DIMENSIONS+2]);

  exahype::reactive::ResilienceStatistics::getInstance().notifyReceivedTask();
  //timestamp
  if (metadata[2 * DIMENSIONS]
      < static_cast<exahype::solvers::ADERDGSolver*>(solver)->getMinTimeStamp()) {
    // if(true) {
    exahype::reactive::ResilienceStatistics::getInstance().notifyLateTask();
    delete data;
    AllocatedSTPsReceive--;
  }
  else {
    tarch::multicore::Lock lock(
        exahype::solvers::ADERDGSolver::EmergencySemaphore);
    tarch::multicore::jobs::Job *recompJob =
        static_cast<exahype::solvers::ADERDGSolver*>(solver)->grabRecomputeJobForCellDescription(
            (const void*) cellDescription);

    //hack: put the job back and ignore very late result
    if (recompJob != nullptr
        && static_cast<MigratablePredictionJob*>(recompJob)->_predictorTimeStamp
            > metadata[2 * DIMENSIONS]) {
      logDebug("receiveBackHandler",
          "job timestamp "<<static_cast<MigratablePredictionJob*>(recompJob)->_predictorTimeStamp
          <<" metadata[2*DIMENSIONS] "<< metadata[2*DIMENSIONS]);
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->addRecomputeJobForCellDescription(
          recompJob, cellDescription);
      recompJob = nullptr;
    }

    //copy into result buffer as I am responsible for result
    if (recompJob != nullptr) {
      double *lduh = static_cast<double*>(cellDescription->getUpdate());
      double *lQhbnd = static_cast<double*>(cellDescription->getExtrapolatedPredictor());
      double *lFhbnd = static_cast<double*>(cellDescription->getFluctuation());
      double *lGradQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient());

      /*if(static_cast<MigratablePredictionJob*>(recompJob)->_predictorTimeStamp!=metadata[2*DIMENSIONS]) {
       logInfo("receiveBackHandler","job timestamp "<<static_cast<MigratablePredictionJob*>(recompJob)->_predictorTimeStamp
       <<" metadata[2*DIMENSIONS] "<< metadata[2*DIMENSIONS]);
       }*/

      assertion(
          static_cast<MigratablePredictionJob*>(recompJob)->_predictorTimeStamp
              == metadata[2 * DIMENSIONS]);

      std::memcpy(lduh, &data->_lduh[0], data->_lduh.size() * sizeof(double));
      std::memcpy(lQhbnd, &data->_lQhbnd[0],
          data->_lQhbnd.size() * sizeof(double));
      std::memcpy(lFhbnd, &data->_lFhbnd[0],
          data->_lFhbnd.size() * sizeof(double));
#if defined(OffloadingGradQhbnd)
      std::memcpy(lGradQhbnd, &data->_lGradQhbnd[0],
          data->_lGradQhbnd.size() * sizeof(double));
#endif
      exahype::reactive::ResilienceStatistics::getInstance().notifySavedTask();

      delete data;
      AllocatedSTPsReceive--;
      delete recompJob;

      tbb::concurrent_hash_map<const CellDescription*, std::pair<int, int>>::accessor a_cellDescToTagRank;
      //logInfo("receiveBackHandler", " cleaning up cell desc to tag/rank for "<<cellDescription);
      found =
          static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapCellDescToTagRank.find(
              a_cellDescToTagRank, cellDescription);
      assertion(found);
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapCellDescToTagRank.erase(
          a_cellDescToTagRank);
      a_cellDescToTagRank.release();

      logDebug("receiveBackHandler", " applied task outcome");

      cellDescription->setHasCompletedLastStep(true);
    }
    //sb else had done it, resolve emergency
    else {
      if (LastEmergencyCell == cellDescription) {
        VetoEmergency = false;
        LastEmergencyCell = nullptr;
      }

      exahype::reactive::ResilienceStatistics::getInstance().notifyLateTask();
      //delete data; //probably not safe to do here
      // AllocatedSTPsReceive--; //probably not safe to do here, race with job
    }
    lock.free();
  }
#endif
  //logInfo("receiveBackHandler", "remote execution took "<<elapsed<<" s ");

  assertion( NumberOfEnclaveJobs>=0 ); assertion( NumberOfRemoteJobs>=0 );
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveHandlerTaskSharing(
    exahype::solvers::Solver* solver, int tag, int team) {
  logDebug("receiveHandlerTaskSharing","successful receive request");

  tbb::concurrent_hash_map<std::pair<int, int>, MigratablePredictionJobData*>::accessor a_tagRankToData;
  MigratablePredictionJobData *data;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToReplicaData.find(
      a_tagRankToData, std::make_pair(team, tag));
  data = a_tagRankToData->second;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToReplicaData.erase(
      a_tagRankToData);
  a_tagRankToData.release();

#if defined(UseSmartMPI) || defined(OffloadingMetadataPacked)
  data->_metadata.unpackContiguousBuffer();
#endif

  MigratablePredictionJobOutcomeKey key(data->_metadata.getCenter(), data->_metadata.getPredictorTimeStamp(),
                                         data->_metadata.getPredictorTimeStepSize(), data->_metadata.getElement());
  logDebug("receiveHandlerTaskSharing", "team "
      <<exahype::reactive::OffloadingContext::getInstance().getTMPITeamNumber()
      <<" received replica job: "
      <<data->_metadata.to_string());

  double minTimeStampToKeep = static_cast<exahype::solvers::ADERDGSolver*>(solver)->getPreviousMinTimeStamp();

  exahype::reactive::TimeStampAndLimiterTeamHistory::getInstance().trackTimeStepAndLimiterActive(team, key._timestamp, key._timestepSize, data->_metadata._isPotSoftErrorTriggered);
  double lastconsistentTimeStamp, lastconsistentTimeStepSize, lastconsistentEstimatedTimeStepSize;
  exahype::reactive::TimeStampAndLimiterTeamHistory::getInstance().getLastConsistentTimeStepData(lastconsistentTimeStamp, lastconsistentTimeStepSize, lastconsistentEstimatedTimeStepSize);

  if(exahype::reactive::OffloadingContext::getInstance().getResilienceStrategy()
    >= exahype::reactive::OffloadingContext::ResilienceStrategy::TaskSharingResilienceChecks) {
    minTimeStampToKeep = std::min(minTimeStampToKeep, lastconsistentTimeStepSize);
  }

  if(key._timestamp<minTimeStampToKeep)
  {
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
      logDebug("receiveHandlerTaskSharing", "team "
          <<exahype::reactive::OffloadingContext::getInstance().getTMPITeamNumber()
          <<" inserted replica job: "
          <<data2->_metadata.to_string());
          //<<std::setprecision(30)<<"timestep "<<data->_metadata.getPredictorTimeStepSize());
    }
    else {
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_outcomeDatabase.insertOutcome(key, data, DeliveryStatus::Received);
      logDebug("receiveHandlerTaskSharing", "team "
          <<exahype::reactive::OffloadingContext::getInstance().getTMPITeamNumber()
          <<" inserted replica job: "
          <<data->_metadata.to_string());
          //<<std::setprecision(30)<<"timestep "<<data->_metadata.getPredictorTimeStepSize());
    }
    static_cast<exahype::solvers::ADERDGSolver*>(solver)->_allocatedOutcomes.push_back(key);
  }
  exahype::reactive::ResilienceStatistics::getInstance().notifyReceivedTask();

}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::sendHandlerTaskSharing(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {

  logDebug("sendHandlerReplication","successfully completed send to other teams");
  exahype::reactive::ResilienceStatistics::getInstance().notifySentTask();
  tbb::concurrent_hash_map<int, MigratablePredictionJobData*>::accessor a_tagToData;
  bool found =
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToSTPData.find(
          a_tagToData, tag);
  assertion(found);
  MigratablePredictionJobData *data = a_tagToData->second;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToSTPData.erase(
      a_tagToData);
  delete data;
  AllocatedSTPsSend--;
  CompletedSentSTPs++;
  logDebug("sendHandlerReplication"," allocated stps send "<<AllocatedSTPsSend);
  a_tagToData.release();
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::packMetaData(MigratablePredictionJobMetaData *metadata) {
#if defined(UseSmartMPI) || defined(OffloadingMetadataPacked)
  char *tmp = metadata->_contiguousBuffer;
  std::memcpy(tmp, _center, DIMENSIONS*sizeof(double)); tmp+= DIMENSIONS*sizeof(double);
  std::memcpy(tmp, _dx, DIMENSIONS*sizeof(double)); tmp+= DIMENSIONS*sizeof(double);
  std::memcpy(tmp, &_predictorTimeStamp, sizeof(double)); tmp+= sizeof(double);
  std::memcpy(tmp, &_predictorTimeStepSize, sizeof(double)); tmp+= sizeof(double);
  std::memcpy(tmp, &_element, sizeof(int)); tmp+= sizeof(int);
  std::memcpy(tmp, &_originRank, sizeof(int)); tmp+= sizeof(int);
  std::memcpy(tmp, &_isPotSoftErrorTriggered, sizeof(bool)); tmp+= sizeof(bool);
  std::memcpy(tmp, &_isCorrupted, sizeof(bool)); tmp+= sizeof(bool);
#else
  metadata->_center[0] = _center[0];
  metadata->_center[1] = _center[1];
#if DIMENSIONS==3
  metadata->_center[2] = _center[2];
#endif

  metadata->_dx[0] = _dx[0];
  metadata->_dx[1] = _dx[1];
#if DIMENSIONS==3
  metadata->_dx[2] = _dx[2];
#endif

  metadata->_predictorTimeStamp = _predictorTimeStamp;
  metadata->_predictorTimeStepSize = _predictorTimeStepSize;
  metadata->_element = _element;
  metadata->_originRank = _originRank;
  metadata->_isPotSoftErrorTriggered = _isPotSoftErrorTriggered;
  metadata->_isCorrupted = _isCorrupted;
#endif
}

std::string exahype::solvers::ADERDGSolver::MigratablePredictionJob::to_string() const {
  std::string result;

  result =  " cellDescriptionIndex = " + std::to_string(_cellDescriptionsIndex);
  result += " center[0] = " + std::to_string(_center[0]);
  result += " center[1] = " + std::to_string(_center[1]);
  #if DIMENSIONS==3
  result += " center[2] = " + std::to_string(_center[2]);
  #endif

  result += " time stamp = " + std::to_string(_predictorTimeStamp);
  result += " time step = " + std::to_string(_predictorTimeStepSize);
  result += " element = " + std::to_string(_element);
  result += " origin = " + std::to_string(_originRank);
  result += " isPotSoftErrorTriggered = " + std::to_string(_isPotSoftErrorTriggered);
  result += " isCorrupted = " + std::to_string(_isCorrupted);

  return result;
}


exahype::solvers::ADERDGSolver::MigratablePredictionJobData::MigratablePredictionJobData(
    ADERDGSolver& solver) :
      _luh(solver.getDataPerCell()),
      _lduh(solver.getUpdateSize()),
      _lQhbnd(solver.getBndTotalSize()),
      _lFhbnd(solver.getBndFluxTotalSize()),
      _lGradQhbnd(solver.getBndGradQTotalSize()){
  AllocatedSTPs++;
}

exahype::solvers::ADERDGSolver::MigratablePredictionJobData::~MigratablePredictionJobData() {
  AllocatedSTPs--;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJobData::matches(const MigratablePredictionJobData& a) const {
  bool equal = true;
  bool tmp;
  assertion(a._lQhbnd.size()==_lQhbnd.size());
  assertion(a._lFhbnd.size()==_lFhbnd.size());
  assertion(a._lduh.size()==_lduh.size());
  assertion(a._lGradQhbnd.size()==_lGradQhbnd.size());

  tmp = _metadata == a._metadata; equal&=tmp;

  tmp = exahype::reactive::ResilienceTools::getInstance().isAdmissibleNumericalError(a._lQhbnd.data(), _lQhbnd.data(), a._lQhbnd.size()); equal&=tmp;
  if(!tmp) {
    logError("matches", "lQhbnd is not (numerically) equal");
  }
#if defined(OffloadingGradQhbnd)
  tmp = exahype::reactive::ResilienceTools::getInstance().isAdmissibleNumericalError(a._lGradQhbnd.data(), _lGradQhbnd, a._lGradQhbnd.size()); equal&=tmp;
  if(!tmp) {
    logError("matches", "lGradQhbnd is not (numerically) equal");
  }
#endif
  tmp = exahype::reactive::ResilienceTools::getInstance().isAdmissibleNumericalError(a._lFhbnd.data(), _lFhbnd.data(), a._lFhbnd.size()); equal&=tmp;
  if(!tmp) {
    logError("matches", "lFhbnd is not  (numerically) equal");
  }
  tmp = exahype::reactive::ResilienceTools::getInstance().isAdmissibleNumericalError(a._lduh.data(), _lduh.data(), a._lduh.size()); equal&=tmp;
  if(!tmp) {
    logError("matches", "lduh is not  (numerically) equal");
  }
  return equal;
}

exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::MigratablePredictionJobMetaData()
: _predictorTimeStamp(0),
  _predictorTimeStepSize(0),
  _element(0),
  _originRank(-1),
  _isPotSoftErrorTriggered(true),
  _isCorrupted(false),
  _contiguousBuffer(nullptr) {
#if defined(UseSmartMPI) || defined(OffloadingMetadataPacked)
  _contiguousBuffer = (char*) allocate_smartmpi(getMessageLen(),0); //todo: does Alignment help here?
#endif
}

exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::~MigratablePredictionJobMetaData() {
#if defined(UseSmartMPI) || defined(OffloadingMetadataPacked)
  free_smartmpi(_contiguousBuffer);
#endif
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::unpackContiguousBuffer() {
  char *tmp = _contiguousBuffer;
  std::memcpy(_center, tmp, DIMENSIONS*sizeof(double)); tmp+= DIMENSIONS*sizeof(double);
  std::memcpy(_dx, tmp, DIMENSIONS*sizeof(double)); tmp+= DIMENSIONS*sizeof(double);
  std::memcpy(&_predictorTimeStamp, tmp, sizeof(double)); tmp+= sizeof(double);
  std::memcpy(&_predictorTimeStepSize, tmp,  sizeof(double)); tmp+= sizeof(double);
  std::memcpy(&_element, tmp, sizeof(int)); tmp+= sizeof(int);
  std::memcpy(&_originRank, tmp, sizeof(int)); tmp+= sizeof(int);
  std::memcpy(&_isPotSoftErrorTriggered, tmp, sizeof(bool)); tmp+= sizeof(bool);
  std::memcpy(&_isCorrupted, tmp, sizeof(bool)); tmp+= sizeof(bool);
}

std::string exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::to_string() const {
  std::string result;

  result = " center[0] = " + std::to_string(_center[0]);
  result += " center[1] = " + std::to_string(_center[1]);
  #if DIMENSIONS==3
  result += " center[2] = " + std::to_string(_center[2]);
  #endif

  result += " time stamp = " + std::to_string(_predictorTimeStamp);
  result += " time step = " + std::to_string(_predictorTimeStepSize);
  result += " element = " + std::to_string(_element);
  result += " origin = " + std::to_string(_originRank);
  result += " isPotSoftErrorTriggered = " + std::to_string(_isPotSoftErrorTriggered);
  result += " isCorrupted = " + std::to_string(_isCorrupted);

  return result;
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::initDatatype() {
#if defined(UseSmartMPI) || defined(OffloadingMetadataPacked)
  _datatype = MPI_BYTE;
#else
  int entries = 2+6;
  MigratablePredictionJobMetaData dummy;

  int blocklengths[] = {DIMENSIONS, DIMENSIONS, 1, 1, 1, 1, sizeof(bool), sizeof(bool)};
  MPI_Datatype subtypes[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_BYTE, MPI_BYTE};
  MPI_Aint displs[] = {0, 0, 0, 0, 0, 0, 0, 0};

  MPI_Aint base;
  MPI_Get_address(&dummy, &base);

  MPI_Get_address(&(dummy._center), &(displs[0]));
  MPI_Get_address(&(dummy._dx), &(displs[1]));
  MPI_Get_address(&(dummy._predictorTimeStamp), &(displs[2]));
  MPI_Get_address(&(dummy._predictorTimeStepSize), &(displs[3]));
  MPI_Get_address(&(dummy._element), &(displs[4]));
  MPI_Get_address(&(dummy._originRank), &(displs[5]));
  MPI_Get_address(&(dummy._isPotSoftErrorTriggered), &(displs[6]));
  MPI_Get_address(&(dummy._isCorrupted), &(displs[7]));

  for(int i=0; i<entries; i++) {
    displs[i] = displs[i]-base;
    logDebug("initDatatype" , " displ["<<i<<"]="<<displs[i]);
  }

  int ierr = MPI_Type_create_struct(entries, blocklengths, displs, subtypes, &_datatype); assert(ierr==MPI_SUCCESS);
  ierr = MPI_Type_commit(&_datatype);  assert(ierr==MPI_SUCCESS);
#endif
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::shutdownDatatype() {
#if !(defined(UseSmartMPI) || defined(OffloadingMetadataPacked))
  MPI_Type_free(&_datatype);
#endif
}

#endif
//#undef assertion
//#define assertion(expr) 
