//#if defined(DistributedOffloading) && defined(Parallel) && defined(SharedTBB)
#if defined(DistributedOffloading) && defined(Parallel) 

#if defined(ScoreP)
#include "scorep/SCOREP_User.h"
#endif

#if defined(FileTrace)
#include "exahype/offloading/STPStatsTracer.h"
#endif

#include "ADERDGSolver.h"

#include "exahype/offloading/PerformanceMonitor.h"
#include "exahype/offloading/OffloadingAnalyser.h"
#include "exahype/offloading/OffloadingProfiler.h"
#include "exahype/offloading/JobTableStatistics.h"
#include "exahype/offloading/MemoryMonitor.h"
#include "exahype/offloading/NoiseGenerator.h"
#include "../offloading/ResilienceTools.h"

#define MAX_PROGRESS_ITS 100
//#undef assertion
//#define assertion assert

MPI_Datatype exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::_datatype;

exahype::solvers::ADERDGSolver::MigratablePredictionJob::MigratablePredictionJob(
    ADERDGSolver& solver, const int cellDescriptionsIndex, const int element,
    const double predictorTimeStamp, const double predictorTimeStepSize) :
      tarch::multicore::jobs::Job(
          tarch::multicore::jobs::JobType::BackgroundTask, 0,
          getTaskPriorityLocalStealableJob(cellDescriptionsIndex, element,
              predictorTimeStamp)),  //this is a locally executed job
      _solver(solver),
      _cellDescriptionsIndex(cellDescriptionsIndex),
      _element(element),
      _predictorTimeStamp(predictorTimeStamp),
      _predictorTimeStepSize(predictorTimeStepSize),
      _originRank(tarch::parallel::Node::getInstance().getRank()),
      _tag(-1),
      _luh(nullptr),
      _lduh(nullptr),
      _lQhbnd(nullptr),
      _lFhbnd(nullptr),
      _lGradQhbnd(nullptr),
      _isLocalReplica(false),
      _isPotSoftErrorTriggered(0),
	  _currentState(State::INITIAL)
{
  LocalStealableSTPCounter++;
  NumberOfEnclaveJobs++;
  exahype::offloading::JobTableStatistics::getInstance().notifySpawnedTask();
  exahype::offloading::PerformanceMonitor::getInstance().incCurrentTasks();

  auto& cellDescription = getCellDescription(cellDescriptionsIndex, element);

  tarch::la::Vector<DIMENSIONS, double> center = cellDescription.getOffset()+0.5*cellDescription.getSize();
  tarch::la::Vector<DIMENSIONS, double> dx = cellDescription.getSize();

  for (int i = 0; i < DIMENSIONS; i++) {
    _center[i] = center[i];
    _dx[i] = dx[i];
  }
}

exahype::solvers::ADERDGSolver::MigratablePredictionJob::MigratablePredictionJob(
    ADERDGSolver& solver, const int cellDescriptionsIndex, const int element,
    const double predictorTimeStamp, const double predictorTimeStepSize,
    double *luh, double *lduh, double *lQhbnd, double *lFhbnd, double *lGradQhbnd, double *dx,
    double *center, const int originRank, const int tag) :
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
      _luh(luh),
      _lduh(lduh),
      _lQhbnd(lQhbnd),
      _lFhbnd(lFhbnd),
      _lGradQhbnd(lGradQhbnd),
      _isLocalReplica(false),
      _isPotSoftErrorTriggered(0),
	  _currentState(State::INITIAL)
{

  for (int i = 0; i < DIMENSIONS; i++) {
    _center[i] = center[i];
    _dx[i] = dx[i];
  }

  if (_originRank != tarch::parallel::Node::getInstance().getRank()) {
    NumberOfStolenJobs++;
  }
  else
    NumberOfEnclaveJobs++;
  exahype::offloading::PerformanceMonitor::getInstance().incCurrentTasks();
}

exahype::solvers::ADERDGSolver::MigratablePredictionJob::~MigratablePredictionJob() {
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::setTrigger(bool flipped) {
  //todo: need to implement some criterion
  //set to true for now, won't benefit from task sharing then, though
  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
      _element);

  _isPotSoftErrorTriggered =  (exahype::offloading::ResilienceTools::TriggerFlipped && flipped)
                           || (exahype::offloading::ResilienceTools::TriggerLimitedCellsOnly && cellDescription.getIsTroubledInLastStep())
                           ||  exahype::offloading::ResilienceTools::TriggerAllMigratableSTPs;
}

//Caution: Compression is not supported yet!
bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::run(
    bool isCalledOnMaster) {

  bool reschedule;

  switch(_currentState)
  {
    case State::INITIAL:
      reschedule = runExecution(isCalledOnMaster);
      break;
    case State::CHECK_REQUIRED:
      reschedule = tryFindOutcomeAndCheck();
      break;
  }
  return reschedule;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::tryFindOutcomeAndCheck() {

#if !defined(OffloadingUseProgressThread)
  exahype::solvers::ADERDGSolver::progressOffloading(&_solver, false, MAX_PROGRESS_ITS);
#endif

  JobOutcomeStatus status;
  MigratablePredictionJobData *outcome;
  bool found = tryToFindAndExtractEquivalentSharedOutcome(status, &outcome);
  bool reschedule;

  if(found && status==JobOutcomeStatus::received) {
	if(matchesOtherOutcome(outcome)) {

	  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
			      _element);
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
      exahype::offloading::OffloadingAnalyser::getInstance().setTimePerSTP(
          watch.getCalendarTime());
    }
  }
  else
    result = handleExecution(isCalledOnMaster, hasComputed);

#if defined(TaskSharing) && !defined(OffloadingUseProgressThread)
 // exahype::solvers::ADERDGSolver::pollForOutstandingCommunicationRequests(&_solver, false, MAX_PROGRESS_ITS);
  exahype::solvers::ADERDGSolver::progressOffloading(&_solver, false, MAX_PROGRESS_ITS);
#endif

#if defined(GenerateNoise)
    exahype::offloading::NoiseGenerator::getInstance().generateNoiseSTP();
#endif

#ifdef USE_ITAC
  VT_end(event_stp);
#endif
  return result;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::handleLocalExecution(
    bool isRunOnMaster, bool& hasComputed) {

  bool reschedule = false;
  bool needToCompute = true;
  bool needToCheck = false;
  bool copyResult = false;
  bool hasResult = false;

  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
	   _element);

  double *luh = static_cast<double*>(cellDescription.getSolution());
  double *lduh = static_cast<double*>(cellDescription.getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double *lGradQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient());
  double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

#if defined(TaskSharing)
  bool needToShare = (AllocatedSTPsSend
          <= exahype::offloading::PerformanceMonitor::getInstance().getTasksPerTimestep());

  JobOutcomeStatus status;
  MigratablePredictionJobData *outcome = nullptr;

  bool found = tryToFindAndExtractEquivalentSharedOutcome(status, &outcome);

  if(found) {
    switch (status) {
      case JobOutcomeStatus::received:
    	hasResult = true;
#if defined(ResilienceChecks)
        needToCompute = outcome->_metadata._isPotSoftErrorTriggered;
        copyResult = !outcome->_metadata._isPotSoftErrorTriggered;
        needToCheck = outcome->_metadata._isPotSoftErrorTriggered;
        needToShare = outcome->_metadata._isPotSoftErrorTriggered;
#else
    	needToShare = false;
		needToCompute = false;
        copyResult = true;
#endif
        break;
      case JobOutcomeStatus::transit:
#ifdef TaskSharingRescheduleIfInTransit
        reschedule = true;
        needToCompute = false;
#endif
        break;
	  }
  }

  if(copyResult) {
    assertion(outcome!=nullptr);
    std::memcpy(lduh, &outcome->_lduh[0], outcome->_lduh.size() * sizeof(double));
    std::memcpy(lQhbnd, &outcome->_lQhbnd[0], outcome->_lQhbnd.size() * sizeof(double));
    std::memcpy(lFhbnd, &outcome->_lFhbnd[0], outcome->_lFhbnd.size() * sizeof(double));
#if OffloadingGradQhbnd
    std::memcpy(lGradQhbnd, &outcome->_lGradQhbnd[0], outcome->_lGradQhbnd.size() * sizeof(double));
#endif
    exahype::offloading::JobTableStatistics::getInstance().notifySavedTask();
    AllocatedSTPsReceive--;
    delete outcome;
  }
#endif /*TaskSharing*/

  if(needToCompute) {
#if defined(USE_ITAC)
    if(_isLocalReplica)
      VT_begin(event_stp_local_replica);
#endif
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

    exahype::offloading::JobTableStatistics::getInstance().notifyExecutedTask();

    bool hasFlipped = exahype::offloading::ResilienceTools::getInstance().generateBitflipErrorInDoubleIfActive(lduh, _solver.getUpdateSize());
    setTrigger(hasFlipped);

#if defined(TaskSharing)
    needToCheck = needToCheck || _isPotSoftErrorTriggered;
    needToShare = needToShare || _isPotSoftErrorTriggered;
#endif

#if defined(FileTrace)
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    exahype::offloading::STPStatsTracer::getInstance().writeTracingEventRunIterations(duration.count(), iterations, exahype::offloading::STPTraceKey::ADERDGOwnMigratable);
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

#if defined(TaskSharing)
  if(!hasResult) {
    hasResult = tryToFindAndExtractEquivalentSharedOutcome(status, &outcome);
  }

  if(needToCheck) {
    if(hasResult) {
#if defined(ResilienceChecks)
      matchesOtherOutcome(outcome);
#endif
    }
    else {
      _currentState = State::CHECK_REQUIRED;
      reschedule = true;
    }
  }
  else {
	if(hasResult)
	  needToShare = false;
  }

  if(needToShare) {
    SentSTPs++;
    _solver.sendTaskOutcomeToOtherTeams(this);
  }
#endif

  if (!reschedule) {
	cellDescription.setHasCompletedLastStep(true);
    exahype::offloading::PerformanceMonitor::getInstance().decRemainingTasks();
    exahype::offloading::PerformanceMonitor::getInstance().decCurrentTasks();
  }

  return reschedule;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::tryToFindAndExtractEquivalentSharedOutcome(JobOutcomeStatus &status, MigratablePredictionJobData **outcome) {
  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
	      _element);

  tarch::la::Vector<DIMENSIONS, double> center;
  center = cellDescription.getOffset() + 0.5 * cellDescription.getSize();

  JobTableKey key;
  for (int i=0; i < DIMENSIONS; i++)
    key.center[i] = center[i];
  key.timestamp = _predictorTimeStamp;
  key.timestepSize = _predictorTimeStepSize;
  key.element = _element;

  tbb::concurrent_hash_map<JobTableKey, JobTableEntry>::accessor a_jobToData;
  bool found = _solver._jobDatabase.find(a_jobToData, key);

  if(found && a_jobToData->second.status == JobOutcomeStatus::received) {
    *outcome = a_jobToData->second.data;
    status = a_jobToData->second.status;
    _solver._jobDatabase.erase(a_jobToData);

    logInfo("tryToFindAndExtractEquivalentSharedOutcome()",
        "team "<<exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank()
        <<" found STP in received jobs:"
        <<(*outcome)->_metadata.to_string());
  }
  else if(found && a_jobToData->second.status == JobOutcomeStatus::transit) {
	*outcome = nullptr;
	status = a_jobToData->second.status;
  }
  else {
	*outcome = nullptr;
  }
  a_jobToData.release();

  return found;
}

//Todo (Philipp): not used anymore, just left as a reference
bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::handleLocalExecutionOld(
    bool isRunOnMaster, bool& hasComputed) {

  bool result = false;
  bool needToCompute = true;
  bool needToCheck = false;

  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
      _element);

  double *luh = static_cast<double*>(cellDescription.getSolution());
  double *lduh = static_cast<double*>(cellDescription.getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double *lGradQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient());
  double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

#if defined(TaskSharing)
  //check if task outcome has been received already
  tarch::la::Vector<DIMENSIONS, double> center;
  center = cellDescription.getOffset() + 0.5 * cellDescription.getSize();

  JobTableKey key;
  for (int i=0; i < DIMENSIONS; i++)
    key.center[i] = center[i];
  key.timestamp = _predictorTimeStamp;
  key.timestepSize = _predictorTimeStepSize;
  key.element = _element;

#if DIMENSION==3
  logInfo("handleLocalExecution()",
      "team "<<exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank()
      <<" looking for job center[0] = "
      << center[0]
      <<" center[1] = "
      << center[1]
      <<" center[2] = "
      << center[2]
      <<" time stamp = "
      <<_predictorTimeStamp
      <<" time step = "
      <<_predictorTimeStepSize
      <<" enclave jobs "<<NumberOfEnclaveJobs <<" remote jobs "<<NumberOfRemoteJobs
      <<" hash = "<<(size_t) key);
#else
  logInfo("handleLocalExecution()",
        "team "<<exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank()
        <<" looking for job center[0] = "
        << center[0]
        <<" center[1] = "
        << center[1]
        <<" time stamp = "
        <<_predictorTimeStamp
        <<" time step = "
        <<_predictorTimeStepSize
        <<" enclave jobs "<<NumberOfEnclaveJobs <<" remote jobs "<<NumberOfRemoteJobs
        <<" hash = "<<(size_t) key);
#endif

//  exahype::solvers::ADERDGSolver::pollForOutstandingCommunicationRequests(&_solver,MAX_PROGRESS_ITS);
//  exahype::solvers::ADERDGSolver::progressOffloading(&_solver, false,MAX_PROGRESS_ITS);

  tbb::concurrent_hash_map<JobTableKey, JobTableEntry>::accessor a_jobToData;
  bool found = _solver._jobDatabase.find(a_jobToData, key);
  MigratablePredictionJobData *data = nullptr;
  if (found && a_jobToData->second.status == JobOutcomeStatus::received) {
    data = a_jobToData->second.data;
    assertion(data->_metadata._predictorTimeStamp == _predictorTimeStamp);
    assertion(data->_metadata._predictorTimeStepSize == _predictorTimeStepSize);

    logInfo("handleLocalExecution()",
        "team "<<exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank()
        <<" found STP in received jobs:"
        <<data->_metadata.to_string());

    if(!data->_metadata._isPotSoftErrorTriggered) {
      std::memcpy(lduh, &data->_lduh[0], data->_lduh.size() * sizeof(double));
      std::memcpy(lQhbnd, &data->_lQhbnd[0], data->_lQhbnd.size() * sizeof(double));
      std::memcpy(lFhbnd, &data->_lFhbnd[0], data->_lFhbnd.size() * sizeof(double));
#if OffloadingGradQhbnd
      std::memcpy(lGradQhbnd, &data->_lGradQhbnd[0], data->_lGradQhbnd.size() * sizeof(double));
#endif
      exahype::offloading::JobTableStatistics::getInstance().notifySavedTask();

      _solver._jobDatabase.erase(a_jobToData);
      a_jobToData.release();
      AllocatedSTPsReceive--;
      delete data;

      //NumberOfRemoteJobs--;
      cellDescription.setHasCompletedLastStep(true);
      result = false;
      needToCompute = false;
    }
    else {
      exahype::offloading::JobTableStatistics::getInstance().notifyExecutedTask();
      needToCompute = true;
#if defined(ResilienceChecks)
      needToCheck = true;
#endif
    }
  }
  else if (found && a_jobToData->second.status == JobOutcomeStatus::transit) {
    logDebug("handleLocalExecution()", "task is in transit, we may want to wait!");
#ifdef TaskSharingRescheduleIfInTransit
    a_jobToData.release();
    result = true;
    needToCompute = false;
#endif
  }
  else {
    a_jobToData.release();
    exahype::offloading::JobTableStatistics::getInstance().notifyExecutedTask();
  }
#endif

  int iterations = 0;

  if (needToCompute) {
    //TODO: add support for lGradQhbnd
#if defined(USE_ITAC)
    if(_isLocalReplica)
    VT_begin(event_stp_local_replica);
#endif
#if defined FileTrace
    auto start = std::chrono::high_resolution_clock::now();
#endif
    iterations = _solver.fusedSpaceTimePredictorVolumeIntegral(lduh,
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

    /*size_t hash_sum = 0;
    std::hash<double> hash_fn_db;
    for(int z=0;z<_solver.getBndTotalSize();z++) {
      hash_sum ^= hash_fn_db(lQhbnd[z]);
    }
    logDebug("handleLocalExecution()",
          "team "<<exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank()
          <<" cellDescription = "<<_cellDescriptionsIndex
          <<" center[0] = "
          << center[0]
          <<" center[1] = "
          << center[1]
          <<" time stamp = "
          <<_predictorTimeStamp
          <<" time step = "
          <<_predictorTimeStepSize
          <<" computed hash = "
          <<hash_sum);*/

    bool hasFlipped = exahype::offloading::ResilienceTools::getInstance().generateBitflipErrorInDoubleIfActive(lduh, _solver.getUpdateSize());
    setTrigger(hasFlipped);
#if defined(TaskSharing) && defined(ResilienceChecks)
    //needToCheck = needToCheck ||(_isPotSoftErrorTriggered==1) ? true : false;

    if(needToCheck) {
      while(!found) {
        found = _solver._jobDatabase.find(a_jobToData, key);
        if(found) {
          data = a_jobToData->second.data;
          assertion(data->_metadata._predictorTimeStamp == _predictorTimeStamp);
        }
        a_jobToData.release();
        exahype::solvers::ADERDGSolver::progressOffloading(&_solver, false, MAX_PROGRESS_ITS);
      }

      assertion(found);
      assertion(data!=nullptr);

      matchesOtherOutcome(data);
    }
#endif

#if defined(FileTrace)
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

   // exahype::offloading::STPStatsTracer::getInstance().writeTracingEventIteration(iterations, exahype::offloading::STPTraceKey::ADERDGRemoteMigratable);
    exahype::offloading::STPStatsTracer::getInstance().writeTracingEventRunIterations(duration.count(), iterations, exahype::offloading::STPTraceKey::ADERDGOwnMigratable);
    //exahype::offloading::STPStatsTracer::getInstance().writeTracingEventIteration(iterations, exahype::offloading::STPTraceKey::ADERDGOwnMigratable);
#endif

#if defined(USE_ITAC)
    if(_isLocalReplica)
    VT_end(event_stp_local_replica);
#endif
    result = false;
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
  if (hasComputed)
    exahype::offloading::JobTableStatistics::getInstance().notifyExecutedTask();
#endif

#if defined (TaskSharing)
  //check one more time
  found = _solver._jobDatabase.find(a_jobToData, key);
  if (found && a_jobToData->second.status == JobOutcomeStatus::received) {
    MigratablePredictionJobData *data = a_jobToData->second.data;
    exahype::offloading::JobTableStatistics::getInstance().notifyLateTask();

    _solver._jobDatabase.erase(a_jobToData);
    a_jobToData.release();
    AllocatedSTPsReceive--;
    delete data;
  }
  else {
    bool sendTaskOutcome = hasComputed &&
        (AllocatedSTPsSend
            <= exahype::offloading::PerformanceMonitor::getInstance().getTasksPerTimestep());
    // && exahype::offloading::MemoryMonitor::getInstance().getFreeMemMB()>10000;
#if defined(TaskSharingUseHandshake)
    if(sendTaskOutcome) {
      _solver.sendKeyOfReplicatedSTPToOtherTeams(this);
    }
#else
    if (sendTaskOutcome) {
      //  if(AllocatedSTPsSend<=1000) {
      SentSTPs++;
      _solver.sendTaskOutcomeToOtherTeams(this);
    }
#endif
  }
#ifndef OffloadingUseProgressThread
  // if(!isRunOnMaster)
  //   exahype::solvers::ADERDGSolver::progressOffloading(&_solver, isRunOnMaster);
#endif
#endif
  if(hasComputed)
    cellDescription.setHasCompletedLastStep(true);

  if (!result) {
    exahype::offloading::PerformanceMonitor::getInstance().decRemainingTasks();
    exahype::offloading::PerformanceMonitor::getInstance().decCurrentTasks();
  }

  return result;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::handleRemoteExecution( bool& hasComputed) {
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
    
#if DIMENSIONS==3
  logInfo("handleRemoteExecution",
        " processJob: center[0] = "<<_center[0]
      <<" center[1] = "<<_center[1]
      <<" center[2] = "<<_center[2]
      <<" time stamp = "<<_predictorTimeStamp
      <<" origin rank = "<<_originRank);
#else
  logInfo("handleRemoteExecution",
        " processJob: center[0] = "<<_center[0]
      <<" center[1] = "<<_center[1]
      <<" time stamp = "<<_predictorTimeStamp
      <<" origin rank = "<<_originRank);
#endif

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

#if defined(GenerateSoftErrors)
   exahype::offloading::SoftErrorInjector::getInstance().generateBitflipErrorInDoubleIfActive(_lduh, _solver.getUpdateSize());
#endif

#if DIMENSIONS==3
   logInfo("handleRemoteExecution",
        " finished job: center[0] = "<<_center[0]
      <<" center[1] = "<<_center[1]
      <<" center[2] = "<<_center[2]
      <<" time stamp = "<<_predictorTimeStamp
      <<" origin rank = "<<_originRank);
#else
   logInfo("handleRemoteExecution",
        " finished job: center[0] = "<<_center[0]
      <<" center[1] = "<<_center[1]
      <<" time stamp = "<<_predictorTimeStamp
      <<" origin rank = "<<_originRank);
#endif
   
   hasComputed = true;
#if defined(USE_ITAC)
    VT_end(event_stp_remote);
#endif
#if defined(FileTrace)
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 

   // exahype::offloading::STPStatsTracer::getInstance().writeTracingEventIteration(iterations, exahype::offloading::STPTraceKey::ADERDGRemoteMigratable);
    exahype::offloading::STPStatsTracer::getInstance().writeTracingEventRunIterations(duration.count(), iterations, exahype::offloading::STPTraceKey::ADERDGRemoteMigratable);
#endif
  return result;
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::sendBackOutcomeToOrigin() {
#if DIMENSIONS==3
  logInfo("sendBackOutcomeToOrigin",
        " send job outcome: center[0] = "<<_center[0]
      <<" center[1] = "<<_center[1]
      <<" center[2] = "<<_center[2]
      <<" time stamp = "<<_predictorTimeStamp);
#else
  logInfo("sendBackOutcomeToOrigin",
        " send job outcome: center[0] = "<<_center[0]
      <<" center[1] = "<<_center[1]
      <<" time stamp = "<<_predictorTimeStamp);
#endif
    //logInfo("handleLocalExecution()", "postSendBack");
#if defined(UseSmartMPI)
  _solver.mpiSendMigratablePredictionJobOutcomeOffload(
      _lduh,
      _lQhbnd,
      _lFhbnd,
      _lGradQhbnd,
      _originRank,
      _tag,
      exahype::offloading::OffloadingManager::getInstance().getMPICommunicatorMapped()
  );
  MigratablePredictionJob::sendBackHandler(&_solver, _tag, _originRank);
#else
  MPI_Request sendBackRequests[NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME];
  _solver.mpiIsendMigratablePredictionJobOutcome(
       _lduh,
       _lQhbnd,
       _lFhbnd,
       _lGradQhbnd,
       _originRank,
       _tag,
       exahype::offloading::OffloadingManager::getInstance().getMPICommunicatorMapped(),
       sendBackRequests);
  exahype::offloading::OffloadingManager::getInstance().submitRequests(
       sendBackRequests,
       NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME,
       _tag,
       _originRank,
       sendBackHandler,
       exahype::offloading::RequestType::sendBack,
       &_solver);
#endif
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::handleExecution(
    bool isRunOnMaster, bool& hasComputed) {
  int myRank = tarch::parallel::Node::getInstance().getRank();
  bool result = false;

#if defined(OffloadingUseProfiler)
  exahype::offloading::OffloadingProfiler::getInstance().beginComputation();
  double time = -MPI_Wtime();
#endif
  //local treatment if this job belongs to the local rank
  if (_originRank == myRank) {
    result = handleLocalExecution(isRunOnMaster, hasComputed);
    if (!_isLocalReplica) NumberOfEnclaveJobs--;
    assertion( NumberOfEnclaveJobs>=0 );
#ifndef OffloadingUseProgressThread
    if (!isRunOnMaster)
      exahype::solvers::ADERDGSolver::progressOffloading(&_solver,
          isRunOnMaster, std::numeric_limits<int>::max());
#endif
  }
  //remote task, need to execute and send it back
  else {
    result = handleRemoteExecution(hasComputed);
  }
#if defined(OffloadingUseProfiler)
  time += MPI_Wtime();
  exahype::offloading::OffloadingProfiler::getInstance().endComputation(time);
#endif

  //send back
  if (_originRank != myRank) {
    sendBackOutcomeToOrigin();
  }
  return result;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::matchesOtherOutcome(MigratablePredictionJobData *data) {
  logDebug("matchesOtherOutcome", "comparing center[0]="<<_center[0]<<" center[1]="<<_center[1]<<" timestamp "<<_predictorTimeStamp<<" with received task outcome "<<data->_metadata.to_string());

  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
      _element);

  double *luh = static_cast<double*>(cellDescription.getSolution());
  double *lduh = static_cast<double*>(cellDescription.getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double *lGradQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient());
  double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

  bool equal = true;
  bool tmp;
  tmp = exahype::offloading::ResilienceTools::getInstance().isAdmissibleNumericalError(data->_lQhbnd.data(), lQhbnd, data->_lQhbnd.size()); equal&=tmp;
  if(!tmp) {
    logError("matchesOtherOutcome", "lQhbnd is not (numerically) equal for cell "<<"center[0]="<<_center[0]<<" center[1]="<<_center[1]<<" timestamp "<<_predictorTimeStamp);
  }

#if defined(OffloadingGradQhbnd)
  tmp = exahype::offloading::ResilienceTools::getInstance().isAdmissibleNumericalError(data->_lGradQhbnd.data(), lGradQhbnd, data->_lGradQhbnd.size()); equal&=tmp;
  if(!tmp) {
    logError("matchesOtherOutcome", "lGradQhbnd is not (numerically) equal for cell "<<"center[0]="<<_center[0]<<" center[1]="<<_center[1]<<" timestamp "<<_predictorTimeStamp);
  }
#endif
  tmp = exahype::offloading::ResilienceTools::getInstance().isAdmissibleNumericalError(data->_lFhbnd.data(), lFhbnd, data->_lFhbnd.size()); equal&=tmp;
  if(!tmp) {
    logError("matchesOtherOutcome", "lFhbnd is not  (numerically) equal for cell "<<"center[0]="<<_center[0]<<" center[1]="<<_center[1]<<" timestamp "<<_predictorTimeStamp);
  }
  tmp = exahype::offloading::ResilienceTools::getInstance().isAdmissibleNumericalError(data->_lduh.data(), lduh, data->_lduh.size()); equal&=tmp;
  if(!tmp) {
    logError("matchesOtherOutcome", "lduh is not  (numerically) equal for cell "<<"center[0]="<<_center[0]<<" center[1]="<<_center[1]<<" timestamp "<<_predictorTimeStamp);
  }

  if(!equal) {
    logError("matchesOtherOutcome", "soft error detected: "<<data->_metadata.to_string());
    exahype::offloading::JobTableStatistics::getInstance().notifyDetectedError();
  }

  logDebug("matchesOtherOutcome", "checked duplicate executions for soft errors, result = "<<equal);
  exahype::offloading::JobTableStatistics::getInstance().notifyDoubleCheckedTask();

  return equal;
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveHandler(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {
  //logInfo("receiveHandler","successful receive request");

  tbb::concurrent_hash_map<std::pair<int, int>, MigratablePredictionJobData*>::accessor a_tagRankToData;
  MigratablePredictionJobData *data;
  bool found = static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToStolenData.find(
      a_tagRankToData, std::make_pair(remoteRank, tag));
  assertion(found);
  data = a_tagRankToData->second;
  a_tagRankToData.release();

  exahype::offloading::OffloadingAnalyser::getInstance().notifyReceivedSTPJob();
  MigratablePredictionJob *job =
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->createFromData(data,
          remoteRank, tag);
  peano::datatraversal::TaskSet spawnedSet(job);

  logInfo("receiveHandler",
      " received task : "<< data->_metadata.to_string());

  exahype::offloading::OffloadingProfiler::getInstance().notifyReceivedTask(
      remoteRank);
}

#if defined(TaskSharing)
void exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveKeyHandlerTaskSharing(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {
  logDebug("receiveKeyHandlerReplica","successful receive request");

  tbb::concurrent_hash_map<std::pair<int, int>, double*>::accessor a_tagRankToData;
  double *key;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToReplicaKey.find(
      a_tagRankToData, std::make_pair(remoteRank, tag));
  key = a_tagRankToData->second;
  //static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagRankToReplicaKey.erase(a_tagRankToData);
  //a_tagRankToData.release();

#if DIMENSIONS==3
  logDebug("receiveKeyHandlerReplica", "team "
      <<" received replica job key: center[0] = "<<key[0]
      <<" center[1] = "<<key[1]
      <<" center[2] = "<<key[2]
      <<" time stamp = "<<key[2*DIMENSIONS]
      <<" element = "<<(int) key[2*DIMENSIONS+2]);
#else
  logDebug("receiveKeyHandlerReplica", "team "
      <<" received replica job key: center[0] = "<<key[0]
      <<" center[1] = "<<key[1]
      <<" time stamp = "<<key[2*DIMENSIONS]
      <<" element = "<<(int) key[2*DIMENSIONS+2]);
#endif

  static_cast<exahype::solvers::ADERDGSolver*>(solver)->sendRequestForJobAndReceive(
      tag, remoteRank, key);

}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveHandlerTaskSharing(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {
  logDebug("receiveHandlerTaskSharing","successful receive request");

  tbb::concurrent_hash_map<std::pair<int, int>, MigratablePredictionJobData*>::accessor a_tagRankToData;
  MigratablePredictionJobData *data;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToReplicaData.find(
      a_tagRankToData, std::make_pair(remoteRank, tag));
  data = a_tagRankToData->second;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToReplicaData.erase(
      a_tagRankToData);
  a_tagRankToData.release();

#if defined(UseSmartMPI)
  data->_metadata.unpackContiguousBuffer();
#endif

  logInfo("receiveHandlerTaskSharing", "team "
      <<exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank()
      <<" received replica job: "
      <<data->_metadata.to_string());

  JobTableKey key; //{&data->_metadata[0], data->_metadata[2*DIMENSIONS], (int) data->_metadata[2*DIMENSIONS+2] };
  const double * center = data->_metadata.getCenter();
  for (int i = 0; i < DIMENSIONS; i++)
    key.center[i] = center[i];
  key.timestamp = data->_metadata.getPredictorTimeStamp();
  key.timestepSize = data->_metadata.getPredictorTimeStepSize();
  key.element = data->_metadata.getElement();

  //bool criticalMemoryConsumption =  exahype::offloading::MemoryMonitor::getInstance().getFreeMemMB()<1000;

  if (key.timestamp
      < static_cast<exahype::solvers::ADERDGSolver*>(solver)->getMinTimeStamp()) { // || criticalMemoryConsumption) {
    exahype::offloading::JobTableStatistics::getInstance().notifyLateTask();
    delete data;
    AllocatedSTPsReceive--;
  }
  else {
    JobTableEntry entry { data, JobOutcomeStatus::received };
    tbb::concurrent_hash_map<JobTableKey, JobTableEntry>::accessor a_jobToData;
    bool found =
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_jobDatabase.find(
            a_jobToData, key);
    if (found) {
      a_jobToData->second.status = JobOutcomeStatus::received;
      a_jobToData.release();
    }
    else {
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_jobDatabase.insert(
          std::make_pair(key, entry));
    }
    static_cast<exahype::solvers::ADERDGSolver*>(solver)->_allocatedJobs.push_back(
        key);
  }
  exahype::offloading::JobTableStatistics::getInstance().notifyReceivedTask();

}
#endif

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveBackHandler(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {

  // logInfo("receiveBackHandler","successful receiveBack request");
  tbb::concurrent_hash_map<int, CellDescription*>::accessor a_tagToCellDesc;
  bool found =
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToCellDesc.find(
          a_tagToCellDesc, tag);
  assertion(found);
  auto cellDescription = a_tagToCellDesc->second;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToCellDesc.erase(
      a_tagToCellDesc);
  a_tagToCellDesc.release();

#ifndef OffloadingLocalRecompute
  tbb::concurrent_hash_map<const CellDescription*, std::pair<int,int>>::accessor a_cellDescToTagRank;
  found = static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapCellDescToTagRank.find(a_cellDescToTagRank, cellDescription);
  assertion(found);
  // do not erase for local recompute as we need this information later on
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapCellDescToTagRank.erase(a_cellDescToTagRank);
  a_cellDescToTagRank.release();
#endif

  tbb::concurrent_hash_map<int, double>::accessor a_tagToOffloadTime;
  found =
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToOffloadTime.find(
          a_tagToOffloadTime, tag);
  double elapsed = MPI_Wtime() + a_tagToOffloadTime->second;
  assertion(found);
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToOffloadTime.erase(
      a_tagToOffloadTime);
  a_tagToOffloadTime.release();

  NumberOfRemoteJobs--;
  NumberOfEnclaveJobs--;
#ifndef OffloadingLocalRecompute
  cellDescription->setHasCompletedLastStep(true);
#else
  logInfo("receiveBackHandler", "received back STP job");
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

  logInfo("receiveBackHandler",
      " received task outcome: center[0] = "<<metadata[0]
      <<" center[1] = "<<metadata[1]
#if DIMENSIONS==3
      <<" center[2] = "<<metadata[2]
#endif
      <<" time stamp = "<<metadata[2*DIMENSIONS] <<" element = "<<(int) metadata[2*DIMENSIONS+2]);

  exahype::offloading::JobTableStatistics::getInstance().notifyReceivedTask();
  //timestamp
  if (metadata[2 * DIMENSIONS]
      < static_cast<exahype::solvers::ADERDGSolver*>(solver)->getMinTimeStamp()) {
    // if(true) {
    exahype::offloading::JobTableStatistics::getInstance().notifyLateTask();
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
      logInfo("receiveBackHandler",
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
      exahype::offloading::JobTableStatistics::getInstance().notifySavedTask();

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

      logInfo("receiveBackHandler", " applied task outcome");

      cellDescription->setHasCompletedLastStep(true);
    }
    //sb else had done it, resolve emergency
    else {
      if (LastEmergencyCell == cellDescription) {
        VetoEmergency = false;
        LastEmergencyCell = nullptr;
      }

      exahype::offloading::JobTableStatistics::getInstance().notifyLateTask();
      //delete data; //probably not safe to do here
      // AllocatedSTPsReceive--; //probably not safe to do here, race with job
    }
    lock.free();
  }
#endif
  //logInfo("receiveBackHandler", "remote execution took "<<elapsed<<" s ");

  assertion( NumberOfEnclaveJobs>=0 ); assertion( NumberOfRemoteJobs>=0 );
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::sendBackHandler(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {
  //logInfo("sendBackHandler","successful sendBack request");
  tbb::concurrent_hash_map<std::pair<int, int>, MigratablePredictionJobData*>::accessor a_tagRankToData;

  MigratablePredictionJobData *data;
  bool found =
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToStolenData.find(
          a_tagRankToData, std::make_pair(remoteRank, tag));
  assertion(found);
  data = a_tagRankToData->second;
  a_tagRankToData.release();
  delete data;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToStolenData.erase(
      std::make_pair(remoteRank, tag));

  NumberOfStolenJobs--;
  assertion( NumberOfStolenJobs>=0 );
}

#if defined(TaskSharing)
void exahype::solvers::ADERDGSolver::MigratablePredictionJob::sendAckHandlerTaskSharing(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {

  logDebug("sendAckHandlerReplication","successfully completed send ack to other team");

}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::sendKeyHandlerTaskSharing(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {

  logDebug("sendKeyHandlerReplication","successfully completed send key to other teams");
  //exahype::offloading::ReplicationStatistics::getInstance().notifySentTask();
  //tbb::concurrent_hash_map<int, double*>::accessor a_tagToData;
  //bool found = static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagToReplicationSendKey.find(a_tagToData, tag);
  //assert(found);
  //double *data = a_tagToData->second;
  //static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagToReplicationSendKey.erase(a_tagToData);

  //delete data;
  //a_tagToData.release();
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::sendHandlerTaskSharing(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {

  logDebug("sendHandlerReplication","successfully completed send to other teams");
  exahype::offloading::JobTableStatistics::getInstance().notifySentTask();
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
#endif

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::sendHandler(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {
#if defined(PerformanceAnalysisOffloadingDetailed)
  static std::atomic<int> cnt=0;
  cnt++;
#endif
  //logInfo("sendHandler","successful send request");
#if !defined(OffloadNoEarlyReceiveBacks) || defined(OffloadingLocalRecompute)
  ADERDGSolver::receiveBackMigratableJob(tag, remoteRank,
      static_cast<exahype::solvers::ADERDGSolver*>(solver));
#endif
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::packMetaData(MigratablePredictionJobMetaData *metadata) {
#if defined(UseSmartMPI)
  char *tmp = metadata->_contiguousBuffer;
  std::memcpy(tmp, _center, DIMENSIONS*sizeof(double)); tmp+= DIMENSIONS*sizeof(double);
  std::memcpy(tmp, _dx, DIMENSIONS*sizeof(double)); tmp+= DIMENSIONS*sizeof(double);
  std::memcpy(tmp, &_predictorTimeStamp, sizeof(double)); tmp+= sizeof(double);
  std::memcpy(tmp, &_predictorTimeStepSize, sizeof(double)); tmp+= sizeof(double);
  std::memcpy(tmp, &_element, sizeof(int)); tmp+= sizeof(int);
  std::memcpy(tmp, &_isPotSoftErrorTriggered, sizeof(char)); tmp+= sizeof(char);
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
  metadata->_isPotSoftErrorTriggered = _isPotSoftErrorTriggered;
#endif
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

exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::MigratablePredictionJobMetaData()
: _predictorTimeStamp(0),
  _predictorTimeStepSize(0),
  _element(0),
  _isPotSoftErrorTriggered(0),
  _contiguousBuffer(nullptr) {
#if defined(UseSmartMPI)
  _contiguousBuffer = new char[getMessageLen()];
#endif
}

exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::~MigratablePredictionJobMetaData() {
#if defined(UseSmartMPI)
  delete[] _contiguousBuffer;
#endif
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::unpackContiguousBuffer() {
  char *tmp = _contiguousBuffer;
  std::memcpy(_center, tmp, DIMENSIONS*sizeof(double)); tmp+= DIMENSIONS*sizeof(double);
  std::memcpy(_dx, tmp, DIMENSIONS*sizeof(double)); tmp+= DIMENSIONS*sizeof(double);
  std::memcpy(&_predictorTimeStamp, tmp, sizeof(double)); tmp+= sizeof(double);
  std::memcpy(&_predictorTimeStepSize, tmp,  sizeof(double)); tmp+= sizeof(double);
  std::memcpy(&_element, tmp, sizeof(int)); tmp+= sizeof(int);
  std::memcpy(&_isPotSoftErrorTriggered, tmp, sizeof(char)); tmp+= sizeof(unsigned char);
}

std::string exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::to_string() const {
  std::string result;

  result = " center[0]  = " + std::to_string(_center[0]);
  result += " center[1] = " + std::to_string(_center[1]);
  #if DIMENSIONS==3
  result += " center[1] = " + std::to_string(_center[2]);
  #endif

  result += " time stamp = " + std::to_string(_predictorTimeStamp);
  result += " element = " + std::to_string(_element);
  result += " isPotSoftErrorTriggered = " + std::to_string(_isPotSoftErrorTriggered);

  return result;
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::initDatatype() {
#if defined(UseSmartMPI)
  _datatype = MPI_BYTE;
#else
  int entries = 2+4;
  MigratablePredictionJobMetaData dummy;

  int blocklengths[] = {DIMENSIONS, DIMENSIONS, 1, 1, 1, 1};
  MPI_Datatype subtypes[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INTEGER, MPI_CHAR};
  MPI_Aint displs[] = {0, 0, 0, 0, 0, 0};

  MPI_Aint base;
  MPI_Get_address(&dummy, &base);

  MPI_Get_address(&(dummy._center), &(displs[0]));
  MPI_Get_address(&(dummy._dx), &(displs[1]));
  MPI_Get_address(&(dummy._predictorTimeStamp), &(displs[2]));
  MPI_Get_address(&(dummy._predictorTimeStepSize), &(displs[3]));
  MPI_Get_address(&(dummy._element), &(displs[4]));
  MPI_Get_address(&(dummy._isPotSoftErrorTriggered), &(displs[5]));

  for(int i=0; i<entries; i++) {
    displs[i] = displs[i]-base;
    logDebug("initDatatype" , " displ["<<i<<"]="<<displs[i]);
  }

  int ierr = MPI_Type_create_struct(entries, blocklengths, displs, subtypes, &_datatype); assert(ierr==MPI_SUCCESS);
  ierr = MPI_Type_commit(&_datatype);  assert(ierr==MPI_SUCCESS);
#endif
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::shutdownDatatype() {
#if !defined(UseSmartMPI)
  MPI_Type_free(&_datatype);
#endif
}

char* exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::getContiguousBuffer() const {
  return _contiguousBuffer;
}

MPI_Datatype exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::getMPIDatatype() {
  return _datatype;
}

size_t exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::getMessageLen() {
#if defined(UseSmartMPI)
  return (2*DIMENSIONS+2)*sizeof(double)+sizeof(int)+sizeof(unsigned char);
#else
  return 1;
#endif
}

const double * exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::getCenter() const {
  return _center;
}

const double * exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::getDx() const {
  return _dx;
}

double exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::getPredictorTimeStamp() const {
  return _predictorTimeStamp;
}

double exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::getPredictorTimeStepSize() const {
  return _predictorTimeStepSize;
}

int exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::getElement() const {
  return _element;
}

#endif

//#undef assertion
//#define assertion(expr) 
