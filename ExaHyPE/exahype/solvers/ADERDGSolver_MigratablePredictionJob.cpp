#if defined(DistributedOffloading) && defined(Parallel) && defined(SharedTBB)

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
      _isLocalReplica(false) {
  LocalStealableSTPCounter++;
  NumberOfEnclaveJobs++;
  exahype::offloading::JobTableStatistics::getInstance().notifySpawnedTask();
  exahype::offloading::PerformanceMonitor::getInstance().incCurrentTasks();
}

exahype::solvers::ADERDGSolver::MigratablePredictionJob::MigratablePredictionJob(
    ADERDGSolver& solver, const int cellDescriptionsIndex, const int element,
    const double predictorTimeStamp, const double predictorTimeStepSize,
    double *luh, double *lduh, double *lQhbnd, double *lFhbnd, double *dx,
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
      _isLocalReplica(false) {

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

//Caution: Compression is not supported yet!
bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::run(
    bool isCalledOnMaster) {

#ifdef USE_ITAC
  VT_begin(event_stp);
#endif

#if defined ScoreP
  SCOREP_USER_REGION( "exahype::solvers::ADERDGSolver::MigratablePredictionJob::run", SCOREP_USER_REGION_TYPE_FUNCTION )
#endif
  #if defined FileTrace
  auto start = std::chrono::high_resolution_clock::now();
  #endif

  bool result = false;
  bool hasComputed = false;
  int curr = std::atomic_fetch_add(&JobCounter, 1);

#if defined(TaskSharing)
 // exahype::solvers::ADERDGSolver::pollForOutstandingCommunicationRequests(&_solver, false, 10000);
  exahype::solvers::ADERDGSolver::progressOffloading(&_solver, false, 10000);
#endif

  if (curr % 1000 == 0) {
    tarch::timing::Watch watch("exahype::MigratablePredictionJob::", "-", false,
        false);
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

#if defined(GenerateNoise)
    exahype::offloading::NoiseGenerator::getInstance().generateNoiseSTP();
#endif

  #if defined FileTrace 
  auto stop = std::chrono::high_resolution_clock::now(); 
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 

  exahype::offloading::STPStatsTracer::getInstance().writeTracingEventRun(duration.count(), exahype::offloading::STPType::ADERDGOwnMigratable);
//  std::stringstream stream;
//  stream<<exahype::parser::Parser::getSTPTracingOutputDirName()<<"/exahype_solvers_ADERDGSolver_MigratablePredictionJob_run_rank_";
//  int rank = tarch::parallel::Node::getInstance().getRank();
//  stream<<rank<<"_";
//  //this will only work for 2 cores per Rank
//#if defined (SharedTBB)
//  int threadId=tarch::multicore::Core::getInstance().getThreadNum();
//  stream<<threadId;
//#endif
//  stream<<".txt";
//  std::string path=stream.str();
//
//  std::ofstream file;
//  file.open(path,std::fstream::app);
//  file << duration.count() << std::endl;
//  file.close();
  #endif

#ifdef USE_ITAC
  VT_end(event_stp);
#endif
  return result;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::handleLocalExecution(
    bool isRunOnMaster, bool& hasComputed) {
  int myRank = tarch::parallel::Node::getInstance().getRank();
  bool result = false;
  bool needToCompute = true;

  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
      _element);

  double *luh = static_cast<double*>(cellDescription.getSolution());
  double *lduh = static_cast<double*>(cellDescription.getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

#if defined(TaskSharing)
  //check if task outcome has been received already
  tarch::la::Vector<DIMENSIONS, double> center;
  center = cellDescription.getOffset() + 0.5 * cellDescription.getSize();

  JobTableKey key;
  for (int i; i < DIMENSIONS; i++)
    key.center[i] = center[i];
  key.timestamp = _predictorTimeStamp;
  key.element = _element;

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
      <<" enclave jobs "<<NumberOfEnclaveJobs <<" remote jobs "<<NumberOfRemoteJobs
      <<" hash = "<<(size_t) key);

//  exahype::solvers::ADERDGSolver::pollForOutstandingCommunicationRequests(&_solver,10000);
//  exahype::solvers::ADERDGSolver::progressOffloading(&_solver, false,10000);

  tbb::concurrent_hash_map<JobTableKey, JobTableEntry>::accessor a_jobToData;
  bool found = _solver._jobDatabase.find(a_jobToData, key);
  if (found && a_jobToData->second.status == JobOutcomeStatus::received) {
    MigratablePredictionJobData *data = a_jobToData->second.data;
    assert(data->_metadata[2 * DIMENSIONS] == _predictorTimeStamp);
    logInfo("handleLocalExecution()",
        "team "<<exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank()
        <<" found STP in received jobs:"
        <<" center[0] = "<<data->_metadata[0]
        <<" center[1] = "<<data->_metadata[1]
        <<" center[2] = "<<data->_metadata[2]
        <<" time stamp = "<<data->_metadata[2*DIMENSIONS]
        <<" element = "<<(int) data->_metadata[2*DIMENSIONS+2]);
    std::memcpy(luh, &data->_luh[0], data->_luh.size() * sizeof(double));
    std::memcpy(lduh, &data->_lduh[0], data->_lduh.size() * sizeof(double));
    std::memcpy(lQhbnd, &data->_lQhbnd[0], data->_lQhbnd.size() * sizeof(double));
    std::memcpy(lFhbnd, &data->_lFhbnd[0], data->_lFhbnd.size() * sizeof(double));
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
    iterations = _solver.fusedSpaceTimePredictorVolumeIntegral(lduh, lQhbnd, nullptr, lFhbnd,
        luh, cellDescription.getOffset() + 0.5 * cellDescription.getSize(),
        cellDescription.getSize(), _predictorTimeStamp, _predictorTimeStepSize,
        true);
    hasComputed = true;

    cellDescription.setHasCompletedLastStep(true);
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

#if defined(FileTrace)
  if(hasComputed) {
    exahype::offloading::STPStatsTracer::getInstance().writeTracingEventIteration(iterations, exahype::offloading::STPType::ADERDGOwnMigratable);
  }
#endif
  //#if defined(FileTrace)
//    std::stringstream stream;
//    stream.str(std::string());
//    stream<<"./TraceOutput/exahype_solvers_ADERDGSolver_OwnMigratableJob_iterations_rank_";
//    int rank=tarch::parallel::Node::getInstance().getRank();
//    stream<<rank<<"_";
//#if defined(SharedTBB)
//    int threadId=tarch::multicore::Core::getInstance().getThreadNum();
//    stream<<threadId;
//#endif
//    stream<<".txt";
//    std::string path=stream.str();
//
//    std::ofstream file;
//    file.open(path,std::fstream::app);
//    file << (iterations+1) << std::endl;
//    file.close();
//#endif


#if defined (TaskSharing)
  //check one more time
  //tbb::concurrent_hash_map<JobTableKey, StealablePredictionJobData*>::accessor a_jobToData;
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
    bool sendTaskOutcome =
        AllocatedSTPsSend
            <= exahype::offloading::PerformanceMonitor::getInstance().getTasksPerTimestep();
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

  if (!result) {
    exahype::offloading::PerformanceMonitor::getInstance().decRemainingTasks();
    exahype::offloading::PerformanceMonitor::getInstance().decCurrentTasks();
  }

  return result;
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
#if defined(USE_ITAC)
    VT_begin(event_stp_remote);
#endif
    result = false;
    //TODO: support for lGradQhbnd
    int iterations=_solver.fusedSpaceTimePredictorVolumeIntegral(
      _lduh,_lQhbnd,nullptr,_lFhbnd,
      _luh,
      _center,
      _dx,
      _predictorTimeStamp,
      _predictorTimeStepSize,
      true);
    hasComputed = true;
#if defined(USE_ITAC)
    VT_end(event_stp_remote);
#endif
#if defined(FileTrace)
    exahype::offloading::STPStatsTracer::getInstance().writeTracingEventIteration(iterations, exahype::offloading::STPType::ADERDGRemoteMigratable);
#endif
//    std::stringstream stream;
//    stream.str(std::string());
//    stream<<"./TraceOutput/exahype_solvers_ADERDGSolver_AlienMigratableJob_iterations_rank_";
//    int rank=tarch::parallel::Node::getInstance().getRank();
//    stream<<rank<<"_";
//#if defined(SharedTBB)
//    int threadId=tarch::multicore::Core::getInstance().getThreadNum();
//    stream<<threadId;
//#endif
//    stream<<".txt";
//    std::string path=stream.str();
//
//    std::ofstream file;
//    file.open(path,std::fstream::app);
//    file << (iterations+1) << std::endl;
//    file.close();
//#endif
  }
#if defined(OffloadingUseProfiler)
  time += MPI_Wtime();
  exahype::offloading::OffloadingProfiler::getInstance().endComputation(time);
#endif

  //send back
  if (_originRank != myRank) {
    MPI_Request sendBackRequests[4];
    logInfo("handleExecution",
        " send job outcome: center[0] = "<<_center[0]
      <<" center[1] = "<<_center[1]
      <<" center[2] = "<<_center[2]
      <<" time stamp = "<<_predictorTimeStamp);
    //logInfo("handleLocalExecution()", "postSendBack");
    _solver.isendMigratablePredictionJob(_luh, _lduh, _lQhbnd, _lFhbnd,
      _originRank,
      _tag,
      exahype::offloading::OffloadingManager::getInstance().getMPICommunicatorMapped(),
      sendBackRequests);
    exahype::offloading::OffloadingManager::getInstance().submitRequests(
      sendBackRequests,
      4,
      _tag,
      _originRank,
      sendBackHandler,
      exahype::offloading::RequestType::sendBack,
      &_solver);
  }
  return result;
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveHandler(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {
  //logInfo("receiveHandler","successful receive request");

  tbb::concurrent_hash_map<std::pair<int, int>, MigratablePredictionJobData*>::accessor a_tagRankToData;
  MigratablePredictionJobData *data;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToStolenData.find(
      a_tagRankToData, std::make_pair(remoteRank, tag));
  data = a_tagRankToData->second;
  a_tagRankToData.release();

  exahype::offloading::OffloadingAnalyser::getInstance().notifyReceivedSTPJob();
  MigratablePredictionJob *job =
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->createFromData(data,
          remoteRank, tag);
  peano::datatraversal::TaskSet spawnedSet(job);

  logInfo("receiveHandler",
      " received task : center[0] = "<<data->_metadata[0]
      <<" center[1] = "<<data->_metadata[1]
      <<" center[2] = "<<data->_metadata[2]
      <<" time stamp = "<<data->_metadata[2*DIMENSIONS]
      <<" element = "<<(int) data->_metadata[2*DIMENSIONS+2]);

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

  logDebug("receiveKeyHandlerReplica", "team "
      <<" received replica job key: center[0] = "<<key[0]
      <<" center[1] = "<<key[1]
      <<" center[2] = "<<key[2]
      <<" time stamp = "<<key[2*DIMENSIONS]
      <<" element = "<<(int) key[2*DIMENSIONS+2]);

  static_cast<exahype::solvers::ADERDGSolver*>(solver)->sendRequestForJobAndReceive(
      tag, remoteRank, key);

}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveHandlerTaskSharing(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {
  logDebug("receiveHandlerReplica","successful receive request");

  tbb::concurrent_hash_map<std::pair<int, int>, MigratablePredictionJobData*>::accessor a_tagRankToData;
  MigratablePredictionJobData *data;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToReplicaData.find(
      a_tagRankToData, std::make_pair(remoteRank, tag));
  data = a_tagRankToData->second;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToReplicaData.erase(
      a_tagRankToData);
  a_tagRankToData.release();

  logDebug("receiveHandlerReplica", "team "
      <<exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank()
      <<" received replica job: center[0] = "<<data->_metadata[0]
      <<" center[1] = "<<data->_metadata[1]
      <<" center[2] = "<<data->_metadata[2]
      <<" time stamp = "<<data->_metadata[2*DIMENSIONS]
      <<" element = "<<(int) data->_metadata[2*DIMENSIONS+2]);

  JobTableKey key; //{&data->_metadata[0], data->_metadata[2*DIMENSIONS], (int) data->_metadata[2*DIMENSIONS+2] };
  for (int i = 0; i < DIMENSIONS; i++)
    key.center[i] = data->_metadata[i];
  key.timestamp = data->_metadata[2 * DIMENSIONS];
  key.element = data->_metadata[2 * DIMENSIONS + 2];

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
    static_cast<exahype::solvers::ADERDGSolver*>(solver)->_allocatedJobs.push(
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
  found =
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToSTPData.find(
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
      <<" center[2] = "<<metadata[2]
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
      double *luh = static_cast<double*>(cellDescription->getSolution());
      double *lduh = static_cast<double*>(cellDescription->getUpdate());
      double *lQhbnd =
          static_cast<double*>(cellDescription->getExtrapolatedPredictor());
      double *lFhbnd = static_cast<double*>(cellDescription->getFluctuation());

      /*if(static_cast<MigratablePredictionJob*>(recompJob)->_predictorTimeStamp!=metadata[2*DIMENSIONS]) {
       logInfo("receiveBackHandler","job timestamp "<<static_cast<MigratablePredictionJob*>(recompJob)->_predictorTimeStamp
       <<" metadata[2*DIMENSIONS] "<< metadata[2*DIMENSIONS]);
       }*/

      assertion(
          static_cast<MigratablePredictionJob*>(recompJob)->_predictorTimeStamp
              == metadata[2 * DIMENSIONS]);

      std::memcpy(luh, &data->_luh[0], data->_luh.size() * sizeof(double));
      std::memcpy(lduh, &data->_lduh[0], data->_lduh.size() * sizeof(double));
      std::memcpy(lQhbnd, &data->_lQhbnd[0],
          data->_lQhbnd.size() * sizeof(double));
      std::memcpy(lFhbnd, &data->_lFhbnd[0],
          data->_lFhbnd.size() * sizeof(double));
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

exahype::solvers::ADERDGSolver::MigratablePredictionJobData::MigratablePredictionJobData(
    ADERDGSolver& solver) :
        _luh(solver.getDataPerCell()), _lduh(solver.getUpdateSize()),
        _lQhbnd(solver.getBndTotalSize()), _lFhbnd(solver.getBndFluxTotalSize()) {
  AllocatedSTPs++;
}

exahype::solvers::ADERDGSolver::MigratablePredictionJobData::~MigratablePredictionJobData() {
  AllocatedSTPs--;
}

#endif
