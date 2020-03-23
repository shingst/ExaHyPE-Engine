#if defined(DistributedOffloading) && defined(Parallel) && defined(SharedTBB)

#include "ADERDGSolver.h"

#include "exahype/offloading/PerformanceMonitor.h"
#include "exahype/offloading/OffloadingAnalyser.h"
#include "exahype/offloading/OffloadingProfiler.h"
#include "exahype/offloading/ReplicationStatistics.h"
#include "exahype/offloading/MemoryMonitor.h"

exahype::solvers::ADERDGSolver::MigratablePredictionJob::MigratablePredictionJob(
    ADERDGSolver& solver,
    const int cellDescriptionsIndex,
    const int element,
    const double predictorTimeStamp,
    const double predictorTimeStepSize) :
       tarch::multicore::jobs::Job(tarch::multicore::jobs::JobType::BackgroundTask, 0, getTaskPriorityLocalStealableJob(cellDescriptionsIndex, element, predictorTimeStamp)),  //this is a locally executed job
        _solver(solver),
        _cellDescriptionsIndex(cellDescriptionsIndex),
        _element(element),
        _predictorTimeStamp(predictorTimeStamp),
        _predictorTimeStepSize(predictorTimeStepSize),
        _originRank(tarch::parallel::Node::getInstance().getRank()),
        _tag(-1),
        _luh(nullptr),_lduh(nullptr),_lQhbnd(nullptr), _lFhbnd(nullptr)
{
  LocalStealableSTPCounter++;
  NumberOfEnclaveJobs++;
  exahype::offloading::ReplicationStatistics::getInstance().notifySpawnedTask();
  exahype::offloading::PerformanceMonitor::getInstance().incCurrentTasks();
};

exahype::solvers::ADERDGSolver::MigratablePredictionJob::MigratablePredictionJob(
    ADERDGSolver& solver,
    const int cellDescriptionsIndex,
    const int element,
    const double predictorTimeStamp,
    const double predictorTimeStepSize,
    double *luh,
    double *lduh,
    double *lQhbnd,
    double *lFhbnd,
    double *dx,
    double *center,
    const int originRank,
    const int tag) :
       tarch::multicore::jobs::Job(tarch::multicore::jobs::JobType::BackgroundTask, 0, tarch::multicore::DefaultPriority*4), //this is a remotely executed job -> high priority
        _solver(solver),
        _cellDescriptionsIndex(cellDescriptionsIndex),
        _element(element),
        _predictorTimeStamp(predictorTimeStamp),
        _predictorTimeStepSize(predictorTimeStepSize),
        _originRank(originRank),
        _tag(tag),
        _luh(luh),_lduh(lduh),_lQhbnd(lQhbnd), _lFhbnd(lFhbnd) {

  for(int i=0; i<DIMENSIONS; i++) {
    _center[i] = center[i];
    _dx[i] = dx[i];
  }

  if(_originRank!=tarch::parallel::Node::getInstance().getRank()) {
    NumberOfStolenJobs++;
  }
  else
    NumberOfEnclaveJobs++;
  exahype::offloading::PerformanceMonitor::getInstance().incCurrentTasks();
};

exahype::solvers::ADERDGSolver::MigratablePredictionJob::~MigratablePredictionJob() {};

//Caution: Compression is not supported yet!
bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::run( bool isCalledOnMaster ) {

#ifdef USE_ITAC
      VT_begin(event_stp);
#endif
     bool result = false;
     int curr = std::atomic_fetch_add(&JobCounter, 1);

     if(curr%1000==0) {
       tarch::timing::Watch watch("exahype::StealablePredictionJob::", "-", false,false);
       watch.startTimer();
       result = handleExecution(isCalledOnMaster);
       watch.stopTimer();

       logDebug("run()","measured time per STP "<<watch.getCalendarTime());

       exahype::offloading::OffloadingAnalyser::getInstance().setTimePerSTP(watch.getCalendarTime());
     }
     else
       result = handleExecution(isCalledOnMaster);

#ifdef USE_ITAC
      VT_end(event_stp);
#endif
  return result;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::handleLocalExecution(bool isRunOnMaster) {
  int myRank = tarch::parallel::Node::getInstance().getRank();
  bool result = false;

  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,_element);

  double *luh    = static_cast<double*>(cellDescription.getSolution());
  double *lduh   = static_cast<double*>(cellDescription.getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

#if defined(TaskSharing) || defined(OffloadingLocalRecompute)
  //check if task outcome has been received already
  tarch::la::Vector<DIMENSIONS, double> center;
  center = cellDescription.getOffset()+0.5*cellDescription.getSize();

  JobTableKey  key;
  for(int i; i<DIMENSIONS; i++)
     key.center[i] = center[i];
  key.timestamp = _predictorTimeStamp;
  key.element = _element;

  logInfo("handleLocalExecution()", "team "<<exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank()
                                   <<" looking for job center[0] = "<< center[0]
                                   <<" center[1] = "<< center[1]
                                   <<" center[2] = "<< center[2]
                                   <<" time stamp = "<<_predictorTimeStamp
                                   <<" hash = "<<(size_t) key);

  //exahype::solvers::ADERDGSolver::progressOffloading(&_solver, false);
  //exahype::solvers::ADERDGSolver::pollForOutstandingCommunicationRequests(&_solver);

  tbb::concurrent_hash_map<JobTableKey, JobTableEntry>::accessor a_jobToData;
  bool found = _solver._jobDatabase.find(a_jobToData, key);
  if(found && a_jobToData->second.status == ReplicationStatus::received ) {
    MigratablePredictionJobData *data = a_jobToData->second.data;
    assert(data->_metadata[2*DIMENSIONS]==_predictorTimeStamp);
    logDebug("handleLocalExecution()", "team "<<exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank()
                                       <<" found STP in received jobs:"
                                       <<" center[0] = "<<data->_metadata[0]
                                       <<" center[1] = "<<data->_metadata[1]
                                       <<" center[2] = "<<data->_metadata[2]
                                       <<" time stamp = "<<data->_metadata[2*DIMENSIONS]
                                       <<" element = "<<(int) data->_metadata[2*DIMENSIONS+2]);
    std::memcpy(luh, &data->_luh[0], data->_luh.size()*sizeof(double));
    std::memcpy(lduh, &data->_lduh[0], data->_lduh.size()*sizeof(double));
    std::memcpy(lQhbnd, &data->_lQhbnd[0], data->_lQhbnd.size()*sizeof(double));
    std::memcpy(lFhbnd, &data->_lFhbnd[0], data->_lFhbnd.size()*sizeof(double));
    exahype::offloading::ReplicationStatistics::getInstance().notifySavedTask();

    _solver._jobDatabase.erase(a_jobToData);
    a_jobToData.release();
    AllocatedSTPsReceive--;
    delete data;

  cellDescription.setHasCompletedLastStep(true);
  result = false;
  }
  else if (found && a_jobToData->second.status == ReplicationStatus::transit) {
    logDebug("handleLocalExecution()", "task is in transit, we may want to wait!");
#ifdef TaskSharingRescheduleIfInTransit
    a_jobToData.release();
    result = true;
#endif
  }
  else {
    a_jobToData.release();
    logInfo("handleLocalExecution()",   "team "<<exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank()
                                          <<" Data not available, gotta do it on my own!"
                                          <<" center[0] = "<<center[0]
                                          <<" center[1] = "<<center[1]
                                          <<" center[2] = "<<center[2]
                                          <<" time stamp = "<<_predictorTimeStamp);
    exahype::offloading::ReplicationStatistics::getInstance().notifyExecutedTask();
  }
#endif

    //TODO: add support for lGradQhbnd
    _solver.fusedSpaceTimePredictorVolumeIntegral(
      lduh,lQhbnd,nullptr, lFhbnd,
      luh,
      cellDescription.getOffset()+0.5*cellDescription.getSize(),
      cellDescription.getSize(),
      _predictorTimeStamp,
      _predictorTimeStepSize,
      true);
    cellDescription.setHasCompletedLastStep(true);

    result = false;

#if defined (TaskSharing)
    //check one more time
    //tbb::concurrent_hash_map<JobTableKey, StealablePredictionJobData*>::accessor a_jobToData;
    found = _solver._jobDatabase.find(a_jobToData, key);
    if(found && a_jobToData->second.status==ReplicationStatus::received) {
      MigratablePredictionJobData *data = a_jobToData->second.data;
      exahype::offloading::ReplicationStatistics::getInstance().notifyLateTask();

      _solver._jobDatabase.erase(a_jobToData);
      a_jobToData.release();
      AllocatedSTPsReceive--;
      delete data;
    }
    else {
      bool sendTaskOutcome = AllocatedSTPsSend<=exahype::offloading::PerformanceMonitor::getInstance().getTasksPerTimestep();
                          // && exahype::offloading::MemoryMonitor::getInstance().getFreeMemMB()>10000;
#if defined(TaskSharingUseHandshake)
      if(sendTaskOutcome) {
        _solver.sendKeyOfReplicatedSTPToOtherTeams(this);
      }
#else
      if(sendTaskOutcome) {
  //  if(AllocatedSTPsSend<=1000) {
        SentSTPs++;
        _solver.sendFullReplicatedSTPToOtherTeams(this);
      }
#endif
    }
#ifndef OffloadingUseProgressThread
   // if(!isRunOnMaster)
   //   exahype::solvers::ADERDGSolver::progressOffloading(&_solver, isRunOnMaster);
#endif
#endif

  if(!result) {
    exahype::offloading::PerformanceMonitor::getInstance().decRemainingTasks();
    exahype::offloading::PerformanceMonitor::getInstance().decCurrentTasks();
  }

  return result;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::handleExecution(bool isRunOnMaster) {
  int myRank = tarch::parallel::Node::getInstance().getRank();
  bool result = false;


#if defined(OffloadingUseProfiler)
  exahype::offloading::OffloadingProfiler::getInstance().beginComputation();
  double time = -MPI_Wtime();
#endif
  //local treatment if this job belongs to the local rank
  if(_originRank==myRank) {
    result = handleLocalExecution(isRunOnMaster);
    NumberOfEnclaveJobs--;
    assertion( NumberOfEnclaveJobs>=0 );
#ifndef OffloadingUseProgressThread
    if(!isRunOnMaster)
      exahype::solvers::ADERDGSolver::progressOffloading(&_solver, isRunOnMaster, std::numeric_limits<int>::max());
#endif
  }
  //remote task, need to execute and send it back
  else {
  result = false;
    //TODO: support for lGradQhbnd
    _solver.fusedSpaceTimePredictorVolumeIntegral(
      _lduh,_lQhbnd,nullptr,_lFhbnd,
      _luh,
      _center,
      _dx,
      _predictorTimeStamp,
      _predictorTimeStepSize,
      true);
  }
#if defined(OffloadingUseProfiler)
  time += MPI_Wtime();
  exahype::offloading::OffloadingProfiler::getInstance().endComputation(time);
#endif

  //send back
  if(_originRank!=myRank) {
#if defined(OffloadingLocalRecompute)
    //send away task outcome for receiving it into job table
    OffloadEntry entry = {
       -1,
       _cellDescriptionsIndex,
       _element,
       _predictorTimeStamp,
       _predictorTimeStepSize};

    tbb::concurrent_hash_map<std::pair<int, int>, MigratablePredictionJobData*>::accessor a_tagRankToData;
    MigratablePredictionJobData *data;
    _solver._mapTagRankToStolenData.find(a_tagRankToData, std::make_pair(_originRank, _tag));
    data = a_tagRankToData->second;
    a_tagRankToData.release();

    MPI_Request sendBackRequests[5];

    //_solver.packMetadataToBuffer(entry, data->_metadata);

    logInfo("handleExecution", " send job outcome: center[0] = "<<data->_metadata[0]
                                               <<" center[1] = "<<data->_metadata[1]
                                               <<" center[2] = "<<data->_metadata[2]
                                               <<" time stamp = "<<_predictorTimeStamp);

   /* _solver.isendStealablePredictionJob(&_luh[0],
                                 &_lduh[0],
                                 &_lQhbnd[0],
                                 &_lFhbnd[0],
                                 _originRank,
                                 _tag,
                                 exahype::offloading::OffloadingManager::getInstance().getMPICommunicatorMapped(),
                                 &sendBackRequests[0],
                                 &data->_metadata[0]);

    exahype::offloading::OffloadingManager::getInstance().submitRequests(sendBackRequests,
                                                                         5,
                                                                         _tag,
                                                                         _originRank,
                                                                         sendBackHandler,
                                                                         exahype::offloading::RequestType::sendBack,
                                                                         &_solver,
                                                                         false); */

#else
    MPI_Request sendBackRequests[4];
    //logInfo("handleLocalExecution()", "postSendBack");
    _solver.isendStealablePredictionJob(_luh,
                                        _lduh,
                                        _lQhbnd,
                                        _lFhbnd,
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
#endif
  }
  return result;
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveHandler(exahype::solvers::Solver* solver, int tag, int remoteRank) {
  //logInfo("receiveHandler","successful receive request");

  tbb::concurrent_hash_map<std::pair<int, int>, MigratablePredictionJobData*>::accessor a_tagRankToData;
  MigratablePredictionJobData *data;
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagRankToStolenData.find(a_tagRankToData, std::make_pair(remoteRank, tag));
  data = a_tagRankToData->second;
  a_tagRankToData.release();

  exahype::offloading::OffloadingAnalyser::getInstance().notifyReceivedSTPJob();
  MigratablePredictionJob *job= static_cast<exahype::solvers::ADERDGSolver*> (solver)->createFromData(data, remoteRank, tag);
  peano::datatraversal::TaskSet spawnedSet(job);

  exahype::offloading::OffloadingProfiler::getInstance().notifyReceivedTask(remoteRank);
}

#if defined(TaskSharing)
void exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveKeyHandlerReplication(exahype::solvers::Solver* solver, int tag, int remoteRank) {
  logDebug("receiveKeyHandlerReplica","successful receive request");

  tbb::concurrent_hash_map<std::pair<int, int>, double*>::accessor a_tagRankToData;
  double *key;
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagRankToReplicaKey.find(a_tagRankToData, std::make_pair(remoteRank, tag));
  key = a_tagRankToData->second;
  //static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagRankToReplicaKey.erase(a_tagRankToData);
  //a_tagRankToData.release();

  logDebug("receiveKeyHandlerReplica", "team "
                           <<" received replica job key: center[0] = "<<key[0]
                         <<" center[1] = "<<key[1]
               <<" center[2] = "<<key[2]
               <<" time stamp = "<<key[2*DIMENSIONS]
               <<" element = "<<(int) key[2*DIMENSIONS+2]);

  static_cast<exahype::solvers::ADERDGSolver*> (solver)->sendRequestForJobAndReceive(tag, remoteRank, key);

}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveHandlerReplication(exahype::solvers::Solver* solver, int tag, int remoteRank) {
  logDebug("receiveHandlerReplica","successful receive request");

  tbb::concurrent_hash_map<std::pair<int, int>, MigratablePredictionJobData*>::accessor a_tagRankToData;
  MigratablePredictionJobData *data;
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagRankToReplicaData.find(a_tagRankToData, std::make_pair(remoteRank, tag));
  data = a_tagRankToData->second;
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagRankToReplicaData.erase(a_tagRankToData);
  a_tagRankToData.release();

  logDebug("receiveHandlerReplica", "team "
                                    <<exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank()
                        <<" received replica job: center[0] = "<<data->_metadata[0]
            <<" center[1] = "<<data->_metadata[1]
            <<" center[2] = "<<data->_metadata[2]
            <<" time stamp = "<<data->_metadata[2*DIMENSIONS]
            <<" element = "<<(int) data->_metadata[2*DIMENSIONS+2]);

  JobTableKey key; //{&data->_metadata[0], data->_metadata[2*DIMENSIONS], (int) data->_metadata[2*DIMENSIONS+2] };
  for(int i=0; i<DIMENSIONS; i++)
    key.center[i] = data->_metadata[i];
  key.timestamp = data->_metadata[2*DIMENSIONS];
  key.element = data->_metadata[2*DIMENSIONS+2];

  //bool criticalMemoryConsumption =  exahype::offloading::MemoryMonitor::getInstance().getFreeMemMB()<1000;

  if(key.timestamp<static_cast<exahype::solvers::ADERDGSolver*> (solver)->getMinTimeStamp()) {// || criticalMemoryConsumption) {
    exahype::offloading::ReplicationStatistics::getInstance().notifyLateTask();
    delete data;
    AllocatedSTPsReceive--;
  }
  else {
    JobTableEntry entry {data, ReplicationStatus::received};
    tbb::concurrent_hash_map<JobTableKey, JobTableEntry>::accessor a_jobToData;
    bool found = static_cast<exahype::solvers::ADERDGSolver*> (solver)->_jobDatabase.find(a_jobToData, key);
    if (found) {
      a_jobToData->second.status = ReplicationStatus::received;
      a_jobToData.release();
    }
    else{
      static_cast<exahype::solvers::ADERDGSolver*> (solver)->_jobDatabase.insert(std::make_pair(key,entry));
    }
    static_cast<exahype::solvers::ADERDGSolver*> (solver)->_allocatedJobs.push(key);
  }
  exahype::offloading::ReplicationStatistics::getInstance().notifyReceivedTask();

}
#endif

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveBackHandler(exahype::solvers::Solver* solver, int tag, int remoteRank) {

 // logInfo("receiveBackHandler","successful receiveBack request");
  tbb::concurrent_hash_map<int, CellDescription*>::accessor a_tagToCellDesc;
  bool found = static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagToCellDesc.find(a_tagToCellDesc, tag);
  assertion(found);
  auto cellDescription = a_tagToCellDesc->second;
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagToCellDesc.erase(a_tagToCellDesc);
  a_tagToCellDesc.release();

  tbb::concurrent_hash_map<const CellDescription*, std::pair<int,int>>::accessor a_cellDescToTagRank;
  found =  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapCellDescToTagRank.find(a_cellDescToTagRank, cellDescription);
  assertion(found);
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapCellDescToTagRank.erase(a_cellDescToTagRank);
  a_cellDescToTagRank.release();

#ifndef OffloadingLocalRecompute
  cellDescription->setHasCompletedLastStep(true);
  NumberOfEnclaveJobs--;
  NumberOfRemoteJobs--;
#else
  // Todo: insert into job table!
#endif

  tbb::concurrent_hash_map<int, double>::accessor a_tagToOffloadTime;
  found = static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagToOffloadTime.find(a_tagToOffloadTime, tag);
  double elapsed = MPI_Wtime() + a_tagToOffloadTime->second;
  assertion(found);
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagToOffloadTime.erase(a_tagToOffloadTime);
  a_tagToOffloadTime.release();
 //logInfo("receiveBackHandler", "remote execution took "<<elapsed<<" s ");

  assertion( NumberOfEnclaveJobs>=0 );
  assertion( NumberOfRemoteJobs>=0 );
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::sendBackHandler(exahype::solvers::Solver* solver, int tag, int remoteRank) {
  //logInfo("sendBackHandler","successful sendBack request");
  tbb::concurrent_hash_map<std::pair<int, int>, MigratablePredictionJobData*>::accessor a_tagRankToData;

  MigratablePredictionJobData *data;
  bool found = static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagRankToStolenData.find(a_tagRankToData, std::make_pair(remoteRank, tag));
  assertion(found);
  data = a_tagRankToData->second;
  a_tagRankToData.release();
  delete data;
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagRankToStolenData.erase(std::make_pair(remoteRank, tag));

  NumberOfStolenJobs--;
  assertion( NumberOfStolenJobs>=0 );
}

#if defined(TaskSharing)
void exahype::solvers::ADERDGSolver::MigratablePredictionJob::sendAckHandlerReplication(
    exahype::solvers::Solver* solver,
    int tag,
    int remoteRank)
{

  logDebug("sendAckHandlerReplication","successfully completed send ack to other team");

}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::sendKeyHandlerReplication(
    exahype::solvers::Solver* solver,
    int tag,
    int remoteRank)
{

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

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::sendHandlerReplication(
    exahype::solvers::Solver* solver,
    int tag,
    int remoteRank)
{

  logDebug("sendHandlerReplication","successfully completed send to other teams");
  exahype::offloading::ReplicationStatistics::getInstance().notifySentTask();
  tbb::concurrent_hash_map<int, MigratablePredictionJobData*>::accessor a_tagToData;
  bool found = static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagToReplicationSendData.find(a_tagToData, tag);
  assertion(found);
  MigratablePredictionJobData *data = a_tagToData->second;
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagToReplicationSendData.erase(a_tagToData);
  delete data;
  AllocatedSTPsSend--;
  CompletedSentSTPs++;
  logDebug("sendHandlerReplication"," allocated stps send "<<AllocatedSTPsSend);
  a_tagToData.release();
}
#endif

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::sendHandler(exahype::solvers::Solver* solver, int tag, int remoteRank) {
#if defined(PerformanceAnalysisOffloadingDetailed)
  static std::atomic<int> cnt=0;
  cnt++;
#endif
  //logInfo("sendHandler","successful send request");
  
#if defined(OffloadingLocalRecompute)
  tbb::concurrent_hash_map<int, MigratablePredictionJobData*>::accessor a_tagToSTPData;
  MigratablePredictionJobData *data;
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagToSTPData.find(a_tagToSTPData, tag);
  assertion(found);
  data = a_tagToSTPData->second;
  a_tagToSTPData.release();
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagToSTPData.erase(tag);
  delete data;

#else
  tbb::concurrent_hash_map<int, double*>::accessor a_tagToMeta;
  double *metaData;
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagToMetaData.find(a_tagToMeta, tag);
  metaData = a_tagToMeta->second;
  a_tagToMeta.release();
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagToMetaData.erase(tag);
  delete[] metaData;

#if !defined(OffloadingNoEarlyReceiveBacks)
  logDebug("sendHandler","posting Irecv back early");
  
  MPI_Comm commMapped = exahype::offloading::OffloadingManager::getInstance().getMPICommunicatorMapped();
  tbb::concurrent_hash_map<int, CellDescription*>::accessor a_tagToCellDesc;
  bool found = static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagToCellDesc.find(a_tagToCellDesc, tag);
  assert(found);
  auto cellDescription = a_tagToCellDesc->second;
  a_tagToCellDesc.release();
  double *luh    = static_cast<double*>(cellDescription->getSolution());
  double *lduh   = static_cast<double*>(cellDescription->getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription->getExtrapolatedPredictor());
  double *lFhbnd = static_cast<double*>(cellDescription->getFluctuation());

  MPI_Request recvRequests[4];
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->irecvStealablePredictionJob(
      luh, lduh, lQhbnd,
      lFhbnd, remoteRank, tag, commMapped, recvRequests);
  exahype::offloading::OffloadingManager::getInstance().submitRequests(
      recvRequests, 4, tag, remoteRank,
      exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveBackHandler,
      exahype::offloading::RequestType::receiveBack, solver, false);
#endif  

#endif

}

exahype::solvers::ADERDGSolver::MigratablePredictionJobData::MigratablePredictionJobData( ADERDGSolver& solver ) :
  _luh(solver.getDataPerCell()),
  _lduh(solver.getUpdateSize()),
  _lQhbnd(solver.getBndTotalSize()),
  _lFhbnd(solver.getBndFluxTotalSize())
{
  AllocatedSTPs++;
};

exahype::solvers::ADERDGSolver::MigratablePredictionJobData::~MigratablePredictionJobData() {
  AllocatedSTPs--;
};

#endif
