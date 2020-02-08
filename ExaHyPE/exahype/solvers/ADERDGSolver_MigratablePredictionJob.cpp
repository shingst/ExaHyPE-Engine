#if defined(DistributedOffloading) && defined(Parallel) && defined(SharedTBB)

#if defined(ScoreP)
#include "scorep/SCOREP_User.h"
#endif

#if defined(FileTrace)
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <unistd.h>
#include <sys/time.h>
#include <sstream>
#include "tarch/parallel/Node.h"
#include "tarch/multicore/Core.h"
#endif

#include "ADERDGSolver.h"

#include "exahype/offloading/PerformanceMonitor.h"
#include "exahype/offloading/OffloadingAnalyser.h"
#include "exahype/offloading/OffloadingProfiler.h"
#include "exahype/offloading/ReplicationStatistics.h"

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
#if defined ScoreP
	SCOREP_USER_REGION( "exahype::solvers::ADERDGSolver::MigratablePredictionJob::run", SCOREP_USER_REGION_TYPE_FUNCTION )
#endif
     #if defined FileTrace
     auto start = std::chrono::high_resolution_clock::now();
     #endif

     bool result = false;
     int curr = std::atomic_fetch_add(&JobCounter, 1);

     if(curr%1000==0) {
       tarch::timing::Watch watch("exahype::StealablePredictionJob::", "-", false,false);
       watch.startTimer();
       result = handleExecution();
       watch.stopTimer();

       logDebug("run()","measured time per STP "<<watch.getCalendarTime());

       exahype::offloading::OffloadingAnalyser::getInstance().setTimePerSTP(watch.getCalendarTime());
     }
     else
       result = handleExecution();

     #if defined FileTrace 
     auto stop = std::chrono::high_resolution_clock::now(); 
     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 

     std::stringstream stream;
     stream<<"./TraceOutput/exahype_solvers_ADERDGSolver_MigratablePredictionJob_run_rank_";
     int rank=tarch::parallel::Node::getInstance().getRank();
     stream<<rank<<"_";
     //this will only work for 2 cores per Rank
     int threadId=tarch::multicore::Core::getInstance().getThreadNum();
     stream<<threadId<<".txt";
     std::string path=stream.str();

     std::ofstream file;
     file.open(path,std::fstream::app);
     file << duration.count() << std::endl;
     file.close();
     #endif

#ifdef USE_ITAC
      VT_end(event_stp);
#endif
  return result;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::handleLocalExecution() {
  int myRank = tarch::parallel::Node::getInstance().getRank();
  bool result = false;

  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,_element);

  double *luh    = static_cast<double*>(cellDescription.getSolution());
  double *lduh   = static_cast<double*>(cellDescription.getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

#if defined(TaskSharing)
  //check if task outcome has been received already
  tarch::la::Vector<DIMENSIONS, double> center;
  center = cellDescription.getOffset()+0.5*cellDescription.getSize();

  JobTableKey  key;
  for(int i; i<DIMENSIONS; i++)
     key.center[i] = center[i];
  key.timestamp = _predictorTimeStamp;
  key.element = _element;

  logDebug("handleLocalExecution()", "team "<<exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank()
  	                          	      <<" looking for replica center[0] = "<< center[0]
					                  <<" center[1] = "<< center[1]
			                          <<" center[2] = "<< center[2]
					                  <<" time stamp = "<<_predictorTimeStamp
	                                  <<" hash = "<<(size_t) key);

  exahype::solvers::ADERDGSolver::progressOffloading(&_solver, false);
  //exahype::solvers::ADERDGSolver::pollForOutstandingCommunicationRequests(&_solver);

  tbb::concurrent_hash_map<JobTableKey, JobTableEntry>::accessor a_jobToData;
  bool found = _solver._jobDatabase.find(a_jobToData, key);
  if(found && a_jobToData->second.status == ReplicationStatus::received ) {
  	StealablePredictionJobData *data = a_jobToData->second.data;
  	assert(data->_metadata[2*DIMENSIONS]==_predictorTimeStamp);
  	logDebug("handleLocalExecution()", "team "<<exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank()
  			                           <<" found STP in received replica jobs:"
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
    logDebug("handleLocalExecution()",   "team "<<exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank()
                                          <<" Data not available, gotta do it on my own!"
                                          <<" center[0] = "<<center[0]
                                          <<" center[1] = "<<center[1]
                                          <<" center[2] = "<<center[2]
                                          <<" time stamp = "<<_predictorTimeStamp);
    exahype::offloading::ReplicationStatistics::getInstance().notifyExecutedTask();
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
      StealablePredictionJobData *data = a_jobToData->second.data;
      exahype::offloading::ReplicationStatistics::getInstance().notifyLateTask();

      _solver._jobDatabase.erase(a_jobToData);
      a_jobToData.release();
      AllocatedSTPsReceive--;
      delete data;
    }
    else {
#if defined(TaskSharingUseHandshake)
      if(AllocatedSTPsSend<=exahype::offloading::PerformanceMonitor::getInstance().getTasksPerTimestep()) {
        _solver.sendKeyOfReplicatedSTPToOtherTeams(this);
      }
#else
      if(AllocatedSTPsSend<=exahype::offloading::PerformanceMonitor::getInstance().getTasksPerTimestep()) {
  //  if(AllocatedSTPsSend<=1000) {
        SentSTPs++;
        _solver.sendFullReplicatedSTPToOtherTeams(this);
      }
#endif
    }
  }
#endif

  if(!result) {
    exahype::offloading::PerformanceMonitor::getInstance().decRemainingTasks();
    exahype::offloading::PerformanceMonitor::getInstance().decCurrentTasks();
  }
  return result;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::handleExecution() {
  int myRank = tarch::parallel::Node::getInstance().getRank();
  bool result = false;

#if defined(OffloadingUseProfiler)
  exahype::offloading::OffloadingProfiler::getInstance().beginComputation();
  double time = -MPI_Wtime();
#endif
  //local treatment if this job belongs to the local rank
  if(_originRank==myRank) {
    result = handleLocalExecution();
    NumberOfEnclaveJobs--;
    assertion( NumberOfEnclaveJobs>=0 );
  }
  //remote task, need to execute and send it back
  else {
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
#if defined(FileTrace)
  std::stringstream stream;
  stream.str(std::string());
  stream<<"./TraceOutput/exahype_solvers_ADERDGSolver_OwnMigratableJob_iterations_rank_";
  int rank=tarch::parallel::Node::getInstance().getRank();
  stream<<rank<<"_";
  //this will only work for 2 cores per Rank
  int threadId=tarch::multicore::Core::getInstance().getThreadNum();
  stream<<threadId<<".txt";
  std::string path=stream.str();

  std::ofstream file;
  file.open(path,std::fstream::app);
  file << iterations << std::endl;
  file.close();
#endif
  }
#if defined(OffloadingUseProfiler)
  time += MPI_Wtime();
  exahype::offloading::OffloadingProfiler::getInstance().endComputation(time);
#endif

  //send back
  if(_originRank!=myRank) {
    MPI_Request sendBackRequests[4];
    //logInfo("handleLocalExecution()", "postSendBack");
    //TODO would be nice to have the amount of picard iterations returned here
    _solver.isendStealablePredictionJob(
					_luh,
         	                        _lduh,
                                        _lQhbnd,
                                        _lFhbnd,
                                        _originRank,
                                        _tag,
                                        exahype::offloading::OffloadingManager::getInstance().getMPICommunicatorMapped(),
                                        sendBackRequests);
    exahype::offloading::OffloadingManager::getInstance().submitRequests(
    		                        sendBackRequests,4,
					_tag,
					_originRank,
					sendBackHandler,
					exahype::offloading::RequestType::sendBack,
					&_solver);
  /*#if defined(FileTrace)
  std::stringstream stream;
  stream.str(std::string());
  stream<<"./TraceOutput/exahype_solvers_ADERDGSolver_ForeignMigratableJob_iterations_rank_";
  int rank=tarch::parallel::Node::getInstance().getRank();
  stream<<rank<<"_";
  //this will only work for 2 cores per Rank
  int threadId=tarch::multicore::Core::getInstance().getThreadNum();
  stream<<threadId<<".txt";
  std::string path=stream.str();

  std::ofstream file;
  file.open(path,std::fstream::app);
  file << iterations << std::endl;
  file.close();
  #endif*/
  }
  return result;
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveHandler(exahype::solvers::Solver* solver, int tag, int remoteRank) {
  //logInfo("receiveHandler","successful receive request");

  tbb::concurrent_hash_map<std::pair<int, int>, StealablePredictionJobData*>::accessor a_tagRankToData;
  StealablePredictionJobData *data;
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

  tbb::concurrent_hash_map<std::pair<int, int>, StealablePredictionJobData*>::accessor a_tagRankToData;
  StealablePredictionJobData *data;
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

  if(key.timestamp<static_cast<exahype::solvers::ADERDGSolver*> (solver)->getMinTimeStamp()) {
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
  cellDescription->setHasCompletedLastStep(true);

  tbb::concurrent_hash_map<const CellDescription*, std::pair<int,int>>::accessor a_cellDescToTagRank;
  found =  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapCellDescToTagRank.find(a_cellDescToTagRank, cellDescription);
  assertion(found);
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapCellDescToTagRank.erase(a_cellDescToTagRank);
  a_cellDescToTagRank.release();

  NumberOfEnclaveJobs--;
  NumberOfRemoteJobs--;

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
  tbb::concurrent_hash_map<std::pair<int, int>, StealablePredictionJobData*>::accessor a_tagRankToData;

  StealablePredictionJobData *data;
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
  tbb::concurrent_hash_map<int, StealablePredictionJobData*>::accessor a_tagToData;
  bool found = static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagToReplicationSendData.find(a_tagToData, tag);
  assert(found);
  StealablePredictionJobData *data = a_tagToData->second;
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

}

exahype::solvers::ADERDGSolver::StealablePredictionJobData::StealablePredictionJobData( ADERDGSolver& solver ) :
  _luh(solver.getDataPerCell()),
  _lduh(solver.getUpdateSize()),
  _lQhbnd(solver.getBndTotalSize()),
  _lFhbnd(solver.getBndFluxTotalSize())
{
  AllocatedSTPs++;
};

exahype::solvers::ADERDGSolver::StealablePredictionJobData::~StealablePredictionJobData() {
  AllocatedSTPs--;
};

#endif
