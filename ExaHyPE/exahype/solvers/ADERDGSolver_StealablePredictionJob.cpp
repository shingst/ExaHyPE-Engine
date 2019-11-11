#ifdef DistributedOffloading

#include "ADERDGSolver.h"

#include "exahype/offloading/PerformanceMonitor.h"
#include "exahype/offloading/OffloadingAnalyser.h"
#include "exahype/offloading/OffloadingProfiler.h"
#include "exahype/offloading/ReplicationStatistics.h"

exahype::solvers::ADERDGSolver::StealablePredictionJob::StealablePredictionJob(
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

exahype::solvers::ADERDGSolver::StealablePredictionJob::StealablePredictionJob(
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

exahype::solvers::ADERDGSolver::StealablePredictionJob::~StealablePredictionJob() {};

//Caution: Compression and restriction are not supported yet!
bool exahype::solvers::ADERDGSolver::StealablePredictionJob::run( bool isCalledOnMaster ) {

#ifdef USE_ITAC
      VT_begin(event_stp);
#endif
     bool result = false;
     int curr = std::atomic_fetch_add(&JobCounter, 1);

     if(curr%1000==0) {
       tarch::timing::Watch watch("exahype::StealablePredictionJob::", "-", false,false);
       watch.startTimer();
       result = handleLocalExecution();
       watch.stopTimer();

       logDebug("run()","measured time per STP "<<watch.getCalendarTime());

       exahype::offloading::OffloadingAnalyser::getInstance().setTimePerSTP(watch.getCalendarTime());
     }
     else
       result = handleLocalExecution();

#ifdef USE_ITAC
      VT_end(event_stp);
#endif
  return result;
}

bool exahype::solvers::ADERDGSolver::StealablePredictionJob::handleLocalExecution() {
  int myRank = tarch::parallel::Node::getInstance().getRank();
  bool result = false;

#if defined(PerformanceAnalysisOffloading)
  tarch::timing::Watch watch("exahype::offloading::", "-", false,false);
  watch.startTimer();
#endif
#if defined(OffloadingUseProfiler)
  exahype::offloading::OffloadingProfiler::getInstance().beginComputation();
  double time = -MPI_Wtime();
#endif

  if(_originRank==myRank) {
    CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,_element);

    double *luh    = static_cast<double*>(cellDescription.getSolution());
    double *lduh   = static_cast<double*>(cellDescription.getUpdate());
    double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
    double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

    logInfo("handleLocalExecution()", "team "<<exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank()<<" cellDescriptionIndex "<<_cellDescriptionsIndex<<" element: "<<_element<<" time stamp: "<<_predictorTimeStamp);

#if defined(ReplicationSaving)
    tarch::la::Vector<DIMENSIONS, double> center;
    center = cellDescription.getOffset()+0.5*cellDescription.getSize();
    JobTableKey  key;
    for(int i; i<DIMENSIONS; i++)
       key.center[i] = center[i];
    key.timestamp = _predictorTimeStamp;
    key.element = _element;
    logInfo("handleLocalExecution()", "team "<<exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank()
    	                          	  <<" looking for replica center[0] = "<< center[0]
					                  <<" center[1] = "<< center[1]
			                          <<" center[2] = "<< center[2]
					                  <<" time stamp = "<<_predictorTimeStamp
	                                  <<" hash = "<<(size_t) key);
    tbb::concurrent_hash_map<JobTableKey, JobTableEntry>::accessor a_jobToData;
    bool found = _solver._jobDatabase.find(a_jobToData, key);
    if(found && a_jobToData->second.status == ReplicationStatus::received ) {
    	StealablePredictionJobData *data = a_jobToData->second.data;
    	assert(data->_metadata[2*DIMENSIONS]==_predictorTimeStamp);
    	logInfo("handleLocalExecution()", "team "<<exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank()
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
    }
    else {
    	if (found && a_jobToData->second.status == ReplicationStatus::transit) {
    	    	logInfo("handleLocalExecution()", "task is in transit, we may want to wait!");
#ifdef ReplicationSavingRescheduleIfInTransit
                a_jobToData.release();
    	        return true;
#endif
    	 }
        a_jobToData.release();
    	logInfo("handleLocalExecution()",   "team "<<exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank()
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

#if defined (ReplicationSaving)
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
#if defined(ReplicationSavingUseHandshake)
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

#endif
    }
    exahype::offloading::PerformanceMonitor::getInstance().decRemainingTasks();
  }
  else {
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
  exahype::offloading::PerformanceMonitor::getInstance().decCurrentTasks();
#if defined(OffloadingUseProfiler)
  time += MPI_Wtime();
  exahype::offloading::OffloadingProfiler::getInstance().endComputation(time);
#endif
#if defined(PerformanceAnalysisOffloading)
  watch.stopTimer();
  if(watch.getCalendarTime() >= 0.0) {
    logInfo(
        "localCompute()",
        "remaining "<<NumberOfEnclaveJobs+NumberOfSkeletonJobs-NumberOfRemoteJobs<<" "<<
        "time=" << std::fixed <<
        watch.getCalendarTime() <<
        ", cpu time=" <<
        watch.getCPUTime()
    );
  }
#endif

  if(_originRank!=myRank) {
#if defined(PerformanceAnalysisOffloading)
    watch.startTimer();
#endif
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
    exahype::offloading::OffloadingManager::getInstance().submitRequests(sendBackRequests, 4, _tag, _originRank, sendBackHandler, exahype::offloading::RequestType::sendBack, &_solver);

#if defined(PerformanceAnalysisOffloading)
    watch.stopTimer();
    if(watch.getCalendarTime() >= 0.0) {
      logInfo(
          "remoteReturnSend()",
          "time=" << std::fixed <<
          watch.getCalendarTime() <<
          ", cpu time=" <<
          watch.getCPUTime()
      );
    }
#endif
  }
  else {
#if defined(PerformanceAnalysisOffloading)
    watch.startTimer();
#endif
    NumberOfEnclaveJobs--;
    assertion( NumberOfEnclaveJobs>=0 );
  }
  return result;
}

void exahype::solvers::ADERDGSolver::StealablePredictionJob::receiveHandler(exahype::solvers::Solver* solver, int tag, int remoteRank) {
#if defined(PerformanceAnalysisOffloadingDetailed)
  logInfo("receiveHandler","successful receive request");
#endif
  //logInfo("receiveHandler","successful receive request");

  tbb::concurrent_hash_map<std::pair<int, int>, StealablePredictionJobData*>::accessor a_tagRankToData;
  StealablePredictionJobData *data;
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagRankToStolenData.find(a_tagRankToData, std::make_pair(remoteRank, tag));
  data = a_tagRankToData->second;
  a_tagRankToData.release();

  exahype::offloading::OffloadingAnalyser::getInstance().notifyReceivedSTPJob();
  StealablePredictionJob *job= static_cast<exahype::solvers::ADERDGSolver*> (solver)->createFromData(data, remoteRank, tag);
  peano::datatraversal::TaskSet spawnedSet(job);

  exahype::offloading::OffloadingProfiler::getInstance().notifyReceivedTask(remoteRank);
}

#if defined(ReplicationSaving)
void exahype::solvers::ADERDGSolver::StealablePredictionJob::receiveKeyHandlerReplication(exahype::solvers::Solver* solver, int tag, int remoteRank) {
#if defined(PerformanceAnalysisOffloadingDetailed)
  logInfo("receiveKeyHandlerReplica","successful receive request");
#endif
  logInfo("receiveKeyHandlerReplica","successful receive request");

  tbb::concurrent_hash_map<std::pair<int, int>, double*>::accessor a_tagRankToData;
  double *key;
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagRankToReplicaKey.find(a_tagRankToData, std::make_pair(remoteRank, tag));
  key = a_tagRankToData->second;
  //static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagRankToReplicaKey.erase(a_tagRankToData);
  //a_tagRankToData.release();

  logInfo("receiveKeyHandlerReplica", "team "
		                           <<" received replica job key: center[0] = "<<key[0]
							       <<" center[1] = "<<key[1]
								   <<" center[2] = "<<key[2]
								   <<" time stamp = "<<key[2*DIMENSIONS]
								   <<" element = "<<(int) key[2*DIMENSIONS+2]);

  static_cast<exahype::solvers::ADERDGSolver*> (solver)->sendRequestForJobAndReceive(tag, remoteRank, key);

}

void exahype::solvers::ADERDGSolver::StealablePredictionJob::receiveHandlerReplication(exahype::solvers::Solver* solver, int tag, int remoteRank) {
#if defined(PerformanceAnalysisOffloadingDetailed)
  logInfo("receiveHandlerReplica","successful receive request");
#endif
  logInfo("receiveHandlerReplica","successful receive request");

  tbb::concurrent_hash_map<std::pair<int, int>, StealablePredictionJobData*>::accessor a_tagRankToData;
  StealablePredictionJobData *data;
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagRankToReplicaData.find(a_tagRankToData, std::make_pair(remoteRank, tag));
  data = a_tagRankToData->second;
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagRankToReplicaData.erase(a_tagRankToData);
  a_tagRankToData.release();

  logInfo("receiveHandlerReplica", "team "<<exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank()
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

void exahype::solvers::ADERDGSolver::StealablePredictionJob::receiveBackHandler(exahype::solvers::Solver* solver, int tag, int remoteRank) {
#if defined(PerformanceAnalysisOffloadingDetailed)
  static std::atomic<int> cnt=0;
  cnt++;
  logInfo("receiveBackHandler","successful receiveBack request, cnt "<<cnt);
#endif

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

void exahype::solvers::ADERDGSolver::StealablePredictionJob::sendBackHandler(exahype::solvers::Solver* solver, int tag, int remoteRank) {
#if defined(PerformanceAnalysisOffloadingDetailed)
  logInfo("sendBackHandler","successful sendBack request");
#endif
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

#if defined(ReplicationSaving)
void exahype::solvers::ADERDGSolver::StealablePredictionJob::sendAckHandlerReplication(
		exahype::solvers::Solver* solver,
		int tag,
		int remoteRank)
{

  logInfo("sendAckHandlerReplication","successfully completed send ack to other team");

}

void exahype::solvers::ADERDGSolver::StealablePredictionJob::sendKeyHandlerReplication(
		exahype::solvers::Solver* solver,
		int tag,
		int remoteRank)
{

  logInfo("sendKeyHandlerReplication","successfully completed send key to other teams");
  //exahype::offloading::ReplicationStatistics::getInstance().notifySentTask();
  //tbb::concurrent_hash_map<int, double*>::accessor a_tagToData;
  //bool found = static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagToReplicationSendKey.find(a_tagToData, tag);
  //assert(found);
  //double *data = a_tagToData->second;
  //static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagToReplicationSendKey.erase(a_tagToData);

  //delete data;
  //a_tagToData.release();
}

void exahype::solvers::ADERDGSolver::StealablePredictionJob::sendHandlerReplication(
		exahype::solvers::Solver* solver,
		int tag,
		int remoteRank)
{

  logInfo("sendHandlerReplication","successfully completed send to other teams");
  exahype::offloading::ReplicationStatistics::getInstance().notifySentTask();
  tbb::concurrent_hash_map<int, StealablePredictionJobData*>::accessor a_tagToData;
  bool found = static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagToReplicationSendData.find(a_tagToData, tag);
  assert(found);
  StealablePredictionJobData *data = a_tagToData->second;
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapTagToReplicationSendData.erase(a_tagToData);
  delete data;
  AllocatedSTPsSend--;
  CompletedSentSTPs++;
  logInfo("sendHandlerReplication"," allocated stps send "<<AllocatedSTPsSend);
  a_tagToData.release();
}
#endif

void exahype::solvers::ADERDGSolver::StealablePredictionJob::sendHandler(exahype::solvers::Solver* solver, int tag, int remoteRank) {
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
      exahype::solvers::ADERDGSolver::StealablePredictionJob::receiveBackHandler,
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
