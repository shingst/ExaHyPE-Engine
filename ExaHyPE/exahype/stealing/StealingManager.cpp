#if  defined(SharedTBB)  && defined(Parallel) && defined(DistributedStealing)

#include "StealingManager.h"

#include <unordered_set>
#include <vector>
#include <string>
#include <algorithm>

#include "tarch/multicore/Lock.h"
#include "tarch/parallel/Node.h"

#include "exahype/stealing/StealingProfiler.h"
#include "exahype/stealing/StaticDistributor.h"
#include "exahype/stealing/DynamicDistributor.h"
#include "exahype/stealing/DiffusiveDistributor.h"
#include "exahype/stealing/AggressiveDistributor.h"
#include "exahype/stealing/AggressiveCCPDistributor.h"
#include "exahype/stealing/AggressiveHybridDistributor.h"
#include "exahype/stealing/PerformanceMonitor.h"

#ifdef USE_ITAC
#include "VT.h"

static int event_handling;
static int event_progress_send;
static int event_progress_receive;
static int event_progress_sendBack;
static int event_progress_receiveBack;
#endif

tarch::logging::Log exahype::stealing::StealingManager::_log( "exahype::stealing::stealingManager" );

exahype::stealing::StealingManager::StealingManager() :
    _nextRequestId(0),
    _nextGroupId(0),
    _runningAndReceivingBack(false),
    _stealingComm(MPI_COMM_NULL),
    _stealingCommMapped(MPI_COMM_NULL),
    _emergencyTriggered(false),
    _numProgressJobs(0),
    _hasNotifiedSendCompleted(false)
    //_numProgressSendJobs(0),
    //_numProgressReceiveJobs(0),
    //_numProgressReceiveBackJobs(0)
{
#ifdef USE_ITAC
  static const char *event_name_handle = "handleRequest";
  int ierr = VT_funcdef(event_name_handle, VT_NOCLASS, &event_handling); assertion(ierr==0);
  static const char *event_name_send = "progressSends";
  ierr = VT_funcdef(event_name_send, VT_NOCLASS, &event_progress_send); assertion(ierr==0);
  static const char *event_name_receive = "progressReceives";
  ierr = VT_funcdef(event_name_receive, VT_NOCLASS, &event_progress_receive); assertion(ierr==0);
  static const char *event_name_sendB = "progressSendBacks";
  ierr = VT_funcdef(event_name_sendB, VT_NOCLASS, &event_progress_sendBack); assertion(ierr==0);
  static const char *event_name_receiveB = "progressReceiveBacks";
  ierr = VT_funcdef(event_name_receiveB, VT_NOCLASS, &event_progress_receiveBack); assertion(ierr==0);
#endif
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  _emergencyHeatMap = new double[nnodes];
  std::fill(&_emergencyHeatMap[0], &_emergencyHeatMap[nnodes], 0);
}

exahype::stealing::StealingManager::~StealingManager()
{
  delete[] _emergencyHeatMap;
}

bool exahype::stealing::StealingManager::getRunningAndReceivingBack() {
	return _runningAndReceivingBack;
}

void exahype::stealing::StealingManager::setRunningAndReceivingBack() {
	_runningAndReceivingBack = true;
}

void exahype::stealing::StealingManager::resetRunningAndReceivingBack() {
	_runningAndReceivingBack = false;
}

void exahype::stealing::StealingManager::createMPICommunicator() {
  int ierr = MPI_Comm_dup(MPI_COMM_WORLD, &_stealingComm);
  assertion(ierr==MPI_SUCCESS);
  ierr = MPI_Comm_dup(MPI_COMM_WORLD, &_stealingCommMapped);
  assertion(ierr==MPI_SUCCESS);
}

void exahype::stealing::StealingManager::destroyMPICommunicator() {
  int ierr = MPI_Comm_free( &_stealingComm);
  assertion(ierr==MPI_SUCCESS);
  ierr = MPI_Comm_free(&_stealingCommMapped);
  assertion(ierr==MPI_SUCCESS);
}


MPI_Comm exahype::stealing::StealingManager::getMPICommunicator() {
  return _stealingComm;
}

MPI_Comm exahype::stealing::StealingManager::getMPICommunicatorMapped() {
  return _stealingCommMapped;
}

exahype::stealing::StealingManager& exahype::stealing::StealingManager::getInstance() {
  static StealingManager stealingManager;
  return stealingManager;
}

int exahype::stealing::StealingManager::requestTypeToMap( RequestType requestType ) {
  return static_cast<int> (requestType);
}

int exahype::stealing::StealingManager::getStealingTag() {
  static std::atomic<int> counter = 1; //0 is reserved for status
  return counter.fetch_add(1);
}

void exahype::stealing::StealingManager::submitRequests(
    MPI_Request *requests,
	int nRequests,
	int tag,
	int remoteRank,
    std::function<void(exahype::solvers::Solver*, int, int)> handler,
    RequestType type,
	exahype::solvers::Solver *solver,
	bool block ) {
    
  //static std::atomic<int> submitted[4];
  int finished = -1;
  MPI_Testall(nRequests, requests, &finished, MPI_STATUSES_IGNORE);
  if(finished) {
     handler(solver, tag, remoteRank);
     return;
  }

  if(block) {
    MPI_Waitall(nRequests, requests, MPI_STATUSES_IGNORE);
    handler(solver, tag, remoteRank);
    return;
  }

  int mapId = requestTypeToMap(type);
  
  //submitted[mapId]++;
  //logInfo("stealingManager","submitted["<<mapId<<"]:"<<submitted[mapId]);

  // assign group id for this request group
  int groupId = _nextGroupId++;
  // insert metadata into maps
  std::pair<int, std::function<void(exahype::solvers::Solver*, int, int)>> handlerElem(groupId, handler);
  std::pair<int, exahype::solvers::Solver*> solverElem(groupId, solver);
  std::pair<int, int> outstandingElem(groupId, nRequests);
  std::pair<int, int> remoteRankElem(groupId, remoteRank);
  std::pair<int, int> remoteTagElem(groupId, tag);

  _handlers[mapId].insert(handlerElem);
  _solvers[mapId].insert(solverElem);
  _outstandingReqsForGroup[mapId].insert(outstandingElem);
  _remoteRanksForGroup[mapId].insert(remoteRankElem);
  _remoteTagsForGroup[mapId].insert(remoteTagElem);

  // push requests into queue
  for(int i=0; i<nRequests; i++) {
    //TODO: avoid overflow!
    int id=_nextRequestId++;

    std::pair<int, MPI_Request> reqElem(id, requests[i]);
    std::pair<int, int> reqGroupElem(id, groupId);
    //logInfo("stealingManager", "inserted groupid "<<groupId<<" for req "<<id);

    _idToRequest[mapId].insert(reqElem);
    _requestToGroup[mapId].insert(reqGroupElem);
    _requests[mapId].push(id);
  }

#ifdef StealingUseProgressTask
  if(_numProgressJobs==0 && type==RequestType::send) {
    //logInfo("submitRequests()", "spawning progress job (high priority)");
    _numProgressJobs++;
    ProgressJob *job = new ProgressJob();
    peano::datatraversal::TaskSet spawnedSet( job);
  }
 /* if(type==RequestType::send && _numProgressSendJobs==0) {
    //logInfo("submitRequests()", "spawning progress send job (high priority)");
    _numProgressSendJobs++;
    ProgressSendJob *job = new ProgressSendJob();
    peano::datatraversal::TaskSet spawnedSet( job, peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible);
  }
  
  if(type==RequestType::receive && _numProgressReceiveJobs==0) {
    logInfo("submitRequests()", "spawning progress receive job (high priority)");
    _numProgressReceiveJobs++;
    ProgressReceiveJob *job = new ProgressReceiveJob();
    peano::datatraversal::TaskSet spawnedSet( job, peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible);
  }

  if(type==RequestType::receiveBack && _numProgressReceiveBackJobs==0) {
    logInfo("submitRequests()", "spawning progress receive back job (high priority)");
    _numProgressReceiveBackJobs++;
    ProgressReceiveBackJob *job = new ProgressReceiveBackJob();
    peano::datatraversal::TaskSet spawnedSet( job, peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible);
  }*/
#endif
}

void exahype::stealing::StealingManager::createRequestArray(
    RequestType type,
    std::vector<MPI_Request> &requests,
    std::unordered_map<int, int> &map) {

  int mapId = requestTypeToMap(type);

  bool gotOne = true;
  int j = 0;

  while(gotOne) {
    int req_id;
    gotOne = _requests[mapId].try_pop(req_id);
    if(gotOne) {
      tbb::concurrent_hash_map<int, MPI_Request>::accessor a_requests;
      bool found = _idToRequest[mapId].find(a_requests, req_id);
      assertion(found);
      MPI_Request request = a_requests->second;
      a_requests.release();
      requests.push_back(request);
      map.insert(std::pair<int, int>(j, req_id));
      j++;
    }
  }
}

bool exahype::stealing::StealingManager::hasOutstandingRequestOfType(RequestType requestType) {
  return (!_requests[requestTypeToMap(requestType)].empty() || !_currentOutstandingRequests[requestTypeToMap(requestType)].size()==0);
}

void exahype::stealing::StealingManager::progressRequests() {

  if(hasOutstandingRequestOfType(RequestType::send)) {
#ifdef USE_ITAC
    VT_begin(event_progress_send);
#endif
    progressRequestsOfType(RequestType::send);
#ifdef USE_ITAC
    VT_end(event_progress_send);
#endif
  }
  else if (hasOutstandingRequestOfType(RequestType::receive)) {
#ifdef USE_ITAC
    VT_begin(event_progress_receive);
#endif
    progressRequestsOfType(RequestType::receive);
#ifdef USE_ITAC
    VT_end(event_progress_receive);
#endif
  }
  else if (hasOutstandingRequestOfType(RequestType::receiveBack)) {
#ifdef USE_ITAC
    VT_begin(event_progress_receiveBack);
#endif
    progressRequestsOfType(RequestType::receiveBack);
#ifdef USE_ITAC
    VT_end(event_progress_receiveBack);
#endif
  }
  else if (hasOutstandingRequestOfType(RequestType::sendBack)) {
#ifdef USE_ITAC
    VT_begin(event_progress_sendBack);
#endif
    progressRequestsOfType(RequestType::sendBack);
#ifdef USE_ITAC
    VT_end(event_progress_sendBack);
#endif
  }
}

void exahype::stealing::StealingManager::progressAnyRequests() {

  if(hasOutstandingRequestOfType(RequestType::send)) {
#ifdef USE_ITAC
    VT_begin(event_progress_send);
#endif
    progressRequestsOfType(RequestType::send);
#ifdef USE_ITAC
    VT_end(event_progress_send);
#endif
  }
  if (hasOutstandingRequestOfType(RequestType::receive)) {
#ifdef USE_ITAC
    VT_begin(event_progress_receive);
#endif
    progressRequestsOfType(RequestType::receive);
#ifdef USE_ITAC
    VT_end(event_progress_receive);
#endif
  }
  if (hasOutstandingRequestOfType(RequestType::receiveBack)) {
#ifdef USE_ITAC
    VT_begin(event_progress_receiveBack);
#endif
    progressRequestsOfType(RequestType::receiveBack);
#ifdef USE_ITAC
    VT_end(event_progress_receiveBack);
#endif
  }
  if (hasOutstandingRequestOfType(RequestType::sendBack)) {
#ifdef USE_ITAC
    VT_begin(event_progress_sendBack);
#endif
    progressRequestsOfType(RequestType::sendBack);
#ifdef USE_ITAC
    VT_end(event_progress_sendBack);
#endif
  }
}

bool exahype::stealing::StealingManager::progressRequestsOfType( RequestType type ) {

  // First, we ensure here that only one thread at a time progresses stealing
  // this attempts to avoid multithreaded MPI problems
  tarch::multicore::Lock lock(_progressSemaphore, false);
  bool canRun = lock.try_lock();
  if(!canRun) {
#if defined(PerformanceAnalysisStealingDetailed)
    watch.stopTimer();
    if(watch.getCalendarTime() >= 0.0) {
      logInfo(
          "progressStealing() ",
          "couldn't run "<<
          "time=" << std::fixed <<
          watch.getCalendarTime() <<
          ", cpu time=" <<
          watch.getCPUTime()
      );
    }
#endif
    return false;
  }

  int mapId = requestTypeToMap(type);

  tbb::concurrent_hash_map<int, std::function<void(exahype::solvers::Solver*, int, int)>>::accessor a_handler;
  tbb::concurrent_hash_map<int, exahype::solvers::Solver*>::accessor                                a_solver;
  tbb::concurrent_hash_map<MPI_Request, int>::accessor                                              a_groupId;
  tbb::concurrent_hash_map<int, int>::accessor                                                      a_outstandingGroup;
  tbb::concurrent_hash_map<int, int>::accessor                                                      a_remoteRank;
  tbb::concurrent_hash_map<int, int>::accessor                                                      a_remoteTag;

  //std::vector<MPI_Request>     outstandingRequests;
  //std::unordered_map<int, int> vecIdToReqId;

  int nRequests = _currentOutstandingRequests[mapId].size();
  if(nRequests==0) {
    //logInfo("progressRequestsOfType()", "begin create req array");
    createRequestArray( type, _currentOutstandingRequests[mapId], _currentOutstandingVecIdxToReqid[mapId] );
    nRequests = _currentOutstandingRequests[mapId].size();
    //logInfo("progressRequestsOfType()", "end create req array");
  }

  std::vector<MPI_Status> stats(nRequests);
  std::vector<int> arrOfIndices(nRequests);
  int outcount = 0;

  double time = -MPI_Wtime();
  exahype::stealing::StealingProfiler::getInstance().beginCommunication();


// For DEBUGGING
//static std::atomic<int> finished_cnt[4];
//  logInfo("stealingManager", "testsome of "<<nRequests<< " of type "<<mapId);
//  std::vector<MPI_Request> copyRequests(outstandingRequests);

//	for(int i=0;i<nRequests;i++) {
//	  MPI_Request search = copyRequests[i];
//	  for(int j=0;j<nRequests;j++) {
//      if(i!=j && copyRequests[j]==search && search!=MPI_REQUEST_NULL) {
//        logInfo("stealingManager", "found duplicate request: i "<<i<<" j "<<j<<" request i: "<<copyRequests[i]<<" request j: "<<copyRequests[j]);
//        assertion(false);
//  	}
//	 }
// }
  /*if(type==RequestType::send)
     logInfo("progressRequestsOfType()", "progressing sends");
  if(type==RequestType::receive)
     logInfo("progressRequestsOfType()", "progressing receives");
  if(type==RequestType::sendBack)
     logInfo("progressRequestsOfType()", "progressing send back");
  if(type==RequestType::receiveBack)
     logInfo("progressRequestsOfType()", "progressing receive back");*/
  //TODO: keine Statusse
  int ierr = MPI_Testsome(nRequests,&_currentOutstandingRequests[mapId][0], &outcount, &arrOfIndices[0], MPI_STATUSES_IGNORE);

//  if(ierr != MPI_SUCCESS) {
//    for(int i=0;i<nRequests;i++) {
//      int ierrstatus = stats[i].MPI_ERROR;
//      if(ierrstatus!=MPI_SUCCESS) {
//        logInfo("stealingManager", "error "<<ierrstatus<<" for request "<<vecIdToReqId[i]<< " source "<<stats[i].MPI_SOURCE<<" tag "<<stats[i].MPI_TAG);
//      }
//      char err_buffer[MPI_MAX_ERROR_STRING];
//      int resultlen = 0;
//      if(ierrstatus!=MPI_SUCCESS) {
//        MPI_Error_string(ierrstatus,err_buffer,&resultlen);
//        fprintf(stderr,err_buffer);
//      }
//    }
//    MPI_Abort(MPI_COMM_WORLD, ierr); /* abort*/
//  }
  time += MPI_Wtime();

  if(outcount>0)
    exahype::stealing::StealingProfiler::getInstance().endCommunication(true, time);
  else
    exahype::stealing::StealingProfiler::getInstance().endCommunication(false, time);

  time = -MPI_Wtime();
  exahype::stealing::StealingProfiler::getInstance().beginHandling();
  bool found=false;
  //handle finished requests
  for(int i=0; i<outcount; i++) {
    int reqIdx = arrOfIndices[i];
    int reqId = _currentOutstandingVecIdxToReqid[mapId][reqIdx];
    assertion(_currentOutstandingRequests[mapId][reqIdx]==MPI_REQUEST_NULL);

    _idToRequest[mapId].erase(reqId);
    int groupId;
    found = _requestToGroup[mapId].find(a_groupId, reqId);

    assertion(found);
    groupId = a_groupId->second;
    _requestToGroup[mapId].erase(a_groupId);
    a_groupId.release();

    found = _outstandingReqsForGroup[mapId].find(a_outstandingGroup, groupId);
    assertion(found);
    a_outstandingGroup->second--;
    bool finished = a_outstandingGroup->second==0;
    a_outstandingGroup.release();

    if(finished) {
      std::function<void(exahype::solvers::Solver*, int ,int)> handler;
      found = _handlers[mapId].find(a_handler, groupId);
      assertion(found);
      handler=a_handler->second;
      a_handler.release();
      exahype::solvers::Solver *solver;
      found=_solvers[mapId].find(a_solver, groupId);
      assertion(found);
      solver = a_solver->second;
      a_solver.release();
      int remoteRank = -1;
      found= _remoteRanksForGroup[mapId].find(a_remoteRank, groupId);
      remoteRank = a_remoteRank->second;
      assertion(found);
      a_remoteRank.release();
      int remoteTag = -1;
      found= _remoteTagsForGroup[mapId].find(a_remoteTag, groupId);
      remoteTag = a_remoteTag->second;
      assertion(found);
      a_remoteTag.release();

      handler(solver, remoteTag, remoteRank);
      //RequestHandlerJob *requestHandlerJob = new RequestHandlerJob(handler, solver, remoteTag, remoteRank);
      //peano::datatraversal::TaskSet spawnedSet( requestHandlerJob, peano::datatraversal::TaskSet::TaskType::Background);

      _handlers[mapId].erase(groupId);
      _outstandingReqsForGroup[mapId].erase(groupId);
      _solvers[mapId].erase(groupId);
      _remoteRanksForGroup[mapId].erase(groupId);
      _remoteTagsForGroup[mapId].erase(groupId);
    }
  }

  // push back all unfinished requests
  //for(int i=0; i<nRequests; i++) {
  //  if(outstandingRequests[i]!=MPI_REQUEST_NULL) {
  //    _requests[mapId].push(vecIdToReqId[i]);
  //  }
  //}
  bool allFinished = true;
  for(int i=0; i<nRequests; i++) { 
     if(_currentOutstandingRequests[mapId][i]!=MPI_REQUEST_NULL)  {
       allFinished = false;
       break;
     }
  }
  if(allFinished) {
    _currentOutstandingRequests[mapId].clear();
    _currentOutstandingVecIdxToReqid[mapId].clear();
  }

  time += MPI_Wtime();
  exahype::stealing::StealingProfiler::getInstance().endHandling(time);

  lock.free();  

  return true;
}

void exahype::stealing::StealingManager::triggerVictimFlag() {
  _isVictim = true;
}

void exahype::stealing::StealingManager::resetVictimFlag() {
  _isVictim = false;
}

bool exahype::stealing::StealingManager::isVictim() {
  return _isVictim;
}

void exahype::stealing::StealingManager::triggerEmergencyForRank(int rank) {
//  if(!_emergencyTriggered) {
//    logInfo("triggerEmergency()","emergency event triggered");
//    _emergencyTriggered = true;
//  }
#ifdef StealingStrategyAggressive
  exahype::stealing::AggressiveDistributor::getInstance().handleEmergencyOnRank(rank);
#elif StealingStrategyAggressiveCCP
  exahype::stealing::AggressiveCCPDistributor::getInstance().handleEmergencyOnRank(rank);
#elif StealingStrategyAggressiveHybrid
  exahype::stealing::AggressiveHybridDistributor::getInstance().handleEmergencyOnRank(rank);
#elif StealingStrategyDiffusive
  exahype::stealing::DiffusiveDistributor::getInstance().handleEmergencyOnRank(rank);
#endif
  _emergencyHeatMap[rank]++;
  exahype::stealing::PerformanceMonitor::getInstance().submitBlacklistValueForRank(_emergencyHeatMap[rank], rank);
  logInfo("triggerEmergencyForRank()","blacklist value for rank "<<rank<<":"<<_emergencyHeatMap[rank]);
}

void exahype::stealing::StealingManager::decreaseHeat() {
  logInfo("decreaseHeat()","decrease heat of emergency heat map");
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  for(int i=0; i<nnodes;i++) {
    _emergencyHeatMap[i]*= 0.9;
    if(_emergencyHeatMap[i]>0)
      exahype::stealing::PerformanceMonitor::getInstance().submitBlacklistValueForRank(_emergencyHeatMap[i], i);
    //if(_emergencyHeatMap[i]>0.5) {
    //  logInfo("decreaseHeat()","blacklist value for rank "<<i<<":"<<_emergencyHeatMap[i]);
    //}
  }
} 

bool exahype::stealing::StealingManager::isBlacklisted(int rank) { 
  //return _emergencyHeatMap[rank]>0.5;
  const double* globalHeatMap = exahype::stealing::PerformanceMonitor::getInstance().getBlacklistSnapshot();
  return globalHeatMap[rank]>0.5;
}

void exahype::stealing::StealingManager::printBlacklist() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  const double* globalHeatMap = exahype::stealing::PerformanceMonitor::getInstance().getBlacklistSnapshot();

  for(int i=0; i<nnodes; i++) {
    if(globalHeatMap[i]>0)
      logInfo("printBlacklist", "blacklist value for rank "<<i<<":"<<globalHeatMap[i]);
  }
}

bool exahype::stealing::StealingManager::isEmergencyTriggered() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  return !std::all_of(&_emergencyHeatMap[0], &_emergencyHeatMap[nnodes], [](double d){return d<0.5;});
}

bool exahype::stealing::StealingManager::isEmergencyTriggeredOnRank(int rank) {
  return !_emergencyHeatMap[rank]<0.5;
}

//void exahype::stealing::StealingManager::resetEmergency() {
//  logInfo("resetEmergency()","emergency flag reset");
//  _emergencyTriggered = false;
//}


bool exahype::stealing::StealingManager::selectVictimRank(int& victim) {
#if defined(StealingStrategyStaticHardcoded)
    return exahype::stealing::StaticDistributor::getInstance().selectVictimRank(victim);
#elif defined(StealingStrategyDiffusive)
    return exahype::stealing::DiffusiveDistributor::getInstance().selectVictimRank(victim);
#elif defined(StealingStrategyAggressive)
    return exahype::stealing::AggressiveDistributor::getInstance().selectVictimRank(victim);
#elif defined(StealingStrategyAggressiveCCP)
    return exahype::stealing::AggressiveCCPDistributor::getInstance().selectVictimRank(victim);
#elif defined(StealingStrategyAggressiveHybrid)
    return exahype::stealing::AggressiveHybridDistributor::getInstance().selectVictimRank(victim);
#else
  double remainingLoadRatio = static_cast<double> (exahype::stealing::PerformanceMonitor::getInstance().getRemainingLocalLoad())
		  	  	  	  	  	  /
							  exahype::stealing::PerformanceMonitor::getInstance().getLocalLoadPerTimestep();
  // this is currently hardcoded: the goal is to refrain from giving tasks away if there is not enough work left
  // for overlap of communication and computation
  if(remainingLoadRatio>0.1) {
#if defined(StealingStrategyStatic) 
    return exahype::stealing::StaticDistributor::getInstance().selectVictimRank(victim);
#elif defined(StealingStrategyDynamic)
    return exahype::stealing::DynamicDistributor::getInstance().selectVictimRank(victim);
#elif defined(StealingStrategyHybrid)
    bool staticDistribution = exahype::stealing::StaticDistributor::getInstance().selectVictimRank(victim);
    if(staticDistribution) return true;
    else {
      return exahype::stealing::DynamicDistributor::getInstance().selectVictimRank(victim);
    }
    return false;
#else
# error "Wrong stealing strategy specified!"
#endif
  }
  else {
     logInfo("stealingManager", "could not select victim remaining load ratio "<<remainingLoadRatio);
    return false;
  }
#endif
}

#ifdef StealingUseProgressTask
void exahype::stealing::StealingManager::resetHasNotifiedSendCompleted() {
  _hasNotifiedSendCompleted = false;
}

void exahype::stealing::StealingManager::notifySendCompleted(int rank) {
  char send = 1;
  MPI_Send(&send, 1, MPI_CHAR, rank, 0, getMPICommunicator());
  //logInfo("notifySendCompleted()","sent status message to "<<rank);
}

void exahype::stealing::StealingManager::receiveCompleted(int rank) {
  char receive = 0;
  MPI_Recv(&receive, 1, MPI_CHAR, rank, 0, getMPICommunicator(), MPI_STATUS_IGNORE);
  //logInfo("receiveCompleted()","received status message from "<<rank);
}

void exahype::stealing::StealingManager::notifyAllVictimsSendCompletedIfNotNotified() {
  if(!_hasNotifiedSendCompleted) {
    _hasNotifiedSendCompleted = true;
    std::vector<int> victimRanks;
#if defined(StealingStrategyAggressiveHybrid) 
    exahype::stealing::AggressiveHybridDistributor::getInstance().getAllVictimRanks(victimRanks);
#elif defined(StealingStrategyStaticHardcoded)
    exahype::stealing::StaticDistributor::getInstance().getAllVictimRanks(victimRanks);
#endif
    for(auto victim : victimRanks) 
      notifySendCompleted(victim);
  }
}
#endif

exahype::stealing::StealingManager::RequestHandlerJob::RequestHandlerJob(
    std::function<void(exahype::solvers::Solver*, int, int)> handleRequest,
    exahype::solvers::Solver* solver,
    int tag,
    int remoteRank) :
        _handleRequest(handleRequest),
		_solver(solver),
		_tag(tag),
		_remoteRank(remoteRank)
{};

bool exahype::stealing::StealingManager::RequestHandlerJob::operator()() {
#ifdef USE_ITAC
  VT_begin(event_handling);
#endif
    _handleRequest(_solver, _tag, _remoteRank);
#ifdef USE_ITAC
  VT_end(event_handling);
#endif
  return false;
}


#ifdef StealingUseProgressTask
exahype::stealing::StealingManager::ProgressJob::ProgressJob() :
   tarch::multicore::jobs::Job(tarch::multicore::jobs::JobType::BackgroundTask, 0, tarch::multicore::DefaultPriority*8)
 {}

bool exahype::stealing::StealingManager::ProgressJob::run() {
   int flag;
   //logInfo("submitRequests()", "executing progress job (high priority)");

   int mapId = StealingManager::requestTypeToMap(RequestType::send);
//   while(StealingManager::getInstance()._requests[0].unsafe_size()>0 || StealingManager::getInstance()._currentOutstandingRequests[0].size()>0
//        || StealingManager::getInstance()._requests[1].unsafe_size()>0 || StealingManager::getInstance()._currentOutstandingRequests[1].size()>0
//        || StealingManager::getInstance()._requests[2].unsafe_size()>0 || StealingManager::getInstance()._currentOutstandingRequests[2].size()>0
//        || StealingManager::getInstance()._requests[3].unsafe_size()>0 || StealingManager::getInstance()._currentOutstandingRequests[3].size()>0
//   ) 
   while(StealingManager::getInstance()._requests[mapId].unsafe_size()>0 || StealingManager::getInstance()._currentOutstandingRequests[mapId].size()>0)
   {
     //getInstance().progressAnyRequests();
     getInstance().progressRequestsOfType(RequestType::send);
     MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, StealingManager::getInstance()._stealingComm, &flag, MPI_STATUS_IGNORE);
     MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, StealingManager::getInstance()._stealingCommMapped, &flag, MPI_STATUS_IGNORE);
   }
   //logInfo("submitRequests()", "terminated progress job (high priority)");



   StealingManager::getInstance()._numProgressJobs--;
   return false;
};

exahype::stealing::StealingManager::ProgressSendJob::ProgressSendJob() {}

bool exahype::stealing::StealingManager::ProgressSendJob::operator()() {
   getInstance().progressRequestsOfType(RequestType::send);
   int mapId = StealingManager::requestTypeToMap(RequestType::send);
//   bool reschedule=StealingManager::getInstance()._requests[mapId].unsafe_size()>0;
 
   while(StealingManager::getInstance()._requests[mapId].unsafe_size()>0 || StealingManager::getInstance()._currentOutstandingRequests[mapId].size()>0) {
     getInstance().progressRequestsOfType(RequestType::send);
   }
 
//   if(!reschedule)
   StealingManager::getInstance()._numProgressSendJobs--;
   return false;
};

exahype::stealing::StealingManager::ProgressReceiveJob::ProgressReceiveJob() {}

bool exahype::stealing::StealingManager::ProgressReceiveJob::operator()() {
   getInstance().progressRequestsOfType(RequestType::receive);
   int mapId = StealingManager::requestTypeToMap(RequestType::receive);

   while(StealingManager::getInstance()._requests[mapId].unsafe_size()>0 || StealingManager::getInstance()._currentOutstandingRequests[mapId].size()>0) {
     getInstance().progressRequestsOfType(RequestType::receive);
   }

   StealingManager::getInstance()._numProgressReceiveJobs--;
   return false;
};

exahype::stealing::StealingManager::ProgressReceiveBackJob::ProgressReceiveBackJob() {}

bool exahype::stealing::StealingManager::ProgressReceiveBackJob::operator()() {
   getInstance().progressRequestsOfType(RequestType::receiveBack);
   int mapId = StealingManager::requestTypeToMap(RequestType::receiveBack);

   while(StealingManager::getInstance()._requests[mapId].unsafe_size()>0 || StealingManager::getInstance()._currentOutstandingRequests[mapId].size()>0) {
     getInstance().progressRequestsOfType(RequestType::receiveBack);
   }

   StealingManager::getInstance()._numProgressReceiveBackJobs--;
   return false;
};
#endif

#endif
