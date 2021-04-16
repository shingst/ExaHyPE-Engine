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

#if  defined(Parallel)

#include "../reactive/RequestManager.h"

#include <unordered_set>
#include <vector>
#include <string>
#include <algorithm>

#include "tarch/multicore/Lock.h"
#include "tarch/multicore/Core.h"
#include "tarch/parallel/Node.h"

#include "../reactive/OffloadingContext.h"
#include "../reactive/OffloadingProfiler.h"


#if defined(UseSmartMPI)
#include "mpi_offloading.h"
#endif

#ifdef USE_ITAC
#include "VT.h"

static int event_handling;
static int event_progress_send;
static int event_progress_receive;
static int event_progress_sendBack;
static int event_progress_receiveBack;
#endif

#if defined(USE_TMPI)
#include "teaMPI.h"
#endif

#include <mpi.h>

//#undef assertion
//#define assertion assert

tarch::logging::Log exahype::reactive::RequestManager::_log( "exahype::reactive::RequestManager" );

exahype::reactive::RequestManager* exahype::reactive::RequestManager::_static_managers[MAX_THREADS];

exahype::reactive::RequestManager::RequestManager(int threadId) :
    _threadId(threadId),
    _nextRequestId(0),
    _nextGroupId(0),
    _numProgressJobs(0) {

  int ierr;
#ifdef USE_ITAC
  static const char *event_name_handle = "handleRequest";
  ierr = VT_funcdef(event_name_handle, VT_NOCLASS, &event_handling); assertion(ierr==0);
  static const char *event_name_send = "progressSends";
  ierr = VT_funcdef(event_name_send, VT_NOCLASS, &event_progress_send); assertion(ierr==0);
  static const char *event_name_receive = "progressReceives";
  ierr = VT_funcdef(event_name_receive, VT_NOCLASS, &event_progress_receive); assertion(ierr==0);
  static const char *event_name_sendB = "progressSendBacks";
  ierr = VT_funcdef(event_name_sendB, VT_NOCLASS, &event_progress_sendBack); assertion(ierr==0);
  static const char *event_name_receiveB = "progressReceiveBacks";
  ierr = VT_funcdef(event_name_receiveB, VT_NOCLASS, &event_progress_receiveBack); assertion(ierr==0);
#endif
  MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes()
      * exahype::reactive::OffloadingContext::getInstance().getTMPINumTeams();

  _postedSendsPerRank = new std::atomic<int>[nnodes];
  std::fill(&_postedSendsPerRank[0], &_postedSendsPerRank[nnodes], 0);
  _postedReceivesPerRank = new std::atomic<int>[nnodes];
  std::fill(&_postedReceivesPerRank[0], &_postedReceivesPerRank[nnodes], 0);
  _postedSendBacksPerRank = new std::atomic<int>[nnodes];
  std::fill(&_postedSendBacksPerRank[0], &_postedSendBacksPerRank[nnodes], 0);
  _postedReceiveBacksPerRank = new std::atomic<int>[nnodes];
  std::fill(&_postedReceiveBacksPerRank[0], &_postedReceiveBacksPerRank[nnodes], 0);
  _postedSendOutcomesPerRank = new std::atomic<int>[nnodes];
  std::fill(&_postedSendOutcomesPerRank[0], &_postedSendOutcomesPerRank[nnodes], 0);
  _postedReceiveOutcomesPerRank = new std::atomic<int>[nnodes];
  std::fill(&_postedReceiveOutcomesPerRank[0], &_postedReceiveOutcomesPerRank[nnodes], 0);

}

exahype::reactive::RequestManager::~RequestManager() {
  delete[] _postedSendsPerRank;
  delete[] _postedReceivesPerRank;
  delete[] _postedSendBacksPerRank;
  delete[] _postedReceiveBacksPerRank;
  delete[] _postedSendOutcomesPerRank;
  delete[] _postedReceiveOutcomesPerRank;
}

void exahype::reactive::RequestManager::resetPostedRequests() {
  logDebug("resetPostedRequests","resetting posted requests statistics");
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  std::fill(&_postedSendsPerRank[0], &_postedSendsPerRank[nnodes], 0);
  std::fill(&_postedReceivesPerRank[0], &_postedReceivesPerRank[nnodes], 0);
  std::fill(&_postedSendBacksPerRank[0], &_postedSendBacksPerRank[nnodes], 0);
  std::fill(&_postedReceiveBacksPerRank[0], &_postedReceiveBacksPerRank[nnodes], 0);
  std::fill(&_postedReceiveOutcomesPerRank[0], &_postedReceiveOutcomesPerRank[nnodes], 0);
  std::fill(&_postedSendOutcomesPerRank[0], &_postedSendOutcomesPerRank[nnodes], 0);

}

void exahype::reactive::RequestManager::printPostedRequests() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  for(int i=0; i<nnodes; i++) {
    std::string str="for rank "+std::to_string(i)+":";
    if(_postedSendsPerRank[i]>0)
    str = str + "posted sends: "+std::to_string(_postedSendsPerRank[i]);
    if(_postedReceivesPerRank[i]>0)
    str = str + "posted receives: "+std::to_string(_postedReceivesPerRank[i]);
    if(_postedReceiveBacksPerRank[i]>0)
    str = str + "posted receiveBacks: "+std::to_string(_postedReceiveBacksPerRank[i]);
    if(_postedSendBacksPerRank[i]>0)
    str = str + "posted sendBacks: "+std::to_string(_postedSendBacksPerRank[i]);
    logDebug("printPostedRequests()", str);
  }
}

int exahype::reactive::RequestManager::getNumberOfOutstandingRequests( RequestType type ) {
  int mapId = requestTypeToMsgQueueIdx(type);
  return _outstandingReqsForGroup[mapId].size() + _outstandingRequests[mapId].unsafe_size();
}

exahype::reactive::RequestManager& exahype::reactive::RequestManager::getInstance() {
  //static RequestManager RequestManager;
  int threadID; // tarch::multicore::Core::getInstance().getThreadNum();
#if defined(UseMPIThreadSplit)
  threadID = tarch::multicore::Core::getInstance().getThreadNum();
#else
  threadID = 0;
#endif

  //logInfo("getInstance()", "getting offloading manager for thread id "<<threadID);

  if(threadID>=MAX_THREADS) {
	  logError("getInstance()","The application is using too many threads (>48), we need to exit...");
	  assertion(false);
	  MPI_Abort(MPI_COMM_WORLD, -1);
  }

  //Todo: need to deallocate somewhere
  if(_static_managers[threadID]==nullptr) {
    _static_managers[threadID] = new RequestManager(threadID);
  }
  return *_static_managers[threadID];
}

int exahype::reactive::RequestManager::requestTypeToMsgQueueIdx( RequestType requestType ) {
  return static_cast<int> (requestType);
}

void exahype::reactive::RequestManager::submitRequests(
    MPI_Request *requests,
    int nRequests,
    int tag,
    int remoteRank,
    std::function<void(exahype::solvers::Solver*, int, int)> handler,
    RequestType type,
    exahype::solvers::Solver *solver,
    bool block ) {
  
  assertion(remoteRank>=0);

  int ierr;
  //bug only appears when using scorep
  #ifdef ScoreP
  for(int i=0; i<nRequests; i++) {
	  assertion(requests[i]!=MPI_REQUEST_NULL);
  }
  #endif

  switch(type) {
    case RequestType::send:
      _postedSendsPerRank[remoteRank]++; break;
    case RequestType::receive:
      _postedReceivesPerRank[remoteRank]++; break;
    case RequestType::sendBack:
      _postedSendBacksPerRank[remoteRank]++; break;
    case RequestType::receiveBack:
      _postedReceiveBacksPerRank[remoteRank]++; break;
    case RequestType::receiveOutcome:
      _postedReceiveOutcomesPerRank[remoteRank]++; break;
    case RequestType::sendOutcome:
      _postedSendOutcomesPerRank[remoteRank]++; break;
  } 

  int finished = -1;

/*  for(int i=0;i<nRequests; i++) {
    assertion(requests[i]!=MPI_REQUEST_NULL); 
    int ierr = MPI_Test(&requests[i], &finished, MPI_STATUS_IGNORE); 
    if(ierr!=MPI_SUCCESS) {
      char err_buffer[MPI_MAX_ERROR_STRING];
      int resultlen = 0;
      MPI_Error_string(ierr,err_buffer,&resultlen);
      fprintf(stderr,err_buffer);
      fprintf(stderr, "request id %d\n", i);
    }

    assertion(ierr==MPI_SUCCESS);
    finished = -1;
  } */ 

  MPI_CHECK("submitRequests", MPI_Testall(nRequests, requests, &finished, MPI_STATUSES_IGNORE));

  assertion(ierr==MPI_SUCCESS);
  if(finished) {
    handler(solver, tag, remoteRank);
    return;
  }
   //check if ok after test
  #ifdef ScoreP
  for(int i=0; i<nRequests; i++) {
        assertion(requests[i]!=MPI_REQUEST_NULL);
  }
  #endif

  if(block) {
    MPI_CHECK("submitRequests", MPI_Waitall(nRequests, requests, MPI_STATUSES_IGNORE));
    assertion(ierr==MPI_SUCCESS);
    handler(solver, tag, remoteRank);
    return;
  }
   //check if ok after wait
  #ifdef ScoreP
  for(int i=0; i<nRequests; i++) {
        assertion(requests[i]!=MPI_REQUEST_NULL);
  }
  #endif

  int mapId = requestTypeToMsgQueueIdx(type);

  //submitted[mapId]++;
  logDebug("submitRequests","submitted["<<mapId<<"]:"<<nRequests);

  // assign group id for this request group
  int groupId = getNextGroupId();
  // insert metadata into maps
  std::pair<int, std::function<void(exahype::solvers::Solver*, int, int)>> handlerElem(groupId, handler);
  std::pair<int, exahype::solvers::Solver*> solverElem(groupId, solver);
  std::pair<int, int> outstandingElem(groupId, nRequests);
  std::pair<int, int> remoteRankElem(groupId, remoteRank);
  std::pair<int, int> remoteTagElem(groupId, tag);

  _handlers[mapId].insert(handlerElem);
  _solvers[mapId].insert(solverElem);
  _outstandingReqsForGroup[mapId].insert(outstandingElem);
  _groupIdToRank[mapId].insert(remoteRankElem);
  _groupIdToTag[mapId].insert(remoteTagElem);

   //check if ok before pushing to queue
  #ifdef ScoreP
  for(int i=0; i<nRequests; i++) {
        assertion(requests[i]!=MPI_REQUEST_NULL);
  }
  #endif

  // push requests into queue
  for(int i=0; i<nRequests; i++) {
    //TODO: avoid overflow!
    int id = getNextRequestId();

    std::pair<int, MPI_Request> reqElem(id, requests[i]);
    std::pair<int, int> reqGroupElem(id, groupId);
    //logInfo("RequestManager", "inserted groupid "<<groupId<<" for req "<<id);

    _reqIdToReqHandle[mapId].insert(reqElem);
    _reqIdToGroup[mapId].insert(reqGroupElem);
    _outstandingRequests[mapId].push(id);
  }

#ifdef OffloadingUseProgressTask
  if(_numProgressJobs==0 && type==RequestType::send) {
    //logInfo("submitRequests()", "spawning progress job (high priority)");
    _numProgressJobs++;
    ProgressJob *job = new ProgressJob();
    peano::datatraversal::TaskSet spawnedSet( job);
  }
#endif
}

void exahype::reactive::RequestManager::createRequestArray(
    RequestType type,
    std::vector<MPI_Request> &requests,
    std::unordered_map<int, int> &map,
    int limit) {

  int mapId = requestTypeToMsgQueueIdx(type);

  bool gotOne = true;
  int j = 0;

  while(gotOne && j<limit) {
    int req_id;
    gotOne = _outstandingRequests[mapId].try_pop(req_id);
    if(gotOne) {
      tbb::concurrent_hash_map<int, MPI_Request>::accessor a_requests;
      bool found = _reqIdToReqHandle[mapId].find(a_requests, req_id);
      if(!found)
        logError("createRequestArray", "didn't find MPI request.");
      assertion(found);
      MPI_Request request = a_requests->second;
      a_requests.release();
      requests.push_back(request);
      map.insert(std::pair<int, int>(j, req_id));
      j++;
    }
  }
}

#if defined(DirtyCleanUp)
void exahype::reactive::RequestManager::cancelOutstandingRequests() {
  std::vector<RequestType> types = {RequestType::send,
      RequestType::sendBack,
      RequestType::receive,
      RequestType::receiveBack,
      RequestType::sendOutcome,
      RequestType::receiveOutcome
  };
  tarch::multicore::Lock lock(_progressSemaphore, false);
  lock.lock();

  for(auto type : types) {
    while(hasOutstandingRequestOfType(type)) {
      int i = requestTypeToMsgQueueIdx(type);
      logInfo("cancelOutstandingRequests ", " cancelling i "<<i<<" size "<<_activeRequests[i].size());
      for(int j=0; j<_activeRequests[i].size(); j++) {
        if(_activeRequests[i][j]==MPI_REQUEST_NULL) continue;      

        int ierr = MPI_Cancel(&_activeRequests[i][j]);
        
        if(ierr!=MPI_SUCCESS) {
           char err_buffer[MPI_MAX_ERROR_STRING];
           int resultlen = 0;
         
           MPI_Error_string(ierr,err_buffer,&resultlen);
          fprintf(stderr,err_buffer);
        }
        assertion(ierr==MPI_SUCCESS);
        ierr = MPI_Request_free(&_activeRequests[i][j]);
        assertion(ierr==MPI_SUCCESS);
      }
      _activeRequests[i].clear();
      createRequestArray( type, _activeRequests[i], _internalIdsOfActiveRequests[i] );
    }
  }
  lock.free();
}
#endif

bool exahype::reactive::RequestManager::hasOutstandingRequestOfType(RequestType requestType) {
  logDebug("hasOutstandingRequestOfType",
           " type "
           <<int(requestType)
           <<" outstanding "
           <<_outstandingRequests[requestTypeToMsgQueueIdx(requestType)].unsafe_size()
           <<" active "<<_activeRequests[requestTypeToMsgQueueIdx(requestType)].size() );
  return (!_outstandingRequests[requestTypeToMsgQueueIdx(requestType)].empty() || !_activeRequests[requestTypeToMsgQueueIdx(requestType)].size()==0);
}

void exahype::reactive::RequestManager::progressRequests() {
  static double lastOutputTimeStamp = 0;

  if(lastOutputTimeStamp==0 || (MPI_Wtime()-lastOutputTimeStamp)>10) {
    lastOutputTimeStamp = MPI_Wtime();
    printPostedRequests();
    logDebug("progressRequests()", "there are "<<getNumberOfOutstandingRequests(RequestType::send)<< " send requests remaining "
        <<","<<getNumberOfOutstandingRequests(RequestType::receive)<<" receive requests remaining"
        <<","<<getNumberOfOutstandingRequests(RequestType::sendBack)<<" sendBack requests remaining"
        <<","<<getNumberOfOutstandingRequests(RequestType::receiveBack)<<" receiveBack requests remaining"
        <<","<<getNumberOfOutstandingRequests(RequestType::sendOutcome)<<" replica send requests remaining"
        <<","<<getNumberOfOutstandingRequests(RequestType::receiveOutcome)<<" replica receive requests remaining");
  }

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
  //Todo: add ITAC events if needed
  if (hasOutstandingRequestOfType(RequestType::receiveOutcome)) {
    progressRequestsOfType(RequestType::receiveOutcome);
  }
  if (hasOutstandingRequestOfType(RequestType::sendOutcome)) {
    progressRequestsOfType(RequestType::sendOutcome);
  }
}

void exahype::reactive::RequestManager::progressAnyRequests() {

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
  //always progress both sends and receives for replica tasks!
  //todo: add ITAC events if needed
  if (hasOutstandingRequestOfType(RequestType::receiveOutcome)) {
    progressRequestsOfType(RequestType::receiveOutcome);
  }
  if (hasOutstandingRequestOfType(RequestType::sendOutcome)) {
    progressRequestsOfType(RequestType::sendOutcome);
  }
}

bool exahype::reactive::RequestManager::progressReceiveBackRequests() {
  return progressRequestsOfType( RequestType::receiveBack );
}

bool exahype::reactive::RequestManager::progressRequestsOfType( RequestType type ) {
  // First, we ensure here that only one thread at a time progresses offloading
  // this attempts to avoid multithreaded MPI problems
  tarch::multicore::Lock lock(_progressSemaphore, false);
  bool canRun = lock.tryLock();
  if(!canRun) {
#if defined(PerformanceAnalysisOffloadingDetailed)
    watch.stopTimer();
    if(watch.getCalendarTime() >= 0.0) {
      logDebug(
          "progressRequestsOfType",
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

  int mapId = requestTypeToMsgQueueIdx(type);

  tbb::concurrent_hash_map<int, std::function<void(exahype::solvers::Solver*, int, int)>>::accessor a_handler;
  tbb::concurrent_hash_map<int, exahype::solvers::Solver*>::accessor a_solver;
  tbb::concurrent_hash_map<int, int>::accessor a_groupId;
  tbb::concurrent_hash_map<int, int>::accessor a_outstandingGroup;
  tbb::concurrent_hash_map<int, int>::accessor a_remoteRank;
  tbb::concurrent_hash_map<int, int>::accessor a_remoteTag;

  //std::vector<MPI_Request>     outstandingRequests;
  //std::unordered_map<int, int> vecIdToReqId;

  int nRequests = _activeRequests[mapId].size();
  if(nRequests==0) {
    //logInfo("progressRequestsOfType()", "begin create req array");
    if(type==RequestType::receiveBack || type==RequestType::sendOutcome)
      createRequestArray( type, _activeRequests[mapId], _internalIdsOfActiveRequests[mapId] );
      //createRequestArray( type, _activeRequests[mapId], _internalIdsOfActiveRequests[mapId], 10 );
    else {
      createRequestArray( type, _activeRequests[mapId], _internalIdsOfActiveRequests[mapId] );
    }
    nRequests = _activeRequests[mapId].size();
    //logInfo("progressRequestsOfType()", "end create req array");
  }

  std::vector<MPI_Status> stats(nRequests);
  std::vector<int> arrOfIndices(nRequests);
  int outcount = 0;

  double time = -MPI_Wtime();
  exahype::reactive::OffloadingProfiler::getInstance().beginCommunication();

  logDebug("progressRequestsOfType()"," type: "<<mapId<< " nreq "<<nRequests);

  int ierr = MPI_Testsome(nRequests,&_activeRequests[mapId][0], &outcount, &arrOfIndices[0], MPI_STATUSES_IGNORE);

  if(ierr != MPI_SUCCESS) {
    for(int i=0;i<nRequests;i++) {
      int ierrstatus = stats[i].MPI_ERROR;
      if(ierrstatus!=MPI_SUCCESS) {
        logError("progressRequestsOfType()", "error "<<ierrstatus<<" for request "<<_internalIdsOfActiveRequests[mapId][i]<< " source "<<stats[i].MPI_SOURCE<<" tag "<<stats[i].MPI_TAG);
      }
      char err_buffer[MPI_MAX_ERROR_STRING];
      int resultlen = 0;
      if(ierrstatus!=MPI_SUCCESS) {
        MPI_Error_string(ierrstatus, err_buffer, &resultlen);
        fprintf(stderr, "%s\n", err_buffer);
      }
    }
    MPI_Abort(MPI_COMM_WORLD, ierr); /* abort*/
  }
  assertion(ierr==MPI_SUCCESS);
  time += MPI_Wtime();

  if(outcount>0)
    exahype::reactive::OffloadingProfiler::getInstance().endCommunication(true, time);
  else
    exahype::reactive::OffloadingProfiler::getInstance().endCommunication(false, time);

  logDebug("progressRequestsOfType()"," finished type: "<<mapId<< " nreq "<<outcount);

  time = -MPI_Wtime();
  exahype::reactive::OffloadingProfiler::getInstance().beginHandling();
  bool found=false;
  //handle finished requests
  for(int i=0; i<outcount; i++) {
    int reqIdx = arrOfIndices[i];
    int reqId = _internalIdsOfActiveRequests[mapId][reqIdx];
    assertion(_activeRequests[mapId][reqIdx]==MPI_REQUEST_NULL);
    _internalIdsOfActiveRequests[mapId].erase(reqIdx);

    _reqIdToReqHandle[mapId].erase(reqId);
    int groupId;
    found = _reqIdToGroup[mapId].find(a_groupId, reqId);
    if(!found)
      logError("progressRequestsOfType()", "Didn't find MPI request...");

    assertion(found);
    groupId = a_groupId->second;
    _reqIdToGroup[mapId].erase(a_groupId);
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
      found= _groupIdToRank[mapId].find(a_remoteRank, groupId);
      remoteRank = a_remoteRank->second;
      assertion(found);
      a_remoteRank.release();
      int remoteTag = -1;
      found= _groupIdToTag[mapId].find(a_remoteTag, groupId);
      remoteTag = a_remoteTag->second;
      assertion(found);
      a_remoteTag.release();

      handler(solver, remoteTag, remoteRank);
      //RequestHandlerJob *requestHandlerJob = new RequestHandlerJob(handler, solver, remoteTag, remoteRank);
      //peano::datatraversal::TaskSet spawnedSet( requestHandlerJob, peano::datatraversal::TaskSet::TaskType::Background);

      _handlers[mapId].erase(groupId);
      _outstandingReqsForGroup[mapId].erase(groupId);
      _solvers[mapId].erase(groupId);
      _groupIdToRank[mapId].erase(groupId);
      _groupIdToTag[mapId].erase(groupId);
    }
    //Currently deactivated per default because it makes things slower in some cases (additional work on master
    //thread?)
    //May be useful in the future, if we can move work away from the master thread or in ExaHyPE2
#if defined(RequestManagerGrabNewRequests)
    //try to grab and add a new request to active requests (thanks to Joseph Schuchart for the suggestion)
    int mapId = requestTypeToMsgQueueIdx(type);

    int new_req_id;
    bool gotOne = _outstandingRequests[mapId].try_pop(new_req_id);

    if(gotOne) {
      tbb::concurrent_hash_map<int, MPI_Request>::accessor a_requests;
      bool found = _reqIdToReqHandle[mapId].find(a_requests, new_req_id);
      assertion(found);
      if(!found)
        logError("progressRequestsOfType", "didn't find MPI request.");
      if(found) {
        MPI_Request request = a_requests->second;
        a_requests.release();
        _activeRequests[mapId][reqIdx] = request; //replace finished request
        _internalIdsOfActiveRequests[mapId].insert(std::pair<int,int>(reqIdx, new_req_id));
        logDebug("progressRequestsOfType","Inserted "<<request<<"  at reqIdx "<<reqIdx<<"  internal id "<<new_req_id);
      }
    }
#endif
  }

  bool allFinished = true;
  for(int i=0; i<nRequests; i++) {
    if(_activeRequests[mapId][i]!=MPI_REQUEST_NULL) {
      allFinished = false;
      break;
    }
  }
  if(allFinished) {
    _activeRequests[mapId].clear();
    _internalIdsOfActiveRequests[mapId].clear();
  }

  time += MPI_Wtime();
  exahype::reactive::OffloadingProfiler::getInstance().endHandling(time);

  lock.free();

  return outcount>0;
}

exahype::reactive::RequestManager::RequestHandlerJob::RequestHandlerJob(
    std::function<void(exahype::solvers::Solver*, int, int)> handleRequest,
    exahype::solvers::Solver* solver,
    int tag,
    int remoteRank) :
  _handleRequest(handleRequest),
  _solver(solver),
  _tag(tag),
  _remoteRank(remoteRank)
{}

bool exahype::reactive::RequestManager::RequestHandlerJob::operator()() {
#ifdef USE_ITAC
  VT_begin(event_handling);
#endif
  _handleRequest(_solver, _tag, _remoteRank);
#ifdef USE_ITAC
  VT_end(event_handling);
#endif
  return false;
}

//todo: I think these can be removed as no longer used
exahype::reactive::RequestManager::ProgressJob::ProgressJob() :
tarch::multicore::jobs::Job(tarch::multicore::jobs::JobType::BackgroundTask, 0, tarch::multicore::DefaultPriority*8)
{}

bool exahype::reactive::RequestManager::ProgressJob::run( bool calledFromMaster ) {

  if(calledFromMaster) return false;

  int flag;
  //logInfo("submitRequests()", "executing progress job (high priority)");

  int mapId = RequestManager::requestTypeToMsgQueueIdx(RequestType::send);
//   while(RequestManager::getInstance()._requests[0].unsafe_size()>0 || RequestManager::getInstance()._currentOutstandingRequests[0].size()>0
//        || RequestManager::getInstance()._requests[1].unsafe_size()>0 || RequestManager::getInstance()._currentOutstandingRequests[1].size()>0
//        || RequestManager::getInstance()._requests[2].unsafe_size()>0 || RequestManager::getInstance()._currentOutstandingRequests[2].size()>0
//        || RequestManager::getInstance()._requests[3].unsafe_size()>0 || RequestManager::getInstance()._currentOutstandingRequests[3].size()>0
//   ) 
  while(RequestManager::getInstance()._outstandingRequests[mapId].unsafe_size()>0 || RequestManager::getInstance()._activeRequests[mapId].size()>0)
  {
    //getInstance().progressAnyRequests();
    getInstance().progressRequestsOfType(RequestType::send);
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, OffloadingContext::getInstance().getMPICommunicator(), &flag, MPI_STATUS_IGNORE);
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, OffloadingContext::getInstance().getMPICommunicatorMapped(), &flag, MPI_STATUS_IGNORE);
  }
  //logInfo("submitRequests()", "terminated progress job (high priority)");

  RequestManager::getInstance()._numProgressJobs--;
  return false;
}

exahype::reactive::RequestManager::ProgressSendJob::ProgressSendJob() {}

bool exahype::reactive::RequestManager::ProgressSendJob::operator()() {
  getInstance().progressRequestsOfType(RequestType::send);
  int mapId = RequestManager::requestTypeToMsgQueueIdx(RequestType::send);
//   bool reschedule=RequestManager::getInstance()._requests[mapId].unsafe_size()>0;

  while(RequestManager::getInstance()._outstandingRequests[mapId].unsafe_size()>0 || RequestManager::getInstance()._activeRequests[mapId].size()>0) {
    getInstance().progressRequestsOfType(RequestType::send);
  }

//   if(!reschedule)
  RequestManager::getInstance()._numProgressSendJobs--;
  return false;
}

exahype::reactive::RequestManager::ProgressReceiveJob::ProgressReceiveJob() {}

bool exahype::reactive::RequestManager::ProgressReceiveJob::operator()() {
  getInstance().progressRequestsOfType(RequestType::receive);
  int mapId = RequestManager::requestTypeToMsgQueueIdx(RequestType::receive);

  while(RequestManager::getInstance()._outstandingRequests[mapId].unsafe_size()>0 || RequestManager::getInstance()._activeRequests[mapId].size()>0) {
    getInstance().progressRequestsOfType(RequestType::receive);
  }

  RequestManager::getInstance()._numProgressReceiveJobs--;
  return false;
}

exahype::reactive::RequestManager::ProgressReceiveBackJob::ProgressReceiveBackJob() {}

bool exahype::reactive::RequestManager::ProgressReceiveBackJob::operator()() {
  getInstance().progressRequestsOfType(RequestType::receiveBack);
  int mapId = RequestManager::requestTypeToMsgQueueIdx(RequestType::receiveBack);

  while(RequestManager::getInstance()._outstandingRequests[mapId].unsafe_size()>0 || RequestManager::getInstance()._activeRequests[mapId].size()>0) {
    getInstance().progressRequestsOfType(RequestType::receiveBack);
  }

  RequestManager::getInstance()._numProgressReceiveBackJobs--;
  return false;
}

#endif

//#undef assertion
//#define assertion(expr) 
