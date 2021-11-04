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

#if defined(Parallel)

#include "exahype/reactive/ReactiveContext.h"

#include <unordered_set>
#include <vector>
#include <string>
#include <algorithm>
#include <mpi.h>

#include "tarch/multicore/Lock.h"
#include "tarch/multicore/Core.h"
#include "tarch/parallel/Node.h"

#include "exahype/reactive/OffloadingProfiler.h"
#include "exahype/reactive/StaticDistributor.h"
#include "exahype/reactive/AggressiveHybridDistributor.h"
#include "exahype/reactive/PerformanceMonitor.h"
#include "exahype/solvers/LimitingADERDGSolver.h"

#ifndef MPI_CHECK
#ifndef Asserts
#define MPI_CHECK(func, x) do { \
  ierr = (x); \
  if (ierr != MPI_SUCCESS) { \
    logError(#func, "Runtime error:"<<#x<<" returned "<<ierr<<" at " << __FILE__<< ":"<< __LINE__); \
  } \
} while (0)
#else
#define MPI_CHECK(func, x) do { \
  ierr = (x); \
  } while (0)
#endif
#endif

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

//#undef assertion
//#define assertion assert

tarch::logging::Log exahype::reactive::ReactiveContext::_log( "exahype::reactive::OffloadingManager" );

int exahype::reactive::ReactiveContext::_numTeams(-1);
int exahype::reactive::ReactiveContext::_interTeamRank(-1);

exahype::reactive::ReactiveContext* exahype::reactive::ReactiveContext::_static_managers[MAX_THREADS];
MPI_Comm exahype::reactive::ReactiveContext::_offloadingComms[MAX_THREADS];
MPI_Comm exahype::reactive::ReactiveContext::_offloadingCommsMapped[MAX_THREADS];

MPI_Comm exahype::reactive::ReactiveContext::_interTeamComms[MAX_THREADS];

exahype::reactive::ReactiveContext::OffloadingStrategy exahype::reactive::ReactiveContext::ChosenOffloadingStrategy(exahype::reactive::ReactiveContext::OffloadingStrategy::None);
exahype::reactive::ReactiveContext::ResilienceStrategy exahype::reactive::ReactiveContext::ChosenResilienceStrategy(exahype::reactive::ReactiveContext::ResilienceStrategy::None);
bool exahype::reactive::ReactiveContext::SaveRedundantComputations(false);
bool exahype::reactive::ReactiveContext::MakeSkeletonsShareable(false);

exahype::reactive::ReactiveContext::ReactiveContext(int threadId) :
    _threadId(threadId),
    _emergencyTriggered(false),
    _hasNotifiedSendCompleted(false) {

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

  //todo: should not be implicitly invoked
  initializeCommunicatorsAndTeamMetadata();

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes()*_numTeams;
  _localBlacklist = new double[nnodes];
  std::fill(&_localBlacklist[0], &_localBlacklist[nnodes], 0);

  int *tag_ptr;
  int flag = 0;
  MPI_CHECK("OffloadingManager()" , MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ptr, &flag));
  _maxTag = *tag_ptr;
  assertion(ierr==MPI_SUCCESS);

  if(!flag) {
    logWarning("OffloadingManager()"," maximum allowed MPI tag could not be determined. Offloading may leave the space of allowed tags for longer runs...");
  }

}

exahype::reactive::ReactiveContext::~ReactiveContext() {
  destroy();
  delete[] _localBlacklist;
}

void exahype::reactive::ReactiveContext::setOffloadingStrategy(OffloadingStrategy strategy) {
  ChosenOffloadingStrategy = strategy;
}

exahype::reactive::ReactiveContext::OffloadingStrategy exahype::reactive::ReactiveContext::getOffloadingStrategy() {
  return ChosenOffloadingStrategy;
}

void exahype::reactive::ReactiveContext::setResilienceStrategy(ResilienceStrategy strategy) {
  ChosenResilienceStrategy = strategy;
}

exahype::reactive::ReactiveContext::ResilienceStrategy exahype::reactive::ReactiveContext::getResilienceStrategy() {
  return ChosenResilienceStrategy;
}

void exahype::reactive::ReactiveContext::setSaveRedundantComputations(bool saveRedundantComputations){
  SaveRedundantComputations = saveRedundantComputations;
}

bool exahype::reactive::ReactiveContext::getSaveRedundantComputations() {
  return SaveRedundantComputations;
}

void exahype::reactive::ReactiveContext::setMakeSkeletonsShareable(bool makeSkeletonsShareable){
  MakeSkeletonsShareable = makeSkeletonsShareable;
}

bool exahype::reactive::ReactiveContext::getMakeSkeletonsShareable() {
  return MakeSkeletonsShareable;
}

bool exahype::reactive::ReactiveContext::isEnabled() {
  return (ChosenOffloadingStrategy!=OffloadingStrategy::None) || (ChosenResilienceStrategy!=ResilienceStrategy::None);
}

bool exahype::reactive::ReactiveContext::usesOffloading() {
  return ChosenOffloadingStrategy!=OffloadingStrategy::None;
}

void exahype::reactive::ReactiveContext::initializeCommunicatorsAndTeamMetadata() {
  static bool initialized = false;
  //Todo:  this collective routine should be exposed to the runner instead of being implicitly invoked by the constructor
  //#if defined(SharedTBB)
  if(!initialized)
    createMPICommunicators();
  //#endif
  initialized = true;
}

void exahype::reactive::ReactiveContext::destroy() {
#if defined(SharedTBB)
  destroyMPICommunicators();
#endif
}

MPI_Comm exahype::reactive::ReactiveContext::getTMPIInterTeamCommunicatorData() {
  return _interTeamComms[_threadId];
}

void exahype::reactive::ReactiveContext::setTMPINumTeams(int nteams) {
  _numTeams = nteams;
}

unsigned int exahype::reactive::ReactiveContext::getTMPINumTeams() {
  return _numTeams;
}

void exahype::reactive::ReactiveContext::setTMPITeamNumber(int interTeamRank) {
  _interTeamRank = interTeamRank;
}

int exahype::reactive::ReactiveContext::getTMPITeamNumber(){
  return _interTeamRank;
}

void exahype::reactive::ReactiveContext::createMPICommunicators() {
  int rank;
  int ierr;
  MPI_CHECK("createMPICommunicators", MPI_Comm_rank(MPI_COMM_WORLD, &rank));
  MPI_Comm interTeamComm = MPI_COMM_NULL;

#if defined(USE_TMPI)
  int worldRank = TMPI_GetWorldRank();
  MPI_CHECK("createMPICommunicators", MPI_Comm_split(MPI_COMM_WORLD, rank, worldRank, &interTeamComm));

  int nteams = TMPI_GetInterTeamCommSize();
   //MPI_Comm interTeamComm = TMPI_GetInterTeamComm();

  int rankInterComm;
  MPI_Comm_rank(interTeamComm, &rankInterComm);
  int team = TMPI_GetTeamNumber();

  _numTeams  = nteams;
  _interTeamRank = rankInterComm;

  //exahype::reactive::OffloadingManager::getInstance().setTMPIInterTeamCommunicators(interTeamComm, interTeamCommDupKey, interTeamCommDupAck);

  logInfo("createMPICommunicators()", " teams: "<<nteams<<", rank in team "
                                                <<team<<" : "<<rank<<", team rank in intercomm: "<<rankInterComm);

#else
  _numTeams = 1;
  _interTeamRank = 0;

  MPI_CHECK("createMPICommunicators", MPI_Comm_rank(MPI_COMM_WORLD, &rank));
  MPI_CHECK("createMPICommunicators", MPI_Comm_split(MPI_COMM_WORLD, rank, rank, &interTeamComm));
#endif

  for(int i=0; i<MAX_THREADS;i++) {
    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Info_set(info, "thread_id", std::to_string(i).c_str());

    MPI_CHECK("createMPICommunicators", MPI_Comm_dup(MPI_COMM_WORLD, &_offloadingComms[i]));
    assertion(ierr==MPI_SUCCESS);
    MPI_CHECK("createMPICommunicators", MPI_Comm_set_info(_offloadingComms[i], info));
    assertion(ierr==MPI_SUCCESS);

    MPI_CHECK("createMPICommunicators", MPI_Comm_dup(MPI_COMM_WORLD, &_offloadingCommsMapped[i]));
    assertion(ierr==MPI_SUCCESS);
    MPI_CHECK("createMPICommunicators", MPI_Comm_set_info(_offloadingCommsMapped[i], info));
    assertion(ierr==MPI_SUCCESS);

    MPI_CHECK("createMPICommunicators", MPI_Comm_dup(interTeamComm, &_interTeamComms[i])); assertion(ierr==MPI_SUCCESS);
    MPI_CHECK("createMPICommunicators", MPI_Comm_set_info(_interTeamComms[i], info));
    assertion(ierr==MPI_SUCCESS);

    MPI_Info_free(&info);
  }
  MPI_CHECK("createMPICommunicators", MPI_Comm_free(&interTeamComm));
}

void exahype::reactive::ReactiveContext::destroyMPICommunicators() {
  int ierr;
  for(int i=0; i<MAX_THREADS;i++) {
    MPI_CHECK("destroyMPICommunicators", MPI_Comm_free(&_offloadingComms[i])); assertion(ierr==MPI_SUCCESS);
    MPI_CHECK("destroyMPICommunicators", MPI_Comm_free(&_offloadingCommsMapped[i])); assertion(ierr==MPI_SUCCESS);
    MPI_CHECK("destroyMPICommunicators", MPI_Comm_free(&_interTeamComms[i])); assertion(ierr==MPI_SUCCESS);
  }
}

MPI_Comm exahype::reactive::ReactiveContext::getMPICommunicator() {
  return _offloadingComms[_threadId];
}

MPI_Comm exahype::reactive::ReactiveContext::getMPICommunicatorMapped() {
  return _offloadingCommsMapped[_threadId];
}

exahype::reactive::ReactiveContext& exahype::reactive::ReactiveContext::getInstance() {
  int threadID;
#if defined(UseMPIThreadSplit)
  threadID = tarch::multicore::Core::getInstance().getThreadNum();
#else
  threadID = 0;
#endif

  if(threadID>=MAX_THREADS) {
	  logError("getInstance()","The application is using too many threads (>48), we need to exit... You may want to increase MAX_THREADS.");
	  assertion(false);
	  MPI_Abort(MPI_COMM_WORLD, -1);
  }

  if(_static_managers[threadID]==nullptr) {
    _static_managers[threadID] = new ReactiveContext(threadID);
  }
  return *_static_managers[threadID];
}


int exahype::reactive::ReactiveContext::getOffloadingTag() {
  static std::atomic<int> next_tag(1); //0 is reserved for status
retry:
  int val =  next_tag.load();
  int val_o = val; 

  if(val>_maxTag-1) {
    logWarning("getOffloadingTag","MPI tag rollover for reactive communication!");
    val = 1;
  }

  int update = val + 1;

  if(!std::atomic_compare_exchange_strong(&next_tag, &val_o, update)) {
    goto retry;
  }
  return val;
}

void exahype::reactive::ReactiveContext::triggerVictimFlag() {
  _isVictim = true;
}

void exahype::reactive::ReactiveContext::resetVictimFlag() {
  _isVictim = false;
}

bool exahype::reactive::ReactiveContext::isVictim() {
  return _isVictim;
}

void exahype::reactive::ReactiveContext::triggerEmergencyForRank(int rank) {
  switch(exahype::reactive::ReactiveContext::getInstance().getOffloadingStrategy()){
    case ReactiveContext::OffloadingStrategy::AggressiveHybrid:
      exahype::reactive::AggressiveHybridDistributor::getInstance().handleEmergencyOnRank(rank); break;
    default:
      //do nothing, no emergencies for other offloading strategies so far
      break;
  }

  _localBlacklist[rank]++;
  exahype::reactive::PerformanceMonitor::getInstance().submitBlacklistValueForRank(_localBlacklist[rank], rank);
}

void exahype::reactive::ReactiveContext::recoverBlacklistedRanks() {
  logDebug("recoverBlacklistedRanks()","decrease heat of emergency heat map");
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  for(int i=0; i<nnodes;i++) {
    _localBlacklist[i]*= 0.9;
    if(_localBlacklist[i]>0)
    exahype::reactive::PerformanceMonitor::getInstance().submitBlacklistValueForRank(_localBlacklist[i], i);
  }
}

bool exahype::reactive::ReactiveContext::isBlacklisted(int rank) {
  const double* globalHeatMap = exahype::reactive::PerformanceMonitor::getInstance().getBlacklistGlobalSnapshot();
  return (globalHeatMap[rank]>0.5) || (_localBlacklist[rank]>0.5);
}

void exahype::reactive::ReactiveContext::printBlacklist() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  const double* globalHeatMap = exahype::reactive::PerformanceMonitor::getInstance().getBlacklistGlobalSnapshot();

  for(int i=0; i<nnodes; i++) {
    if(globalHeatMap[i]>0 || _localBlacklist[i]>0)
    logInfo("printBlacklist", "blacklist value for rank "<<i<<":"<<globalHeatMap[i]);
  }
}

bool exahype::reactive::ReactiveContext::isEmergencyTriggered() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  return !std::all_of(&_localBlacklist[0], &_localBlacklist[nnodes], [](double d) {return d<0.5;});
}

bool exahype::reactive::ReactiveContext::isEmergencyTriggeredOnRank(int rank) {
  return !(_localBlacklist[rank]<0.5);
}

bool exahype::reactive::ReactiveContext::selectVictimRank(int& victim, bool& last) {
  last = false;
  if(ChosenOffloadingStrategy==OffloadingStrategy::StaticHardcoded)
    return exahype::reactive::StaticDistributor::getInstance().selectVictimRank(victim);
  if(ChosenOffloadingStrategy==OffloadingStrategy::AggressiveHybrid)
    return exahype::reactive::AggressiveHybridDistributor::getInstance().selectVictimRank(victim, last);

  return false;
}

#endif

//#undef assertion
//#define assertion(expr) 
