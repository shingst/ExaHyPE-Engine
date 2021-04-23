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

#include "OffloadingContext.h"

#include <unordered_set>
#include <vector>
#include <string>
#include <algorithm>

#include "tarch/multicore/Lock.h"
#include "tarch/multicore/Core.h"
#include "tarch/parallel/Node.h"

#include "../reactive/OffloadingProfiler.h"
#include "../reactive/StaticDistributor.h"
#include "../reactive/DynamicDistributor.h"
#include "../reactive/DiffusiveDistributor.h"
#include "../reactive/AggressiveDistributor.h"
#include "../reactive/AggressiveCCPDistributor.h"
#include "../reactive/AggressiveHybridDistributor.h"
#include "../reactive/PerformanceMonitor.h"
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

#include <mpi.h>

//#undef assertion
//#define assertion assert

tarch::logging::Log exahype::reactive::OffloadingContext::_log( "exahype::reactive::OffloadingManager" );

int exahype::reactive::OffloadingContext::_numTeams(-1);
int exahype::reactive::OffloadingContext::_interTeamRank(-1);

exahype::reactive::OffloadingContext* exahype::reactive::OffloadingContext::_static_managers[MAX_THREADS];
MPI_Comm exahype::reactive::OffloadingContext::_offloadingComms[MAX_THREADS];
MPI_Comm exahype::reactive::OffloadingContext::_offloadingCommsMapped[MAX_THREADS];

MPI_Comm exahype::reactive::OffloadingContext::_interTeamComms[MAX_THREADS];
MPI_Comm exahype::reactive::OffloadingContext::_interTeamCommsKey[MAX_THREADS];
MPI_Comm exahype::reactive::OffloadingContext::_interTeamCommsAck[MAX_THREADS];

exahype::reactive::OffloadingContext::OffloadingStrategy exahype::reactive::OffloadingContext::_offloadingStrategy(exahype::reactive::OffloadingContext::OffloadingStrategy::None);
exahype::reactive::OffloadingContext::ResilienceStrategy exahype::reactive::OffloadingContext::_resilienceStrategy(exahype::reactive::OffloadingContext::ResilienceStrategy::None);

exahype::reactive::OffloadingContext::OffloadingContext(int threadId) :
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

exahype::reactive::OffloadingContext::~OffloadingContext() {
  destroy();
  delete[] _localBlacklist;
}

void exahype::reactive::OffloadingContext::setOffloadingStrategy(OffloadingStrategy strategy) {
  _offloadingStrategy = strategy;
}

exahype::reactive::OffloadingContext::OffloadingStrategy exahype::reactive::OffloadingContext::getOffloadingStrategy() {
  return _offloadingStrategy;
}

void exahype::reactive::OffloadingContext::setResilienceStrategy(ResilienceStrategy strategy) {
  _resilienceStrategy = strategy;
}

exahype::reactive::OffloadingContext::ResilienceStrategy exahype::reactive::OffloadingContext::getResilienceStrategy() {
  return _resilienceStrategy;
}

bool exahype::reactive::OffloadingContext::isEnabled() {
  return _offloadingStrategy!=OffloadingStrategy::None || _resilienceStrategy!=ResilienceStrategy::None;
}

void exahype::reactive::OffloadingContext::initializeCommunicatorsAndTeamMetadata() {
  //exahype::reactive::OffloadingManager::getInstance().createMPICommunicator();
 static bool initialized = false;
//Todo:  this collective routine should be exposed to the runner instead of being implicitly invoked by the constructor 
// can cause deadlocks if only a single thread calls the first getInstance() for getting an OffloadingContext
//#if defined(SharedTBB)
 if(!initialized)
  createMPICommunicators();
//#endif

 initialized = true;
}

void exahype::reactive::OffloadingContext::destroy() {
#if defined(SharedTBB)
  destroyMPICommunicators();
#endif
}

MPI_Comm exahype::reactive::OffloadingContext::getTMPIInterTeamCommunicatorData() {
  return _interTeamComms[_threadId];
}

MPI_Comm exahype::reactive::OffloadingContext::getTMPIInterTeamCommunicatorKey() {
  return _interTeamCommsKey[_threadId];
}

MPI_Comm exahype::reactive::OffloadingContext::getTMPIInterTeamCommunicatorAck() {
  return _interTeamCommsAck[_threadId];
}

void exahype::reactive::OffloadingContext::setTMPINumTeams(int nteams) {
  _numTeams = nteams;
}

int exahype::reactive::OffloadingContext::getTMPINumTeams() {
  return _numTeams;
}

void exahype::reactive::OffloadingContext::setTMPIInterTeamRank(int interTeamRank) {
  _interTeamRank = interTeamRank;
}

int exahype::reactive::OffloadingContext::getTMPIInterTeamRank(){
  return _interTeamRank;
}

void exahype::reactive::OffloadingContext::createMPICommunicators() {
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

    MPI_CHECK("createMPICommunicators", MPI_Comm_dup(interTeamComm, &_interTeamCommsKey[i])); assertion(ierr==MPI_SUCCESS);
    MPI_CHECK("createMPICommunicators", MPI_Comm_set_info(_interTeamCommsKey[i], info));
    assertion(ierr==MPI_SUCCESS);

    MPI_CHECK("createMPICommunicators", MPI_Comm_dup(interTeamComm, &_interTeamCommsAck[i])); assertion(ierr==MPI_SUCCESS);
    MPI_CHECK("createMPICommunicators", MPI_Comm_set_info(_interTeamCommsAck[i], info));
    assertion(ierr==MPI_SUCCESS);

    MPI_Info_free(&info);
  }
  MPI_CHECK("createMPICommunicators", MPI_Comm_free(&interTeamComm));
}

void exahype::reactive::OffloadingContext::destroyMPICommunicators() {
  int ierr;
  for(int i=0; i<MAX_THREADS;i++) {
    MPI_CHECK("destroyMPICommunicators", MPI_Comm_free(&_offloadingComms[i])); assertion(ierr==MPI_SUCCESS);
    MPI_CHECK("destroyMPICommunicators", MPI_Comm_free(&_offloadingCommsMapped[i])); assertion(ierr==MPI_SUCCESS);

    MPI_CHECK("destroyMPICommunicators", MPI_Comm_free(&_interTeamComms[i])); assertion(ierr==MPI_SUCCESS);
    MPI_CHECK("destroyMPICommunicators", MPI_Comm_free(&_interTeamCommsKey[i])); assertion(ierr==MPI_SUCCESS);
    MPI_CHECK("destroyMPICommunicators", MPI_Comm_free(&_interTeamCommsAck[i])); assertion(ierr==MPI_SUCCESS);
  }
}

MPI_Comm exahype::reactive::OffloadingContext::getMPICommunicator() {
  return _offloadingComms[_threadId];
}

MPI_Comm exahype::reactive::OffloadingContext::getMPICommunicatorMapped() {
  return _offloadingCommsMapped[_threadId];
}

exahype::reactive::OffloadingContext& exahype::reactive::OffloadingContext::getInstance() {
  //static OffloadingManager offloadingManager;
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

  if(_static_managers[threadID]==nullptr) {
    _static_managers[threadID] = new OffloadingContext(threadID);
  }
  return *_static_managers[threadID];
}


int exahype::reactive::OffloadingContext::getOffloadingTag() {
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

void exahype::reactive::OffloadingContext::triggerVictimFlag() {
  _isVictim = true;
}

void exahype::reactive::OffloadingContext::resetVictimFlag() {
  _isVictim = false;
}

bool exahype::reactive::OffloadingContext::isVictim() {
  return _isVictim;
}

void exahype::reactive::OffloadingContext::triggerEmergencyForRank(int rank) {
//  if(!_emergencyTriggered) {
//    logInfo("triggerEmergency()","emergency event triggered");
//    _emergencyTriggered = true;
//  }
  switch(exahype::reactive::OffloadingContext::getInstance().getOffloadingStrategy()){
    case OffloadingContext::OffloadingStrategy::Aggressive:
      exahype::reactive::AggressiveDistributor::getInstance().handleEmergencyOnRank(rank); break;
    case OffloadingContext::OffloadingStrategy::AggressiveCCP:
      exahype::reactive::AggressiveCCPDistributor::getInstance().handleEmergencyOnRank(rank); break;
    case OffloadingContext::OffloadingStrategy::AggressiveHybrid:
      exahype::reactive::AggressiveHybridDistributor::getInstance().handleEmergencyOnRank(rank); break;
    case OffloadingContext::OffloadingStrategy::Diffusive:
      exahype::reactive::DiffusiveDistributor::getInstance().handleEmergencyOnRank(rank); break;
    default:
      //do nothing, not implemented yet
      break;
  }

  _localBlacklist[rank]++;
  exahype::reactive::PerformanceMonitor::getInstance().submitBlacklistValueForRank(_localBlacklist[rank], rank);
  //logInfo("triggerEmergencyForRank()","blacklist value for rank "<<rank<<":"<<_localBlacklist[rank]
	//                                 <<" NumberOfRemoteJobs"<<  exahype::solvers::ADERDGSolver::NumberOfRemoteJobs
	//								                 <<" NumberOfEnclaveJobs"<<  exahype::solvers::ADERDGSolver::NumberOfEnclaveJobs);
}

void exahype::reactive::OffloadingContext::recoverBlacklistedRanks() {
  logDebug("decreaseHeat()","decrease heat of emergency heat map");
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  for(int i=0; i<nnodes;i++) {
    _localBlacklist[i]*= 0.9;
    if(_localBlacklist[i]>0)
    exahype::reactive::PerformanceMonitor::getInstance().submitBlacklistValueForRank(_localBlacklist[i], i);
    //if(_emergencyHeatMap[i]>0.5) {
    //  logInfo("decreaseHeat()","blacklist value for rank "<<i<<":"<<_emergencyHeatMap[i]);
    //}
  }
}

bool exahype::reactive::OffloadingContext::isBlacklisted(int rank) {
  //return _emergencyHeatMap[rank]>0.5;
  const double* globalHeatMap = exahype::reactive::PerformanceMonitor::getInstance().getBlacklistSnapshot();
  return (globalHeatMap[rank]>0.5) || (_localBlacklist[rank]>0.5);
}

void exahype::reactive::OffloadingContext::printBlacklist() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  const double* globalHeatMap = exahype::reactive::PerformanceMonitor::getInstance().getBlacklistSnapshot();

  for(int i=0; i<nnodes; i++) {
    if(globalHeatMap[i]>0 || _localBlacklist[i]>0)
    logInfo("printBlacklist", "blacklist value for rank "<<i<<":"<<globalHeatMap[i]);
  }
}

bool exahype::reactive::OffloadingContext::isEmergencyTriggered() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  return !std::all_of(&_localBlacklist[0], &_localBlacklist[nnodes], [](double d) {return d<0.5;});
}

bool exahype::reactive::OffloadingContext::isEmergencyTriggeredOnRank(int rank) {
  return !(_localBlacklist[rank]<0.5);
}

bool exahype::reactive::OffloadingContext::selectVictimRank(int& victim, bool& last) {
  last = false;
  if(_offloadingStrategy==OffloadingStrategy::StaticHardcoded)
    return exahype::reactive::StaticDistributor::getInstance().selectVictimRank(victim);
  if(_offloadingStrategy==OffloadingStrategy::Diffusive)
    return exahype::reactive::DiffusiveDistributor::getInstance().selectVictimRank(victim);
  if(_offloadingStrategy==OffloadingStrategy::Aggressive)
    return exahype::reactive::AggressiveDistributor::getInstance().selectVictimRank(victim);
  if(_offloadingStrategy==OffloadingStrategy::AggressiveCCP)
    return exahype::reactive::AggressiveCCPDistributor::getInstance().selectVictimRank(victim);
  if(_offloadingStrategy==OffloadingStrategy::AggressiveHybrid)
    return exahype::reactive::AggressiveHybridDistributor::getInstance().selectVictimRank(victim, last);

  return false;
/*  double remainingLoadRatio = static_cast<double> (exahype::reactive::PerformanceMonitor::getInstance().getRemainingTasks())
  /
  exahype::reactive::PerformanceMonitor::getInstance().getTasksPerTimestep();
  // this is currently hardcoded: the goal is to refrain from giving tasks away if there is not enough work left
  // for overlap of communication and computation
  if(remainingLoadRatio>0.1) {
#if defined(OffloadingStrategyStatic)
    return exahype::reactive::StaticDistributor::getInstance().selectVictimRank(victim);
#elif defined(OffloadingStrategyDynamic)
    return exahype::reactive::DynamicDistributor::getInstance().selectVictimRank(victim);
#elif defined(OffloadingStrategyHybrid)
    bool staticDistribution = exahype::reactive::StaticDistributor::getInstance().selectVictimRank(victim);
    if(staticDistribution) return true;
    else {
      return exahype::reactive::DynamicDistributor::getInstance().selectVictimRank(victim);
    }
    return false;
#else
    return false;
#endif
  }
  else {
    logDebug("offloadingManager", "could not select victim remaining load ratio "<<remainingLoadRatio);
    return false;
  }*/
}

#ifdef OffloadingUseProgressTask
void exahype::reactive::OffloadingContext::resetHasNotifiedSendCompleted() {
  _hasNotifiedSendCompleted = false;
  logDebug("resetHasNotifiedSendCompleted","resetting flag");
}

void exahype::reactive::OffloadingContext::notifySendCompleted(int rank) {
  char send = 1;
  MPI_Send(&send, 1, MPI_CHAR, rank, 0, getMPICommunicator());
  //logInfo("notifySendCompleted()","sent status message to "<<rank);
}

void exahype::reactive::OffloadingContext::receiveCompleted(int rank, int rail) {
  char receive = 0;
#if defined (UseSmartMPI)
  MPI_Status_Offload stat;
  MPI_Recv_offload(&receive, 1, MPI_CHAR, rank, 0, getMPICommunicator(), &stat, rail);
#else
  MPI_Recv(&receive, 1, MPI_CHAR, rank, 0, getMPICommunicator(), MPI_STATUS_IGNORE);
#endif
  //logInfo("receiveCompleted()","received status message from "<<rank);
}

void exahype::reactive::OffloadingContext::notifyAllVictimsSendCompletedIfNotNotified() {
  if(!_hasNotifiedSendCompleted) {
    //logInfo("notifyAllVictimsSendCompleted","notifying that last job was sent to victims");
    _hasNotifiedSendCompleted = true;
    std::vector<int> victimRanks;
    if(_offloadingStrategy==OffloadingStrategy::AggressiveHybrid)
      exahype::reactive::AggressiveHybridDistributor::getInstance().getAllVictimRanks(victimRanks);
    if(_offloadingStrategy==OffloadingStrategy::Static)
      exahype::reactive::StaticDistributor::getInstance().getAllVictimRanks(victimRanks);

    for(auto victim : victimRanks)
      notifySendCompleted(victim);
  }
}
#endif

#endif

//#undef assertion
//#define assertion(expr) 
