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

tarch::logging::Log exahype::reactive::ReactiveContext::_log( "exahype::reactive::ReactiveContext" );

int exahype::reactive::ReactiveContext::NumTeams(-1);
int exahype::reactive::ReactiveContext::Team(-1);

int exahype::reactive::ReactiveContext::MaxSupportedTag(-1);

exahype::reactive::ReactiveContext* exahype::reactive::ReactiveContext::StaticManagers[MAX_THREADS];
MPI_Comm exahype::reactive::ReactiveContext::OffloadingComms[MAX_THREADS];
MPI_Comm exahype::reactive::ReactiveContext::OffloadingCommsMapped[MAX_THREADS];

MPI_Comm exahype::reactive::ReactiveContext::InterTeamComms[MAX_THREADS];

exahype::reactive::ReactiveContext::OffloadingStrategy exahype::reactive::ReactiveContext::ChosenOffloadingStrategy(exahype::reactive::ReactiveContext::OffloadingStrategy::None);
exahype::reactive::ReactiveContext::ResilienceStrategy exahype::reactive::ReactiveContext::ChosenResilienceStrategy(exahype::reactive::ReactiveContext::ResilienceStrategy::None);

bool exahype::reactive::ReactiveContext::SaveRedundantComputations(false);
bool exahype::reactive::ReactiveContext::MakeSkeletonsShareable(false);
double exahype::reactive::ReactiveContext::ResilienceChecksTimeout(300);


std::atomic<bool> exahype::reactive::ReactiveContext::IsVictim(false);

exahype::reactive::ReactiveContext::ReactiveContext(int threadId) :
    _threadId(threadId),
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

  initialize();

  int *tag_ptr;
  int flag = 0;
  MPI_CHECK("OffloadingManager()", MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ptr, &flag));
  MaxSupportedTag = *tag_ptr;
  assertion(ierr==MPI_SUCCESS);

  if(!flag) {
    logWarning("OffloadingManager()"," maximum allowed MPI tag could not be determined. Offloading may leave the space of allowed tags for longer runs...");
  }
}

exahype::reactive::ReactiveContext::~ReactiveContext() {
  destroy();
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

double exahype::reactive::ReactiveContext::getResilienceChecksTimeout() {
  return ResilienceChecksTimeout;
}

void exahype::reactive::ReactiveContext::setResilienceChecksTimeout(double timeout) {
  ResilienceChecksTimeout = timeout;
}

bool exahype::reactive::ReactiveContext::isReactivityEnabled() {
  return (ChosenOffloadingStrategy!=OffloadingStrategy::None) || (ChosenResilienceStrategy!=ResilienceStrategy::None);
}

bool exahype::reactive::ReactiveContext::isReactiveOffloadingEnabled() {
  return ChosenOffloadingStrategy!=OffloadingStrategy::None;
}

void exahype::reactive::ReactiveContext::triggerVictimFlag() {
  IsVictim = true;
}

void exahype::reactive::ReactiveContext::resetVictimFlag() {
  IsVictim = false;
}

bool exahype::reactive::ReactiveContext::isVictim() {
  return IsVictim;
}


void exahype::reactive::ReactiveContext::initialize() {
  static bool initialized = false;
  //Todo:  this collective routine should probably be exposed to the runner instead of being implicitly invoked by the constructor
  //#if defined(SharedTBB)
  if(!initialized)
    createMPICommunicators();
  //#endif
  initialized = true;
}

void exahype::reactive::ReactiveContext::destroy() {
  //#if defined(SharedTBB)
  destroyMPICommunicators();
  //#endif
}

MPI_Comm exahype::reactive::ReactiveContext::getTMPIInterTeamCommunicatorData() const {
  return InterTeamComms[_threadId];
}

void exahype::reactive::ReactiveContext::setTMPINumTeams(int nteams) {
  NumTeams = nteams;
}

unsigned int exahype::reactive::ReactiveContext::getTMPINumTeams() {
  return NumTeams;
}

void exahype::reactive::ReactiveContext::setTMPITeamNumber(int interTeamRank) {
  Team = interTeamRank;
}

int exahype::reactive::ReactiveContext::getTMPITeamNumber() {
  return Team;
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

  NumTeams  = nteams;
  Team = rankInterComm;

  logInfo("createMPICommunicators()", " teams: "<<nteams<<", rank in team "
                                                <<team<<" : "<<rank<<", team rank in intercomm: "<<rankInterComm);

#else
  NumTeams = 1;
  Team = 0;

  MPI_CHECK("createMPICommunicators", MPI_Comm_rank(MPI_COMM_WORLD, &rank));
  MPI_CHECK("createMPICommunicators", MPI_Comm_split(MPI_COMM_WORLD, rank, rank, &interTeamComm));
#endif

  for(int i=0; i<MAX_THREADS;i++) {
    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Info_set(info, "thread_id", std::to_string(i).c_str());

    MPI_CHECK("createMPICommunicators", MPI_Comm_dup(MPI_COMM_WORLD, &OffloadingComms[i]));
    assertion(ierr==MPI_SUCCESS);
    MPI_CHECK("createMPICommunicators", MPI_Comm_set_info(OffloadingComms[i], info));
    assertion(ierr==MPI_SUCCESS);

    MPI_CHECK("createMPICommunicators", MPI_Comm_dup(MPI_COMM_WORLD, &OffloadingCommsMapped[i]));
    assertion(ierr==MPI_SUCCESS);
    MPI_CHECK("createMPICommunicators", MPI_Comm_set_info(OffloadingCommsMapped[i], info));
    assertion(ierr==MPI_SUCCESS);

    MPI_CHECK("createMPICommunicators", MPI_Comm_dup(interTeamComm, &InterTeamComms[i])); assertion(ierr==MPI_SUCCESS);
    MPI_CHECK("createMPICommunicators", MPI_Comm_set_info(InterTeamComms[i], info));
    assertion(ierr==MPI_SUCCESS);

    MPI_Info_free(&info);
  }
  MPI_CHECK("createMPICommunicators", MPI_Comm_free(&interTeamComm));
}

void exahype::reactive::ReactiveContext::destroyMPICommunicators() {
  int ierr;
  for(int i=0; i<MAX_THREADS;i++) {
    MPI_CHECK("destroyMPICommunicators", MPI_Comm_free(&OffloadingComms[i])); assertion(ierr==MPI_SUCCESS);
    MPI_CHECK("destroyMPICommunicators", MPI_Comm_free(&OffloadingCommsMapped[i])); assertion(ierr==MPI_SUCCESS);
    MPI_CHECK("destroyMPICommunicators", MPI_Comm_free(&InterTeamComms[i])); assertion(ierr==MPI_SUCCESS);
  }
}

MPI_Comm exahype::reactive::ReactiveContext::getMPICommunicator() const {
  return OffloadingComms[_threadId];
}

MPI_Comm exahype::reactive::ReactiveContext::getMPICommunicatorMapped() const {
  return OffloadingCommsMapped[_threadId];
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

  if(StaticManagers[threadID]==nullptr) {
    StaticManagers[threadID] = new ReactiveContext(threadID);
  }
  return *StaticManagers[threadID];
}


int exahype::reactive::ReactiveContext::getNextMPITag() {
  static std::atomic<int> next_tag(1); //0 is reserved for status
retry:
  int val =  next_tag.load();
  int val_o = val; 

  if(val>MaxSupportedTag-1) {
    logWarning("getOffloadingTag","MPI tag rollover for reactive communication!");
    val = 1;
  }

  int update = val + 1;

  if(!std::atomic_compare_exchange_strong(&next_tag, &val_o, update)) {
    goto retry;
  }
  return val;
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
