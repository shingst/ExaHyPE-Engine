#include "SHMController.h"
#include "SHMLockTask.h"


#include <iostream>
#include <thread>
#include <sstream>


tbb::task_group_context  shminvade::SHMController::InvasiveTaskGroupContext;


shminvade::SHMController::SHMController():
  _switchedOn( true ) {
  #if SHM_INVADE_DEBUG>=1
  std::cout << SHM_DEBUG_PREFIX <<  "Create SHMController, maxCores=" << getMaxAvailableCores() <<
    " (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
  #endif

  InvasiveTaskGroupContext.set_priority( tbb::priority_low );

  registerNewThread(
    (pid_t) syscall (__NR_gettid),
    ThreadType::Master
  );
}


shminvade::SHMController::~SHMController() {
  #if SHM_INVADE_DEBUG>=1
  std::cout << SHM_DEBUG_PREFIX <<  "Destroy SHMController (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
  #endif

  shutdown();
}


shminvade::SHMController&  shminvade::SHMController::getInstance() {
  static SHMController singleton;
  return singleton;
}


void shminvade::SHMController::switchOn() {
  _switchedOn = true;
}


void shminvade::SHMController::switchOff() {
  _switchedOn = false;
}


int shminvade::SHMController::getMaxAvailableCores() const {
  return std::thread::hardware_concurrency();
}


int shminvade::SHMController::getFreeCores() const {
  return getMaxAvailableCores() - getBookedCores();
}


int shminvade::SHMController::getBookedCores() const {
  int result = 1;

  for (auto p: _threads) {
    if ( getThreadType(p.first)==ThreadType::ExclusivelyOwned) {
      result++;
    }
  }

  return result;
}


shminvade::SHMController::ThreadType shminvade::SHMController::getThreadType( pid_t pid ) const {
  // switch to non-const to ensure that noone else is modifying the entry at the very moment
  ThreadTable::const_accessor a;
  _threads.find(a,pid);
  return a->second->type;
}


void shminvade::SHMController::shutdown() {
  _switchedOn = false;

  #if SHM_INVADE_DEBUG>=2
  std::cout << SHM_DEBUG_PREFIX <<  "start to instruct all threads to shut down (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
  #endif

  for (auto p: _threads) {
    ThreadTable::accessor a;
    _threads.find(a,p.first);
    ThreadState::Mutex::scoped_lock lock( a->second->mutex );
    a->second->type = ThreadType::Shutdown;
  }

  __TBB_Yield();
  #ifdef SHM_MIN_SLEEP
  sleep(SHM_MIN_SLEEP);
  #endif

  int totalNumberOfLockTasks = 0;
  for (auto p: _threads) {
    ThreadTable::accessor a;
    _threads.find(a,p.first);
    ThreadState::Mutex::scoped_lock lock( a->second->mutex );
    totalNumberOfLockTasks += a->second->numberOfExistingLockTasks;
  }
  if (totalNumberOfLockTasks>0) {
    #if SHM_INVADE_DEBUG>=1
    std::cout << SHM_DEBUG_PREFIX <<  "wait for all lock threads to terminate as lock tasks seem to be alive (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
    #endif
    sleep(SHM_MAX_SLEEP);
  }

  #if SHM_INVADE_DEBUG>=2
  std::cout << SHM_DEBUG_PREFIX <<  "assume all threads have terminated (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
  #endif
}


bool shminvade::SHMController::tryToBookThread( pid_t pid ) {
  if (!_switchedOn) return false;

  bool result = false;
  ThreadTable::accessor a;
  _threads.find(a,pid);

  ThreadState::Mutex::scoped_lock  lock( a->second->mutex);
  if ( a->second->type==SHMController::ThreadType::NotOwned ) {
    a->second->type = SHMController::ThreadType::ExclusivelyOwned;
    result = true;
    #if SHM_INVADE_DEBUG>=4
    std::cout << SHM_DEBUG_PREFIX <<  "invade thread " << pid << " (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
    #endif
  }

  lock.release();
  a.release();

  return result;
}


void shminvade::SHMController::retreatFromThread( pid_t pid ) {
  ThreadTable::accessor a;
  _threads.find(a,pid);

  ThreadState::Mutex::scoped_lock  lock( a->second->mutex);

  if ( a->second->type==ThreadType::ExclusivelyOwned ) {
    a->second->type = ThreadType::NotOwned;
    #if SHM_INVADE_DEBUG>=4
    std::cout << SHM_DEBUG_PREFIX <<  "retreat from thread " << pid << " (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
    #endif
  }

  if (
    a->second->type!=ThreadType::Shutdown
    &&
    a->second->numberOfExistingLockTasks==0
  ) {
    a->second->numberOfExistingLockTasks++;
    tbb::task &t = *new(tbb::task::allocate_root(InvasiveTaskGroupContext)) SHMLockTask(pid);
    tbb::task::enqueue(t);
    #if SHM_INVADE_DEBUG>=4
    std::cout << SHM_DEBUG_PREFIX <<  "issue new lock task for thread " << pid << " (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
    #endif
  }

  lock.release();
  a.release();
}


void shminvade::SHMController::retreatFromAllCores() {
  for (auto p: _threads) {
    if (getThreadType(p.first)==ThreadType::ExclusivelyOwned) {
      retreatFromThread(p.first);
    }
  }
}


void shminvade::SHMController::registerNewThread(int threadId, ThreadType initialType) {
  ThreadState* newThread = new ThreadState(initialType);
  _threads.insert( std::pair<pid_t, ThreadState* >(threadId,newThread) );

  #if SHM_INVADE_DEBUG>=1
  std::cout << SHM_DEBUG_PREFIX <<  "register new thread " << threadId << " (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
  #endif

  retreatFromThread(threadId);
}


std::string shminvade::SHMController::toString( ThreadState state ) {
  std::ostringstream msg;
  msg << "(";
  switch (state.type) {
    case ThreadType::Master:
      msg << "master";
      break;
    case ThreadType::ExclusivelyOwned:
      msg << "exclusively-owned";
      break;
    case ThreadType::NotOwned:
      msg << "not-owned";
      break;
    case ThreadType::Shutdown:
      msg << "shutdown";
      break;
  }
  msg << ",no-of-existing-lock-tasks" << state.numberOfExistingLockTasks << ")";
  return msg.str();
}
