#include "SHMController.h"
#include "SHMLockTask.h"


#include <iostream>
#include <thread>
#include <sstream>
#include <cassert>


tbb::task_group_context  shminvade::SHMController::InvasiveTaskGroupContext;


shminvade::SHMController::SHMController():
  _pinningObserver(),
  _globalThreadCountControl(tbb::global_control::max_allowed_parallelism,std::thread::hardware_concurrency()),
  _switchedOn( true ) {
  _pinningObserver.observe(true);

  InvasiveTaskGroupContext.set_priority( tbb::priority_low );

  init( true, 1, 1 );
}


shminvade::SHMController::~SHMController() {
  #if SHM_INVADE_DEBUG>=1
  std::cout << getSHMDebugPrefix() <<  "Destroy SHMController (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
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


int shminvade::SHMController::getMaxAvailableCores(bool useHyperthreading) const {
  return useHyperthreading ? std::thread::hardware_concurrency() : std::thread::hardware_concurrency()/2;
}


int shminvade::SHMController::getFreeCores(bool useHyperthreading) const {
  return getMaxAvailableCores(useHyperthreading) - getBookedCores();
}


int shminvade::SHMController::getBookedCores() const {
  int result = 1;

  for (auto p: _cores) {
    if ( getThreadType(p.first)==ThreadType::Owned) {
      result++;
    }
  }

  return result;
}


shminvade::SHMController::ThreadType shminvade::SHMController::getThreadType( int core ) const {
  // switch to non-const to ensure that noone else is modifying the entry at the very moment
  ThreadTable::const_accessor a;
  _cores.find(a,core);
  return a->second->type;
}


void shminvade::SHMController::shutdown() {
  _switchedOn = false;

  #if SHM_INVADE_DEBUG>=2
  std::cout << getSHMDebugPrefix() <<  "start to instruct all threads to shut down (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
  #endif

  for (auto p: _cores) {
    ThreadTable::accessor a;
    _cores.find(a,p.first);
    ThreadState::Mutex::scoped_lock lock( a->second->mutex );
    a->second->type = ThreadType::Shutdown;
  }

  __TBB_Yield();
  #ifdef SHM_MIN_SLEEP
  sleep(SHM_MIN_SLEEP);
  #endif

  int totalNumberOfLockTasks = 0;
  for (auto p: _cores) {
    ThreadTable::accessor a;
    _cores.find(a,p.first);
    ThreadState::Mutex::scoped_lock lock( a->second->mutex );
    totalNumberOfLockTasks += a->second->numberOfExistingLockTasks;
  }
  if (totalNumberOfLockTasks>0) {
    #if SHM_INVADE_DEBUG>=1
    std::cout << getSHMDebugPrefix() <<  "wait for all lock threads to terminate as lock tasks seem to be alive (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
    #endif
    sleep(SHM_MAX_SLEEP);
  }

  #if SHM_INVADE_DEBUG>=2
  std::cout << getSHMDebugPrefix() <<  "Assume all threads have terminated (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
  #endif
}


bool shminvade::SHMController::tryToBookCore( int core ) {
  if (!_switchedOn) return false;

  bool result = false;
  ThreadTable::accessor a;
  _cores.find(a,core);

  ThreadState::Mutex::scoped_lock  lock( a->second->mutex);
  if ( a->second->type==SHMController::ThreadType::NotOwned ) {
    a->second->type = SHMController::ThreadType::Owned;
    result = true;
    #if SHM_INVADE_DEBUG>=4
    std::cout << getSHMDebugPrefix() <<  "Invade core " << core << " (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
    #endif
  }

  lock.release();
  a.release();

  return result;
}


void shminvade::SHMController::retreat( int core ) {
  ThreadTable::accessor a;
  _cores.find(a,core);

  ThreadState::Mutex::scoped_lock  lock( a->second->mutex);

  if ( a->second->type==ThreadType::Owned ) {
    a->second->type = ThreadType::NotOwned;
    #if SHM_INVADE_DEBUG>=4
    std::cout << getSHMDebugPrefix() <<  "Retreat from core " << core << " (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
    #endif
  }

  if (
    a->second->type!=ThreadType::Shutdown
    &&
    a->second->numberOfExistingLockTasks==0
  ) {
    a->second->numberOfExistingLockTasks++;
    tbb::task &t = *new(tbb::task::allocate_root(InvasiveTaskGroupContext)) SHMLockTask(core);
    tbb::task::enqueue(t);
    #if SHM_INVADE_DEBUG>=4
    std::cout << getSHMDebugPrefix() <<  "Issue new lock task for core " << core << " (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
    #endif
  }

  lock.release();
  a.release();
}


void shminvade::SHMController::retreatFromAllCores() {
  for (auto p: _cores) {
    if (getThreadType(p.first)==ThreadType::Owned) {
      retreat(p.first);
    }
  }
}


void shminvade::SHMController::registerNewCore(int core, ThreadType initialType) {
  ThreadState* newThread = new ThreadState(initialType);
  assert( initialType!=ThreadType::Owned );

  if (_cores.count(core)==0) {
    _cores.insert( std::pair<int, ThreadState* >(core,newThread) );
    if ( initialType==ThreadType::NotOwned ) {
      retreat(core);
    }
  }
  else {
    ThreadTable::accessor a;
	_cores.find(a,core);
    if (a->second->type==ThreadType::Master and initialType!=ThreadType::Master) {
      a->second->type = initialType;
      retreat(core);
    }
    else {
      a->second->type = initialType;
    }
  }

  #if SHM_INVADE_DEBUG>=1
  std::cout << getSHMDebugPrefix() <<  "Register new core " << core << " as " << newThread->toString() << " (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
  #endif
}


std::string shminvade::SHMController::ThreadState::toString() const {
  std::ostringstream msg;
  msg << "(";
  switch (type) {
    case ThreadType::Master:
      msg << "master";
      break;
    case ThreadType::Owned:
      msg << "owned";
      break;
    case ThreadType::NotOwned:
      msg << "not-owned";
      break;
    case ThreadType::Shutdown:
      msg << "shutdown";
      break;
  }
  msg << ",no-of-existing-lock-tasks=" << numberOfExistingLockTasks << ")";
  return msg.str();
}


void shminvade::SHMController::init( bool useHyperthreading, int ranksPerNode, int rank ) {
  const int localRankNumber = (rank % ranksPerNode);
  const int coresPerRank    = getMaxAvailableCores(useHyperthreading) /  ranksPerNode;
  const int masterCore      = localRankNumber * coresPerRank + coresPerRank/2;

  if (coresPerRank<1) {
    std::cerr << getSHMDebugPrefix() <<  "Init called with " << useHyperthreading << " hyperthreading, "
    		  << ranksPerNode << " ranks per node on rank " << rank << " whichc yields " << coresPerRank << " cores per rank" << std::endl;
  }

  #if SHM_INVADE_DEBUG>=1
  std::cout << getSHMDebugPrefix() <<  "Init called with " << useHyperthreading << " hyperthreading, "
  		    << ranksPerNode << " ranks per node on rank " << rank << " which yields " << coresPerRank << " cores per rank"
			<< " (line:" << __LINE__ << ",file:" << __FILE__ << ")"
			<< std::endl;
  std::cout << getSHMDebugPrefix() <<  "Create " << getMaxAvailableCores(true) << " threads " << std::endl;
  #endif

  _globalThreadCountControl = tbb::global_control(tbb::global_control::max_allowed_parallelism,std::thread::hardware_concurrency());

  _cores.clear();

  for (int i=0; i<getMaxAvailableCores(true); i++ ) {
	if (i==masterCore) {
      registerNewCore(i,ThreadType::Master);
	}
	else if (i<getMaxAvailableCores(useHyperthreading)) {
      registerNewCore(i,ThreadType::NotOwned );
	}
	else {
      registerNewCore(i,ThreadType::Shutdown );
	}
  }
}
