#include "SHMLockTask.h"
#include "SHMController.h"


#include <iostream>
#include <assert.h>


shminvade::SHMLockTask::SHMLockTask(int core, int sleepTime):
  _core( core ),
  _sleepTime( sleepTime ) {
  // TBB complains and clarifies that affinities are ignored for enqueued
  // tasks
  // set_affinity( _pid_t );
}


void shminvade::SHMLockTask::reenqueue() {
  if (SHMController::getInstance().getThreadType(_core) == SHMController::ThreadType::Shutdown) {
    #if SHM_INVADE_DEBUG>=8
    std::cout << SHM_DEBUG_PREFIX <<  "Controller seems to be down already, so stop lock task too (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
    #endif
  }
  else {
    tbb::task &t = *new(tbb::task::allocate_root(SHMController::InvasiveTaskGroupContext)) SHMLockTask(
      _core,std::min( SHM_MAX_SLEEP,_sleepTime+1)
    );
    tbb::task::enqueue(t);
  }
}


void shminvade::SHMLockTask::terminate() {
  #if SHM_INVADE_DEBUG>=8
  std::cout << SHM_DEBUG_PREFIX <<  "Task for core " << _core << " terminates now (line:" << __LINE__ << ",file: " << __FILE__ << ")" << std::endl;
  #endif

  SHMController::ThreadTable::accessor a;
  SHMController::getInstance()._cores.find(a,_core);
  SHMController::ThreadState::Mutex::scoped_lock lock( a->second->mutex );
//  assert( a->second->numberOfExistingLockTasks>0);
  if (a->second->numberOfExistingLockTasks>0) {
    a->second->numberOfExistingLockTasks--;
  }
  lock.release();
  a.release();
}



tbb::task* shminvade::SHMLockTask::execute() {
  const int   currentCore     = sched_getcpu();

  SHMController::ThreadType state = SHMController::getInstance().getThreadType(_core);
  switch (state) {
    case SHMController::ThreadType::Master:
      // there should be no lock tasks on the master so let this one die
      #if SHM_INVADE_DEBUG>=1
      std::cout << SHM_DEBUG_PREFIX <<  "Core " << _core << " is the master. Retreat immediately as lock task found on core " << currentCore << " (line:" << __LINE__ << ",file: " << __FILE__ << ")" << std::endl;
      #endif
      terminate();
      return nullptr;
    case SHMController::ThreadType::Owned:
      // we own it so let this one die
      #if SHM_INVADE_DEBUG>=8
      std::cout << SHM_DEBUG_PREFIX <<  "Core " << _core << " is owned by process. Lock task found on core " << currentCore << ". Terminate (line:" << __LINE__ << ",file: " << __FILE__ << ")" << std::endl;
      #endif
      terminate();
      return nullptr;
      case SHMController::ThreadType::NotOwned:
        if ( currentCore!=_core ) {
          #if SHM_INVADE_DEBUG>=8
    	  std::cout << SHM_DEBUG_PREFIX <<  "Lock task for core " << _core <<
            " has been invoked on core " << currentCore <<
    	    " so re-enqueue (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
    	  #endif
    	  reenqueue();
    	  return nullptr;
    	}
        else {
          __TBB_Yield();
          if (_sleepTime<SHM_MAX_SLEEP) {
            _sleepTime *= 2;
          }
          #if SHM_INVADE_DEBUG>=2
          std::cout << SHM_DEBUG_PREFIX <<  "Thread on core " << _core << " should not be used. Make lock task sleep for " << _sleepTime << "s before we reenqueue (line:" << __LINE__ << ",file: " << __FILE__ << ")" << std::endl;
          #endif
          sleep(_sleepTime);
          reenqueue();
          return nullptr;
        }
      case SHMController::ThreadType::Shutdown:
        #if SHM_INVADE_DEBUG>=4
        std::cout << SHM_DEBUG_PREFIX <<  "Core " << _core << " is marked to shut down. Lock task found on core " << currentCore << ". Terminate (line:" << __LINE__ << ",file: " << __FILE__ << ")" << std::endl;
        #endif
        // we should die anyway
        terminate();
        return nullptr;
  }
  return nullptr;
}
