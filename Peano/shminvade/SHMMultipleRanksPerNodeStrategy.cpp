#include "SHMMultipleRanksPerNodeStrategy.h"
#include "SHMController.h"
#include "SHMSharedMemoryBetweenTasks.h"
#include "SHMMacros.h"


#include <iostream>


shminvade::SHMMultipleRanksPerNodeStrategy::SHMMultipleRanksPerNodeStrategy() {
}


shminvade::SHMMultipleRanksPerNodeStrategy::~SHMMultipleRanksPerNodeStrategy() {
}


std::set<pid_t> shminvade::SHMMultipleRanksPerNodeStrategy::invadeThreads(int wantedNumberOfThreads) {
  std::set<pid_t> bookedCores;

  for (auto p: SHMController::getInstance()._threads) {
    if ( SHMController::getInstance().tryToBookThread(p.first) ) {
      const bool success = SHMSharedMemoryBetweenTasks::getInstance().tryToBookThreadForProcess(p.first);
      if (success) {
        bookedCores.insert(p.first);
        wantedNumberOfThreads--;
        if (wantedNumberOfThreads==0) break;
      }
      else {
        SHMController::getInstance().retreatFromThread(p.first);
      }
    }
  }

  #if SHM_INVADE_DEBUG>=4
  if (!bookedCores.empty()) {
    std::cout << SHM_DEBUG_PREFIX <<  "invaded " << bookedCores.size() << " thread(s) in total with " << wantedNumberOfThreads << " open requests (line:" << __LINE__  << ",file: " << __FILE__ << ")" << std::endl;
    std::cout << SHM_DEBUG_PREFIX <<  "known thread-process association: " << SHMSharedMemoryBetweenTasks::getInstance().getThreadProcessAssociation() << std::endl;
  }
  #endif

  return bookedCores;
}


void shminvade::SHMMultipleRanksPerNodeStrategy::cleanUp() {
  SHMSharedMemoryBetweenTasks::getInstance().cleanUp();
}


void shminvade::SHMMultipleRanksPerNodeStrategy::retreat(const std::set<pid_t>& threadIds) {
  for (auto p: threadIds) {
    SHMSharedMemoryBetweenTasks::getInstance().freeThread(p);
    SHMController::getInstance().retreatFromThread(p);
  }
  #if SHM_INVADE_DEBUG>=4
  std::cout << SHM_DEBUG_PREFIX <<  "retreated from " << threadIds.size() << " thread(s) (line:" << __LINE__  << ",file: " << __FILE__ << ")" << std::endl;
  std::cout << SHM_DEBUG_PREFIX <<  "known thread-process association: " << SHMSharedMemoryBetweenTasks::getInstance().getThreadProcessAssociation() << std::endl;
  #endif
}
