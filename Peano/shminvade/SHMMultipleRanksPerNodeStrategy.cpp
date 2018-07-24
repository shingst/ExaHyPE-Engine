#include "SHMMultipleRanksPerNodeStrategy.h"
#include "SHMController.h"
#include "SHMSharedMemoryBetweenTasks.h"
#include "SHMMacros.h"


#include <iostream>


shminvade::SHMMultipleRanksPerNodeStrategy::SHMMultipleRanksPerNodeStrategy() {
  #if SHM_INVADE_DEBUG>=4
  std::cout << SHM_DEBUG_PREFIX <<  "created SHMMultipleRanksPerNodeStrategy" << std::endl;
  #endif
}


shminvade::SHMMultipleRanksPerNodeStrategy::~SHMMultipleRanksPerNodeStrategy() {
}


std::set<int> shminvade::SHMMultipleRanksPerNodeStrategy::invade(int wantedNumberOfCores) {
  std::set<int> bookedCores;

  for (auto p: SHMController::getInstance()._cores) {
    if ( SHMController::getInstance().tryToBookCore(p.first) ){
      const bool success = SHMSharedMemoryBetweenTasks::getInstance().tryToBookCoreForProcess(p.first);
      if (success) {
        bookedCores.insert(p.first);
        wantedNumberOfCores--;
        if (wantedNumberOfCores==0) break;
      }
      else {
        SHMController::getInstance().retreat(p.first);
      }
    }
    else if (
      p.second->type==SHMController::ThreadType::Master
	  and
	  !SHMSharedMemoryBetweenTasks::getInstance().isBooked(p.first)
	) {
      #if SHM_INVADE_DEBUG>=1
      std::cout << SHM_DEBUG_PREFIX <<  "Core " << p.first << " is own master but is not flagged as booked (line:" << __LINE__  << ",file: " << __FILE__ << ")" << std::endl;
      #endif
      const bool success = SHMSharedMemoryBetweenTasks::getInstance().tryToBookCoreForProcess(p.first);
      if (!success) {
        #if SHM_INVADE_DEBUG>=2
        std::cout << SHM_DEBUG_PREFIX <<  "WARNING: Core " << p.first << " is own master but we failed to book it (line:" << __LINE__  << ",file: " << __FILE__ << ")" << std::endl;
        #endif
      }
    }
  }

  #if SHM_INVADE_DEBUG>=4
  if (!bookedCores.empty()) {
    std::cout << SHM_DEBUG_PREFIX <<  "Invaded " << bookedCores.size() << " thread(s) in total with " << wantedNumberOfCores << " open requests (line:" << __LINE__  << ",file: " << __FILE__ << ")" << std::endl;
    std::cout << SHM_DEBUG_PREFIX <<  "Known core-process association: " << SHMSharedMemoryBetweenTasks::getInstance().getCoreProcessAssociation() << std::endl;
  }
  #endif

  return bookedCores;
}


void shminvade::SHMMultipleRanksPerNodeStrategy::cleanUp() {
  SHMSharedMemoryBetweenTasks::getInstance().cleanUp();
}


void shminvade::SHMMultipleRanksPerNodeStrategy::retreat(const std::set<int>& cores) {
  for (auto p: cores) {
    SHMSharedMemoryBetweenTasks::getInstance().freeCore(p);
    SHMController::getInstance().retreat(p);
  }
  #if SHM_INVADE_DEBUG>=4
  std::cout << SHM_DEBUG_PREFIX <<  "retreated from " << cores.size() << " core(s) (line:" << __LINE__  << ",file: " << __FILE__ << ")" << std::endl;
  std::cout << SHM_DEBUG_PREFIX <<  "known core-process association: " << SHMSharedMemoryBetweenTasks::getInstance().getCoreProcessAssociation() << std::endl;
  #endif
}
