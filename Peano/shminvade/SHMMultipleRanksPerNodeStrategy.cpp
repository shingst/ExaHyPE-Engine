#include "SHMMultipleRanksPerNodeStrategy.h"
#include "SHMController.h"
#include "SHMSharedMemoryBetweenTasks.h"
#include "SHMMacros.h"


#include <iostream>


shminvade::SHMMultipleRanksPerNodeStrategy::SHMMultipleRanksPerNodeStrategy() {
  #if SHM_INVADE_DEBUG>=4
  std::cout << getSHMDebugPrefix() <<  "created SHMMultipleRanksPerNodeStrategy" << std::endl;
  #endif
}


shminvade::SHMMultipleRanksPerNodeStrategy::~SHMMultipleRanksPerNodeStrategy() {
}


std::set<int> shminvade::SHMMultipleRanksPerNodeStrategy::invade(int wantedNumberOfCores) {
  std::set<int> bookedCores;

  const int masterCore = SHMController::getInstance().getMasterCore();
  // lazy registration of master core in shared mem
  if ( !SHMSharedMemoryBetweenTasks::getInstance().isBooked(masterCore) ) {
    #if SHM_INVADE_DEBUG>=1
    std::cout << getSHMDebugPrefix() <<  "Core " << masterCore << " is own master but is not flagged as booked (line:" << __LINE__  << ",file: " << __FILE__ << ")" << std::endl;
    #endif
    const bool success = SHMSharedMemoryBetweenTasks::getInstance().tryToBookCoreForProcess(masterCore);
    if (!success) {
      #if SHM_INVADE_DEBUG>=2
      std::cout << getSHMDebugPrefix() <<  "WARNING: Core " << masterCore << " is own master but we failed to book it (line:" << __LINE__  << ",file: " << __FILE__ << ")" << std::endl;
      #endif
    }
  }

  const int numberOfCores = SHMController::getInstance()._cores.size();

  // Alternating search around master
  if (wantedNumberOfCores>0) {
    for (int i=1; i<numberOfCores/2+1; ) {
  	const int core = (masterCore + i + numberOfCores) % numberOfCores;
      if ( SHMController::getInstance().tryToBookCore(core) ){
        const bool success = SHMSharedMemoryBetweenTasks::getInstance().tryToBookCoreForProcess(core);
        if (success) {
          bookedCores.insert(core);
          wantedNumberOfCores--;
          if (wantedNumberOfCores==0) break;
        }
        else {
          i = -numberOfCores;
        }
      }

      i = i>0 ? -i : -i+1;
    }
  }

/*
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
  }

*/
  #if SHM_INVADE_DEBUG>=4
  if (!bookedCores.empty()) {
    std::cout << getSHMDebugPrefix() <<  "Invaded " << bookedCores.size() << " thread(s) in total with " << wantedNumberOfCores << " open requests (line:" << __LINE__  << ",file: " << __FILE__ << ")" << std::endl;
    std::cout << getSHMDebugPrefix() <<  "Known core-process association: " << SHMSharedMemoryBetweenTasks::getInstance().getCoreProcessAssociation() << std::endl;
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
  std::cout << getSHMDebugPrefix() <<  "retreated from " << cores.size() << " core(s) (line:" << __LINE__  << ",file: " << __FILE__ << ")" << std::endl;
  std::cout << getSHMDebugPrefix() <<  "known core-process association: " << SHMSharedMemoryBetweenTasks::getInstance().getCoreProcessAssociation() << std::endl;
  #endif
}
