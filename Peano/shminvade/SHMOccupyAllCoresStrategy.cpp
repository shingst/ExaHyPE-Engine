#include "SHMOccupyAllCoresStrategy.h"
#include "SHMController.h"
#include "SHMMacros.h"


#include <iostream>


shminvade::SHMOccupyAllCoresStrategy::SHMOccupyAllCoresStrategy() {
}


shminvade::SHMOccupyAllCoresStrategy::~SHMOccupyAllCoresStrategy() {
}


std::set<pid_t> shminvade::SHMOccupyAllCoresStrategy::invadeThreads(int wantedNumberOfThreads) {
  std::set<pid_t> bookedCores;

  for (auto p: SHMController::getInstance()._threads) {
    if (SHMController::getInstance().tryToBookThread(p.first)) {
      bookedCores.insert(p.first);
      wantedNumberOfThreads--;
      if (wantedNumberOfThreads==0) break;
    }
  }

  return bookedCores;
}


void shminvade::SHMOccupyAllCoresStrategy::cleanUp() {
}


void shminvade::SHMOccupyAllCoresStrategy::retreat(const std::set<pid_t>& threadIds) {
  for (auto p: threadIds) {
    SHMController::getInstance().retreatFromThread(p);
  }
}
