#include "SHMOccupyAllCoresStrategy.h"
#include "SHMController.h"
#include "SHMMacros.h"


#include <iostream>


shminvade::SHMOccupyAllCoresStrategy::SHMOccupyAllCoresStrategy() {
  #if SHM_INVADE_DEBUG>=4
  std::cout << SHM_DEBUG_PREFIX <<  "created SHMOccupyAllCoresStrategy" << std::endl;
  #endif
}


shminvade::SHMOccupyAllCoresStrategy::~SHMOccupyAllCoresStrategy() {
}


std::set<int> shminvade::SHMOccupyAllCoresStrategy::invade(int wantedNumberOfCores) {
  std::set<int> bookedCores;

  for (auto p: SHMController::getInstance()._cores) {
    if (
	   SHMController::getInstance().tryToBookCore(p.first)
	) {
      bookedCores.insert(p.first);
      wantedNumberOfCores--;
      if (wantedNumberOfCores==0) break;
    }
  }

  return bookedCores;
}


void shminvade::SHMOccupyAllCoresStrategy::cleanUp() {
}


void shminvade::SHMOccupyAllCoresStrategy::retreat(const std::set<int>& cores) {
  for (auto p: cores) {
    SHMController::getInstance().retreat(p);
  }
}
