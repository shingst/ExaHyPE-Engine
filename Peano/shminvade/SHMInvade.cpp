#include "SHMInvade.h"
#include "SHMStrategy.h"
#include "SHMController.h"


#include <assert.h>


shminvade::SHMInvade::SHMInvade(int cores):
  _occupiedCores() {
  assert(cores>0 || cores==MaxCores);

  if (cores==MaxCores) {
    cores = SHMController::getInstance().getMaxAvailableCores(true);
  }

  _occupiedCores = SHMStrategy::getInstance().invade(cores);
}


shminvade::SHMInvade::~SHMInvade() {
  retreat();
}


void shminvade::SHMInvade::retreat() {
  if (!_occupiedCores.empty()) {
    SHMStrategy::getInstance().retreat(_occupiedCores);
  }
  _occupiedCores.clear();
}
