#include "SHMInvade.h"
#include "SHMStrategy.h"
#include "SHMController.h"


#include <assert.h>


shminvade::SHMInvade::SHMInvade(int threads):
  _occupiedThreads() {
  assert(threads>0 || threads==MaxThreads);

  if (threads==MaxThreads) {
    threads = SHMController::getInstance().getMaxAvailableCores();
  }

  _occupiedThreads = SHMStrategy::getInstance().invadeThreads(threads);
}


shminvade::SHMInvade::~SHMInvade() {
  retreat();
}


void shminvade::SHMInvade::retreat() {
  if (!_occupiedThreads.empty()) {
    SHMStrategy::getInstance().retreat(_occupiedThreads);
  }
  _occupiedThreads.clear();
}
