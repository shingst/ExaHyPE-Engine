#include "SHMStrategy.h"
#include "SHMOccupyAllCoresStrategy.h"
#include "SHMController.h"

#include <assert.h>


shminvade::SHMStrategy*  shminvade::SHMStrategy::_activeStrategy( new SHMOccupyAllCoresStrategy() );


shminvade::SHMStrategy::~SHMStrategy() {
}


shminvade::SHMStrategy& shminvade::SHMStrategy::getInstance() {
  assert( _activeStrategy!=nullptr );
  return *_activeStrategy;
}


void shminvade::SHMStrategy::setStrategy(SHMStrategy* strategy) {
  assert( strategy!=nullptr );
  delete _activeStrategy;
  SHMController::getInstance().retreatFromAllCores();
  _activeStrategy = strategy;
}

