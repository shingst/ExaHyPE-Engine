#if defined(SharedTBB)  && defined(Parallel) && defined(DistributedStealing)

#include "exahype/stealing/StealingProgressService.h"
#include "tarch/services/ServiceFactory.h"
#include "tarch/multicore/Jobs.h"
registerService(exahype::stealing::StealingProgressService);


tarch::logging::Log exahype::stealing::StealingProgressService::_log("exahype::stealing::StealingProgressService");

exahype::stealing::StealingProgressService::StealingProgressService() : _isSet(false), _solver(nullptr) {};

exahype::stealing::StealingProgressService::~StealingProgressService() {};

void exahype::stealing::StealingProgressService::receiveDanglingMessages() {
  if(_isSet) {
    exahype::solvers::ADERDGSolver::progressStealing(_solver);
  }
}

exahype::stealing::StealingProgressService& exahype::stealing::StealingProgressService::getInstance() {
  static StealingProgressService service;
  return service;
}

void exahype::stealing::StealingProgressService::setSolver(exahype::solvers::ADERDGSolver *solver) {
  _solver = solver;
  _isSet = true;
}

#endif

