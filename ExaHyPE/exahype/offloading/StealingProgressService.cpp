/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#if defined(SharedTBB)  && defined(Parallel) && defined(DistributedStealing)

#include "exahype/offloading/StealingProgressService.h"
#include "tarch/services/ServiceFactory.h"
#include "tarch/multicore/Jobs.h"
registerService(exahype::offloading::StealingProgressService);

tarch::logging::Log exahype::offloading::StealingProgressService::_log("exahype::stealing::StealingProgressService");

exahype::offloading::StealingProgressService::StealingProgressService()
: _isSet(false), _solver(nullptr) {};

exahype::offloading::StealingProgressService::~StealingProgressService() {};

void exahype::offloading::StealingProgressService::receiveDanglingMessages() {
  if(_isSet) {
    exahype::solvers::ADERDGSolver::setMaxNumberOfIprobesInProgressStealing(1);
    exahype::solvers::ADERDGSolver::progressStealing(_solver);
    exahype::solvers::ADERDGSolver::setMaxNumberOfIprobesInProgressStealing(std::numeric_limits<int>::max());
  }
}

exahype::offloading::StealingProgressService& exahype::offloading::StealingProgressService::getInstance() {
  static StealingProgressService service;
  return service;
}

void exahype::offloading::StealingProgressService::setSolver(exahype::solvers::ADERDGSolver *solver) {
  _solver = solver;
  _isSet = true;
}

#endif

