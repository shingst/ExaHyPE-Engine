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

#if defined(SharedTBB)  && defined(Parallel) && defined(DistributedOffloading)

#include "exahype/offloading/OffloadingProgressService.h"
#include "tarch/services/ServiceFactory.h"
#include "tarch/multicore/Jobs.h"
registerService(exahype::offloading::OffloadingProgressService);

tarch::logging::Log exahype::offloading::OffloadingProgressService::_log("exahype::offloading::OffloadingProgressService");

exahype::offloading::OffloadingProgressService::OffloadingProgressService()
: _isSet(false), _isEnabled(false), _solver(nullptr) {};

exahype::offloading::OffloadingProgressService::~OffloadingProgressService() {};

void exahype::offloading::OffloadingProgressService::enable() {
  _isEnabled = true;
}

void exahype::offloading::OffloadingProgressService::receiveDanglingMessages() {
  if(_isSet && _isEnabled) {
    exahype::solvers::ADERDGSolver::setMaxNumberOfIprobesInProgressOffloading(1);
    // ToDo (Philipp): pass number of iterations through progress engine directly
    exahype::solvers::ADERDGSolver::progressOffloading(_solver);
    exahype::solvers::ADERDGSolver::setMaxNumberOfIprobesInProgressOffloading(std::numeric_limits<int>::max());
  }
}

exahype::offloading::OffloadingProgressService& exahype::offloading::OffloadingProgressService::getInstance() {
  static OffloadingProgressService service;
  return service;
}

void exahype::offloading::OffloadingProgressService::setSolver(exahype::solvers::ADERDGSolver *solver) {
  _solver = solver;
  _isSet = true;
}

#endif

