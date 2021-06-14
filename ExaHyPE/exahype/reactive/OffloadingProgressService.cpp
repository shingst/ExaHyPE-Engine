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

#if defined(Parallel) && defined(SharedTBB)

#include "../reactive/OffloadingProgressService.h"
#include "tarch/services/ServiceFactory.h"
#include "tarch/multicore/Jobs.h"

//#ifndef TaskSharing

#ifdef MPIWaitsProgressOffloading
#ifndef OffloadingUseProgressThread
registerService(exahype::reactive::OffloadingProgressService)
#endif
#endif

tarch::logging::Log exahype::reactive::OffloadingProgressService::_log("exahype::reactive::OffloadingProgressService");

exahype::reactive::OffloadingProgressService::OffloadingProgressService()
: _solver(nullptr),  _isSet(false), _isEnabled(false){}

exahype::reactive::OffloadingProgressService::~OffloadingProgressService() {}

void exahype::reactive::OffloadingProgressService::enable() {
  _isEnabled = true;
}

void exahype::reactive::OffloadingProgressService::receiveDanglingMessages() {
  if(_isSet && _isEnabled) {
    exahype::solvers::ADERDGSolver::progressOffloading(_solver, true, 1);
  }
}

exahype::reactive::OffloadingProgressService& exahype::reactive::OffloadingProgressService::getInstance() {
  static OffloadingProgressService service;
  return service;
}

void exahype::reactive::OffloadingProgressService::setSolver(exahype::solvers::ADERDGSolver *solver) {
  _solver = solver;
  _isSet = true;
}
#endif

