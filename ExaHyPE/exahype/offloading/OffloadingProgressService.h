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

#if !defined(_EXAHYPE_STEALINGPROGRESSSERVICE_H_) && defined(SharedTBB) && defined(Parallel) && defined(DistributedOffloading)
#define _EXAHYPE_STEALINGPROGRESSSERVICE_H_

#include "tarch/logging/Log.h"
#include "tarch/services/Service.h"

#include "exahype/solvers/ADERDGSolver.h"

namespace exahype {
  namespace offloading {
    class OffloadingProgressService;
  }
}

/**
 * A peano service to make progress on outstanding offloading-related
 * communication.
 *
 * This service make progress on outstanding offloading communication
 * whenever the master thread calls receiveDangling().
 * This ensures that task offloading communication actually progresses in
 * the background. A single solver is attached for now, on which
 * the offloading progress method is invoked. In the future, this may
 * need to be extended to support multiple solvers.
 */
class exahype::offloading::OffloadingProgressService : public tarch::services::Service {

  private:
  /**
   *  Solver on which to make offloading progress.
   */
  exahype::solvers::ADERDGSolver* _solver;
  /**
   * Flag indicating if a solver has been registered.
   */
  bool _isSet;

  bool _isEnabled;
  /**
   * The log device of this class.
   */
  static tarch::logging::Log _log;

  public:
  OffloadingProgressService();
  static OffloadingProgressService& getInstance();
 
  virtual void enable();

  virtual ~OffloadingProgressService();
  /**
   *  Progress method invoked by Peano
   */
  virtual void receiveDanglingMessages();
  /**
   * Registers a solver on which offloading progress is made.
   */
  void setSolver(exahype::solvers::ADERDGSolver *solver);
};

#endif

