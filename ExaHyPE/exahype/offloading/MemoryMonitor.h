
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


#ifndef EXAHYPE_OFFLOADING_MEMORYMONITOR_H_
#define EXAHYPE_OFFLOADING_MEMORYMONITOR_H_

#include "tarch/logging/Log.h"
#include "tarch/services/Service.h"

#include "exahype/solvers/ADERDGSolver.h"

namespace exahype {
  namespace offloading {
    class MemoryMonitor;
  }
}


class exahype::offloading::MemoryMonitor : public tarch::services::Service {

  private:

  /**
   * The log device of this class.
   */
  static tarch::logging::Log _log;

  public:
  MemoryMonitor();
  static MemoryMonitor& getInstance();

  static size_t getFreeMemMB();

  virtual ~MemoryMonitor();
  /**
   *  Progress method invoked by Peano
   */
  virtual void receiveDanglingMessages();

};


#endif /* EXAHYPE_OFFLOADING_MEMORYMONITOR_H_ */
