
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

#include <chrono>
#include <ctime>

namespace exahype {
  namespace offloading {
    class MemoryMonitor;
    class MemoryMeasurement;
  }
}


class exahype::offloading::MemoryMeasurement {
  private:
    unsigned long _elapsed;
    size_t _freeMem;
    size_t _usedMem;
  public:
    MemoryMeasurement(unsigned long, std::size_t freeMem, std::size_t usedMem );

    std::string to_string();
};

class exahype::offloading::MemoryMonitor : public tarch::services::Service {

  private:

  /**
   * The log device of this class.
   */
  static tarch::logging::Log _log;

  std::vector<MemoryMeasurement> _measurements;

  std::chrono::system_clock::time_point _lastMeasurementTimestamp, _start;

  std::string _output_dir;

  public:
  MemoryMonitor();
  static MemoryMonitor& getInstance();

  static size_t getFreeMemMB();

  virtual ~MemoryMonitor();
  /**
   *  Progress method invoked by Peano
   */
  virtual void receiveDanglingMessages();

  virtual void dumpMemoryUsage();

  virtual void setOutputDir(std::string output_dir);
};


#endif /* EXAHYPE_OFFLOADING_MEMORYMONITOR_H_ */
