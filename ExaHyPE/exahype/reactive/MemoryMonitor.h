
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
  namespace reactive {
    class MemoryMonitor;
    class MemorySample;
  }
}

/**
 * Represents a single sample point of a memory consumption snapshot.
 */
class exahype::reactive::MemorySample {
  private:
    /**
     * Elapsed wall clock time since the start of the application.
     */
    unsigned long _elapsedSeconds;
    size_t _freeMemMB;
    size_t _usedMemMB;
  public:
    MemorySample(unsigned long, std::size_t freeMem, std::size_t usedMem );

    std::string to_string();
};

/**
 * The memory monitor is a service invoked by Peano which repeatedly
 * samples ExaHyPE's memory usage at runtime. The samples
 * can be dumped into an output text file for postprocessing.
 * Memory monitoring is deactivated per default.
 */
class exahype::reactive::MemoryMonitor : public tarch::services::Service {

  private:

  /**
   * The log device of this class.
   */
  static tarch::logging::Log _log;

  /**
   * Vector of memory state samples.
   */
  std::vector<MemorySample> _samples;

  /**
   * Timestamp at which the last sample was taken.
   */
  std::chrono::system_clock::time_point _timestampLastSample;

  /**
   * Timestamp at which the monitoring started.
   */
  std::chrono::system_clock::time_point _timestampStart;

  /**
   * Directory path where the output should be stored.
   */
  std::string _output_dir;

  public:

  MemoryMonitor();
  static MemoryMonitor& getInstance();

  /**
   * Reads the amount of free memory in MB from the linux filesystem.
   */
  static size_t getFreeMemMB();

  virtual ~MemoryMonitor();
  /**
   *  Progress method invoked by Peano
   */
  virtual void receiveDanglingMessages();

  /**
   * Dumps memory usage to a file.
   */
  virtual void dumpMemoryUsage();

  virtual void setOutputDir(std::string output_dir);
};


#endif /* EXAHYPE_OFFLOADING_MEMORYMONITOR_H_ */
