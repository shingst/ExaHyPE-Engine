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

#ifndef EXAHYPE_EXAHYPE_OFFLOADING_REPLICATIONSTATISTICS_H_
#define EXAHYPE_EXAHYPE_OFFLOADING_REPLICATIONSTATISTICS_H_

#include <atomic>
#include "tarch/logging/Log.h"

namespace exahype {
namespace reactive {

/**
 * Gathers statistics about job execution and job migration at runtime (only with local recompute or task sharing).
 */
class JobTableStatistics {

private:

  static tarch::logging::Log _log;

  std::atomic<int> _spawnedTasks;
  std::atomic<int> _executedTasks;
  std::atomic<int> _savedTasks;
  std::atomic<int> _receivedTasks;
  std::atomic<int> _sentTasks;
  std::atomic<int> _sentKeys;
  std::atomic<int> _receivedKeys;
  std::atomic<int> _declinedTasks;
  std::atomic<int> _lateTasks;
  std::atomic<int> _recomputedTasks;
  std::atomic<int> _doubleCheckedTasks;
  std::atomic<int> _softErrorsDetected;
  std::atomic<int> _healedTasks;

  JobTableStatistics();
  virtual ~JobTableStatistics();
public:

  static JobTableStatistics& getInstance();

  void printStatistics();

  void notifyLateTask();
  void notifyDeclinedTask();
  void notifyReceivedTask();
  void notifySentTask();
  void notifySentKey();
  void notifyReceivedKey();
  void notifySavedTask();
  void notifySpawnedTask();
  void notifyExecutedTask();
  void notifyRecomputedTask();
  void notifyDoubleCheckedTask();
  void notifyDetectedError();
  void notifyHealedTask();
};

} /* namespace offloading */
} /* namespace exahype */

#endif /* EXAHYPE_EXAHYPE_OFFLOADING_REPLICATIONSTATISTICS_H_ */
