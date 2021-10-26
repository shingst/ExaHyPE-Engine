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

#ifndef EXAHYPE_OFFLOADING_STPSTATSTRACER_H_
#define EXAHYPE_OFFLOADING_STPSTATSTRACER_H_

#include <string>
#include <vector>
#include "tarch/logging/Log.h"
#include "tarch/multicore/Lock.h"
#include "tarch/multicore/BooleanSemaphore.h"

namespace exahype {
namespace reactive {

enum STPTraceKey {ADERDGPrediction = 0, ADERDGOwnMigratable  = 1, ADERDGRemoteMigratable = 2, LimitingFusedTimeStep = 3};


class STPStatsTracer {

private:
	std::string _outputDir;
  static tarch::logging::Log  _log;

  std::vector<unsigned long long>  _iterations[4];
  std::vector<unsigned long long>  _elapsed[4];

  tarch::multicore::BooleanSemaphore* _semaphores[4];
  std::vector<tarch::multicore::Lock> _locks[4];

  int _dumpInterval;
  int _dumpCnt;

  bool isActive(int timestep);

public:
	STPStatsTracer();
	virtual ~STPStatsTracer();

	void dumpAndResetTraceIfActive();

  void writeTracingEventIteration(unsigned int iterations, STPTraceKey type);
  void writeTracingEventRun(unsigned int elapsed, STPTraceKey type);
  void writeTracingEventRunIterations(unsigned int elapsed, unsigned int iterations, STPTraceKey type);

  void writeTracingEventIterationDetailed(unsigned int iterations, STPTraceKey type);
  void writeTracingEventRunDetailed(unsigned int elapsed, STPTraceKey type);
  void writeTracingEventRunIterationsDetailed(unsigned int elapsed, unsigned int iterations, STPTraceKey type);

  void setDumpInterval(int interval);
  void setOutputDir(std::string directory);
  static STPStatsTracer& getInstance();
};

} /* namespace offloading */
} /* namespace exahype */

#endif /* EXAHYPE_OFFLOADING_STPSTATSTRACER_H_ */
