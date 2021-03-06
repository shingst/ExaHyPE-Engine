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

#ifndef EXAHYPE_EXAHYPE_REACTIVE_NOISEGENERATIONSTRATEGYROUNDROBIN_H_
#define EXAHYPE_EXAHYPE_REACTIVE_NOISEGENERATIONSTRATEGYROUNDROBIN_H_

#if defined(Parallel)

#include "exahype/reactive/NoiseGenerationStrategy.h"
#include "tarch/logging/Log.h"

namespace exahype {
namespace reactive {

/**
 * Implements a noise generation strategy where noisy ranks
 * are selected in a round-robin fashion.
 */

class NoiseGenerationStrategyRoundRobin : public NoiseGenerationStrategy {

public:
  static tarch::logging::Log     _log;
  NoiseGenerationStrategyRoundRobin();
  NoiseGenerationStrategyRoundRobin(int stepsBetweenDisturbance, double waitFractionTimestep);
  virtual ~NoiseGenerationStrategyRoundRobin();

  virtual void generateNoiseIfActive(const int rank);

private:
  /**
   * The number k of steps after which the next noisy rank is selected
   */
  int _stepsBetweenDisturbance;
  // factor determines the strength of the noise
  double _waitFractionTimestep;
};

} /* namespace offloading */
} /* namespace exahype */

#endif

#endif /* EXAHYPE_EXAHYPE_REACTIVE_NOISEGENERATIONSTRATEGYROUNDROBIN_H_ */
