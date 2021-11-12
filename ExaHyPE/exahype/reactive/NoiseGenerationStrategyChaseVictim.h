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

#ifndef EXAHYPE_EXAHYPE_REACTIVE_NOISEGENERATIONSTRATEGYCHASEVICTIM_H_
#define EXAHYPE_EXAHYPE_REACTIVE_NOISEGENERATIONSTRATEGYCHASEVICTIM_H_

#include "exahype/reactive/NoiseGenerationStrategy.h"
#include "tarch/logging/Log.h"

namespace exahype {
namespace reactive {

/**
 * Implements a noise generation strategy where only offloading victim ranks are disturbed.
 */
class NoiseGenerationStrategyChaseVictim : public NoiseGenerationStrategy {

public:
  static tarch::logging::Log     _log;

  NoiseGenerationStrategyChaseVictim();
  NoiseGenerationStrategyChaseVictim(double factor, double baseNoise);

  virtual ~NoiseGenerationStrategyChaseVictim();

  virtual void generateNoiseIfActive(const int myRank);
private:
  double _scalingFactor;
  double _baseNoise;
};

} /* namespace offloading */
} /* namespace exahype */

#endif /* EXAHYPE_EXAHYPE_REACTIVE_NOISEGENERATIONSTRATEGYCHASEVICTIM_H_ */
