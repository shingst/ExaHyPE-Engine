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

#ifndef EXAHYPE_EXAHYPE_REACTIVE_NOISEGENERATIONSTRATEGY_H_
#define EXAHYPE_EXAHYPE_REACTIVE_NOISEGENERATIONSTRATEGY_H_

#include <chrono>

namespace exahype {
namespace reactive {

/**
 * Generic strategy pattern for noise generation. A specific noise
 * generation strategy needs to implement this interface. The
 * noise generator context class will call the implemented
 * noise generation methods during the simulation.
 */

class NoiseGenerationStrategy {
public:
  /**
   * Disturbs the calling core according to a noise generation strategy. To be called once per time step.
   * @param rank The current rank
   */
	virtual void generateNoiseIfActive(const int rank) = 0;

	virtual ~NoiseGenerationStrategy(){};
};

} /* namespace offloading */
} /* namespace exahype */

#endif /* EXAHYPE_EXAHYPE_REACTIVE_NOISEGENERATIONSTRATEGY_H_ */
