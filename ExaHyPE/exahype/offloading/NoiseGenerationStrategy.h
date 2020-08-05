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

#ifndef EXAHYPE_EXAHYPE_OFFLOADING_NOISEGENERATIONSTRATEGY_H_
#define EXAHYPE_EXAHYPE_OFFLOADING_NOISEGENERATIONSTRATEGY_H_

#include <chrono>

namespace exahype {
namespace offloading {

/**
 * Generic strategy pattern for noise generation. A specific noise
 * generation strategy needs to implement this interface. The
 * noise generator context class will call the contained
 * noise generation methods during the simulation.
 */

class NoiseGenerationStrategy {
public:
    // generates noise at the beginning of a time step
	virtual void generateNoise(int rank, std::chrono::system_clock::time_point timestamp ) = 0;

	//generates noise when executing STPs
	virtual void generateNoiseSTP(int rank, std::chrono::system_clock::time_point timestamp ) = 0;

	virtual ~NoiseGenerationStrategy(){};
};

} /* namespace offloading */
} /* namespace exahype */

#endif /* EXAHYPE_EXAHYPE_OFFLOADING_NOISEGENERATIONSTRATEGY_H_ */
