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

#include "NoiseGenerationStrategyRoundRobin.h"

#include <unistd.h>

namespace exahype {
namespace offloading {

NoiseGenerationStrategyRoundRobin::NoiseGenerationStrategyRoundRobin() {
	// TODO Auto-generated constructor stub

}

NoiseGenerationStrategyRoundRobin::~NoiseGenerationStrategyRoundRobin() {
	// TODO Auto-generated destructor stub
}

void NoiseGenerationStrategyRoundRobin::generateNoise(int rank, std::chrono::system_clock::time_point timestamp) {
  usleep(10000);
}

} /* namespace offloading */
} /* namespace exahype */
