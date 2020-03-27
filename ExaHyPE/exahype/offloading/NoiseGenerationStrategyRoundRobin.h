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

#ifndef EXAHYPE_EXAHYPE_OFFLOADING_NOISEGENERATIONSTRATEGYROUNDROBIN_H_
#define EXAHYPE_EXAHYPE_OFFLOADING_NOISEGENERATIONSTRATEGYROUNDROBIN_H_

#include "NoiseGenerationStrategy.h"
#include "tarch/logging/Log.h"

namespace exahype {
namespace offloading {

class NoiseGenerationStrategyRoundRobin : public NoiseGenerationStrategy {
public:
  static tarch::logging::Log     _log;
  NoiseGenerationStrategyRoundRobin();
  virtual ~NoiseGenerationStrategyRoundRobin();

  virtual void generateNoise(int rank, std::chrono::system_clock::time_point timestamp);
};

} /* namespace offloading */
} /* namespace exahype */

#endif /* EXAHYPE_EXAHYPE_OFFLOADING_NOISEGENERATIONSTRATEGYROUNDROBIN_H_ */
