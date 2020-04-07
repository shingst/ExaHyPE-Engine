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

#ifndef EXAHYPE_EXAHYPE_OFFLOADING_NOISEGENERATIONSTRATEGYCHASEVICTIM_H_
#define EXAHYPE_EXAHYPE_OFFLOADING_NOISEGENERATIONSTRATEGYCHASEVICTIM_H_

#include "NoiseGenerationStrategy.h"
#include "tarch/logging/Log.h"

namespace exahype {
namespace offloading {

class NoiseGenerationStrategyChaseVictim : public NoiseGenerationStrategy {
public:
  static tarch::logging::Log     _log;
  NoiseGenerationStrategyChaseVictim();
  NoiseGenerationStrategyChaseVictim(int frequency, double factor);
  virtual ~NoiseGenerationStrategyChaseVictim();

  virtual void generateNoise(int rank, std::chrono::system_clock::time_point timestamp);

  virtual void generateNoiseSTP(int rank, std::chrono::system_clock::time_point timestamp);
private:
  int _frequency;
  double _factor;
};

} /* namespace offloading */
} /* namespace exahype */

#endif /* EXAHYPE_EXAHYPE_OFFLOADING_NOISEGENERATIONSTRATEGYCHASEVICTIM_H_ */
