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

#ifndef EXAHYPE_EXAHYPE_OFFLOADING_NOISEGENERATOR_H_
#define EXAHYPE_EXAHYPE_OFFLOADING_NOISEGENERATOR_H_

#include "NoiseGenerationStrategy.h"

namespace exahype {
 namespace offloading {
   class NoiseGenerator;
 }
}

class exahype::offloading::NoiseGenerator {
public:
  NoiseGenerator();
  virtual ~NoiseGenerator();

  static NoiseGenerator& getInstance();

  void generateNoise();
  void setStrategy(NoiseGenerationStrategy* strategy);
private:
  NoiseGenerationStrategy *_strategy;
};

#endif /* EXAHYPE_EXAHYPE_OFFLOADING_NOISEGENERATOR_H_ */
