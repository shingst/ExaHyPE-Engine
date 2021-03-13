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

#ifndef EXAHYPE_EXAHYPE_REACTIVE_NOISEGENERATOR_H_
#define EXAHYPE_EXAHYPE_REACTIVE_NOISEGENERATOR_H_

#include "../reactive/NoiseGenerationStrategy.h"

namespace exahype {
 namespace reactive {
   class NoiseGenerator;
 }
}

/**
 * Implements a simple noise generator that is used to disturb the application at runtime.
 * It provides the context class in the strategy pattern where the specific
 * noise generation strategy is provided by the NoiseGenerationStrategy interface.
 */
class exahype::reactive::NoiseGenerator {
public:
  NoiseGenerator();
  virtual ~NoiseGenerator();

  static NoiseGenerator& getInstance();

  void generateNoise();
  void generateNoiseSTP();

  void setStrategy(NoiseGenerationStrategy* strategy);
private:
  NoiseGenerationStrategy *_strategy;
};

#endif /* EXAHYPE_EXAHYPE_REACTIVE_NOISEGENERATOR_H_ */
