/**;
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

#include "exahype/reactive/NoiseGenerator.h"

#include <unistd.h>

#include "tarch/parallel/Node.h"
#include <chrono>
#include "exahype/reactive/NoiseGenerationStrategy.h"
#include "exahype/reactive/NoiseGenerationStrategyRoundRobin.h"

#if defined(Parallel)

exahype::reactive::NoiseGenerator::NoiseGenerator() {
  _strategy = new NoiseGenerationStrategyRoundRobin();
}

exahype::reactive::NoiseGenerator::~NoiseGenerator() {
  delete _strategy;
}

exahype::reactive::NoiseGenerator& exahype::reactive::NoiseGenerator::getInstance() {
  static NoiseGenerator instance;
  return instance;
}

void exahype::reactive::NoiseGenerator::generateNoiseIfActive() {
  //usleep(1000);
  if(_strategy!=nullptr)
	  _strategy->generateNoiseIfActive(tarch::parallel::Node::getInstance().getRank());
}

void exahype::reactive::NoiseGenerator::setStrategy(NoiseGenerationStrategy *strategy) {
  delete _strategy;
  _strategy = strategy;
}

#endif
