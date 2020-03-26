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

#include "NoiseGenerator.h"
#include "NoiseGenerationStrategy.h"
#include "NoiseGenerationStrategyRoundRobin.h"

#include <unistd.h>

#include "tarch/parallel/Node.h"
#include <chrono>

exahype::offloading::NoiseGenerator::NoiseGenerator() {
  _strategy = new NoiseGenerationStrategyRoundRobin();
}

exahype::offloading::NoiseGenerator::~NoiseGenerator() {
  delete _strategy;
}

exahype::offloading::NoiseGenerator& exahype::offloading::NoiseGenerator::getInstance() {
  static NoiseGenerator instance;
  return instance;
}

void exahype::offloading::NoiseGenerator::generateNoise() {
  //usleep(1000);
  if(_strategy!=nullptr)
	_strategy->generateNoise(tarch::parallel::Node::getInstance().getRank(), std::chrono::system_clock::now());
}

void exahype::offloading::NoiseGenerator::setStrategy(NoiseGenerationStrategy *strategy) {
  delete _strategy;
  _strategy = strategy;
}
