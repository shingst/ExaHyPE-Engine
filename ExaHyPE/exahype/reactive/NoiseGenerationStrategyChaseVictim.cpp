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

#include "exahype/reactive/NoiseGenerationStrategyChaseVictim.h"

#include <unistd.h>
#include <cstdlib>
#include <string>

#include "tarch/parallel/Node.h"

#include "exahype/reactive/AggressiveHybridDistributor.h"
#include "exahype/reactive/OffloadingAnalyser.h"
#include "exahype/reactive/ReactiveContext.h"


namespace exahype {
namespace reactive {

tarch::logging::Log  exahype::reactive::NoiseGenerationStrategyChaseVictim::_log( "exahype::reactive::NoiseGenerationStrategyChaseVictim" );

NoiseGenerationStrategyChaseVictim::NoiseGenerationStrategyChaseVictim() : _scalingFactor(0.5), _baseNoise(1){
}

NoiseGenerationStrategyChaseVictim::NoiseGenerationStrategyChaseVictim(double scalingFactor, double baseNoise)
 : _scalingFactor(scalingFactor), _baseNoise(baseNoise) {
   assertion(scalingFactor>=0);
}

NoiseGenerationStrategyChaseVictim::~NoiseGenerationStrategyChaseVictim() {
}

void NoiseGenerationStrategyChaseVictim::generateNoiseIfActive(const int myRank ) {
//does not make sense without parallelism
#if defined(Parallel)
  pid_t pid = getpid();
  static const int stepsBetweenDisturbance = 10;
  static int stepsRemainingBeforeDisturbance = stepsBetweenDisturbance;

  stepsRemainingBeforeDisturbance--;
  bool isTimeToDisturb = (stepsRemainingBeforeDisturbance==0);

  if (ReactiveContext::getInstance().isVictim()
    && isTimeToDisturb) {

    double timeToWait = _baseNoise*_scalingFactor;
    std::string call = " kill -STOP "+std::to_string(pid)+" ; sleep "+std::to_string(timeToWait)+"; kill -CONT "+std::to_string(pid);

    logDebug("generateNoise()", "running cmd "<<call<<std::endl);
    int ret = std::system( call.c_str() );

    if(ret==-1) {
      logError("generateNoise", "Sleep command could not be called, exiting...");
      exit(-1);
    }
  }

  if(isTimeToDisturb) {
    stepsRemainingBeforeDisturbance = stepsBetweenDisturbance;
  }
#endif
}

} /* namespace offloading */
} /* namespace exahype */
