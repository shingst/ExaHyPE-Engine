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

#include "exahype/reactive/NoiseGenerationStrategyRoundRobin.h"

#include <unistd.h>
#include <chrono>
#include <cstdlib>
#include <string>

#include "../../../Peano/tarch/parallel/Node.h"
#include "exahype/reactive/OffloadingAnalyser.h"

#if defined(Parallel)
namespace exahype {
namespace reactive {

tarch::logging::Log  exahype::reactive::NoiseGenerationStrategyRoundRobin::_log( "exahype::reactive::NoiseGenerationStrategyRoundRobin" );


NoiseGenerationStrategyRoundRobin::NoiseGenerationStrategyRoundRobin() : _stepsBetweenDisturbance(1), _waitFractionTimestep(0.5){
}

NoiseGenerationStrategyRoundRobin::NoiseGenerationStrategyRoundRobin(int stepsBetweenDisturbance, double waitFraction)
 : _stepsBetweenDisturbance(stepsBetweenDisturbance), _waitFractionTimestep(waitFraction){
}

NoiseGenerationStrategyRoundRobin::~NoiseGenerationStrategyRoundRobin() {
}

void NoiseGenerationStrategyRoundRobin::generateNoiseIfActive(const int myRank) {
  pid_t pid = getpid();
  static int cnt = 0;
  static int nextRankToDisturb = 0;
  
  int nranks = tarch::parallel::Node::getInstance().getNumberOfNodes();

  if(nextRankToDisturb==myRank) {
    double timePerTimeStep = exahype::reactive::OffloadingAnalyser::getInstance().getTimePerTimeStep();
    double timeToWait = timePerTimeStep*_waitFractionTimestep;

    std::string call = " kill -STOP "+std::to_string(pid)+" ; sleep "+std::to_string(timeToWait)+"; kill -CONT "+std::to_string(pid);

    logDebug("generateNoise()", "running cmd "<<call<<std::endl);
    int ret = std::system( call.c_str() );

    if(ret==-1) {
      logError("generateNoise", "Sleep command could not be called, exiting...");
      exit(-1);
    }
  }
  cnt = cnt + 1;
  if(cnt==_stepsBetweenDisturbance) {
    nextRankToDisturb = (nextRankToDisturb+1) % nranks;
    cnt = 0;
  }
}

} /* namespace offloading */

} /* namespace exahype */

#endif
