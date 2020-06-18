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
#include "OffloadingAnalyser.h"
#include "tarch/parallel/Node.h"

#include <unistd.h>
#include <cstdlib>
#include <string>

namespace exahype {
namespace offloading {

tarch::logging::Log  exahype::offloading::NoiseGenerationStrategyRoundRobin::_log( "exahype::offloading::NoiseGenerationStrategyRoundRobin" );


NoiseGenerationStrategyRoundRobin::NoiseGenerationStrategyRoundRobin() : _frequency(1), _factor(0.5){
}

NoiseGenerationStrategyRoundRobin::NoiseGenerationStrategyRoundRobin(int frequency, double factor)
 : _frequency(frequency), _factor(factor){
}

NoiseGenerationStrategyRoundRobin::~NoiseGenerationStrategyRoundRobin() {
}

void NoiseGenerationStrategyRoundRobin::generateNoise(int rank, std::chrono::system_clock::time_point timestamp) {
  pid_t pid = getpid();
  static int cnt = 0;
  static int phase_cnt = 0;
  
  int nranks = tarch::parallel::Node::getInstance().getNumberOfNodes();

  if(phase_cnt==rank) {
    double timePerTimeStep = exahype::offloading::OffloadingAnalyser::getInstance().getTimePerTimeStep();
    double timeToWait = timePerTimeStep*_factor;

    std::string call = " kill -STOP "+std::to_string(pid)+" ; sleep "+std::to_string(timeToWait)+"; kill -CONT "+std::to_string(pid);

    logInfo("generateNoise()", "running cmd "<<call<<std::endl);
    std::system( call.c_str() );
  }
  cnt = cnt + 1;
  if(cnt==_frequency) {
    phase_cnt = (phase_cnt+1) % nranks;
    cnt = 0;
  }
}


void NoiseGenerationStrategyRoundRobin::generateNoiseSTP(int rank, std::chrono::system_clock::time_point timestamp) {}
  //Todo: nothing here yet
} /* namespace offloading */

} /* namespace exahype */
