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


NoiseGenerationStrategyRoundRobin::NoiseGenerationStrategyRoundRobin() {
	// TODO Auto-generated constructor stub

}

NoiseGenerationStrategyRoundRobin::~NoiseGenerationStrategyRoundRobin() {
	// TODO Auto-generated destructor stub
}

void NoiseGenerationStrategyRoundRobin::generateNoise(int rank, std::chrono::system_clock::time_point timestamp) {
  pid_t pid = getpid();
  static int cnt = 0;
  
  int nranks = tarch::parallel::Node::getInstance().getNumberOfNodes();

  if(cnt==rank) {
    double timePerTimeStep = exahype::offloading::OffloadingAnalyser::getInstance().getTimePerTimeStep();
    double timeToWait = timePerTimeStep/2;

    std::string call = " kill -STOP "+std::to_string(pid)+" ; sleep "+std::to_string(timeToWait)+"; kill -CONT "+std::to_string(pid);

    logInfo("generateNoise()", "running cmd "<<call<<std::endl);
    std::system( call.c_str() );
  }
  cnt = (cnt+1)%nranks;
  //usleep(10000000);
}

} /* namespace offloading */
} /* namespace exahype */
