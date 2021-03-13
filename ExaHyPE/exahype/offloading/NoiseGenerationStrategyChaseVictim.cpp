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

#include "NoiseGenerationStrategyChaseVictim.h"
#include "OffloadingAnalyser.h"
#include "tarch/parallel/Node.h"
#include "AggressiveHybridDistributor.h"
#include "OffloadingManager.h"

#include <unistd.h>
#include <cstdlib>
#include <string>

namespace exahype {
namespace offloading {

tarch::logging::Log  exahype::offloading::NoiseGenerationStrategyChaseVictim::_log( "exahype::offloading::NoiseGenerationStrategyChaseVictim" );


NoiseGenerationStrategyChaseVictim::NoiseGenerationStrategyChaseVictim() : _factor(0.5), _baseNoise(1){
}

NoiseGenerationStrategyChaseVictim::NoiseGenerationStrategyChaseVictim(double factor, double baseNoise)
 : _factor(factor), _baseNoise(baseNoise) {
}

NoiseGenerationStrategyChaseVictim::~NoiseGenerationStrategyChaseVictim() {
}

void NoiseGenerationStrategyChaseVictim::generateNoise(int rank, std::chrono::system_clock::time_point timestamp) {
  pid_t pid = getpid();
  static int cnt = 0;


  static bool triggeredVictim = false;

  static const int kRecover = 10;
  static int kLastTriggered = -1;

  if(triggeredVictim && (cnt-kLastTriggered)==kRecover){
    triggeredVictim = false;
    kLastTriggered = -1;
  }

  if (OffloadingManager::getInstance().isVictim() && !triggeredVictim) {
    //double timePerTimeStep = exahype::offloading::OffloadingAnalyser::getInstance().getTimePerTimeStep();
    double timeToWait = _baseNoise*_factor;
    std::string call = " kill -STOP "+std::to_string(pid)+" ; sleep "+std::to_string(timeToWait)+"; kill -CONT "+std::to_string(pid);

    logDebug("generateNoise()", "running cmd "<<call<<std::endl);
    int ret = std::system( call.c_str() );

    if(ret==-1) {
      logError("generateNoise", "Sleep command could not be called, exiting...");
      exit(-1);
    }

    if(OffloadingManager::getInstance().isVictim()) {
      triggeredVictim = true;
      kLastTriggered = cnt;
    }
  }

  cnt++;
}

} /* namespace offloading */
} /* namespace exahype */
