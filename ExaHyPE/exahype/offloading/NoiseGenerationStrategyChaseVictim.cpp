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


NoiseGenerationStrategyChaseVictim::NoiseGenerationStrategyChaseVictim() : _frequency(1), _factor(0.5){
}

NoiseGenerationStrategyChaseVictim::NoiseGenerationStrategyChaseVictim(int frequency, double factor)
 : _frequency(frequency), _factor(factor){
}

NoiseGenerationStrategyChaseVictim::~NoiseGenerationStrategyChaseVictim() {
}

void NoiseGenerationStrategyChaseVictim::generateNoise(int rank, std::chrono::system_clock::time_point timestamp) {
  pid_t pid = getpid();
  static int cnt = 0;
  static int phase_cnt = 0;
  
  /*int nranks = tarch::parallel::Node::getInstance().getNumberOfNodes();

  int current_critical_rank = AggressiveHybridDistributor::getInstance().getCurrentCriticalRank();
  int current_optimal_victim = AggressiveHybridDistributor::getInstance().getCurrentOptimalVictim();

 // int noise_victim = (current_critical_rank!=-1 && current_optimal_victim!=-1) ? current_optimal_victim : phase_cnt;

  if((phase_cnt==rank-1 && rank!=current_critical_rank)  || OffloadingManager::getInstance().isVictim()) {
    double timePerTimeStep = exahype::offloading::OffloadingAnalyser::getInstance().getTimePerTimeStep();
    double timeToWait = timePerTimeStep*_factor;

    std::string call = " kill -STOP "+std::to_string(pid)+" ; sleep "+std::to_string(timeToWait)+"; kill -CONT "+std::to_string(pid);

    logInfo("generateNoise()", "running cmd "<<call<<std::endl);
    std::system( call.c_str() );
  }
  cnt = cnt + 1;
  if(cnt==_frequency) {
	  phase_cnt = (phase_cnt+1) % (nranks-1);
	  cnt = 0;
  }*/

  if(rank==1) {
    double timePerTimeStep = exahype::offloading::OffloadingAnalyser::getInstance().getTimePerTimeStep();
    double timeToWait = timePerTimeStep*_factor;

    std::string call = " kill -STOP "+std::to_string(pid)+" ; sleep "+std::to_string(timeToWait)+"; kill -CONT "+std::to_string(pid);

    logInfo("generateNoise()", "running cmd "<<call<<std::endl);
    std::system( call.c_str() );
  }

  //usleep(10000000);
}

void NoiseGenerationStrategyChaseVictim::generateNoiseSTP(int rank, std::chrono::system_clock::time_point timestamp) {
  //pid_t pid = getpid();
  //static int cnt = 0;
  //static int phase_cnt = 0;

  int nranks = tarch::parallel::Node::getInstance().getNumberOfNodes();

  int current_critical_rank = AggressiveHybridDistributor::getInstance().getCurrentCriticalRank();
  int current_optimal_victim = AggressiveHybridDistributor::getInstance().getCurrentOptimalVictim();

 // int noise_victim = (current_critical_rank!=-1 && current_optimal_victim!=-1) ? current_optimal_victim : phase_cnt;

  if( OffloadingManager::getInstance().isVictim()) {
    double timePerSTP = exahype::offloading::OffloadingAnalyser::getInstance().getTimePerSTP();
    double timeToWait = timePerSTP*_factor*1e06;

    logInfo("generateNoiseSTP()", "sleeping "<<timeToWait<<std::endl);
    //usleep(timeToWait*1e6);

    while( std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now()-timestamp).count()<timeToWait) {
    }

    logInfo("generateNoiseSTP()", "slept "<<std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now()-timestamp).count()<<std::endl);
    //std::string call = " kill -STOP "+std::to_string(pid)+" ; sleep "+std::to_string(timeToWait)+"; kill -CONT "+std::to_string(pid);

    //logInfo("generateNoise()", "running cmd "<<call<<std::endl);
    //std::system( call.c_str() );
  }

  //usleep(10000000);
}

} /* namespace offloading */
} /* namespace exahype */
