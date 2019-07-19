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

#if  defined(SharedTBB)  && defined(Parallel) && defined(DistributedOffloading)
#include "exahype/offloading/DiffusiveDistributor.h"

#include <algorithm>
#include <numeric>

#include "tarch/multicore/Lock.h"
#include "tarch/parallel/Node.h"

#include "exahype/offloading/OffloadingProfiler.h"
#include "exahype/offloading/PerformanceMonitor.h"
#include "tarch/multicore/Core.h"
#include "tarch/multicore/Jobs.h"

tarch::logging::Log exahype::offloading::DiffusiveDistributor::_log( "exahype::offloading::DiffusiveDistributor" );

exahype::offloading::DiffusiveDistributor::DiffusiveDistributor() :
  _zeroThreshold(10*2000)
{
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  _tasksToOffload          = new int[nnodes];
  _remainingTasksToOffload = new std::atomic<int>[nnodes];

  std::fill( &_remainingTasksToOffload[0], &_remainingTasksToOffload[nnodes], 0);
  std::fill( &_tasksToOffload[0], &_tasksToOffload[nnodes], 0);
}

exahype::offloading::DiffusiveDistributor::~DiffusiveDistributor() {
  delete[] _tasksToOffload;
  delete[] _remainingTasksToOffload;
}

void exahype::offloading::DiffusiveDistributor::updateZeroThreshold(int threshold) {
  _zeroThreshold = threshold;
}

void exahype::offloading::DiffusiveDistributor::updateLoadDistribution() {

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  double *waitingTimesSnapshot = new double[nnodes*nnodes];
  const double* currentWaitingTimesSnapshot = exahype::offloading::PerformanceMonitor::getInstance().getWaitingTimesSnapshot();
  std::copy(&currentWaitingTimesSnapshot[0], &currentWaitingTimesSnapshot[nnodes*nnodes], waitingTimesSnapshot);

  //waitingTimesSnapshot[myRank] = waitingTime;

  bool isVictim = exahype::offloading::OffloadingManager::getInstance().isVictim();
  bool emergencyTriggered = exahype::offloading::OffloadingManager::getInstance().isEmergencyTriggered();

  logInfo("updateLoadDistribution()", " isVictim: "<<isVictim<<" emergency event: "<<emergencyTriggered);

  //determine who is the fastest rank which is not blacklisted and who is a critical rank
  int *waitingRanks = new int[nnodes];
  bool *isWaitingForSomeone = new bool[nnodes];
  std::fill(waitingRanks, waitingRanks+nnodes, 0);

  double currentLongestWaitTime = -1;
  int currentOptimalVictim = -1;
  int k = 0;

  for(int i=0; i<nnodes; i++) {
    bool waitingForSomeone = false;
    for(int j=0; j<nnodes; j++) {
      //logInfo("updateLoadDistribution()","rank "<<i<<" waiting for "<<waitingTimesSnapshot[k+j]<<" for rank "<<j);
      if(waitingTimesSnapshot[k+j]>currentLongestWaitTime && !exahype::offloading::OffloadingManager::getInstance().isBlacklisted(i)) {
        currentLongestWaitTime = waitingTimesSnapshot[k+j];
        currentOptimalVictim = i;
      }
      if(waitingTimesSnapshot[k+j]>_zeroThreshold) {
        waitingRanks[j]++;   
        waitingForSomeone = true;
      }
    }
    isWaitingForSomeone[i]= waitingForSomeone;
    k+= nnodes;
  }

  int criticalRank = -1;
  for(int i=0; i<nnodes; i++) {
    if(!isWaitingForSomeone[i] && waitingRanks[i]>0) {
      criticalRank = i; 
      break;
    }
  }
 
  logInfo("updateLoadDistribution()", "optimal victim: "<<currentOptimalVictim<<" critical rank:"<<criticalRank);

  /*if(criticalRank==currentOptimalVictim) {
    k= 0;
    for(int i=0; i<nnodes; i++) {
      for(int j=0; j<nnodes; j++) {
      logInfo("updateLoadDistribution()","critical==optimal; rank "<<i<<" waiting for "<<waitingTimesSnapshot[k+j]<<" for rank "<<j);  
      }
      k+= nnodes;
    }
  }*/

  if(myRank == criticalRank && criticalRank!=currentOptimalVictim) {
    if(!isVictim) {
      _tasksToOffload[currentOptimalVictim]++;
      logInfo("updateLoadDistribution()", "I am a critical rank, increment, send "<<_tasksToOffload[currentOptimalVictim]<<" to rank "<<currentOptimalVictim );
    }
    /*else if(emergencyTriggered){
      logInfo("updateLoadDistribution()", "I was a critical rank, but emergency event happened."); //TODO: emergency handling in separate method
      for(int i=0; i<nnodes; i++) {
        if(i!=myRank && _tasksToOffload[i]>0) {
          _tasksToOffload[i]--;
          logInfo("updateLoadDistribution()", "decrement, send "<<_tasksToOffload[i]<<" to rank "<<i );
        }
      }
      //exahype::offloading::OffloadingManager::getInstance().resetEmergency();
    } */  
  }
  else if(_tasksToOffload[criticalRank]>0) {
    _tasksToOffload[criticalRank]--;
    logInfo("updateLoadDistribution()", "decrement, send "<<_tasksToOffload[criticalRank]<<" to rank "<<criticalRank );
  }

  for(int i=0; i<nnodes; i++) 
    _remainingTasksToOffload[i] = _tasksToOffload[i];

  delete[] isWaitingForSomeone;
  delete[] waitingRanks;
  delete[] waitingTimesSnapshot;

}

void exahype::offloading::DiffusiveDistributor::handleEmergencyOnRank(int rank) {
   logInfo("handleEmergencyOnRank()", "Emergency event happened for rank "<<rank);
   _tasksToOffload[rank]--;
   logInfo("handleEmergencyOnRank()", "decrement, send "<<_tasksToOffload[rank]<<" to rank "<<rank); 
}

exahype::offloading::DiffusiveDistributor& exahype::offloading::DiffusiveDistributor::getInstance() {
  static DiffusiveDistributor diffusiveDist;
  return diffusiveDist;
}

bool exahype::offloading::DiffusiveDistributor::selectVictimRank(int& victim) {

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();
  victim = myRank;

  static std::atomic<int> rank_cnt = 0;

  int l_rank = rank_cnt;

  for(int i=0; i<nnodes; i++) {
    if(l_rank!=myRank && _remainingTasksToOffload[l_rank].fetch_sub(1)>0) {
      victim = l_rank;
      l_rank = (l_rank + 1)%nnodes;
      break;
    }
    else
    _remainingTasksToOffload[l_rank]=0;
    l_rank = (l_rank + 1)%nnodes;
  }
  rank_cnt=l_rank;

  //if(victim!=myRank)
  // logInfo("selectVictimRank", "chose victim "<<victim<<" _remainingTasksToOffload "<<_remainingTasksToOffload[victim]);
  
  return victim != myRank;
}
#endif
