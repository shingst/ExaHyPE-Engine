#if  defined(SharedTBB)  && defined(Parallel) && defined(DistributedStealing)
#include "DiffusiveDistributor.h"

#include <algorithm>
#include <numeric>

#include "tarch/multicore/Lock.h"
#include "tarch/parallel/Node.h"

#include "exahype/stealing/StealingProfiler.h"
#include "exahype/stealing/PerformanceMonitor.h"
#include "tarch/multicore/Core.h"
#include "tarch/multicore/Jobs.h"

tarch::logging::Log exahype::stealing::DiffusiveDistributor::_log( "exahype::stealing::DiffusiveDistributor" );

exahype::stealing::DiffusiveDistributor::DiffusiveDistributor() :
  _zeroThreshold(10*2000)
{
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  _tasksToOffload          = new int[nnodes];
  _remainingTasksToOffload = new std::atomic<int>[nnodes];

  std::fill( &_remainingTasksToOffload[0], &_remainingTasksToOffload[nnodes], 0);
  std::fill( &_tasksToOffload[0], &_tasksToOffload[nnodes], 0);
}

exahype::stealing::DiffusiveDistributor::~DiffusiveDistributor() {
  delete[] _tasksToOffload;
  delete[] _remainingTasksToOffload;
}

void exahype::stealing::DiffusiveDistributor::updateZeroThreshold(int threshold) {
  _zeroThreshold = threshold;
}

void exahype::stealing::DiffusiveDistributor::updateLoadDistribution(int currentLoad) {

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  int *loadSnapshot = new int[nnodes];
  const int* currentLoadSnapshot = exahype::stealing::PerformanceMonitor::getInstance().getCurrentLoadSnapshot();
  std::copy(&currentLoadSnapshot[0], &currentLoadSnapshot[nnodes], loadSnapshot);

  loadSnapshot[myRank] = currentLoad;

  bool isVictim = exahype::stealing::StealingManager::getInstance().isVictim();
  bool emergencyTriggered = exahype::stealing::StealingManager::getInstance().isEmergencyTriggered();

  logInfo("updateLoadDistribution()", "current maximum wait time "<<currentLoad
                                      <<" isVictim: "<<isVictim<<" emergency event: "<<emergencyTriggered);

  //determine who is the fastest rank which is not blacklisted
  //int fastestRank = std::distance(&loadSnapshot[0], std::max_element(&loadSnapshot[0], &loadSnapshot[nnodes]));
  int currentLongestWaitTime = -1;
  int currentFastestRank = -1;
  for(int i=0; i<nnodes; i++) {
    if(loadSnapshot[i]>currentLongestWaitTime && !exahype::stealing::StealingManager::getInstance().isBlacklisted(i)) {
      currentLongestWaitTime = loadSnapshot[i];
      currentFastestRank = i;
    }
  }

  //determine who is slowest
  int slowestRank = std::distance(&loadSnapshot[0], std::min_element(&loadSnapshot[0], &loadSnapshot[nnodes]));
  
  logInfo("updateLoadDistribution()", "fastest: "<<currentFastestRank<<" slowest:"<<slowestRank);

  if(myRank == slowestRank && slowestRank!=currentFastestRank) {
    if(!isVictim && *std::min_element(&loadSnapshot[0], &loadSnapshot[nnodes])<_zeroThreshold) {
      _tasksToOffload[currentFastestRank]++;
      logInfo("updateLoadDistribution()", "I am a critical rank, increment, send "<<_tasksToOffload[currentFastestRank]<<" to rank "<<currentFastestRank );
    }
    /*else if(emergencyTriggered){
      logInfo("updateLoadDistribution()", "I was a critical rank, but emergency event happened."); //TODO: emergency handling in separate method
      for(int i=0; i<nnodes; i++) {
        if(i!=myRank && _tasksToOffload[i]>0) {
          _tasksToOffload[i]--;
          logInfo("updateLoadDistribution()", "decrement, send "<<_tasksToOffload[i]<<" to rank "<<i );
        }
      }
      //exahype::stealing::StealingManager::getInstance().resetEmergency();
    } */  
  }
  else if(_tasksToOffload[slowestRank]>0) {
    _tasksToOffload[slowestRank]--;
    logInfo("updateLoadDistribution()", "decrement, send "<<_tasksToOffload[slowestRank]<<" to rank "<<slowestRank );
  }

  for(int i=0; i<nnodes; i++) 
    _remainingTasksToOffload[i] = _tasksToOffload[i];

  delete[] loadSnapshot;

}

void exahype::stealing::DiffusiveDistributor::handleEmergencyOnRank(int rank) {
   logInfo("handleEmergencyOnRank()", "Emergency event happened for rank "<<rank);
   _tasksToOffload[rank]--;
   logInfo("handleEmergencyOnRank()", "decrement, send "<<_tasksToOffload[rank]<<" to rank "<<rank); 
}

exahype::stealing::DiffusiveDistributor& exahype::stealing::DiffusiveDistributor::getInstance() {
  static DiffusiveDistributor diffusiveDist;
  return diffusiveDist;
}

bool exahype::stealing::DiffusiveDistributor::selectVictimRank(int& victim) {

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
