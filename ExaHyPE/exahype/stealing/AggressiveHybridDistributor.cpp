#if  defined(SharedTBB)  && defined(Parallel) && defined(DistributedStealing)
#include "AggressiveHybridDistributor.h"

#include <algorithm>
#include <numeric>
#include <cmath>

#include "tarch/multicore/Lock.h"
#include "tarch/parallel/Node.h"
#include "tarch/timing/Watch.h"

#include "exahype/stealing/StealingProfiler.h"
#include "exahype/stealing/PerformanceMonitor.h"
#include "exahype/stealing/StealingAnalyser.h"
#include "tarch/multicore/Core.h"
#include "tarch/multicore/tbb/Jobs.h"

tarch::logging::Log exahype::stealing::AggressiveHybridDistributor::_log( "exahype::stealing::AggressiveHybridDistributor" );

exahype::stealing::AggressiveHybridDistributor::AggressiveHybridDistributor() :
  _isEnabled(false),
  _temperature(0.5),
  _totalTasksOffloaded(0),
  _oldTotalTasksOffloaded(0),
  _useCCP(true)
{

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  _initialLoadPerRank      = new int[nnodes];
  _newLoadDistribution     = new int[nnodes];
  _idealTasksToOffload     = new int[nnodes];
  _tasksToOffload          = new int[nnodes];
  _remainingTasksToOffload = new std::atomic<int>[nnodes];
  _emergenciesPerRank      = new int[nnodes];
  _notOffloaded            = new int[nnodes];
 
  //for(int i=0; i<nnodes;i++) {
  //  _consumersPerRank[i] = std::max(1, tarch::multicore::Core::getInstance().getNumberOfThreads()-1);
 //   logInfo("AggressiveHybridDistributor()","weight "<<_consumersPerRank[i]<<" for rank "<<i);
  //}
 // _consumersPerRank[myRank] = std::max(1, tarch::multicore::Core::getInstance().getNumberOfThreads());
 // _consumersPerRank[0]     = std::max(1, tarch::multicore::Core::getInstance().getNumberOfThreads()-1);
 
  std::fill( &_remainingTasksToOffload[0], &_remainingTasksToOffload[nnodes], 0);
  std::fill( &_tasksToOffload[0], &_tasksToOffload[nnodes], 0);
  std::fill( &_idealTasksToOffload[0], &_idealTasksToOffload[nnodes], 0);
  std::fill( &_newLoadDistribution[0], &_newLoadDistribution[nnodes], 0);
  std::fill( &_initialLoadPerRank[0], &_initialLoadPerRank[nnodes], 0);
  std::fill( &_notOffloaded[0], &_notOffloaded[nnodes], 0);
  std::fill( &_emergenciesPerRank[0], &_emergenciesPerRank[nnodes], 0);
}

exahype::stealing::AggressiveHybridDistributor::~AggressiveHybridDistributor() {
  delete[] _emergenciesPerRank;
  delete[] _notOffloaded;
  delete[] _initialLoadPerRank;
  delete[] _newLoadDistribution;
  delete[] _tasksToOffload;
  delete[] _remainingTasksToOffload;
  delete[] _idealTasksToOffload;
}

void exahype::stealing::AggressiveHybridDistributor::enable() {
  _isEnabled = true;
}

void exahype::stealing::AggressiveHybridDistributor::disable() {
  _isEnabled = false;
}

void exahype::stealing::AggressiveHybridDistributor::computeIdealLoadDistribution(int enclaveCells, int skeletonCells) {

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  int *newLoadDist = new int[nnodes];
  std::fill(&newLoadDist[0], &newLoadDist[nnodes], 0);

  int totalCells   = enclaveCells + skeletonCells;
  MPI_Allgather(&totalCells, 1, MPI_INTEGER, _initialLoadPerRank, 1, MPI_INTEGER, MPI_COMM_WORLD);

#if defined(STEALING_USE_MASTER)
  int input_r=0, input_l=0;
  int output_r=0, output_l=0;
#else
  int input_r=1, input_l=0;
  int output_r=1, output_l=0;
#endif

  int total_l=0;
  total_l = std::accumulate(&_initialLoadPerRank[0], &_initialLoadPerRank[nnodes], total_l);

#if defined(STEALING_USE_MASTER)
  int avg_l = total_l / nnodes;
#else
  int avg_l = total_l / (nnodes-1);
  newLoadDist[0] = 0;
  _idealTasksToOffload[0] = 0;
#endif

  input_l = _initialLoadPerRank[input_r];
  output_l = _initialLoadPerRank[output_r];

  while(output_r<nnodes) {
    int target_load_out = avg_l;
    int target_load_in = avg_l;

    while(output_l<target_load_out) {
      int diff_l = target_load_out-output_l;

      if(output_r==input_r) {
        input_r++;
        assert(input_r<=nnodes);
        input_l = _initialLoadPerRank[input_r];
        continue;
      }

      int moveable = input_l-target_load_in;
      if(moveable>0) {
        int inc_l = std::min( diff_l, moveable );
        output_l += inc_l;
        input_l -= inc_l;
        //logInfo("performance monitor", " moving "<<inc_l<<" from rank "<<input_r<<" to rank "<<output_r);
        newLoadDist[output_r] = output_l;
        newLoadDist[input_r]  = input_l;

        if(input_r==myRank) {
          logInfo("computeIdeal","inc_l="<<inc_l);
          _idealTasksToOffload[output_r] = inc_l;
          stealing::StealingProfiler::getInstance().notifyTargetOffloadedTask(inc_l, output_r);
          //_tasksToOffload[output_r]= std::min(inc_l,1);
          //stealing::StealingProfiler::getInstance().notifyTargetOffloadedTask(std::min(inc_l,1), output_r);
        }
      }

      if(input_l <=target_load_in ) {
        input_r++;
        if(input_r<nnodes) {
          input_l = _initialLoadPerRank[input_r];
          target_load_in = avg_l;
        }
      }
    }
    output_r++;
    if(output_r<nnodes)
    output_l = _initialLoadPerRank[output_r];
  }

  std::string str="ideal load distribution ";
  for(int i=0;i<nnodes;i++) str=str+" , "+std::to_string(newLoadDist[i]);
  logInfo("aggressive distributor", str);
  str="ideal tasks to offload ";
  for(int i=0;i<nnodes;i++) str=str+" , "+std::to_string(_idealTasksToOffload[i]);
  logInfo("aggressive distributor", str);

  delete[] newLoadDist;
}

exahype::stealing::AggressiveHybridDistributor& exahype::stealing::AggressiveHybridDistributor::getInstance() {
  static AggressiveHybridDistributor aggressiveDist;
  return aggressiveDist;
}

void exahype::stealing::AggressiveHybridDistributor::printOffloadingStatistics() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  for(int i=0; i<nnodes; i++) {
    if(i==myRank)
      continue;
    logInfo("printOffloadingStatistics()", "target tasks to rank "<<i<<" ntasks "<<_tasksToOffload[i]<<" not offloaded "<<_notOffloaded[i]); 
    _notOffloaded[i]=0;
  }
  logInfo("printOffloadingStatistics()", "temperature value "<<_temperature );

}

void exahype::stealing::AggressiveHybridDistributor::resetRemainingTasksToOffload() {

  //logInfo("resetRemainingTasksToOffload()", "resetting remaining tasks to offload");

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  for(int i=0; i<nnodes; i++) {
    if(i==myRank)
      continue;
    _remainingTasksToOffload[i] = _tasksToOffload[i];
    _newLoadDistribution[i] = _tasksToOffload[i] + _initialLoadPerRank[i];
    //logInfo("resetRemainingTasksToOffload()", "to rank "<<i<<" ntasks "<<_tasksToOffload[i]<<" new load "<<_newLoadDistribution[i]); 
  }
}

void exahype::stealing::AggressiveHybridDistributor::handleEmergencyOnRank(int rank) {
  _emergenciesPerRank[rank]++;
  logInfo("handleEmergencyOnRank()","emergencies for rank:"<<_emergenciesPerRank[rank]);
}

void exahype::stealing::AggressiveHybridDistributor::updateLoadDistribution() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  if(!_isEnabled) {
    return;
  }
  logInfo("updateLoadDistribution()","total offloaded: "<<_totalTasksOffloaded<<" previous: "<<_oldTotalTasksOffloaded);

  if(_totalTasksOffloaded>0) {
    if(_totalTasksOffloaded-_oldTotalTasksOffloaded>0)
      _temperature = std::min(1.1, _temperature*1.1);
    else if(_totalTasksOffloaded-_oldTotalTasksOffloaded<0)
      _temperature = std::max(0.1, _temperature*0.9);
  }

  _oldTotalTasksOffloaded = _totalTasksOffloaded;

  if(_useCCP) {
    updateLoadDistributionCCP();
  }
  else
    updateLoadDistributionDiffusive();

  _totalTasksOffloaded = 0;
  for(int i=0; i<nnodes; i++) {
#ifdef DistributedStealingDisable
          assert(_tasksToOffload[i]==0);
#endif
    _tasksToOffload[i] = _tasksToOffload[i]-_emergenciesPerRank[i];
    _emergenciesPerRank[i] = 0;
    _totalTasksOffloaded += _tasksToOffload[i];

#ifdef DistributedStealingDisable
          assert(_tasksToOffload[i]==0);
#endif
  }

  resetRemainingTasksToOffload();

}

void exahype::stealing::AggressiveHybridDistributor::updateLoadDistributionCCP() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  int *waitingTimesSnapshot = new int[nnodes*nnodes];
 // const int* currentWaitingTimesSnapshot = exahype::stealing::PerformanceMonitor::getInstance().getWaitingTimesSnapshot();
  const int* currentWaitingTimesSnapshot = exahype::stealing::StealingAnalyser::getInstance().getFilteredWaitingTimesSnapshot();
  std::copy(&currentWaitingTimesSnapshot[0], &currentWaitingTimesSnapshot[nnodes*nnodes], waitingTimesSnapshot);

  //waitingTimesSnapshot[myRank] = waitingTime;

  //determine who is the fastest rank which is not blacklisted and who is a critical rank
  int *waitingRanks = new int[nnodes];
  bool *isWaitingForSomeone = new bool[nnodes];
  std::fill(waitingRanks, waitingRanks+nnodes, 0);

  int k = 0;

  for(int i=0; i<nnodes; i++) {
    bool waitingForSomeone = false;
    for(int j=0; j<nnodes; j++) {
      if(waitingTimesSnapshot[k+j]>0)
        logInfo("updateLoadDistribution()","rank "<<i<<" waiting for "<<waitingTimesSnapshot[k+j]<<" for rank "<<j);
      if(waitingTimesSnapshot[k+j]> exahype::stealing::StealingAnalyser::getInstance().getZeroThreshold()) {
        waitingRanks[j]++;   
        waitingForSomeone = true;
      }
    }
    isWaitingForSomeone[i]= waitingForSomeone;
    k+= nnodes;
  }

  int currentCriticalRank = -1;
  int currentLongestWaitTimeCritical = -1;

  k = 0;
  for(int i=0; i<nnodes; i++) {
    for(int j=0; j<nnodes; j++) {
      if(waitingTimesSnapshot[k+j]> exahype::stealing::StealingAnalyser::getInstance().getZeroThreshold()) {
         if(waitingTimesSnapshot[k+j]>currentLongestWaitTimeCritical
           && !isWaitingForSomeone[j]
           && waitingRanks[j]>0) {
             currentCriticalRank = j;
             currentLongestWaitTimeCritical = waitingTimesSnapshot[k+j];
         }
      }
    }
    k+= nnodes;
  }

  logInfo("updateLoadDistribution()", " current critical rank:"<<currentCriticalRank);
  
  bool isVictim = exahype::stealing::StealingManager::getInstance().isVictim();
  if(myRank == currentCriticalRank && !isVictim) {
     _useCCP = false;
     for(int i=0; i<nnodes; i++) {
       if(_idealTasksToOffload[i]>0) {
         //we have a potential victim rank
         if(!exahype::stealing::StealingManager::getInstance().isBlacklisted(i)) {
          //logInfo("updateLoadDistribution", "tasks for victim "<<i<<" before recomputation:"<<_tasksToOffload[i]<< " ideal: "<<_idealTasksToOffload[i]<< " temp "<<_temperature);
          //logInfo("updateLoadDistribution", "first : "<<(1.0-_temperature)*_tasksToOffload[i]<<" second:"<<_temperature*_idealTasksToOffload[i]);
          _tasksToOffload[i] = std::ceil(std::max((1.0-_temperature), 0.0)*_tasksToOffload[i] + _temperature*_idealTasksToOffload[i]);
#ifdef DistributedStealingDisable
          _tasksToOffload[i] = 0;
#endif
          //logInfo("updateLoadDistribution", "tasks for victim "<<i<<" after recomputation:"<<_tasksToOffload[i]);
         }
       }
     }
  }
  else if(_tasksToOffload[currentCriticalRank]>0) {
    _tasksToOffload[currentCriticalRank]--;
  }
 
  delete[] isWaitingForSomeone;
  delete[] waitingRanks;
  delete[] waitingTimesSnapshot;

}

void exahype::stealing::AggressiveHybridDistributor::updateLoadDistributionDiffusive() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  int *waitingTimesSnapshot = new int[nnodes*nnodes];
  const int* currentWaitingTimesSnapshot = exahype::stealing::StealingAnalyser::getInstance().getFilteredWaitingTimesSnapshot();
  std::copy(&currentWaitingTimesSnapshot[0], &currentWaitingTimesSnapshot[nnodes*nnodes], waitingTimesSnapshot);

  //determine who is the fastest rank which is not blacklisted and who is a critical rank
  int *waitingRanks = new int[nnodes];
  bool *isWaitingForSomeone = new bool[nnodes];
  std::fill(waitingRanks, waitingRanks+nnodes, 0);

  int currentLongestWaitTimeVictim = -1;
  int currentOptimalVictim = -1;
  int currentLongestWaitTimeCritical = -1;
  int currentCriticalRank = -1;
  int k = 0;

  for(int i=0; i<nnodes; i++) {
    bool waitingForSomeone = false;
    for(int j=0; j<nnodes; j++) {
      if(waitingTimesSnapshot[k+j]>0)
        logInfo("updateLoadDistribution()","rank "<<i<<" waiting for "<<waitingTimesSnapshot[k+j]<<" for rank "<<j);
//      if(waitingTimesSnapshot[k+j]>currentLongestWaitTime && !exahype::stealing::StealingManager::getInstance().isBlacklisted(i)) {
//        currentLongestWaitTime = waitingTimesSnapshot[k+j];
//        currentOptimalVictim = i;
//      }
      if(waitingTimesSnapshot[k+j]> exahype::stealing::StealingAnalyser::getInstance().getZeroThreshold()) {
        waitingRanks[j]++;   
        waitingForSomeone = true;
      }
    }
    isWaitingForSomeone[i]= waitingForSomeone;
    k+= nnodes;
  }

  k = 0;
  for(int i=0; i<nnodes; i++) {
    for(int j=0; j<nnodes; j++) {
      if(waitingTimesSnapshot[k+j]> exahype::stealing::StealingAnalyser::getInstance().getZeroThreshold()) {
         if(waitingTimesSnapshot[k+j]>currentLongestWaitTimeVictim 
           && !exahype::stealing::StealingManager::getInstance().isBlacklisted(i)
           && waitingRanks[i]==0) {
             currentOptimalVictim = i;
             currentLongestWaitTimeVictim = waitingTimesSnapshot[k+j];
         }
         if(waitingTimesSnapshot[k+j]>currentLongestWaitTimeCritical
           && !isWaitingForSomeone[j]
           && waitingRanks[j]>0) {
             currentCriticalRank = j;
             currentLongestWaitTimeCritical = waitingTimesSnapshot[k+j];
         } 
      }
    }
    k+= nnodes;
  }

  logInfo("updateLoadDistributionDiffusive()", "optimal victim: "<<currentOptimalVictim<<" critical rank:"<<currentCriticalRank);

  bool isVictim = exahype::stealing::StealingManager::getInstance().isVictim();
  if(myRank == currentCriticalRank && currentCriticalRank!=currentOptimalVictim) {
    if(!isVictim) {
      int currentTasksCritical = _initialLoadPerRank[currentCriticalRank];
      for(int i=0; i<nnodes; i++) {
        currentTasksCritical -= _tasksToOffload[i];
      }
      int currentTasksOptimal = _initialLoadPerRank[currentOptimalVictim]+_tasksToOffload[currentOptimalVictim];

      _tasksToOffload[currentOptimalVictim] += _temperature*0.5*(currentTasksCritical-currentTasksOptimal);
     }
  }
  else if(_tasksToOffload[currentCriticalRank]>0) {
    _tasksToOffload[currentCriticalRank]--;
    //logInfo("updateLoadDistribution()", "decrement, send "<<_tasksToOffload[criticalRank]<<" to rank "<<criticalRank );
  }

  resetRemainingTasksToOffload();

  delete[] isWaitingForSomeone;
  delete[] waitingRanks;
  delete[] waitingTimesSnapshot;
}

bool exahype::stealing::AggressiveHybridDistributor::selectVictimRank(int& victim) {

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
 
  //logInfo("selectVictimRank", "chose victim "<<victim<<" _remainingTasksToOffload "<<_remainingTasksToOffload[victim]);

  int threshold = 1+std::max(1, tarch::multicore::Core::getInstance().getNumberOfThreads()-1)*tarch::multicore::jobs::internal::_minimalNumberOfJobsPerConsumerRun;
  threshold = std::max(threshold, 20);


//  logInfo("selectVictimRank","waiting "<<tarch::multicore::jobs::getNumberOfWaitingBackgroundJobs()<<" criterion "<<threshold);
 
  if(tarch::multicore::jobs::getNumberOfWaitingBackgroundJobs()<
        threshold) {
    //logInfo("selectVictimRank", "number of running consumers: "<<tarch::multicore::jobs::internal::_numberOfRunningJobConsumerTasks.load()<<" max running "<<tarch::multicore::Core::getInstance().getNumberOfThreads()-1);
    _notOffloaded[victim]++;
    victim = myRank;
  }
  //if(victim!=myRank)
  // logInfo("selectVictimRank", "chose victim "<<victim<<" _remainingTasksToOffload "<<_remainingTasksToOffload[victim]);
  
  return victim != myRank;
}
#endif
