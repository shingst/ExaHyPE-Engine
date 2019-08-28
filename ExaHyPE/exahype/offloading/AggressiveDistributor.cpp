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
#include "exahype/offloading/AggressiveDistributor.h"

#include <algorithm>
#include <numeric>

#include "tarch/multicore/Lock.h"
#include "tarch/parallel/Node.h"
#include "tarch/timing/Watch.h"

#include "exahype/offloading/OffloadingProfiler.h"
#include "exahype/offloading/PerformanceMonitor.h"
#include "tarch/multicore/Core.h"
#include "tarch/multicore/tbb/Jobs.h"

tarch::logging::Log exahype::offloading::AggressiveDistributor::_log( "exahype::offloading::AggressiveDistributor" );

exahype::offloading::AggressiveDistributor::AggressiveDistributor() :
  _zeroThreshold(2000*10),
  _isEnabled(false) {

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  _initialLoadPerRank      = new int[nnodes];
  _newLoadDistribution     = new int[nnodes];
  _idealTasksToOffload     = new int[nnodes];
  _tasksToOffload          = new int[nnodes];
  _remainingTasksToOffload = new std::atomic<int>[nnodes];

  _consumersPerRank        = new int[nnodes];
  _notOffloaded            = new int[nnodes];
 
  for(int i=0; i<nnodes;i++) {
    _consumersPerRank[i] = std::max(1, tarch::multicore::Core::getInstance().getNumberOfThreads()-1);
 //   logInfo("AggressiveDistributor()","weight "<<_consumersPerRank[i]<<" for rank "<<i);
  }
 // _consumersPerRank[myRank] = std::max(1, tarch::multicore::Core::getInstance().getNumberOfThreads());
 // _consumersPerRank[0]     = std::max(1, tarch::multicore::Core::getInstance().getNumberOfThreads()-1);
 
  std::fill( &_remainingTasksToOffload[0], &_remainingTasksToOffload[nnodes], 0);
  std::fill( &_tasksToOffload[0], &_tasksToOffload[nnodes], 0);
  std::fill( &_idealTasksToOffload[0], &_idealTasksToOffload[nnodes], 0);
  std::fill( &_newLoadDistribution[0], &_newLoadDistribution[nnodes], 0);
  std::fill( &_initialLoadPerRank[0], &_initialLoadPerRank[nnodes], 0);
  std::fill( &_notOffloaded[0], &_notOffloaded[nnodes], 0);
}

exahype::offloading::AggressiveDistributor::~AggressiveDistributor() {
  delete[] _notOffloaded;
  delete[] _initialLoadPerRank;
  delete[] _newLoadDistribution;
  delete[] _tasksToOffload;
  delete[] _remainingTasksToOffload;
  delete[] _consumersPerRank;
  delete[] _idealTasksToOffload;
}

void exahype::offloading::AggressiveDistributor::enable() {
  _isEnabled = true;
}

void exahype::offloading::AggressiveDistributor::disable() {
  _isEnabled = false;
}

void exahype::offloading::AggressiveDistributor::computeIdealLoadDistribution(int enclaveCells, int skeletonCells) {

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  int *newLoadDist = new int[nnodes];

  int totalCells   = enclaveCells + skeletonCells;
  MPI_Allgather(&totalCells, 1, MPI_INTEGER, _initialLoadPerRank, 1, MPI_INTEGER, MPI_COMM_WORLD);

  int input_r=0, input_l=0;
  int output_r=0, output_l=0;

  int total_l=0;
  total_l = std::accumulate(&_initialLoadPerRank[0], &_initialLoadPerRank[nnodes], total_l);

  int total_consumers = 0;
  total_consumers = std::accumulate(&_consumersPerRank[0], &_consumersPerRank[nnodes], total_consumers);


  int avg_l_per_consumer = 0;
  avg_l_per_consumer = total_l / total_consumers;

  input_l = _initialLoadPerRank[input_r];
  output_l = _initialLoadPerRank[output_r];

  while(output_r<nnodes) {
    int target_load_out = _consumersPerRank[output_r] * avg_l_per_consumer;
    int target_load_in = _consumersPerRank[input_r] * avg_l_per_consumer;

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
          _idealTasksToOffload[output_r] = inc_l;
          offloading::OffloadingProfiler::getInstance().notifyTargetOffloadedTask(inc_l, output_r);
          //_tasksToOffload[output_r]= std::min(inc_l,1);
          //offloading::OffloadingProfiler::getInstance().notifyTargetOffloadedTask(std::min(inc_l,1), output_r);
        }
      }

      if(input_l <=target_load_in ) {
        input_r++;
        if(input_r<nnodes) {
          input_l = _initialLoadPerRank[input_r];
          target_load_in = _consumersPerRank[input_r] * avg_l_per_consumer;
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
  str="tasks to offload ";
  for(int i=0;i<nnodes;i++) str=str+" , "+std::to_string(_idealTasksToOffload[i]);
  logInfo("aggressive distributor", str);

  delete[] newLoadDist;
}

exahype::offloading::AggressiveDistributor& exahype::offloading::AggressiveDistributor::getInstance() {
  static AggressiveDistributor aggressiveDist;
  return aggressiveDist;
}

void exahype::offloading::AggressiveDistributor::printOffloadingStatistics() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  for(int i=0; i<nnodes; i++) {
    if(i==myRank)
      continue;
    logInfo("printOffloadingStatistics()", "target tasks to rank "<<i<<" ntasks "<<_tasksToOffload[i]<<" not offloaded "<<_notOffloaded[i]); 
    _notOffloaded[i]=0;
  }
}

void exahype::offloading::AggressiveDistributor::resetRemainingTasksToOffload() {

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

void exahype::offloading::AggressiveDistributor::handleEmergencyOnRank(int rank) {
  _tasksToOffload[rank]--;
  logInfo("handleEmergencyOnRank()","decrement for rank:"<<rank<<" tasks to offload "<<_tasksToOffload[rank]);
  //_remainingTasksToOffload[rank] = _tasksToOffload[rank];
}

void exahype::offloading::AggressiveDistributor::updateLoadDistribution() {

  if(!_isEnabled) {
    return;
  }

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  double *waitingTimesSnapshot = new double[nnodes*nnodes];
  const double* currentWaitingTimesSnapshot = exahype::offloading::PerformanceMonitor::getInstance().getWaitingTimesSnapshot();
  std::copy(&currentWaitingTimesSnapshot[0], &currentWaitingTimesSnapshot[nnodes*nnodes], waitingTimesSnapshot);

  //waitingTimesSnapshot[myRank] = waitingTime;

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
      if(waitingTimesSnapshot[k+j]>0)
        logInfo("updateLoadDistribution()","rank "<<i<<" waiting for "<<waitingTimesSnapshot[k+j]<<" for rank "<<j);
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

  bool isVictim = exahype::offloading::OffloadingManager::getInstance().isVictim();
  if(myRank == criticalRank && criticalRank!=currentOptimalVictim) {
    if(!isVictim) {
      int currentTasksCritical = _initialLoadPerRank[criticalRank];
      for(int i=0; i<nnodes; i++) {
        currentTasksCritical -= _tasksToOffload[i];
      }
      int currentTasksOptimal = _initialLoadPerRank[currentOptimalVictim]+_tasksToOffload[currentOptimalVictim];

      _tasksToOffload[currentOptimalVictim] += 0.5*(currentTasksCritical-currentTasksOptimal);
      //logInfo("updateLoadDistribution()", "I am a critical rank, increment, send "<<_tasksToOffload[currentOptimalVictim]<<" to rank "<<currentOptimalVictim );
    }
  }
  else if(_tasksToOffload[criticalRank]>0) {
    _tasksToOffload[criticalRank]--;
    //logInfo("updateLoadDistribution()", "decrement, send "<<_tasksToOffload[criticalRank]<<" to rank "<<criticalRank );
  }

//  if(myRank == slowestRank) {
//    if(!isVictim && !emergencyTriggered && *std::min_element(&waitingTimesSnapshot[0], &waitingTimesSnapshot[nnodes])<_zeroThreshold) {
//      for(int i=0; i<nnodes; i++) {
//        if(i!=myRank) { 
//          _tasksToOffload[i] = 0.5* (_tasksToOffload[i] + _idealTasksToOffload[i]);
//          logInfo("updateLoadDistribution()", "I am a critical rank, increment, send "<<_tasksToOffload[i]<<" to rank "<<i );          
//        }
//      }
//    }
//    else if(emergencyTriggered){
//      //TODO: maybe step back a bit
//      //exahype::offloading::OffloadingManager::getInstance().resetEmergency();
//    }   
//  }

  resetRemainingTasksToOffload();

  delete[] isWaitingForSomeone;
  delete[] waitingRanks;
  delete[] waitingTimesSnapshot;

}

bool exahype::offloading::AggressiveDistributor::selectVictimRank(int& victim) {

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();
  victim = myRank;

  static std::atomic<int> rank_cnt (0);

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
