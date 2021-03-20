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

#if  defined(SharedTBB)  && defined(Parallel)
#include "../reactive/DynamicDistributor.h"

#include <algorithm>
#include <numeric>

#include "tarch/multicore/Lock.h"
#include "tarch/parallel/Node.h"
#include "tarch/timing/Watch.h"

#include "../reactive/OffloadingProfiler.h"
#include "tarch/multicore/Core.h"

tarch::logging::Log exahype::reactive::DynamicDistributor::_log( "exahype::reactive::DynamicDistributor" );

exahype::reactive::DynamicDistributor::DynamicDistributor() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  _tasksToOffload = new int[nnodes];
  _consumersPerRank = new int[nnodes];
  _remainingTasksToOffload = new std::atomic<int>[nnodes];

  std::fill( &_consumersPerRank[0], &_consumersPerRank[nnodes], 0);
  for(int i=0; i<nnodes;i++) {
    _consumersPerRank[i] = std::max(1, tarch::multicore::Core::getInstance().getNumberOfThreads()-1);
 //   logInfo("AggressiveDistributor()","weight "<<_consumersPerRank[i]<<" for rank "<<i);
  }

  std::fill( &_tasksToOffload[0], &_tasksToOffload[nnodes], 0);
  std::fill( &_remainingTasksToOffload[0], &_remainingTasksToOffload[nnodes], 0);
}

exahype::reactive::DynamicDistributor& exahype::reactive::DynamicDistributor::getInstance() {
  static DynamicDistributor dynamicDist;
  return dynamicDist;
}

exahype::reactive::DynamicDistributor::~DynamicDistributor() {
  delete[] _tasksToOffload;
  delete[] _consumersPerRank;
  delete[] _remainingTasksToOffload;

}

void exahype::reactive::DynamicDistributor::computeNewLoadDistribution(int *currentLoadSnapshot) {
  int input_r=0, input_l=0;
  int output_r=0, output_l=0;

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  //int * newLoadDist = new int[nnodes];

  int total_l=0;
  total_l = std::accumulate(&currentLoadSnapshot[0], &currentLoadSnapshot[nnodes], total_l);
  assertion(total_l>=0);

  int total_consumers = 0;
  total_consumers = std::accumulate(&_consumersPerRank[0], &_consumersPerRank[nnodes], total_consumers);

  int avg_l_per_consumer = 0;
  avg_l_per_consumer = total_l / total_consumers;

  if(avg_l_per_consumer == 0) return;
  assertion(avg_l_per_consumer>=0);

  input_l = currentLoadSnapshot[input_r];
  output_l = currentLoadSnapshot[output_r];
  assertion(input_l>=0);
  assertion(output_l>=0);

  while(output_r<nnodes) {
    if(input_r==nnodes) {
      break;
    }

    assertion(output_r<nnodes);
    assertion(input_r<nnodes);
    int target_load_out = _consumersPerRank[output_r] * avg_l_per_consumer;
    int target_load_in  = _consumersPerRank[input_r] * avg_l_per_consumer;
    assertion(target_load_out>=0);
    assertion(target_load_in>=0);

    while(output_l<target_load_out) {

      assertion(output_l>=0);
      int diff_l = target_load_out-output_l;

      if(output_r==input_r) {
        input_r++;
        assertion(input_r<nnodes);
        input_l = currentLoadSnapshot[input_r];
        continue;
      }

      int moveable = input_l-target_load_in;

      if(moveable>0) {
        int inc_l = std::min( diff_l, moveable );
        output_l += inc_l;
        input_l -= inc_l;
        //logInfo("performance monitor", " moving "<<inc_l<<" from rank "<<input_r<<" to rank "<<output_r);
        //newLoadDist[output_r]= output_l;
        //newLoadDist[input_r]= input_l;

        if(input_r==myRank) {
          //_tasksToOffload[output_r]= inc_l;
          //offloading::OffloadingProfiler::getInstance().notifyTargetOffloadedTask(inc_l, output_r);
          _tasksToOffload[output_r]= std::min(inc_l,1);
          reactive::OffloadingProfiler::getInstance().notifyTargetOffloadedTask(std::min(inc_l,1), output_r);
        }
      }

      if(input_l <=target_load_in ) {
        input_r++;
        if(input_r<nnodes) {
          input_l = currentLoadSnapshot[input_r];
          target_load_in = _consumersPerRank[input_r] * avg_l_per_consumer;
        }
      }
    }
    output_r++;
    if(output_r<nnodes)
    output_l = currentLoadSnapshot[output_r];
  }

//  std::string str="new (intended!) load distribution ";
//  for(int i=0;i<nnodes;i++) str=str+" , "+std::to_string(newLoadDist[i]);
//  logInfo("dynamic distributor", str);
//  str="tasks to offload ";
//  for(int i=0;i<nnodes;i++) str=str+" , "+std::to_string(_tasksToOffload[i]);
//  logInfo("dynamic distributor", str);

  for(int i=0; i<nnodes; i++) {
    _remainingTasksToOffload[i] = _tasksToOffload[i];
  }

//  delete[] newLoadDist;
}

bool exahype::reactive::DynamicDistributor::selectVictimRank(int& victim) {

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();
  victim = myRank;

  static std::atomic<int> rank_cnt (0);

  int l_rank = rank_cnt;

  for(int i=0; i<nnodes; i++) {
    if(l_rank!=myRank
       &&
	   _remainingTasksToOffload[l_rank].fetch_sub(1)>0) {
      victim = l_rank;
      l_rank = (l_rank + 1)%nnodes;
      break;
    }
    else
    _tasksToOffload[l_rank]=0;
    l_rank = (l_rank + 1)%nnodes;
  }
  rank_cnt=l_rank;

//  logInfo("performance monitor", "chose victim "<<victim);
  return victim!=myRank;
}

#endif
