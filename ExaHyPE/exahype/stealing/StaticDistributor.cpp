#if  defined(SharedTBB)  && defined(Parallel)
#include "StaticDistributor.h"

#include <algorithm>
#include <numeric>

#include "tarch/multicore/Lock.h"
#include "tarch/parallel/Node.h"
#include "tarch/timing/Watch.h"

#include "exahype/stealing/StealingProfiler.h"
#include "exahype/stealing/PerformanceMonitor.h"
#include "tarch/multicore/Core.h"
#include "tarch/multicore/Jobs.h"

tarch::logging::Log exahype::stealing::StaticDistributor::_log( "exahype::stealing::StaticDistributor" );

exahype::stealing::StaticDistributor::StaticDistributor() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  _tasksToOffload          = new int[nnodes];
  _remainingTasksToOffload = new std::atomic<int>[nnodes];

  _consumersPerRank        = new int[nnodes];
  _consumersPerRank[0]     = tarch::multicore::jobs::Job::getMaxNumberOfRunningBackgroundThreads()-1;
  for(int i=1; i<nnodes;i++) {
    _consumersPerRank[i] = tarch::multicore::jobs::Job::getMaxNumberOfRunningBackgroundThreads();
    logInfo("StaticDistributor()","weight "<<_consumersPerRank[i]<<" for rank "<<i);
  }

  std::fill( &_remainingTasksToOffload[0], &_remainingTasksToOffload[nnodes], 0);
  std::fill( &_tasksToOffload[0], &_tasksToOffload[nnodes], 0);
}

exahype::stealing::StaticDistributor::~StaticDistributor() {
  delete[] _tasksToOffload;
  delete[] _remainingTasksToOffload;
  delete[] _consumersPerRank;
}

#ifdef StealingStrategyStaticHardcoded
void exahype::stealing::StaticDistributor::loadDistributionFromFile(const std::string& filename) {

  logInfo("loadDistributionFromFile()", "loading from file "<<filename);
  int myRank = tarch::parallel::Node::getInstance().getRank();
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  std::ifstream file(filename);
  std::string str = "";

  if(!file.is_open()) {
    logError("loadDistributionFromFile()", "could not open file, no tasks will be offloaded");
  }

  while(std::getline(file,str)) {
    std::string rightString = str;

    try { 
      std::string source_rank = rightString.substr(0, rightString.find(","));
      rightString             = rightString.substr( source_rank.size()+1 );
      std::string dest_rank   = rightString.substr(0, rightString.find(","));
      rightString             = rightString.substr( dest_rank.size()+1 );
      std::string num_tasks   = rightString;

      int i_source_rank = stoi(source_rank);
      int i_dest_rank   = stoi(dest_rank);
      int i_num_tasks   = stoi(num_tasks);
    
      //logInfo("loadDistributionFromFile()", "found rule "<<i_source_rank<< ","<<i_dest_rank<<","<<i_num_tasks);

      if(i_num_tasks<0 || i_source_rank<0 || i_source_rank>= nnodes || i_dest_rank<0 || i_dest_rank>= nnodes) {
        logError("loadDistributionFromFile()", "found flawed stealing rule "
        										<<i_source_rank<< ","
												<<i_dest_rank<<","
												<<i_num_tasks);
        logError("loadDistributionFromFile()", "ignoring.. ");
        continue;
      }

      if(i_source_rank==myRank) {
         _tasksToOffload[i_dest_rank] = i_num_tasks;
      }
    }
    catch(std::out_of_range& exception) {
      logError("loadDistributionFromFile    ()", 
               "failed to parse hardcoded stealing distribution file" << filename
			   << " with error in " << exception.what());
      logError("loadDistributionFromFile()",
               "flawed string: " <<str);
    }
  }
}
#endif

void exahype::stealing::StaticDistributor::computeNewLoadDistribution(int enclaveCells, int skeletonCells) {

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  int *loadPerRank = new int[nnodes];
  //int *newLoadDist = new int[nnodes];

  int totalCells   = enclaveCells + skeletonCells;
  MPI_Allgather(&totalCells, 1, MPI_INTEGER, loadPerRank, 1, MPI_INTEGER, MPI_COMM_WORLD);

  int input_r=0, input_l=0;
  int output_r=0, output_l=0;

  int total_l=0;
  total_l = std::accumulate(&loadPerRank[0], &loadPerRank[nnodes], total_l);

  int total_consumers = 0;
  total_consumers = std::accumulate(&_consumersPerRank[0], &_consumersPerRank[nnodes], total_consumers);

  int avg_l_per_consumer = 0;
  avg_l_per_consumer = total_l / total_consumers;

  input_l = loadPerRank[input_r];
  output_l = loadPerRank[output_r];

  while(output_r<nnodes) {
    int target_load_out = _consumersPerRank[output_r] * avg_l_per_consumer;
    int target_load_in = _consumersPerRank[input_r] * avg_l_per_consumer;

    while(output_l<target_load_out) {
      int diff_l = target_load_out-output_l;

      if(output_r==input_r) {
        input_r++;
        assert(input_r<=nnodes);
        input_l = loadPerRank[input_r];
        continue;
      }

      int moveable = input_l-target_load_in;
      if(moveable>0) {
        int inc_l = std::min( diff_l, moveable );
        output_l += inc_l;
        input_l -= inc_l;
        //logInfo("performance monitor", " moving "<<inc_l<<" from rank "<<input_r<<" to rank "<<output_r);
        //newLoadDist[output_r] = output_l;
        //newLoadDist[input_r]  = input_l;

        if(input_r==myRank) {
          _tasksToOffload[output_r] = inc_l;
          stealing::StealingProfiler::getInstance().notifyTargetOffloadedTask(inc_l, output_r);
          //_tasksToOffload[output_r]= std::min(inc_l,1);
          //stealing::StealingProfiler::getInstance().notifyTargetOffloadedTask(std::min(inc_l,1), output_r);
        }
      }

      if(input_l <=target_load_in ) {
        input_r++;
        if(input_r<nnodes) {
          input_l = loadPerRank[input_r];
          target_load_in = _consumersPerRank[input_r] * avg_l_per_consumer;
        }
      }
    }
    output_r++;
    if(output_r<nnodes)
    output_l = loadPerRank[output_r];
  }

//  std::string str="new (intended!) load distribution ";
//  for(int i=0;i<nnodes;i++) str=str+" , "+std::to_string(newLoadDist[i]);
//  logInfo("static distributor", str);
//  str="tasks to offload ";
//  for(int i=0;i<nnodes;i++) str=str+" , "+std::to_string(_tasksToOffload[i]);
//  logInfo("static distributor", str);

  resetRemainingTasksToOffload();

  //delete[] newLoadDist;
  delete[] loadPerRank;
}

exahype::stealing::StaticDistributor& exahype::stealing::StaticDistributor::getInstance() {
  static StaticDistributor staticDist;
  return staticDist;
}

void exahype::stealing::StaticDistributor::resetRemainingTasksToOffload() {

  logInfo("resetRemainingTasksToOffload()", "resetting remaining tasks to offload");

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  for(int i=0; i<nnodes; i++) {
    _remainingTasksToOffload[i] = _tasksToOffload[i];
    logInfo("resetRemainingTasksToOffload()", "to rank "<<i<<" ntasks "<<_tasksToOffload[i]);
    
  }
}

bool exahype::stealing::StaticDistributor::selectVictimRank(int& victim) {

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

  if(victim!=myRank)
   logInfo("selectVictimRank", "chose victim "<<victim<<" _remainingTasksToOffload "<<_remainingTasksToOffload[victim]);
  
  return victim != myRank;
}
#endif
