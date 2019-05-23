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
#include "exahype/solvers/ADERDGSolver.h"
#include "tarch/multicore/Core.h"
#include "tarch/multicore/tbb/Jobs.h"

tarch::logging::Log exahype::stealing::AggressiveHybridDistributor::_log( "exahype::stealing::AggressiveHybridDistributor" );

exahype::stealing::AggressiveHybridDistributor::AggressiveHybridDistributor() :
  _isEnabled(false),
  _temperatureCCP(1),
  _temperatureDiffusion(0.5),
  _totalTasksOffloaded(0),
  _oldTotalTasksOffloaded(0),
  _incrementCurrent(0),
  _incrementPrevious(0)
{

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  _initialLoadPerRank      = new int[nnodes];
  _idealTasksToOffload     = new int[nnodes];
  _tasksToOffload          = new int[nnodes];
  _optimalTasksPerRank     = new int[nnodes];
  _remainingTasksToOffload = new std::atomic<int>[nnodes];
  _emergenciesPerRank      = new int[nnodes];
  _notOffloaded            = new int[nnodes];
 
  std::fill( &_remainingTasksToOffload[0], &_remainingTasksToOffload[nnodes], 0);
  std::fill( &_tasksToOffload[0], &_tasksToOffload[nnodes], 0);
  std::fill( &_idealTasksToOffload[0], &_idealTasksToOffload[nnodes], 0);
  std::fill( &_initialLoadPerRank[0], &_initialLoadPerRank[nnodes], 0);
  std::fill( &_notOffloaded[0], &_notOffloaded[nnodes], 0);
  std::fill( &_emergenciesPerRank[0], &_emergenciesPerRank[nnodes], 0);
  std::fill( &_optimalTasksPerRank[0], &_optimalTasksPerRank[nnodes], 0);
}

exahype::stealing::AggressiveHybridDistributor::~AggressiveHybridDistributor() {
  delete[] _optimalTasksPerRank;
  delete[] _emergenciesPerRank;
  delete[] _notOffloaded;
  delete[] _initialLoadPerRank;
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
  
  if(nnodes<=1)
    return;

  int *newLoadDist = new int[nnodes];
  std::fill(&newLoadDist[0], &newLoadDist[nnodes], 0);

  int totalCells   = enclaveCells + skeletonCells;
  MPI_Allgather(&totalCells, 1, MPI_INTEGER, _initialLoadPerRank, 1, MPI_INTEGER, MPI_COMM_WORLD);

#if defined(StealingUseMaster)
  int input_r=0, input_l=0;
  int output_r=0, output_l=0;
#else
  int input_r=1, input_l=0;
  int output_r=1, output_l=0;
#endif

  int total_l=0;
  total_l = std::accumulate(&_initialLoadPerRank[0], &_initialLoadPerRank[nnodes], total_l);

#if defined(StealingUseMaster)
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
  logInfo("computeIdealLoadDistribution()", str);
  str="ideal tasks to offload ";
  for(int i=0;i<nnodes;i++) str=str+" , "+std::to_string(_idealTasksToOffload[i]);
  logInfo("computeIdealLoadDistribution()", str);

  delete[] newLoadDist;
}

int exahype::stealing::AggressiveHybridDistributor::determineCriticalRank() {

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();  

  const double* waitingTimesSnapshot = exahype::stealing::StealingAnalyser::getInstance().getFilteredWaitingTimesSnapshot();
 
  int *waitingRanks = new int[nnodes];
  bool *isWaitingForSomeone = new bool[nnodes];
  std::fill(waitingRanks, waitingRanks+nnodes, 0);

  double currentLongestWaitTimeCritical = -1;
  int currentCriticalRank = -1;
  int k = 0;

  for(int i=0; i<nnodes; i++) {
    bool waitingForSomeone = false;
    for(int j=0; j<nnodes; j++) {
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
  
  delete[] isWaitingForSomeone;
  delete[] waitingRanks;

  return currentCriticalRank;
}

void exahype::stealing::AggressiveHybridDistributor::determineOptimalVictim(int &optimalVictim, double &waitingTimeOptimalVictim) {
  
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();  

  const double* waitingTimesSnapshot = exahype::stealing::StealingAnalyser::getInstance().getFilteredWaitingTimesSnapshot();

  int *waitingRanks = new int[nnodes];
  bool *isWaitingForSomeone = new bool[nnodes];
  std::fill(waitingRanks, waitingRanks+nnodes, 0);

  double currentLongestWaitTimeVictim = -1;
  int currentOptimalVictim = -1;

  int k = 0;
  for(int i=0; i<nnodes; i++) {
    bool waitingForSomeone = false;
    for(int j=0; j<nnodes; j++) {
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
      //     && i!=0) { //exclude rank 0 as optimal victim
             currentOptimalVictim = i;
             currentLongestWaitTimeVictim = waitingTimesSnapshot[k+j];
         }
      }
    }
    k+= nnodes;
  }

#ifndef StealingUseMaster
  if(currentOptimalVictim!=0) {
#endif
  optimalVictim = currentOptimalVictim;
  waitingTimeOptimalVictim = currentLongestWaitTimeVictim;
#ifndef StealingUseMaster
  }
#endif

  delete[] isWaitingForSomeone;
  delete[] waitingRanks;
}


exahype::stealing::AggressiveHybridDistributor& exahype::stealing::AggressiveHybridDistributor::getInstance() {
  static AggressiveHybridDistributor aggressiveDist;
  return aggressiveDist;
}

void exahype::stealing::AggressiveHybridDistributor::configure(
                   double startTempCCP, double startTempDiffusion,
                   int CCPFrequency, int CCPStepsPerPhase,
                   bool adaptTemperature,
                   double thresholdTempAdaptation) {

  logDebug("configure", "start temp CCP: "<<startTempCCP<<
                       " start temp diffusion: "<<startTempDiffusion<<
                       " CCP frequency: "<<CCPFrequency<<
                       " CCP steps per phase: "<<CCPStepsPerPhase<<
                       " adapt temperature: "<<adaptTemperature<<
                       " adapt temperature threshold:"<<thresholdTempAdaptation );

  _temperatureCCP = startTempCCP;
  _temperatureDiffusion = startTempDiffusion;
  _adaptTemperature = adaptTemperature;
  _CCPFrequency = CCPFrequency;
  _CCPStepsPerPhase = CCPStepsPerPhase;
  _thresholdTempAdaptation = thresholdTempAdaptation;
}

void exahype::stealing::AggressiveHybridDistributor::printOffloadingStatistics() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  for(int i=0; i<nnodes; i++) {
    if(i==myRank)
      continue;
    if(_tasksToOffload[i]>0)
    logDebug("printOffloadingStatistics()", "target tasks to rank "<<i<<" ntasks "<<_tasksToOffload[i]<<" not offloaded "<<_notOffloaded[i]); 
    _notOffloaded[i]=0;
  }
  logDebug("printOffloadingStatistics()", "temperature value CCP "<<_temperatureCCP );
  logDebug("printOffloadingStatistics()", "temperature value diffusion "<<_temperatureDiffusion );
  logDebug("printOffloadingStatistics()", "time per STP  "<< exahype::stealing::StealingAnalyser::getInstance().getTimePerSTP());
}


void exahype::stealing::AggressiveHybridDistributor::resetRemainingTasksToOffload() {

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  for(int i=0; i<nnodes; i++) {
    if(i==myRank)
      continue;
    _remainingTasksToOffload[i] = _tasksToOffload[i];
  }
}

void exahype::stealing::AggressiveHybridDistributor::handleEmergencyOnRank(int rank) {
  _emergenciesPerRank[rank]++;
  logDebug("handleEmergencyOnRank()","emergencies for rank:"<<_emergenciesPerRank[rank]);
}

void exahype::stealing::AggressiveHybridDistributor::updateLoadDistribution() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  static int iterationCounter = 0;
  static int remainingCCPSteps = 0;

  if(!_isEnabled) {
    return;
  }
  logInfo("updateLoadDistribution()","total offloaded: "<<_totalTasksOffloaded<<" previous: "<<_oldTotalTasksOffloaded);
  logInfo("updateLoadDistribution()","increment current "<<_incrementCurrent<<" previous: "<<_incrementPrevious);


  _oldTotalTasksOffloaded = _totalTasksOffloaded;

  std::copy(&_tasksToOffload[0], &_tasksToOffload[nnodes], &_optimalTasksPerRank[0]);
  int *tasksOffloadedPreviously = new int[nnodes];
  std::copy(&_tasksToOffload[0], &_tasksToOffload[nnodes], &tasksOffloadedPreviously[0]);

  double *temperature;

  bool useCCP = false;
  if(_CCPFrequency==0) {
    useCCP = false;
  } 
  else if(remainingCCPSteps>0) {
    useCCP = true;
  }
  else if(iterationCounter%_CCPFrequency==0) {
    useCCP = true;
    remainingCCPSteps = _CCPStepsPerPhase;
  }

  if(useCCP) {
    temperature = &_temperatureCCP;
    updateLoadDistributionCCP();
    remainingCCPSteps--;
  }
  else {
    temperature = &_temperatureDiffusion;
    updateLoadDistributionDiffusive();
  }

  _totalTasksOffloaded = 0;
  for(int i=0; i<nnodes; i++) {
    if(_emergenciesPerRank[i]>0) {   
      _optimalTasksPerRank[i] = 0; 
      _tasksToOffload[i] = std::max( (1-*temperature)* _tasksToOffload[i], 0.0);
      _emergenciesPerRank[i] = 0;
    }
    _totalTasksOffloaded += _tasksToOffload[i];
  }
  //Todo:: handling if we offload to critical rank here? -> currently only in diffusion, CCP does not care about critical ranks anymore

  //compute new increment
  _incrementPrevious = _incrementCurrent;
  _incrementCurrent = 0;
  for(int i=0; i<nnodes; i++) {
    _incrementCurrent += std::abs( _optimalTasksPerRank[i]-tasksOffloadedPreviously[i] );
  }
  delete[] tasksOffloadedPreviously;

  resetRemainingTasksToOffload();

  iterationCounter++;
}

void exahype::stealing::AggressiveHybridDistributor::updateLoadDistributionCCP() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  if(_incrementPrevious>0 && _incrementCurrent>0 && _adaptTemperature) {
    if((_incrementCurrent*1.0f/_incrementPrevious)>_thresholdTempAdaptation)
      _temperatureCCP = std::min(1.1, _temperatureCCP*1.1);
    else 
      _temperatureCCP = std::max(0.1, _temperatureCCP*0.9);
  }
  
  for(int i=0; i<nnodes; i++) {
    if(_idealTasksToOffload[i]>0) {
      //we have a potential victim rank
//     if(!exahype::stealing::StealingManager::getInstance().isBlacklisted(i)) {
      logDebug("updateLoadDistributionCCP", "offloading to "<<i<<" tasks "<<_temperatureCCP*_idealTasksToOffload[i]);
      _optimalTasksPerRank[i] = _idealTasksToOffload[i];
      _tasksToOffload[i] = std::ceil(std::max((1.0-_temperatureCCP), 0.0)*_tasksToOffload[i] + _temperatureCCP*_idealTasksToOffload[i]);
#ifdef DistributedStealingDisable
       _tasksToOffload[i] = 0;
#endif
//       }
     }
   }

}

void exahype::stealing::AggressiveHybridDistributor::updateLoadDistributionDiffusive() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  if(_incrementPrevious>0 && _incrementCurrent>0 && _adaptTemperature) {
    if((_incrementCurrent*1.0f/_incrementPrevious)>_thresholdTempAdaptation)
      _temperatureDiffusion = std::min(1.1, _temperatureDiffusion*1.1);
    else 
      _temperatureDiffusion = std::max(0.1, _temperatureDiffusion*0.9);
  }

  int currentCriticalRank = determineCriticalRank();

  int currentOptimalVictim=-1;
  double currentLongestWaitTimeVictim=-1;
  determineOptimalVictim(currentOptimalVictim, currentLongestWaitTimeVictim);

  logInfo("updateLoadDistributionDiffusive()", "optimal victim: "<<currentOptimalVictim<<" critical rank:"<<currentCriticalRank);

#ifndef DistributedStealingDisable

  bool isVictim = exahype::stealing::StealingManager::getInstance().isVictim();
  if(currentOptimalVictim>=0 && myRank == currentCriticalRank && currentCriticalRank!=currentOptimalVictim) {
    if(!isVictim) {
      int currentTasksCritical = _initialLoadPerRank[currentCriticalRank];
      for(int i=0; i<nnodes; i++) {
        currentTasksCritical -= _tasksToOffload[i];
      }
      int currentTasksOptimal = _initialLoadPerRank[currentOptimalVictim]+_tasksToOffload[currentOptimalVictim];

      //int optimalTasksToOffload = currentLongestWaitTimeVictim/(2*exahype::stealing::StealingAnalyser::getInstance().getTimePerSTP());
      int optimalTasksToOffload = (tarch::multicore::Core::getInstance().getNumberOfThreads()-1) *
                                  currentLongestWaitTimeVictim/(2*exahype::stealing::StealingAnalyser::getInstance().getTimePerSTP());
      logDebug("updateLoadDistributionDiffusive()", "optimal tasks to offload "<<optimalTasksToOffload);

      _optimalTasksPerRank[currentOptimalVictim] = optimalTasksToOffload;
      _tasksToOffload[currentOptimalVictim] = std::max((1-_temperatureDiffusion), 0.0)*_tasksToOffload[currentOptimalVictim] 
                                             + _temperatureDiffusion*optimalTasksToOffload;
     }
  }
  else if(_tasksToOffload[currentCriticalRank]>0) {
    _optimalTasksPerRank[currentCriticalRank] = 0;
    _tasksToOffload[currentCriticalRank] = std::max( (1-_temperatureDiffusion), 0.0)*_tasksToOffload[currentCriticalRank];
  }
#endif
  resetRemainingTasksToOffload();

}

void exahype::stealing::AggressiveHybridDistributor::getAllVictimRanks(std::vector<int> &victims ) {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  for(int i=0; i<nnodes; i++) {
    if(_tasksToOffload[i]>0 && _notOffloaded[i]!=_tasksToOffload[i] ) victims.push_back(i);
  }
}

bool exahype::stealing::AggressiveHybridDistributor::selectVictimRank(int& victim, bool& last) {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();
  victim = myRank;

  static std::atomic<int> rank_cnt = 0;

  int l_rank = rank_cnt;

  for(int i=0; i<nnodes; i++) {
    if(l_rank!=myRank) {
      int lastVal = _remainingTasksToOffload[l_rank].fetch_sub(1);
      if(lastVal>0) {
        victim = l_rank;
        l_rank = (l_rank + 1)%nnodes;
        break;
      }
      else
        _remainingTasksToOffload[l_rank]=0;
    }
    //else
    //  _remainingTasksToOffload[l_rank]=0;
    l_rank = (l_rank + 1)%nnodes;
  }
  rank_cnt=l_rank;

#ifdef StealingUseProgressTask
  if(victim == myRank) 
    last = true;
    //exahype::stealing::StealingManager::getInstance().notifyAllVictimsSendCompletedIfNotNotified();
#endif

  int threshold = 1+std::max(1, tarch::multicore::Core::getInstance().getNumberOfThreads()-1)*tarch::multicore::jobs::internal::_minimalNumberOfJobsPerConsumerRun;
  logDebug("selectVictimRank", "threshold "<<threshold<< " there are "<<exahype::solvers::ADERDGSolver::NumberOfEnclaveJobs-exahype::solvers::ADERDGSolver::NumberOfRemoteJobs<<" jobs "<< " and "<<tarch::multicore::jobs::getNumberOfWaitingBackgroundJobs());
  threshold = std::max(threshold, 20);
  //threshold = 0; //TODO:test

  //logInfo("selectVictimRank","waiting "<<tarch::multicore::jobs::getNumberOfWaitingBackgroundJobs()<<" criterion "<<threshold);
 
  //if(tarch::multicore::jobs::getNumberOfWaitingBackgroundJobs()<
  //      threshold) {
  if(exahype::solvers::ADERDGSolver::NumberOfEnclaveJobs-exahype::solvers::ADERDGSolver::NumberOfRemoteJobs<
        threshold) {
    //logInfo("selectVictimRank", "number of running consumers: "<<tarch::multicore::jobs::internal::_numberOfRunningJobConsumerTasks.load()<<" max running "<<tarch::multicore::Core::getInstance().getNumberOfThreads()-1);
    _notOffloaded[victim]++;
    victim = myRank;
  }
  //if(victim!=myRank)
   logInfo("selectVictimRank", "chose victim "<<victim<<" _remainingTasksToOffload "<<_remainingTasksToOffload[victim]);
  
  return victim != myRank;
}
#endif
