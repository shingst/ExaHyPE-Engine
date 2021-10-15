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

#if defined(Parallel)
#include "exahype/reactive/AggressiveHybridDistributor.h"

#include <algorithm>
#include <numeric>
#include <cmath>

#include "tarch/multicore/Lock.h"
#include "tarch/parallel/Node.h"
#include "tarch/timing/Watch.h"

#include "exahype/reactive/OffloadingProfiler.h"
#include "exahype/reactive/PerformanceMonitor.h"
#include "exahype/reactive/OffloadingAnalyser.h"
#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/LimitingADERDGSolver.h"
#include "tarch/multicore/Core.h"
#include "tarch/multicore/tbb/Jobs.h"

tarch::logging::Log exahype::reactive::AggressiveHybridDistributor::_log( "exahype::reactive::AggressiveHybridDistributor" );

#define MINIMUM_TASKS_IN_QUEUE 20 //a magic constant which controls that there is always a minimum of 20 tasks in the queue

exahype::reactive::AggressiveHybridDistributor::AggressiveHybridDistributor() :
  _temperatureDiffusion(0.5),
  _temperatureCCP(1),
  _thresholdTempAdaptation(1),
  _adaptTemperature(false),
  _CCPFrequency(0),
  _CCPStepsPerPhase(0),
  _totalTasksOffloaded(0),
  _oldTotalTasksOffloaded(0),
  _incrementCurrent(0),
  _incrementPrevious(0),
  _isEnabled(false),
  _currentOptimalVictim(-1),
  _currentCriticalRank(-1),
  _localStarvationThreshold(-1)
{

  int nnodes                = tarch::parallel::Node::getInstance().getNumberOfNodes();

  _initialLoad              = new int[nnodes];
  _idealTasksToOffloadCCP   = new int[nnodes];
  _tasksToOffload           = new int[nnodes];
  _optimalTasks             = new int[nnodes];
  _remainingTasksToOffload  = new std::atomic<int>[nnodes];
  _emergencies              = new int[nnodes];
  _tasksNotOffloaded        = new std::atomic<int>[nnodes];
  _tasksActuallyOffloaded   = new std::atomic<int>[nnodes];
 
  std::fill( &_remainingTasksToOffload[0], &_remainingTasksToOffload[nnodes], 0);
  std::fill( &_tasksToOffload[0], &_tasksToOffload[nnodes], 0);
  std::fill( &_idealTasksToOffloadCCP[0], &_idealTasksToOffloadCCP[nnodes], 0);
  std::fill( &_initialLoad[0], &_initialLoad[nnodes], 0);
  std::fill( &_tasksNotOffloaded[0], &_tasksNotOffloaded[nnodes], 0);
  std::fill( &_emergencies[0], &_emergencies[nnodes], 0);
  std::fill( &_optimalTasks[0], &_optimalTasks[nnodes], 0);
  std::fill( &_tasksActuallyOffloaded[0], &_tasksActuallyOffloaded[nnodes], 0);
}

exahype::reactive::AggressiveHybridDistributor::~AggressiveHybridDistributor() {
  delete[] _tasksActuallyOffloaded;
  delete[] _optimalTasks;
  delete[] _emergencies;
  delete[] _tasksNotOffloaded;
  delete[] _initialLoad;
  delete[] _tasksToOffload;
  delete[] _remainingTasksToOffload;
  delete[] _idealTasksToOffloadCCP;
}

void exahype::reactive::AggressiveHybridDistributor::enable() {
  _isEnabled = true;
}

void exahype::reactive::AggressiveHybridDistributor::disable() {
  _isEnabled = false;
}

void exahype::reactive::AggressiveHybridDistributor::computeIdealLoadDistribution(int enclaveCells, int skeletonCells) {

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();
  
  // nothing to do if there is a single rank only
  if(nnodes<=1)
    return;

  int *newLoadDist = new int[nnodes];
  std::fill(&newLoadDist[0], &newLoadDist[nnodes], 0);

  int totalCells   = enclaveCells + skeletonCells;
  MPI_Allgather(&totalCells, 1, MPI_INT, _initialLoad, 1, MPI_INT, MPI_COMM_WORLD);

#if defined(OffloadingUseMaster)
  int input_r=0, input_l=0;
  int output_r=0, output_l=0;
#else
  int input_r=1, input_l=0;
  int output_r=1, output_l=0;
#endif

  int totalLoad=0;
  totalLoad = std::accumulate(&_initialLoad[0], &_initialLoad[nnodes], totalLoad);

#if defined(OffloadingUseMaster)
  int avg_l = totalLoad / nnodes;
#else
  int avg_l = totalLoad / (nnodes-1);
  newLoadDist[0] = 0;
  _idealTasksToOffloadCCP[0] = 0;
#endif

  input_l = _initialLoad[input_r];
  output_l = _initialLoad[output_r];

  while(output_r<nnodes) {
    int target_load_out = avg_l;
    int target_load_in = avg_l;

    while(output_l<target_load_out) {
      int diff_l = target_load_out-output_l;

      if(output_r==input_r) {
        input_r++;
        assert(input_r<=nnodes);
        input_l = _initialLoad[input_r];
        continue;
      }

      int moveable = input_l-target_load_in;
      if(moveable>0) {
        int inc_l = std::min( diff_l, moveable );
        output_l += inc_l;
        input_l -= inc_l;
        logDebug("computeIdealLoadDistribution()", " moving "<<inc_l<<" tasks from rank "<<input_r<<" to rank "<<output_r);
        newLoadDist[output_r] = output_l;
        newLoadDist[input_r]  = input_l;

        if(input_r==myRank) {
          logInfo("computeIdealLoadDistribution()"," inc_l="<<inc_l);
          _idealTasksToOffloadCCP[output_r] = inc_l;
        }
      }

      if(input_l <=target_load_in ) {
        input_r++;
        if(input_r<nnodes) {
          input_l = _initialLoad[input_r];
          target_load_in = avg_l;
        }
      }
    }
    output_r++;
    if(output_r<nnodes)
    output_l = _initialLoad[output_r];
  }

  std::string str="ideal load distribution ";
  for(int i=0;i<nnodes;i++) str=str+" , "+std::to_string(newLoadDist[i]);
    logDebug("computeIdealLoadDistribution()", str);
  str="ideal tasks to offload ";
  for(int i=0;i<nnodes;i++) str=str+" , "+std::to_string(_idealTasksToOffloadCCP[i]);
    logDebug("computeIdealLoadDistribution()", str);

  delete[] newLoadDist;
}

int exahype::reactive::AggressiveHybridDistributor::determineCriticalRank() {

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();  

  const double* waitingTimesSnapshot =
		 exahype::reactive::OffloadingAnalyser::getInstance().getFilteredWaitingTimesSnapshot();
 
  int *waitingRanks = new int[nnodes];
  bool *isWaitingForRank = new bool[nnodes];
  std::fill(waitingRanks, waitingRanks+nnodes, 0);

  double currentLongestWaitTimeCritical = -1;
  int currentCriticalRank = -1;
  int k = 0;

  for(int i=0; i<nnodes; i++) {
    bool isWaitingForSomeone = false;
    for(int j=0; j<nnodes; j++) {
      if(waitingTimesSnapshot[k+j]> exahype::reactive::OffloadingAnalyser::getInstance().getZeroThreshold()) {
        waitingRanks[j]++;   
        isWaitingForSomeone = true;
      }
    }
    isWaitingForRank[i]= isWaitingForSomeone;
    k+= nnodes;
  }

  k = 0;
  for(int i=0; i<nnodes; i++) {
    for(int j=0; j<nnodes; j++) {
      if(waitingTimesSnapshot[k+j]> exahype::reactive::OffloadingAnalyser::getInstance().getZeroThreshold()) {
         if(waitingTimesSnapshot[k+j]>currentLongestWaitTimeCritical
           && !isWaitingForRank[j]
           && waitingRanks[j]>0) {
             currentCriticalRank = j;
             currentLongestWaitTimeCritical = waitingTimesSnapshot[k+j];
         } 
      }
    }
    k+= nnodes;
  }
  
  delete[] isWaitingForRank;
  delete[] waitingRanks;

  return currentCriticalRank;
}

void exahype::reactive::AggressiveHybridDistributor::determineOptimalVictim(
    int &optimalVictim,
    double &waitingTimeOptimalVictim) {
  
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();  

  const double* waitingTimesSnapshot = exahype::reactive::OffloadingAnalyser::getInstance().getFilteredWaitingTimesSnapshot();

  int *waitingRanks = new int[nnodes];
  bool *isWaitingForSomeone = new bool[nnodes];
  std::fill(waitingRanks, waitingRanks+nnodes, 0);

  double currentLongestWaitTimeVictim = -1;
  int currentOptimalVictim = -1;

  int k = 0;
  for(int i=0; i<nnodes; i++) {
    bool waitingForSomeone = false;
    for(int j=0; j<nnodes; j++) {
      if(waitingTimesSnapshot[k+j]> exahype::reactive::OffloadingAnalyser::getInstance().getZeroThreshold()) {
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
      if(waitingTimesSnapshot[k+j]> exahype::reactive::OffloadingAnalyser::getInstance().getZeroThreshold()) {
         if(waitingTimesSnapshot[k+j]>currentLongestWaitTimeVictim 
           && !exahype::reactive::ReactiveContext::getInstance().isBlacklisted(i)
           && waitingRanks[i]==0) {
      //     && i!=0) { //exclude rank 0 as optimal victim
             currentOptimalVictim = i;
             currentLongestWaitTimeVictim = waitingTimesSnapshot[k+j];
         }
      }
    }
    k+= nnodes;
  }

#ifndef OffloadingUseMaster
  if(currentOptimalVictim!=0) {
#endif
    optimalVictim = currentOptimalVictim;
    waitingTimeOptimalVictim = currentLongestWaitTimeVictim;
#ifndef OffloadingUseMaster
  }
#endif

  delete[] isWaitingForSomeone;
  delete[] waitingRanks;
}


exahype::reactive::AggressiveHybridDistributor& exahype::reactive::AggressiveHybridDistributor::getInstance() {
  static AggressiveHybridDistributor aggressiveDist;
  return aggressiveDist;
}

void exahype::reactive::AggressiveHybridDistributor::configure(
    double startTempCCP,
    double startTempDiffusion,
    int CCPFrequency,
    int CCPStepsPerPhase,
    bool adaptTemperature,
    double thresholdTempAdaptation,
    int localStarvationThreshold) {

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
  _localStarvationThreshold = localStarvationThreshold;
}

void exahype::reactive::AggressiveHybridDistributor::printOffloadingStatistics() {

  if(!_isEnabled) return;

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  for(int i=0; i<nnodes; i++) {
    if(i==myRank)
      continue;
    if(_tasksToOffload[i]>0)
      logDebug("printOffloadingStatistics()", "target tasks to rank "<<i<<" ntasks "<<_tasksToOffload[i]<<" not offloaded "
                                            <<_tasksNotOffloaded[i]<<" actually offloaded "<<_tasksActuallyOffloaded[i]);
    _tasksNotOffloaded[i] = 0;
    _tasksActuallyOffloaded[i] = 0;
  }
  logDebug("printOffloadingStatistics()", "temperature value CCP "<<_temperatureCCP );
  logDebug("printOffloadingStatistics()", "temperature value diffusion "<<_temperatureDiffusion );
  logDebug("printOffloadingStatistics()", "time per STP  "<< exahype::reactive::OffloadingAnalyser::getInstance().getTimePerSTP());
}


void exahype::reactive::AggressiveHybridDistributor::resetRemainingTasksToOffload() {

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();

  for(int i=0; i<nnodes; i++) {
    if(i==myRank)
      continue;
    _remainingTasksToOffload[i] = _tasksToOffload[i];
    reactive::OffloadingProfiler::getInstance().notifyTargetOffloadedTask(_tasksToOffload[i], i);
  }
}

void exahype::reactive::AggressiveHybridDistributor::resetTasksToOffload() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  for(int i=0; i<nnodes; i++) {
    _tasksToOffload[i] = 0;
  }
}

void exahype::reactive::AggressiveHybridDistributor::handleEmergencyOnRank(int rank) {
  _emergencies[rank]++;
  logDebug("handleEmergencyOnRank()","emergencies for rank:"<<_emergencies[rank]);
}

void exahype::reactive::AggressiveHybridDistributor::updateLoadDistribution() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  static int iterationCounter = 0;
  static int remainingCCPSteps = 0;

  if(!_isEnabled) {
    return;
  }
  logDebug("updateLoadDistribution()","total offloaded (target): "<<_totalTasksOffloaded<<" previous: "<<_oldTotalTasksOffloaded);
  logDebug("updateLoadDistribution()","increment current "<<_incrementCurrent<<" previous: "<<_incrementPrevious);

  _oldTotalTasksOffloaded = _totalTasksOffloaded;

  std::copy(&_tasksToOffload[0], &_tasksToOffload[nnodes], &_optimalTasks[0]);
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
    if(_emergencies[i]>0) {   
      _optimalTasks[i] = 0; 
      _tasksToOffload[i] = std::max( (1-*temperature)* _tasksToOffload[i], 0.0);
      _emergencies[i] = 0;
    }
    _totalTasksOffloaded += _tasksToOffload[i];
  }

  //compute new increment
  _incrementPrevious = _incrementCurrent;
  _incrementCurrent = 0;
  for(int i=0; i<nnodes; i++) {
    _incrementCurrent += std::abs( _optimalTasks[i]-tasksOffloadedPreviously[i] );
  }
  delete[] tasksOffloadedPreviously;

  resetRemainingTasksToOffload();

  iterationCounter++;
}

void exahype::reactive::AggressiveHybridDistributor::updateLoadDistributionCCP() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  if(_incrementPrevious>0 && _incrementCurrent>0 && _adaptTemperature) {
    if((_incrementCurrent*1.0f/_incrementPrevious)>_thresholdTempAdaptation)
      _temperatureCCP = std::min(1.1, _temperatureCCP*1.1);
    else 
      _temperatureCCP = std::max(0.1, _temperatureCCP*0.9);
  }
  
  for(int i=0; i<nnodes; i++) {
    if(_idealTasksToOffloadCCP[i]>0) {
       _optimalTasks[i] = _idealTasksToOffloadCCP[i];
      _tasksToOffload[i] = std::ceil(std::max((1.0-_temperatureCCP), 0.0)*_tasksToOffload[i] + _temperatureCCP*_idealTasksToOffloadCCP[i]);
#ifdef DistributedOffloadingDisable
       _tasksToOffload[i] = 0;
#endif
     }
   }

}

void exahype::reactive::AggressiveHybridDistributor::updateLoadDistributionDiffusive() {
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

  logDebug("updateLoadDistributionDiffusive()", "optimal victim: "<<currentOptimalVictim<<" critical rank:"<<currentCriticalRank);

#ifndef DistributedOffloadingDisable
  bool isVictim = exahype::reactive::ReactiveContext::getInstance().isVictim();

  if(currentOptimalVictim>=0
  && myRank == currentCriticalRank
  && currentCriticalRank!=currentOptimalVictim) {
    if(!isVictim) {
      int currentTasksCritical = _initialLoad[currentCriticalRank];
      for(int i=0; i<nnodes; i++) {
        currentTasksCritical -= _tasksToOffload[i];
      }

      int optimalTasksToOffload = (tarch::multicore::Core::getInstance().getNumberOfThreads()-1) *
                                  currentLongestWaitTimeVictim/
								  (2*exahype::reactive::OffloadingAnalyser::getInstance().getTimePerSTP());
      logDebug("updateLoadDistributionDiffusive()", "optimal tasks to offload "<<optimalTasksToOffload);

      _optimalTasks[currentOptimalVictim] = optimalTasksToOffload;
      _tasksToOffload[currentOptimalVictim] = std::max((1-_temperatureDiffusion), 0.0)*_tasksToOffload[currentOptimalVictim] 
                                             + _temperatureDiffusion*optimalTasksToOffload;
     }
  }
  else if(currentCriticalRank>=0 && _tasksToOffload[currentCriticalRank]>0) {
    _optimalTasks[currentCriticalRank] = 0;
    _tasksToOffload[currentCriticalRank] = std::max( (1-_temperatureDiffusion), 0.0)*_tasksToOffload[currentCriticalRank];
  }
#endif

  //store decision
  _currentOptimalVictim = currentOptimalVictim;
  _currentCriticalRank = currentCriticalRank;

  resetRemainingTasksToOffload();

}

void exahype::reactive::AggressiveHybridDistributor::getAllVictimRanks(std::vector<int> &victims ) const {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  for(int i=0; i<nnodes; i++) {
    if(_tasksToOffload[i]>0 && _tasksNotOffloaded[i]!=_tasksToOffload[i] ) victims.push_back(i);
  }
}

int exahype::reactive::AggressiveHybridDistributor::getCurrentOptimalVictim() const {
  return _currentOptimalVictim;
}

int exahype::reactive::AggressiveHybridDistributor::getCurrentCriticalRank() const {
  return _currentCriticalRank;
}

bool exahype::reactive::AggressiveHybridDistributor::selectVictimRank(int& victim, bool& isLastVictim) {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int myRank = tarch::parallel::Node::getInstance().getRank();
  victim = myRank;

  static std::atomic<int> nextRoundRobinRank (0);

  int tmpRank = nextRoundRobinRank;

  for(int i=0; i<nnodes; i++) {
    if(tmpRank!=myRank) {
      int lastVal = _remainingTasksToOffload[tmpRank].fetch_sub(1);
      if(lastVal>0) {
        victim = tmpRank;
        tmpRank = (tmpRank + 1)%nnodes;
        break;
      }
      else
        _remainingTasksToOffload[tmpRank]=0;
    }
    tmpRank = (tmpRank + 1)%nnodes;
  }
  nextRoundRobinRank = tmpRank;

#if defined(SharedTBB)
  int minimumTasksInQueue = (_localStarvationThreshold==-1) ? 1
		         + std::max(1, tarch::multicore::Core::getInstance().getNumberOfThreads()-1)
                 * tarch::multicore::jobs::internal::_minimalNumberOfJobsPerConsumerRun
		 : _localStarvationThreshold;
#else
  int minimumTasksInQueue = MINIMUM_TASKS_IN_QUEUE;
#endif 
            
  minimumTasksInQueue = std::max(minimumTasksInQueue, MINIMUM_TASKS_IN_QUEUE);

  if(exahype::solvers::ADERDGSolver::NumberOfEnclaveJobs-exahype::solvers::ADERDGSolver::NumberOfRemoteJobs //current number of jobs
      < minimumTasksInQueue) {
	  logDebug("selectVictimRank", "threshold "<<minimumTasksInQueue
	  		                       << " there are "<<exahype::solvers::ADERDGSolver::NumberOfEnclaveJobs-exahype::solvers::ADERDGSolver::NumberOfRemoteJobs
	                             <<" jobs "<< " and "<<tarch::multicore::jobs::getNumberOfWaitingBackgroundJobs()<<" waiting background jobs");

    reactive::OffloadingProfiler::getInstance().notifyThresholdFail();
    _tasksNotOffloaded[victim]++;
    victim = myRank;
  }
  else if(victim!=myRank) {
    _tasksActuallyOffloaded[victim]++;
    //task is tracked for profiling in ADERDGSolver
  }

  logDebug("selectVictimRank", "chose victim "<<victim<<" _remainingTasksToOffload "<<_remainingTasksToOffload[victim]);
  
  return victim != myRank;
}
#endif
