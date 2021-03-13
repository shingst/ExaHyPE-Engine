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

#if  defined(Parallel)
#include "../reactive/OffloadingAnalyser.h"

#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"
#include "tarch/logging/CommandLineLogger.h"

#include "tarch/multicore/Lock.h"
#include "tarch/multicore/Core.h"
#include "tarch/multicore/Jobs.h"

#include <thread>
#include <limits>

#include "../reactive/PerformanceMonitor.h"
#include "../reactive/DiffusiveDistributor.h"
#include "../reactive/AggressiveDistributor.h"
#include "../reactive/AggressiveCCPDistributor.h"
#include "../reactive/AggressiveHybridDistributor.h"
#include "../reactive/OffloadingManager.h"

tarch::logging::Log  exahype::reactive::OffloadingAnalyser::_log( "exahype::offloading::OffloadingAnalyser" );

#ifdef USE_ITAC
#include "VT.h"
static int event_waitForWorker = -1;
static const char *event_name_waitForWorker = "waitForWorker";
#endif

exahype::reactive::OffloadingAnalyser::OffloadingAnalyser():
  _isSwitchedOn(true),
  _waitForWorkerDataWatch("exahype::offloading::OffloadingAnalyser", "-", false,false),
  _waitForMasterDataWatch("exahype::offloading::OffloadingAnalyser", "-", false,false),
  _waitForGlobalMasterDataWatch("exahype::offloading::OffloadingAnalyser", "-", false,false),
  _timeStepWatch("exahype::offloading::OffloadingAnalyser", "-", false,false),
  _waitForOtherRank(0),
  _currentZeroThreshold(0.002),
  _iterationCounter(0),
  _currentAccumulatedWorkerTime(0),
  _estimatedWtimeForPendingJobs(0),
  _lateSTPJobs(0)
{
  enable(true);
#ifdef USE_ITAC
  VT_funcdef(event_name_waitForWorker, VT_NOCLASS, &event_waitForWorker ); assertion(ierr==0);
#endif
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  _currentFilteredWaitingTimesSnapshot = new double[nnodes*nnodes];
  std::fill(&_currentFilteredWaitingTimesSnapshot[0], &_currentFilteredWaitingTimesSnapshot[nnodes*nnodes], 0);

  for(int i=0; i<nnodes; i++) 
     //_waitForOtherRank.push_back(tarch::timing::GlidingAverageMeasurement(0.1,16));
     _waitForOtherRank.push_back(tarch::timing::GlidingAverageMeasurement());
}


exahype::reactive::OffloadingAnalyser& exahype::reactive::OffloadingAnalyser::getInstance() {
  static OffloadingAnalyser* analyser = nullptr;

  if(analyser==nullptr) {
    analyser = new OffloadingAnalyser();
  }  
  return *analyser;
}

exahype::reactive::OffloadingAnalyser::~OffloadingAnalyser() {
    delete[] _currentFilteredWaitingTimesSnapshot;
}


void exahype::reactive::OffloadingAnalyser::enable(bool value) {
  _isSwitchedOn=value;
}

const double* exahype::reactive::OffloadingAnalyser::getFilteredWaitingTimesSnapshot() {
  return _currentFilteredWaitingTimesSnapshot;
}

double exahype::reactive::OffloadingAnalyser::getZeroThreshold() {
  return _currentZeroThreshold;
}

void exahype::reactive::OffloadingAnalyser::setTimePerSTP(double timePerSTP) {
  _timePerMigratablePredictionJob.setValue(timePerSTP);
  //logInfo("setTimePerSTP()", "submitted new STP measurement, current time per stp: "<<getTimePerSTP());
}

double exahype::reactive::OffloadingAnalyser::getTimePerSTP() {
  return _timePerMigratablePredictionJob.getValue();
}

void exahype::reactive::OffloadingAnalyser::setTimePerTimeStep(double timePerStep) {
  return _timePerTimeStep.setValue(timePerStep);
}

double exahype::reactive::OffloadingAnalyser::getTimePerTimeStep() {
  return _timePerTimeStep.getValue();
}

void exahype::reactive::OffloadingAnalyser::updateZeroTresholdAndFilteredSnapshot() {
  if(!_isSwitchedOn) return;

//#if !defined(AnalyseWaitingTimes)
  const double* currentWaitingTimesSnapshot = exahype::reactive::PerformanceMonitor::getInstance().getWaitingTimesSnapshot();
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  //check if valid: there is an entry>0 for every rank
  bool isValid = true;
  for(int i=0; i<nnodes; i++) {
    bool hasEntry = false;
    for(int j=0; j<nnodes; j++) {
      double currentWaitingTime = currentWaitingTimesSnapshot[i*nnodes+j];
      logDebug("updateZeroTresholdAndFilteredSnapshot()","rank "<<i<<" waiting for "<<currentWaitingTime<<" for rank "<<j);
      if(currentWaitingTime>0) {
        hasEntry = true;
        break;
      }
    }
    if(!hasEntry) {
       isValid=false;
       break;
    }
  }
  
  if(!isValid)  {
    logDebug("updateZeroThresholdAndFilteredSnapshot","not valid yet!");
    return;
  }

  double min, max;
  min = std::numeric_limits<double>::max();
  max = std::numeric_limits<double>::min();


  for(int i=0; i<nnodes; i++) {
    for(int j=0; j<nnodes; j++) {
      double currentWaitingTime = currentWaitingTimesSnapshot[i*nnodes+j];
      if(currentWaitingTime>0)
        logDebug("updateZeroTresholdAndFilteredSnapshot()","rank "<<i<<" waiting for "<<currentWaitingTime<<" for rank "<<j);
      if(currentWaitingTime>max) {
	 max = currentWaitingTime;
      }
      if(currentWaitingTime>0 && currentWaitingTime<min) {
        min = currentWaitingTime;
      }

      if(currentWaitingTime > _currentZeroThreshold) {
         _currentFilteredWaitingTimesSnapshot[i*nnodes+j] = currentWaitingTime;
      }
      else
        _currentFilteredWaitingTimesSnapshot[i*nnodes+j] = 0;
    }
  }
  if(min < std::numeric_limits<double>::max()) {
    double newThreshold = 0.95*min+0.05*max;
    _currentZeroThreshold = newThreshold;
    logDebug("updateZeroTresholdAndFilteredSnapshot()", " zero threshold set to "<< newThreshold);
  }
//#endif
}

void exahype::reactive::OffloadingAnalyser::beginIteration() {
  static int timestep_cnt = 0;

  if(_iterationCounter%2 !=0) return;

  if (_isSwitchedOn) {
    if(_timeStepWatch.isOn()) {
      _timeStepWatch.stopTimer();
      if(timestep_cnt <=5) {
        setTimePerTimeStep(_timeStepWatch.getCalendarTime());
      }
      timestep_cnt++;
    }
    _timeStepWatch.startTimer();

    exahype::reactive::OffloadingManager::getInstance().resetVictimFlag(); //TODO: correct position here?
    exahype::reactive::OffloadingManager::getInstance().recoverBlacklistedRanks();
  }
}


void exahype::reactive::OffloadingAnalyser::endIteration(double numberOfInnerLeafCells, double numberOfOuterLeafCells, double numberOfInnerCells, double numberOfOuterCells, double numberOfLocalCells, double numberOfLocalVertices) {
  if(!_isSwitchedOn) return;

  exahype::reactive::OffloadingManager::getInstance().printBlacklist();

  if(_iterationCounter%2 !=0) {
    _iterationCounter++;
    return;
  }

  _currentAccumulatedWorkerTime = 0;
 
  for(size_t i=0; i<_waitForOtherRank.size(); i++) {
    if(static_cast<int> (i) != tarch::parallel::Node::getInstance().getRank()) {
      logDebug("endIteration()", "wait for rank "<<i<<_waitForOtherRank[i].toString());
      if(_waitForOtherRank[i].isAccurateValue())
       exahype::reactive::PerformanceMonitor::getInstance().submitWaitingTimeForRank(_waitForOtherRank[i].getValue(), i);
    }     
  }

  updateZeroTresholdAndFilteredSnapshot();
  printWaitingTimes();

  switch(exahype::reactive::OffloadingManager::getInstance().getOffloadingStrategy()){
    case exahype::reactive::OffloadingManager::OffloadingStrategy::Diffusive:
      exahype::reactive::DiffusiveDistributor::getInstance().updateLoadDistribution(); break;
    case exahype::reactive::OffloadingManager::OffloadingStrategy::Aggressive:
      exahype::reactive::AggressiveDistributor::getInstance().printOffloadingStatistics();
      exahype::reactive::AggressiveDistributor::getInstance().updateLoadDistribution();
      break;
    case exahype::reactive::OffloadingManager::OffloadingStrategy::AggressiveCCP:
      exahype::reactive::AggressiveCCPDistributor::getInstance().printOffloadingStatistics();
      exahype::reactive::AggressiveCCPDistributor::getInstance().updateLoadDistribution();
      break;
    case exahype::reactive::OffloadingManager::OffloadingStrategy::AggressiveHybrid:
      exahype::reactive::AggressiveHybridDistributor::getInstance().printOffloadingStatistics();
      exahype::reactive::AggressiveHybridDistributor::getInstance().updateLoadDistribution();
      break;
  }
  exahype::reactive::OffloadingManager::getInstance().printBlacklist();
  //exahype::offloading::OffloadingManager::getInstance().resetPostedRequests();

  _iterationCounter++;
}

void exahype::reactive::OffloadingAnalyser::printWaitingTimes() {
  if(!_isSwitchedOn) return;
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  const double* waitingTimesSnapshot = getFilteredWaitingTimesSnapshot();
  int k = 0;
  for(int i=0; i<nnodes; i++) {
    for(int j=0; j<nnodes; j++) {
      if(waitingTimesSnapshot[k+j]>0)
        logInfo("printWaitingTimes()","rank "<<i<<" waiting for "<<waitingTimesSnapshot[k+j]<<" for rank "<<j);
    }
    k+= nnodes;
  }

}

void exahype::reactive::OffloadingAnalyser::beginToReceiveDataFromWorker() {
  if(!_isSwitchedOn) return;

  if (_isSwitchedOn) {
    _waitForWorkerDataWatch.startTimer();
    int pendingJobs  = tarch::multicore::jobs::getNumberOfWaitingBackgroundJobs();
    _lateSTPJobs = 0;
    _estimatedWtimeForPendingJobs = pendingJobs * getTimePerSTP() / tarch::multicore::Core::getInstance().getNumberOfThreads();
    logDebug("beginToReceiveDataFromWorker()","there are "<<pendingJobs<<" jobs left, estimated time "<<_estimatedWtimeForPendingJobs);

#ifdef USE_ITAC
    VT_begin(event_waitForWorker);
#endif
  }
}


void exahype::reactive::OffloadingAnalyser::endToReceiveDataFromWorker( int fromRank ) {

  if(_iterationCounter%2 !=0) {
     _iterationCounter++; 
     return;
  }

  if (_isSwitchedOn) {
 
#ifdef USE_ITAC
    VT_end(event_waitForWorker);
#endif
    _waitForWorkerDataWatch.stopTimer();
    double estimatedTimeForLateSTPs = _lateSTPJobs * getTimePerSTP()/ tarch::multicore::Core::getInstance().getNumberOfThreads();
    const double elapsedTime = std::max(0.000001,_waitForWorkerDataWatch.getCalendarTime()-_estimatedWtimeForPendingJobs
                                            -estimatedTimeForLateSTPs);

    if(elapsedTime>_currentZeroThreshold) {
      _currentAccumulatedWorkerTime += elapsedTime;
      logDebug("endToReceiveDataFromWorker", "elapsed "<<elapsedTime<<" accumulated "<<_currentAccumulatedWorkerTime);
      _waitForOtherRank[fromRank].setValue(_currentAccumulatedWorkerTime);
    }
    else {
      _waitForOtherRank[fromRank].setValue(elapsedTime);
    }
#if defined(Debug)
    double currentAvg = _waitForOtherRank[fromRank].getValue();
    
    if (tarch::la::greater(elapsedTime,0.0)) {
      logDebug(
        "endToReceiveDataFromWorker()",
        "rank had to wait for worker " << fromRank << " for "<< elapsedTime<<
        " currentAvg "<< currentAvg << "s"
      );
    }
#endif
//#ifdef USE_ITAC
//        VT_begin(event_waitForWorker);
//#endif
    _waitForWorkerDataWatch.startTimer();
  }
}

void exahype::reactive::OffloadingAnalyser::resetMeasurements() {
  if(!_isSwitchedOn) return;
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  for(int i=0; i<nnodes; i++)
     //_waitForOtherRank.push_back(tarch::timing::GlidingAverageMeasurement(0.1,16));
     _waitForOtherRank[i].erase();
}

void exahype::reactive::OffloadingAnalyser::beginToReceiveDataFromMaster(int master) {
  //if (_isSwitchedOn && !_waitForMasterDataWatch.isOn()) {
  //  _waitForMasterDataWatch.startTimer();
  //}
}


void exahype::reactive::OffloadingAnalyser::endToReceiveDataFromMaster(int master) {
  if(!_isSwitchedOn) return;
  if(master==0) 
     endToReceiveDataFromGlobalMaster();
  /*if (_isSwitchedOn && _waitForMasterDataWatch.isOn()) {
    _waitForMasterDataWatch.stopTimer();
    const double elapsedTime = _waitForMasterDataWatch.getCalendarTime();
    //int myMaster = tarch::parallel::NodePool::getInstance().getMasterRank();
    
    //_waitForOtherRank[myMaster].setValue(elapsedTime);

    //double currentAvg = _waitForOtherRank[myMaster].getValue();

    if (tarch::la::greater(elapsedTime,0.0)) {
      logInfo(
        "endToReceiveDataFromMaster()",
        "rank had to wait for master " << master << " for "<< elapsedTime 
      );
    }
  }*/
  
}

void exahype::reactive::OffloadingAnalyser::beginToSendDataToMaster() {
  //if (_isSwitchedOn && !_waitForMasterDataWatch.isOn()) {
  //  _waitForMasterDataWatch.startTimer();
  //}
}


void exahype::reactive::OffloadingAnalyser::endToSendDataToMaster() {
  /*if (_isSwitchedOn && _waitForMasterDataWatch.isOn()) {
    _waitForMasterDataWatch.stopTimer();
    const double elapsedTime = _waitForMasterDataWatch.getCalendarTime();
    int myMaster = tarch::parallel::NodePool::getInstance().getMasterRank();
    
    //_waitForOtherRank[myMaster].setValue(elapsedTime);

    //double currentAvg = _waitForOtherRank[myMaster].getValue();

    if (tarch::la::greater(elapsedTime,0.0)) {
      logInfo(
        "endToSendDataToMaster()",
        "rank had to wait for send to master " << myMaster << " for "<< elapsedTime 
      );
    }
    _waitForMasterDataWatch.startTimer();
  }*/
  if (_isSwitchedOn && !_waitForMasterDataWatch.isOn()) {
    int pendingJobs  = tarch::multicore::jobs::getNumberOfWaitingBackgroundJobs();
    

    _estimatedWtimeForPendingJobs = pendingJobs * getTimePerSTP() / tarch::multicore::Core::getInstance().getNumberOfThreads();
    _lateSTPJobs = 0;
    logDebug("endToSendDataToMaster()","there are "<<pendingJobs<<" left estimated time "<<_estimatedWtimeForPendingJobs); 
    _waitForMasterDataWatch.startTimer();
  }
}

void exahype::reactive::OffloadingAnalyser::beginToSendDataToWorker() {
  if (_isSwitchedOn && !_waitForWorkerDataWatch.isOn()) {
    _waitForWorkerDataWatch.startTimer();
  }
}


void exahype::reactive::OffloadingAnalyser::endToSendDataToWorker(int worker) {
  if (_isSwitchedOn && _waitForWorkerDataWatch.isOn()) {
    _waitForWorkerDataWatch.stopTimer();
    const double elapsedTime = _waitForWorkerDataWatch.getCalendarTime();

    //_waitForOtherRank[worker].setValue(elapsedTime);

    //double currentAvg = _waitForOtherRank[worker].getValue();

    if (tarch::la::greater(elapsedTime,0.0)) {
      logDebug(
        "endToSendDataToWorker()",
        "rank had to wait for send to worker " << worker << " for "<< elapsedTime 
      );
    }
  }
}

void exahype::reactive::OffloadingAnalyser::notifyReceivedSTPJob() {
  _lateSTPJobs++;
}

void exahype::reactive::OffloadingAnalyser::beginToReceiveDataFromGlobalMaster() {
  //if (_isSwitchedOn && !_waitForGlobalMasterDataWatch.isOn()) {
    //_waitForGlobalMasterDataWatch.startTimer();
  //}
}


void exahype::reactive::OffloadingAnalyser::endToReceiveDataFromGlobalMaster() {
  if (_isSwitchedOn && _waitForMasterDataWatch.isOn()) {
    _waitForMasterDataWatch.stopTimer();
    double estimatedTimeForLateSTPs = _lateSTPJobs * getTimePerSTP()/ tarch::multicore::Core::getInstance().getNumberOfThreads();
    logDebug("endToReceiveDataFromGlobalMaster()","estimate for late STPs "<<estimatedTimeForLateSTPs<<"s  estimate for pending jobs "<<_estimatedWtimeForPendingJobs); 
    const double elapsedTime = std::max(0.000001, _waitForMasterDataWatch.getCalendarTime()-_estimatedWtimeForPendingJobs
                                            -estimatedTimeForLateSTPs);
    _estimatedWtimeForPendingJobs = 0;
    _lateSTPJobs = 0;
    
    _waitForOtherRank[0].setValue(elapsedTime); //0 is global master

#if defined(Debug)
    double currentAvg = _waitForOtherRank[0].getValue();

    if (tarch::la::greater(elapsedTime,0.0)) {
      logDebug(
        "endToReceiveDataFromGlobalMaster()",
        "rank had to wait for global master "  << " for "<< elapsedTime <<
        " currentAvg "<< currentAvg << "s"
      );
    }
#endif
  }
}

void exahype::reactive::OffloadingAnalyser::dataWasNotReceivedInBackground( int fromRank, int tag, int cardinality, int pageSize ) {
}


void exahype::reactive::OffloadingAnalyser::beginToReleaseSynchronousHeapData() {
}


void exahype::reactive::OffloadingAnalyser::endToReleaseSynchronousHeapData() {
}


void exahype::reactive::OffloadingAnalyser::beginToPrepareAsynchronousHeapDataExchange() {
}


void exahype::reactive::OffloadingAnalyser::endToPrepareAsynchronousHeapDataExchange() {
}


void exahype::reactive::OffloadingAnalyser::beginReleaseOfJoinData() {
}


void exahype::reactive::OffloadingAnalyser::endReleaseOfJoinData() {
}


void exahype::reactive::OffloadingAnalyser::beginReleaseOfBoundaryData() {
}


void exahype::reactive::OffloadingAnalyser::endReleaseOfBoundaryData() {
}


void exahype::reactive::OffloadingAnalyser::changeConcurrencyLevel(int actualChange, int maxPossibleChange) {
}


void exahype::reactive::OffloadingAnalyser::minuteNumberOfBackgroundTasks(int taskCount) {
}


void exahype::reactive::OffloadingAnalyser::beginProcessingBackgroundJobs() {
}


void exahype::reactive::OffloadingAnalyser::endProcessingBackgroundJobs() {
}


void exahype::reactive::OffloadingAnalyser::enterCentralElementOfEnclosingSpacetree() {
}


void exahype::reactive::OffloadingAnalyser::leaveCentralElementOfEnclosingSpacetree() {
}


void exahype::reactive::OffloadingAnalyser::addWorker(
  int                                 workerRank,
  int                                 level
) {
}


void exahype::reactive::OffloadingAnalyser::removeWorker(
  int                                 workerRank,
  int                                 level
) {
}
#endif

