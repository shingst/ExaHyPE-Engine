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
#include "exahype/reactive/OffloadingAnalyser.h"

#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"
#include "tarch/logging/CommandLineLogger.h"

#include "tarch/multicore/Lock.h"
#include "tarch/multicore/Core.h"
#include "tarch/multicore/Jobs.h"

#include <thread>
#include <limits>

#include "exahype/reactive/OffloadingProfiler.h"

#include "exahype/reactive/PerformanceMonitor.h"
#include "exahype/reactive/AggressiveHybridDistributor.h"
#include "exahype/reactive/ReactiveContext.h"
#include "exahype/reactive/Blacklist.h"

tarch::logging::Log  exahype::reactive::OffloadingAnalyser::_log( "exahype::reactive::OffloadingAnalyser" );

#ifdef USE_ITAC
#include "VT.h"
static int event_waitForWorker = -1;
static const char *event_name_waitForWorker = "waitForWorker";
#endif

exahype::reactive::OffloadingAnalyser::OffloadingAnalyser():
  _isSwitchedOn(true),
  _waitForWorkerDataWatch("exahype::reactive::OffloadingAnalyser", "-", false,false),
  _waitForMasterDataWatch("exahype::reactive::OffloadingAnalyser", "-", false,false),
  _waitForGlobalMasterDataWatch("exahype::reactive::OffloadingAnalyser", "-", false,false),
  _timeStepWatch("exahype::reactive::OffloadingAnalyser", "-", false,false),
  _waitAvgForOtherRank(0),
  _currentZeroThreshold(0.002),
  _iterationCounter(0),
  _currentAccumulatedWorkerTime(0),
  _estWalltimeForPendingJobsAtReceiveForWorker(0),
  _numLateIncomingJobsWhileWaiting(0)
{
  enable(true);
#ifdef USE_ITAC
  VT_funcdef(event_name_waitForWorker, VT_NOCLASS, &event_waitForWorker ); assertion(ierr==0);
#endif
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  _currentFilteredWaitingTimesSnapshot.resize(nnodes*nnodes);
  std::fill(_currentFilteredWaitingTimesSnapshot.begin(), _currentFilteredWaitingTimesSnapshot.end(), 0);

  for(int i=0; i<nnodes; i++) 
     _waitAvgForOtherRank.push_back(tarch::timing::GlidingAverageMeasurement());
}

exahype::reactive::OffloadingAnalyser& exahype::reactive::OffloadingAnalyser::getInstance() {
  static OffloadingAnalyser *analyser = nullptr; //must be a raw pointer since Peano takes ownership and frees up the analyser eventually!
  //fixme: this can cause a (non-critical) resource leak at the end of the application if the analyser is not used by Peano
  //Peano should not necessarily have to take owernship such that we could change analyser to a normal variable here...
  if(analyser==nullptr) {
    analyser = new OffloadingAnalyser();
  }
  return *analyser;
}

void exahype::reactive::OffloadingAnalyser::enable(bool value) {
  _isSwitchedOn=value;
}

const std::vector<double>& exahype::reactive::OffloadingAnalyser::getFilteredWaitingTimesSnapshot() const {
  return _currentFilteredWaitingTimesSnapshot;
}

double exahype::reactive::OffloadingAnalyser::getZeroThreshold() const {
  return _currentZeroThreshold;
}

void exahype::reactive::OffloadingAnalyser::setTimePerSTP(double timePerSTP) {
  _avgTimePerSTP.setValue(timePerSTP);
}

double exahype::reactive::OffloadingAnalyser::getTimePerSTP() const {
  return _avgTimePerSTP.getValue();
}

void exahype::reactive::OffloadingAnalyser::setTimePerTimeStep(double timePerStep) {
  return _avgTimePerTimeStep.setValue(timePerStep);
}

double exahype::reactive::OffloadingAnalyser::getTimePerTimeStep() const {
  return _avgTimePerTimeStep.getValue();
}

bool exahype::reactive::OffloadingAnalyser::isValidSnapshot() const {
  const std::vector<double>& currentWaitingTimesSnapshot = exahype::reactive::PerformanceMonitor::getInstance().getWaitingTimesGlobalSnapshot();
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
  return isValid;
}

void exahype::reactive::OffloadingAnalyser::updateZeroTresholdAndFilteredSnapshot() {
  if(!_isSwitchedOn) return;

  const std::vector<double>& currentWaitingTimesSnapshot = exahype::reactive::PerformanceMonitor::getInstance().getWaitingTimesGlobalSnapshot();
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  //check if valid: there is an entry>0 for every rank
  bool isValid = isValidSnapshot();
  
  if(!isValid)  {
    logDebug("updateZeroThresholdAndFilteredSnapshot","not valid yet!");
    return;
  }

  double min, max;
  min = std::numeric_limits<double>::max();
  max = std::numeric_limits<double>::min();

  //filter out zero waiting times
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

  //update zero threshold
  if(min < std::numeric_limits<double>::max()) {
    double newThreshold = 0.95*min+0.05*max;
    _currentZeroThreshold = newThreshold;
    logDebug("updateZeroTresholdAndFilteredSnapshot()", " zero threshold set to "<< newThreshold);
  }
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

    exahype::reactive::ReactiveContext::getInstance().resetVictimFlag(); //todo: find better position -> fused time step mapping
    exahype::reactive::Blacklist::getInstance().recoverBlacklistedRanks();
  }
}


void exahype::reactive::OffloadingAnalyser::endIteration(double numberOfInnerLeafCells, double numberOfOuterLeafCells, double numberOfInnerCells, double numberOfOuterCells, double numberOfLocalCells, double numberOfLocalVertices) {
  if(!_isSwitchedOn) return;

  if(_iterationCounter%2 !=0) { //2 mesh sweeps for fused time stepping, todo: is that info stored somewhere?
    _iterationCounter++;
    return;
  }

  _currentAccumulatedWorkerTime = 0;
 
  for(size_t i=0; i<_waitAvgForOtherRank.size(); i++) {
    if(static_cast<int> (i) != tarch::parallel::Node::getInstance().getRank()) {
      logDebug("endIteration()", "wait for rank "<<i<<_waitAvgForOtherRank[i].toString());
      if(_waitAvgForOtherRank[i].isAccurateValue())
       exahype::reactive::PerformanceMonitor::getInstance().submitWaitingTimeForRank(_waitAvgForOtherRank[i].getValue(), i);
    }     
  }

  updateZeroTresholdAndFilteredSnapshot();
  //printWaitingTimes();
  //exahype::reactive::LocalBlacklist::getInstance().printBlacklist(); -> move to FusedTimeStep mapping

  switch(exahype::reactive::ReactiveContext::getInstance().getOffloadingStrategy()){
    case exahype::reactive::ReactiveContext::OffloadingStrategy::AggressiveHybrid:
      exahype::reactive::AggressiveHybridDistributor::getInstance().printOffloadingStatistics();
      exahype::reactive::AggressiveHybridDistributor::getInstance().updateLoadDistribution();
      break;
    default:
      //do nothing for now
      break;
  }

  _iterationCounter++;
}

void exahype::reactive::OffloadingAnalyser::printWaitingTimes() const {
  if(!_isSwitchedOn) return;
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  const std::vector<double> waitingTimesSnapshot = getFilteredWaitingTimesSnapshot();
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
    _numLateIncomingJobsWhileWaiting = 0;
    _estWalltimeForPendingJobsAtReceiveForWorker = pendingJobs * getTimePerSTP() / tarch::multicore::Core::getInstance().getNumberOfThreads();
    logDebug("beginToReceiveDataFromWorker()","there are "<<pendingJobs<<" jobs left, estimated time "<<_estWalltimeForPendingJobsAtReceiveForWorker);

#ifdef USE_ITAC
    VT_begin(event_waitForWorker);
#endif
  }
}

void exahype::reactive::OffloadingAnalyser::endToReceiveDataFromWorker( int fromRank ) {

  if (_isSwitchedOn) {
#ifdef USE_ITAC
    VT_end(event_waitForWorker);
#endif
    _waitForWorkerDataWatch.stopTimer();
    double estimatedTimeForLateIncomingJobs = _numLateIncomingJobsWhileWaiting * getTimePerSTP()/ tarch::multicore::Core::getInstance().getNumberOfThreads();
    const double elapsedTime = std::max(0.000001,_waitForWorkerDataWatch.getCalendarTime()-_estWalltimeForPendingJobsAtReceiveForWorker
                                            -estimatedTimeForLateIncomingJobs);

#if defined(OffloadingUseProfiler)
    exahype::reactive::OffloadingProfiler::getInstance().endWaitForWorker(elapsedTime);
#endif

    if(elapsedTime>_currentZeroThreshold) {
      _currentAccumulatedWorkerTime += elapsedTime;
      logDebug("endToReceiveDataFromWorker", "elapsed "<<elapsedTime<<" accumulated "<<_currentAccumulatedWorkerTime);
      _waitAvgForOtherRank[fromRank].setValue(_currentAccumulatedWorkerTime);
    }
    else {
      _waitAvgForOtherRank[fromRank].setValue(elapsedTime);
    }
#if defined(Debug)
    double currentAvg = _waitAvgForOtherRank[fromRank].getValue();
    
    if (tarch::la::greater(elapsedTime,0.0)) {
      logDebug(
        "endToReceiveDataFromWorker()",
        "rank had to wait for worker " << fromRank << " for "<< elapsedTime<<
        " currentAvg "<< currentAvg << "s"
      );
    }
#endif

    _waitForWorkerDataWatch.startTimer();
  }
}

void exahype::reactive::OffloadingAnalyser::resetWaitingTimes() {
  if(!_isSwitchedOn) return;

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  for(int i=0; i<nnodes; i++)
     _waitAvgForOtherRank[i].erase();
}

void exahype::reactive::OffloadingAnalyser::beginToReceiveDataFromMaster(int master) {
  //nop
}


void exahype::reactive::OffloadingAnalyser::endToReceiveDataFromMaster(int master) {
  if(!_isSwitchedOn) return;
  if(master==0) 
     endToReceiveDataFromGlobalMaster();
}

void exahype::reactive::OffloadingAnalyser::beginToSendDataToMaster() {
  //nop
}


void exahype::reactive::OffloadingAnalyser::endToSendDataToMaster() {
  if (_isSwitchedOn && !_waitForMasterDataWatch.isOn()) {
    int pendingJobs  = tarch::multicore::jobs::getNumberOfWaitingBackgroundJobs();
    _estWalltimeForPendingJobsAtReceiveForWorker = pendingJobs * getTimePerSTP() / tarch::multicore::Core::getInstance().getNumberOfThreads();
    _numLateIncomingJobsWhileWaiting = 0;  //reset here when waiting for global master data
    logDebug("endToSendDataToMaster()","there are "<<pendingJobs<<" left estimated time "<<_estWalltimeForPendingJobsAtReceiveForWorker); 
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

    if (tarch::la::greater(elapsedTime,0.0)) {
      logDebug(
        "endToSendDataToWorker()",
        "rank had to wait for send to worker " << worker << " for "<< elapsedTime 
      );
    }
  }
}

void exahype::reactive::OffloadingAnalyser::notifyReceivedSTPJob() {
  _numLateIncomingJobsWhileWaiting++;
}

void exahype::reactive::OffloadingAnalyser::beginToReceiveDataFromGlobalMaster() {
  //nop
}

void exahype::reactive::OffloadingAnalyser::endToReceiveDataFromGlobalMaster() {
  if (_isSwitchedOn && _waitForMasterDataWatch.isOn()) {
    _waitForMasterDataWatch.stopTimer();
    double estimatedTimeForLateIncomingJobs = _numLateIncomingJobsWhileWaiting * getTimePerSTP()/ tarch::multicore::Core::getInstance().getNumberOfThreads();
    logDebug("endToReceiveDataFromGlobalMaster()","estimate for late STPs "<<estimatedTimeForLateIncomingJobs<<"s  estimate for pending jobs "<<_estWalltimeForPendingJobsAtReceiveForWorker); 

    const double elapsedWaitTime = std::max(0.000001, _waitForMasterDataWatch.getCalendarTime()-_estWalltimeForPendingJobsAtReceiveForWorker
                                            -estimatedTimeForLateIncomingJobs);
    _estWalltimeForPendingJobsAtReceiveForWorker = 0;
    _numLateIncomingJobsWhileWaiting = 0;
    
    _waitAvgForOtherRank[0].setValue(elapsedWaitTime); //0 is global master

#if defined(OffloadingUseProfiler)
    exahype::reactive::OffloadingProfiler::getInstance().endWaitForGlobalMaster(elapsedWaitTime);
#endif

#if defined(Debug)
    double currentAvg = _waitAvgForOtherRank[0].getValue();

    if (tarch::la::greater(elapsedWaitTime,0.0)) {
      logDebug(
        "endToReceiveDataFromGlobalMaster()",
        "rank had to wait for global master "  << " for "<< elapsedWaitTime <<
        " currentAvg "<< currentAvg << "s"
      );
    }
#endif
  }
}

//noops
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

