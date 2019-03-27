#if  defined(SharedTBB)  && defined(Parallel) && (defined(DistributedStealing) || defined(AnalyseWaitingTimes))
#include "exahype/stealing/StealingAnalyser.h"

#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"
#include "tarch/logging/CommandLineLogger.h"

#include "tarch/multicore/Lock.h"
#include "tarch/multicore/Core.h"
#include "tarch/multicore/Jobs.h"

#include <thread>
#include <limits>

#include "exahype/stealing/PerformanceMonitor.h"
#include "exahype/stealing/DiffusiveDistributor.h"
#include "exahype/stealing/AggressiveDistributor.h"
#include "exahype/stealing/AggressiveCCPDistributor.h"
#include "exahype/stealing/AggressiveHybridDistributor.h"
#include "exahype/stealing/StealingManager.h"

tarch::logging::Log  exahype::stealing::StealingAnalyser::_log( "exahype::stealing::StealingAnalyser" );

#ifdef USE_ITAC
#include "VT.h"
static int event_waitForWorker = -1;
static const char *event_name_waitForWorker = "waitForWorker";
#endif

exahype::stealing::StealingAnalyser::StealingAnalyser():
  _isSwitchedOn(true),
  _waitForWorkerDataWatch("exahype::stealing::StealingAnalyser", "-", false,false),
  _waitForMasterDataWatch("exahype::stealing::StealingAnalyser", "-", false,false),
  _waitForGlobalMasterDataWatch("exahype::stealing::StealingAnalyser", "-", false,false),
  _waitForOtherRank(0),
  _currentZeroThreshold(0.002),
  _iterationCounter(0),
  _currentAccumulatedWorkerTime(0),
  _estimatedWtimeForPendingJobs(0),
  _lateSTPJobs(0)
{
  enable(true);
#ifdef USE_ITAC
  VT_funcdef(event_name_waitForWorker, VT_NOCLASS, &event_waitForWorker ); assertion(ierr==0)
#endif
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  _currentFilteredWaitingTimesSnapshot = new double[nnodes*nnodes];
  std::fill(&_currentFilteredWaitingTimesSnapshot[0], &_currentFilteredWaitingTimesSnapshot[nnodes*nnodes], 0);

  for(int i=0; i<nnodes; i++) 
     //_waitForOtherRank.push_back(tarch::timing::GlidingAverageMeasurement(0.1,16));
     _waitForOtherRank.push_back(tarch::timing::GlidingAverageMeasurement());
}


exahype::stealing::StealingAnalyser& exahype::stealing::StealingAnalyser::getInstance() {
  static StealingAnalyser* analyser = nullptr;

  if(analyser==nullptr) {
    analyser = new StealingAnalyser();
  }  


  return *analyser;
}

exahype::stealing::StealingAnalyser::~StealingAnalyser() {
    delete[] _currentFilteredWaitingTimesSnapshot;
}


void exahype::stealing::StealingAnalyser::enable(bool value) {
  _isSwitchedOn=value;
}

const double* exahype::stealing::StealingAnalyser::getFilteredWaitingTimesSnapshot() {
  return _currentFilteredWaitingTimesSnapshot;
}

double exahype::stealing::StealingAnalyser::getZeroThreshold() {
  return _currentZeroThreshold;
}

void exahype::stealing::StealingAnalyser::setTimePerSTP(double timePerSTP) {
  _timePerStealablePredictionJob.setValue(timePerSTP);
  //logInfo("setTimePerSTP()", "submitted new STP measurement, current time per stp: "<<getTimePerSTP());
}

double exahype::stealing::StealingAnalyser::getTimePerSTP() {
  return _timePerStealablePredictionJob.getValue();
}

void exahype::stealing::StealingAnalyser::updateZeroTresholdAndFilteredSnapshot() {
#if !defined(AnalyseWaitingTimes)
  const double* currentWaitingTimesSnapshot = exahype::stealing::PerformanceMonitor::getInstance().getWaitingTimesSnapshot();
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  //check if valid: there is an entry>0 for every rank
  bool isValid = true;
  for(int i=0; i<nnodes; i++) {
    bool hasEntry = false;
    for(int j=0; j<nnodes; j++) {
      double currentWaitingTime = currentWaitingTimesSnapshot[i*nnodes+j];
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
  
  if(!isValid) return;

  double min, max;
  min = std::numeric_limits<double>::max();
  max = std::numeric_limits<double>::min();


  for(int i=0; i<nnodes; i++) {
    for(int j=0; j<nnodes; j++) {
      double currentWaitingTime = currentWaitingTimesSnapshot[i*nnodes+j];
      if(currentWaitingTime>0)
        logInfo("updateZeroTresholdAndFilteredSnapshot()","rank "<<i<<" waiting for "<<currentWaitingTime<<" for rank "<<j);
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
    logInfo("updateZeroTresholdAndFilteredSnapshot()", " zero threshold set to "<< newThreshold);
  }
#endif
}

void exahype::stealing::StealingAnalyser::beginIteration() {
#if !defined(AnalyseWaitingTimes)
  if(_iterationCounter%2 !=0) return;

  exahype::stealing::StealingManager::getInstance().resetVictimFlag(); //TODO: correct position here?
  exahype::stealing::StealingManager::getInstance().decreaseHeat();
#endif
}


void exahype::stealing::StealingAnalyser::endIteration() {
#if !defined(AnalyseWaitingTimes)
  if(_iterationCounter%2 !=0) {
     _iterationCounter++; 
     return;
  }
  _currentAccumulatedWorkerTime = 0;
 
  for(int i=0; i<_waitForOtherRank.size(); i++) {
    if(i != tarch::parallel::Node::getInstance().getRank()) {
      logInfo("endIteration()", "wait for rank "<<i<<_waitForOtherRank[i].toString());
      if(_waitForOtherRank[i].isAccurateValue())
       exahype::stealing::PerformanceMonitor::getInstance().submitWaitingTimeForRank(_waitForOtherRank[i].getValue(), i);
    }     
  }

  updateZeroTresholdAndFilteredSnapshot();
  printWaitingTimes();
#if defined(StealingStrategyDiffusive)
  exahype::stealing::DiffusiveDistributor::getInstance().updateLoadDistribution();
#elif defined(StealingStrategyAggressive)
  exahype::stealing::AggressiveDistributor::getInstance().printOffloadingStatistics();
  exahype::stealing::AggressiveDistributor::getInstance().updateLoadDistribution();
#elif defined(StealingStrategyAggressiveCCP)
  exahype::stealing::AggressiveCCPDistributor::getInstance().printOffloadingStatistics();
  exahype::stealing::AggressiveCCPDistributor::getInstance().updateLoadDistribution();
#elif defined(StealingStrategyAggressiveHybrid)
  exahype::stealing::AggressiveHybridDistributor::getInstance().printOffloadingStatistics();
  exahype::stealing::AggressiveHybridDistributor::getInstance().updateLoadDistribution();
#endif
  exahype::stealing::StealingManager::getInstance().printBlacklist();

  _iterationCounter++;
#endif
}

void exahype::stealing::StealingAnalyser::printWaitingTimes() {
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

void exahype::stealing::StealingAnalyser::beginToReceiveDataFromWorker() {
  if (_isSwitchedOn) {
    _waitForWorkerDataWatch.startTimer();
    int pendingJobs  = tarch::multicore::jobs::getNumberOfWaitingBackgroundJobs();
    _lateSTPJobs = 0;
    _estimatedWtimeForPendingJobs = pendingJobs * getTimePerSTP() / tarch::multicore::Core::getInstance().getNumberOfThreads();
    logInfo("beginToReceiveDataFromWorker()","there are "<<pendingJobs<<" jobs left, estimated time "<<_estimatedWtimeForPendingJobs);

#ifdef USE_ITAC
    VT_begin(event_waitForWorker);
#endif
  }
}


void exahype::stealing::StealingAnalyser::endToReceiveDataFromWorker( int fromRank ) {

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
      logInfo("endToReceiveDataFromWorker", "elapsed "<<elapsedTime<<" accumulated "<<_currentAccumulatedWorkerTime);
      _waitForOtherRank[fromRank].setValue(_currentAccumulatedWorkerTime);
    }
    else {
      _waitForOtherRank[fromRank].setValue(elapsedTime);
    }

    double currentAvg = _waitForOtherRank[fromRank].getValue();
    
    if (tarch::la::greater(elapsedTime,0.0)) {
      logInfo(
        "endToReceiveDataFromWorker()",
        "rank had to wait for worker " << fromRank << " for "<< elapsedTime<<
        " currentAvg "<< currentAvg << "s"
      );
    }
//#ifdef USE_ITAC
//        VT_begin(event_waitForWorker);
//#endif
    _waitForWorkerDataWatch.startTimer();
  }
}


void exahype::stealing::StealingAnalyser::beginToReceiveDataFromMaster() {
  if (_isSwitchedOn && !_waitForMasterDataWatch.isOn()) {
    _waitForMasterDataWatch.startTimer();
  }
}


void exahype::stealing::StealingAnalyser::endToReceiveDataFromMaster() {
  if (_isSwitchedOn && _waitForMasterDataWatch.isOn()) {
    _waitForMasterDataWatch.stopTimer();
    const double elapsedTime = _waitForMasterDataWatch.getCalendarTime();
    int myMaster = tarch::parallel::NodePool::getInstance().getMasterRank();
    
    //_waitForOtherRank[myMaster].setValue(elapsedTime);

    //double currentAvg = _waitForOtherRank[myMaster].getValue();

    if (tarch::la::greater(elapsedTime,0.0)) {
      logInfo(
        "endToReceiveDataFromMaster()",
        "rank had to wait for master " << myMaster << " for "<< elapsedTime 
      );
    }
  }
}

void exahype::stealing::StealingAnalyser::beginToSendDataToMaster() {
  //if (_isSwitchedOn && !_waitForMasterDataWatch.isOn()) {
  //  _waitForMasterDataWatch.startTimer();
  //}
}


void exahype::stealing::StealingAnalyser::endToSendDataToMaster() {
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
    logInfo("endToSendDataToMaster()","there are "<<pendingJobs<<" left estimated time "<<_estimatedWtimeForPendingJobs); 
    _waitForMasterDataWatch.startTimer();
  }
}

void exahype::stealing::StealingAnalyser::beginToSendDataToWorker() {
  if (_isSwitchedOn && !_waitForWorkerDataWatch.isOn()) {
    _waitForWorkerDataWatch.startTimer();
  }
}


void exahype::stealing::StealingAnalyser::endToSendDataToWorker(int worker) {
  if (_isSwitchedOn && _waitForWorkerDataWatch.isOn()) {
    _waitForWorkerDataWatch.stopTimer();
    const double elapsedTime = _waitForWorkerDataWatch.getCalendarTime();

    //_waitForOtherRank[worker].setValue(elapsedTime);

    //double currentAvg = _waitForOtherRank[worker].getValue();

    if (tarch::la::greater(elapsedTime,0.0)) {
      logInfo(
        "endToSendDataToWorker()",
        "rank had to wait for send to worker " << worker << " for "<< elapsedTime 
      );
    }
  }
}

void exahype::stealing::StealingAnalyser::notifyReceivedSTPJob() {
  _lateSTPJobs++;
}

void exahype::stealing::StealingAnalyser::beginToReceiveDataFromGlobalMaster() {
  //if (_isSwitchedOn && !_waitForGlobalMasterDataWatch.isOn()) {
    //_waitForGlobalMasterDataWatch.startTimer();
  //}
}


void exahype::stealing::StealingAnalyser::endToReceiveDataFromGlobalMaster() {
  if (_isSwitchedOn && _waitForMasterDataWatch.isOn()) {
    _waitForMasterDataWatch.stopTimer();
    double estimatedTimeForLateSTPs = _lateSTPJobs * getTimePerSTP()/ tarch::multicore::Core::getInstance().getNumberOfThreads();
    logInfo("endToReceiveDataFromGlobalMaster()","estimate for late STPs "<<estimatedTimeForLateSTPs<<"s  estimate for pending jobs "<<_estimatedWtimeForPendingJobs); 
    const double elapsedTime = std::max(0.000001, _waitForMasterDataWatch.getCalendarTime()-_estimatedWtimeForPendingJobs
                                            -estimatedTimeForLateSTPs);
    _estimatedWtimeForPendingJobs = 0;
    _lateSTPJobs = 0;
    int myMaster = 0;
    
    _waitForOtherRank[0].setValue(elapsedTime); //0 is global master

    double currentAvg = _waitForOtherRank[0].getValue();

    if (tarch::la::greater(elapsedTime,0.0)) {
      logInfo(
        "endToReceiveDataFromGlobalMaster()",
        "rank had to wait for master " << myMaster << " for "<< elapsedTime <<
        " currentAvg "<< currentAvg << "s"
      );
    }
  }
}

void exahype::stealing::StealingAnalyser::dataWasNotReceivedInBackground( int fromRank, int tag, int cardinality, int pageSize ) {
}


void exahype::stealing::StealingAnalyser::beginToReleaseSynchronousHeapData() {
}


void exahype::stealing::StealingAnalyser::endToReleaseSynchronousHeapData() {
}


void exahype::stealing::StealingAnalyser::beginToPrepareAsynchronousHeapDataExchange() {
}


void exahype::stealing::StealingAnalyser::endToPrepareAsynchronousHeapDataExchange() {
}


void exahype::stealing::StealingAnalyser::beginReleaseOfJoinData() {
}


void exahype::stealing::StealingAnalyser::endReleaseOfJoinData() {
}


void exahype::stealing::StealingAnalyser::beginReleaseOfBoundaryData() {
}


void exahype::stealing::StealingAnalyser::endReleaseOfBoundaryData() {
}


void exahype::stealing::StealingAnalyser::changeConcurrencyLevel(int actualChange, int maxPossibleChange) {
}


void exahype::stealing::StealingAnalyser::minuteNumberOfBackgroundTasks(int taskCount) {
}


void exahype::stealing::StealingAnalyser::beginProcessingBackgroundJobs() {
}


void exahype::stealing::StealingAnalyser::endProcessingBackgroundJobs() {
}


void exahype::stealing::StealingAnalyser::enterCentralElementOfEnclosingSpacetree() {
}


void exahype::stealing::StealingAnalyser::leaveCentralElementOfEnclosingSpacetree() {
}


void exahype::stealing::StealingAnalyser::addWorker(
  int                                 workerRank,
  int                                 level
) {
}


void exahype::stealing::StealingAnalyser::removeWorker(
  int                                 workerRank,
  int                                 level
) {
}
#endif

