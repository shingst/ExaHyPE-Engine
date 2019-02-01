#if  defined(SharedTBB)  && defined(Parallel) && defined(DistributedStealing)
#include "exahype/stealing/StealingAnalyser.h"

#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"
#include "tarch/logging/CommandLineLogger.h"

#include "tarch/multicore/Lock.h"
#include "tarch/multicore/Core.h"

#include <thread>
#include <limits>

#include "exahype/stealing/PerformanceMonitor.h"
#include "exahype/stealing/DiffusiveDistributor.h"
#include "exahype/stealing/AggressiveDistributor.h"
#include "exahype/stealing/AggressiveCCPDistributor.h"
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
  _waitForOtherRank(tarch::parallel::Node::getInstance().getNumberOfNodes()),
  _currentZeroThreshold(20000),
  _iterationCounter(0)
{
  enable(true);
#ifdef USE_ITAC
  VT_funcdef(event_name_waitForWorker, VT_NOCLASS, &event_waitForWorker ); assertion(ierr==0)
#endif
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  _currentFilteredWaitingTimesSnapshot = new int[nnodes*nnodes];
  std::fill(&_currentFilteredWaitingTimesSnapshot[0], &_currentFilteredWaitingTimesSnapshot[nnodes*nnodes], 0);
}

exahype::stealing::StealingAnalyser& exahype::stealing::StealingAnalyser::getInstance() {
  static StealingAnalyser analyser;
  return analyser;
}

exahype::stealing::StealingAnalyser::~StealingAnalyser() {
  delete[] _currentFilteredWaitingTimesSnapshot;
}


void exahype::stealing::StealingAnalyser::enable(bool value) {
  _isSwitchedOn=value;
}

const int* exahype::stealing::StealingAnalyser::getFilteredWaitingTimesSnapshot() {
  return _currentFilteredWaitingTimesSnapshot;
}

int  exahype::stealing::StealingAnalyser::getZeroThreshold() {
  return _currentZeroThreshold;
}

void exahype::stealing::StealingAnalyser::updateZeroTresholdAndFilteredSnapshot() {
  const int* currentWaitingTimesSnapshot = exahype::stealing::PerformanceMonitor::getInstance().getWaitingTimesSnapshot();
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  
  int min, max;
  min = std::numeric_limits<int>::max();
  max = std::numeric_limits<int>::min();

  for(int i=0; i<nnodes; i++) {
    for(int j=0; j<nnodes; j++) {
      int currentWaitingTime = currentWaitingTimesSnapshot[i*nnodes+j];
      // logInfo("updateZeroTresholdAndFilteredSnapshot()","rank "<<i<<" waiting for "<<currentWaitingTime<<" for rank "<<j);
      if(currentWaitingTime>max) {
	 max = currentWaitingTime;
      }
      else if(currentWaitingTime>0 && currentWaitingTime<min) {
        min = currentWaitingTime;
      }
      _currentFilteredWaitingTimesSnapshot[i*nnodes+j] = (currentWaitingTime < _currentZeroThreshold) ? 0 : currentWaitingTime;
    }
  }
  if(min < std::numeric_limits<int>::max()) {
    int newThreshold = static_cast<int>(0.9*min+0.1*max);
    _currentZeroThreshold = newThreshold;
    logInfo("updateZeroTresholdAndFilteredSnapshot()", " zero treshold set to "<< newThreshold);
  }

}

void exahype::stealing::StealingAnalyser::beginIteration() {
  if(_iterationCounter%2 !=0) return;

  exahype::stealing::StealingManager::getInstance().resetVictimFlag(); //TODO: correct position here?
  exahype::stealing::StealingManager::getInstance().decreaseHeat();
}


void exahype::stealing::StealingAnalyser::endIteration(double numberOfInnerLeafCells, double numberOfOuterLeafCells, double numberOfInnerCells, double numberOfOuterCells, double numberOfLocalCells, double numberOfLocalVertices) {
  if(_iterationCounter%2 !=0) {
     _iterationCounter++; 
     return;
  }

  for(int i=0; i<_waitForOtherRank.size(); i++) {
    if(i != tarch::parallel::Node::getInstance().getRank()) {
      logInfo("endIteration()", "wait for rank "<<i<<_waitForOtherRank[i].toString());
      exahype::stealing::PerformanceMonitor::getInstance().submitWaitingTimeForRank(static_cast<int>(_waitForOtherRank[i].getValue()*1e06), i);
    }     
  }

  updateZeroTresholdAndFilteredSnapshot();
#if defined(StealingStrategyDiffusive)
  exahype::stealing::DiffusiveDistributor::getInstance().updateLoadDistribution();
#elif defined(StealingStrategyAggressive)
  exahype::stealing::AggressiveDistributor::getInstance().printOffloadingStatistics();
  exahype::stealing::AggressiveDistributor::getInstance().updateLoadDistribution();
#elif defined(StealingStrategyAggressiveCCP)
  exahype::stealing::AggressiveCCPDistributor::getInstance().printOffloadingStatistics();
  exahype::stealing::AggressiveCCPDistributor::getInstance().updateLoadDistribution();
#endif

  _iterationCounter++;
}


void exahype::stealing::StealingAnalyser::beginToReceiveDataFromWorker() {
  if (_isSwitchedOn) {
    _waitForWorkerDataWatch.startTimer();
#ifdef USE_ITAC
    VT_begin(event_waitForWorker);
#endif
  }
}


void exahype::stealing::StealingAnalyser::endToReceiveDataFromWorker( int fromRank ) {
  if (_isSwitchedOn) {
#ifdef USE_ITAC
    VT_end(event_waitForWorker);
#endif
    _waitForWorkerDataWatch.stopTimer();
    const double elapsedTime = _waitForWorkerDataWatch.getCalendarTime();

    _waitForOtherRank[fromRank].setValue(elapsedTime);

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
    
    _waitForOtherRank[myMaster].setValue(elapsedTime);

    double currentAvg = _waitForOtherRank[myMaster].getValue();

    if (tarch::la::greater(elapsedTime,0.0)) {
      logInfo(
        "endToReceiveDataFromMaster()",
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

