#if  defined(SharedTBB)  && defined(Parallel) && defined(DistributedStealing)
#include "exahype/stealing/StealingAnalyser.h"

#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"
#include "tarch/logging/CommandLineLogger.h"

#include "tarch/multicore/Lock.h"
#include "tarch/multicore/Core.h"

#include <thread>

#include "exahype/stealing/PerformanceMonitor.h"
#include "exahype/stealing/DiffusiveDistributor.h"

tarch::logging::Log  exahype::stealing::StealingAnalyser::_log( "exahype::stealing::StealingAnalyser" );


exahype::stealing::StealingAnalyser::StealingAnalyser():
  _isSwitchedOn(true),
  _waitForWorkerDataWatch("exahype::stealing::StealingAnalyser", "-", false,false),
  _waitForMasterDataWatch("exahype::stealing::StealingAnalyser", "-", false,false),
  _waitForOtherRank(tarch::parallel::Node::getInstance().getNumberOfNodes()),
  _currentMaxWaitTime(0)
{
  enable(true);
}


exahype::stealing::StealingAnalyser::~StealingAnalyser() {
}


void exahype::stealing::StealingAnalyser::enable(bool value) {
  _isSwitchedOn=value;
}


void exahype::stealing::StealingAnalyser::beginIteration() {
  _currentMaxWaitTime = 0;
  exahype::stealing::DiffusiveDistributor::getInstance().resetVictimFlag(); //TODO: correct position here?
}


void exahype::stealing::StealingAnalyser::endIteration(double numberOfInnerLeafCells, double numberOfOuterLeafCells, double numberOfInnerCells, double numberOfOuterCells, double numberOfLocalCells, double numberOfLocalVertices) {

  for(int i=0; i<_waitForOtherRank.size(); i++) {
    if(i != tarch::parallel::Node::getInstance().getRank())
      logInfo("endIteration()", "wait for rank "<<i<<_waitForOtherRank[i].toString());
  }

  exahype::stealing::PerformanceMonitor::getInstance().setCurrentLoad(static_cast<int>(_currentMaxWaitTime*1e06));
}


void exahype::stealing::StealingAnalyser::beginToReceiveDataFromWorker() {
  if (_isSwitchedOn) {
    _waitForWorkerDataWatch.startTimer();
  }
}


void exahype::stealing::StealingAnalyser::endToReceiveDataFromWorker( int fromRank ) {
  if (_isSwitchedOn) {
    _waitForWorkerDataWatch.stopTimer();
    const double elapsedTime = _waitForWorkerDataWatch.getCalendarTime();

    _waitForOtherRank[fromRank].setValue(elapsedTime);

    double currentAvg = _waitForOtherRank[fromRank].getValue();

    if(_currentMaxWaitTime < currentAvg) {
      _currentMaxWaitTime = currentAvg;
      exahype::stealing::DiffusiveDistributor::getInstance().updateLoadDistribution(static_cast<int>(_currentMaxWaitTime*1e06));
    }
    
    if (tarch::la::greater(elapsedTime,0.0)) {
      logInfo(
        "endToReceiveDataFromWorker()",
        "rank had to wait for worker " << fromRank << " for "<< elapsedTime<<
        " currentAvg "<< currentAvg << "s"<<
        " currentMaxWaitTime "<<_currentMaxWaitTime<<"s"
      );
    }

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
        " currentAvg "<< currentAvg << "s"<<
        " currentMaxWaitTime "<<_currentMaxWaitTime<<"s"
      );
    }

    if(_currentMaxWaitTime< currentAvg) {
      _currentMaxWaitTime = currentAvg;
      exahype::stealing::DiffusiveDistributor::getInstance().updateLoadDistribution(static_cast<int>(_currentMaxWaitTime*1e06));
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

