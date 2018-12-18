#if  defined(SharedTBB)  && defined(Parallel) && defined(DistributedStealing)
#include "PerformanceMonitor.h"

#include <algorithm>
#include <numeric>

#include "tarch/multicore/Lock.h"
#include "tarch/parallel/Node.h"
#include "tarch/timing/Watch.h"

#include "exahype/stealing/StealingProfiler.h"
#include "exahype/stealing/DynamicDistributor.h"
#include "exahype/stealing/StaticDistributor.h"
#include "exahype/stealing/StealingManager.h"


#define TERMINATE_SIGNAL -1

tarch::logging::Log exahype::stealing::PerformanceMonitor::_log( "exahype::stealing::PerformanceMonitor" );

exahype::stealing::PerformanceMonitor::PerformanceMonitor() :
    _isStarted(true),
    _gatherTasksRequest(MPI_REQUEST_NULL),
    _gatherWaitingTimesRequest(MPI_REQUEST_NULL),
    _currentTasks(0),
    _currentTasksSendBuffer(0),
    _remainingTasks(0),
    _tasksPerTimestep(0) {

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  _currentTasksSnapshot = new int[nnodes];
  _currentTasksReceiveBuffer   = new int[nnodes];

  std::fill(_currentTasksSnapshot, _currentTasksSnapshot+nnodes, 0);
  std::fill(_currentTasksReceiveBuffer, _currentTasksReceiveBuffer+nnodes, 0);
}

exahype::stealing::PerformanceMonitor::~PerformanceMonitor() {
  delete[] _currentTasksSnapshot;
  delete[] _currentTasksReceiveBuffer;
}

void exahype::stealing::PerformanceMonitor::setLocalLoadPerTimestep(int load) {
   logInfo("setLocalLoadPerTimestep", "setting local load per timestep to "<<load);
  _tasksPerTimestep = load;
  _remainingTasks = _tasksPerTimestep;
}

int exahype::stealing::PerformanceMonitor::getLocalLoadPerTimestep() {
  return _tasksPerTimestep;
}

int exahype::stealing::PerformanceMonitor::getRemainingLocalLoad() {
  return _remainingTasks;
}

const int* exahype::stealing::PerformanceMonitor::getCurrentLoadSnapshot() {
  return _currentTasksSnapshot;
}

exahype::stealing::PerformanceMonitor& exahype::stealing::PerformanceMonitor::getInstance() {
  static PerformanceMonitor perfMon;
  return perfMon;
}

void exahype::stealing::PerformanceMonitor::stop() {
    _isStarted=false;
}

void exahype::stealing::PerformanceMonitor::setCurrentLoad(int num) {
  logInfo("performance monitor", "setting current load to "<<num);
  int myRank = tarch::parallel::Node::getInstance().getRank();
  tarch::multicore::Lock lock(_semaphore);
  _currentTasksSnapshot[myRank] = num;
  _currentTasks = num;
  lock.free();
}

void exahype::stealing::PerformanceMonitor::incCurrentLoad() {
#ifndef StealingStrategyDiffusive
  assertion(_currentTasks>=0);
  _currentTasks++;
#endif
}

void exahype::stealing::PerformanceMonitor::decCurrentLoad() {
#ifndef StealingStrategyDiffusive
  _currentTasks--;
  assertion(_currentTasks>=0);
#endif
}

void exahype::stealing::PerformanceMonitor::decRemainingLocalLoad() {
#ifndef StealingStrategyDiffusive
  tarch::multicore::Lock lock(_semaphore);
  _remainingTasks--;
  if(_remainingTasks==0) {
    _remainingTasks=_tasksPerTimestep;
  }
  lock.free();
  assertion(_remainingTasks>=0);
#endif
}

void exahype::stealing::PerformanceMonitor::run() {
  progressGather();
}

void exahype::stealing::PerformanceMonitor::progressGather() {
  int myRank    = tarch::parallel::Node::getInstance().getRank();
  int nnodes    = tarch::parallel::Node::getInstance().getNumberOfNodes();

  int completed = 0;
#if defined(PerformanceAnalysisStealing)
  double timeSinceLastGather=0;
  static std::atomic<double> lastGather = 0;
  static std::atomic<int> unsuccessful  = 0;
  static std::atomic<int> successful    = 0;

  tarch::timing::Watch watch("exahype::stealing::", "-", false,false);
  watch.startTimer();

  timeSinceLastGather = lastGather + MPI_Wtime();
#endif

  tarch::multicore::Lock lock(_semaphore);

  if( !isGloballyTerminated() && _gatherTasksRequest!=MPI_REQUEST_NULL) {
	double time = - MPI_Wtime();
    exahype::stealing::StealingProfiler::getInstance().beginCommunication();
    MPI_Test(&_gatherTasksRequest, &completed, MPI_STATUS_IGNORE);
    time += MPI_Wtime();

#if defined(PerformanceAnalysisStealing)
    if(completed) {
      exahype::stealing::StealingProfiler::getInstance().endCommunication(true, time);
      successful++;
    }
    else {
      exahype::stealing::StealingProfiler::getInstance().endCommunication(false, time);
      unsuccessful++;
    }
    if(successful%1000==0 || unsuccessful%10000==0) {
      logInfo("performance monitor", " successful "<<successful<<" unsuccessful "<<unsuccessful<<" ratio "<<1.0f*successful/unsuccessful);
    }
#endif
  }

#if defined(PerformanceAnalysisStealing)
  watch.stopTimer();
  if(watch.getCalendarTime() >= 0.00001) {
    logInfo(
        "performance monitor()",
        "MPI " <<
        "time=" << std::fixed <<
        watch.getCalendarTime() <<
        ", cpu time=" <<
        watch.getCPUTime()
    );
  }
  watch.startTimer();
#endif

  if(completed) {
    stealing::StealingProfiler::getInstance().notifyPerformanceUpdate();
    std::copy(&_currentTasksReceiveBuffer[0], &_currentTasksReceiveBuffer[nnodes], &_currentTasksSnapshot[0]);
#if defined(PerformanceAnalysisStealing)
    std::string str="received new update, current load "+std::to_string(_currentTasks.load());
//    if(timeSinceLastGather>0.001) {
//      str=str+ " took too long: "+std::to_string(timeSinceLastGather);
//      stealing::StealingProfiler::getInstance().notifyLatePerformanceUpdate();
//    }
    for(int i=0;i<nnodes;i++) str=str+" , "+std::to_string(_currentTasksBuffer[i]);
    str+="\n";
    logInfo("performance monitor", str);
#endif
    if(_currentTasks.load()>0 && std::all_of(&_currentTasksSnapshot[0], &_currentTasksSnapshot[nnodes], [](int i) {return i>=0;})) {
      exahype::stealing::DynamicDistributor::getInstance().computeNewLoadDistribution(_currentTasksSnapshot);
      stealing::StealingProfiler::getInstance().notifyStealingDecision();
    }
#if defined(PerformanceAnalysisStealing)
    watch.stopTimer();
    if(watch.getCalendarTime() >= 0.00001) {
      logInfo(
          "performance monitor()",
          "completion" <<
          "time=" << std::fixed <<
          watch.getCalendarTime() <<
          ", cpu time=" <<
          watch.getCPUTime()
      );
    }
    watch.startTimer();
#endif
  }

  if(_gatherTasksRequest==MPI_REQUEST_NULL && !isGloballyTerminated()) {
    postGatherTasks();
#if defined(PerformanceAnalysisStealing)
    lastGather=-MPI_Wtime();
    watch.stopTimer();
    if(watch.getCalendarTime() >= 0.00001) {
      logInfo(
          "performance monitor()",
          "post gather" <<
          "time=" << std::fixed <<
          watch.getCalendarTime() <<
          ", cpu time=" <<
          watch.getCPUTime()
      );
    }
    watch.startTimer();
#endif
  }
  lock.free();
}

void exahype::stealing::PerformanceMonitor::postGatherTasks() {
  int myRank = tarch::parallel::Node::getInstance().getRank();

  _currentTasksSendBuffer = _isStarted ? _currentTasks.load() : TERMINATE_SIGNAL;

  MPI_Iallgather(&_currentTasksSendBuffer, 1, MPI_INTEGER, _currentTasksReceiveBuffer, 1, MPI_INTEGER, exahype::stealing::StealingManager::getInstance().getMPICommunicator(), &_gatherTasksRequest);
}

bool exahype::stealing::PerformanceMonitor::isGloballyTerminated() {
  static bool globalTermination=false;
  if(globalTermination) return true;
  if(_isStarted) return false;

  bool result=true;
  for(int i=0; i<tarch::parallel::Node::getInstance().getNumberOfNodes(); i++)
    result &= _currentTasksReceiveBuffer[i]==TERMINATE_SIGNAL;

  globalTermination=result;
  return globalTermination;
}

#endif
