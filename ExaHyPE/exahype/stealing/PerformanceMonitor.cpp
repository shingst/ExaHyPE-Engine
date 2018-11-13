#if  defined(SharedTBB)  && defined(Parallel)
#include "PerformanceMonitor.h"

#include <algorithm>
#include <numeric>

#include "tarch/multicore/Lock.h"
#include "tarch/parallel/Node.h"
#include "tarch/timing/Watch.h"

#include "exahype/stealing/StealingProfiler.h"
#include "exahype/stealing/DynamicDistributor.h"
#include "exahype/stealing/StaticDistributor.h"

#define TERMINATE_SIGNAL -1

tarch::logging::Log exahype::stealing::PerformanceMonitor::_log( "exahype::stealing::PerformanceMonitor" );

exahype::stealing::PerformanceMonitor::PerformanceMonitor() :
    _isStarted(true),
    _gather_request(MPI_REQUEST_NULL),
    _currentLoadLocal(0),
    _currentLoadLocalBuffer(0),
    _remainingLoadLocal(0),
    _localLoadPerTimestep(0) {

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  _currentLoadSnapshot = new int[nnodes];
  _currentLoadBuffer   = new int[nnodes];

  std::fill(_currentLoadSnapshot, _currentLoadSnapshot+nnodes, 0);
  std::fill(_currentLoadBuffer, _currentLoadBuffer+nnodes, 0);
}

exahype::stealing::PerformanceMonitor::~PerformanceMonitor() {
  delete[] _currentLoadSnapshot;
  delete[] _currentLoadBuffer;
}

void exahype::stealing::PerformanceMonitor::setLocalLoadPerTimestep(int load) {
   logInfo("setLocalLoadPerTimestep", "setting local load per timestep to "<<load);
  _localLoadPerTimestep = load;
  _remainingLoadLocal   = _localLoadPerTimestep;
}

int exahype::stealing::PerformanceMonitor::getLocalLoadPerTimestep() {
  return _localLoadPerTimestep;
}

int exahype::stealing::PerformanceMonitor::getRemainingLocalLoad() {
  return _remainingLoadLocal;
}

exahype::stealing::PerformanceMonitor& exahype::stealing::PerformanceMonitor::getInstance() {
  static PerformanceMonitor perfMon;
  return perfMon;
}

void exahype::stealing::PerformanceMonitor::stop() {
  //logInfo("stop()", "stopping");
  //assertion(_isStarted);
    _isStarted=false;
}

//void exahype::stealing::PerformanceMonitor::setCurrentLoad(int num) {
//  logInfo("performance monitor", "setting current load to "<<num);
//  int myRank = tarch::parallel::Node::getInstance().getRank();
//  tarch::multicore::Lock lock(_semaphore);
//  _currentLoadSnapshot[myRank]=num;
//  _currentLoadLocal=num;
//  lock.free();
//}

void exahype::stealing::PerformanceMonitor::incCurrentLoad() {
  assertion(_currentLoadLocal>=0);
  _currentLoadLocal++;
}

void exahype::stealing::PerformanceMonitor::decCurrentLoad() {
  _currentLoadLocal--;
  assertion(_currentLoadLocal>=0);
}

void exahype::stealing::PerformanceMonitor::decRemainingLocalLoad() {
  tarch::multicore::Lock lock(_semaphore);
  _remainingLoadLocal--;
  if(_remainingLoadLocal==0) {
    _remainingLoadLocal=_localLoadPerTimestep;
  }
  lock.free();
  assertion(_remainingLoadLocal>=0);
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

  if( !isGloballyTerminated() && _gather_request!=MPI_REQUEST_NULL) {
	double time = - MPI_Wtime();
    exahype::stealing::StealingProfiler::getInstance().beginCommunication();
    MPI_Test(&_gather_request, &completed, MPI_STATUS_IGNORE);
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
    std::copy(&_currentLoadBuffer[0], &_currentLoadBuffer[nnodes], &_currentLoadSnapshot[0]);
#if defined(PerformanceAnalysisStealing)
    std::string str="received new update, current load "+std::to_string(_currentLoadLocal.load());
    if(timeSinceLastGather>0.001) {
      str=str+ " took too long: "+std::to_string(timeSinceLastGather);
      stealing::StealingProfiler::getInstance().notifyLatePerformanceUpdate();
    }
    for(int i=0;i<nnodes;i++) str=str+" , "+std::to_string(_currentLoadBuffer[i]);
    str+="\n";
    logInfo("performance monitor", str);
#endif
    if(_currentLoadLocal.load()>0 && std::all_of(&_currentLoadSnapshot[0], &_currentLoadSnapshot[nnodes], [](int i) {return i>=0;})) {
      exahype::stealing::DynamicDistributor::getInstance().computeNewLoadDistribution(_currentLoadSnapshot);
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

  if(_gather_request==MPI_REQUEST_NULL && !isGloballyTerminated()) {
    postGather();
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

void exahype::stealing::PerformanceMonitor::postGather() {
  int myRank = tarch::parallel::Node::getInstance().getRank();

  _currentLoadLocalBuffer = _isStarted ? _currentLoadLocal.load() : TERMINATE_SIGNAL;

  MPI_Iallgather(&_currentLoadLocalBuffer, 1, MPI_INTEGER, _currentLoadBuffer, 1, MPI_INTEGER, MPI_COMM_WORLD, &_gather_request);
}

bool exahype::stealing::PerformanceMonitor::isGloballyTerminated() {
  static bool globalTermination=false;
  if(globalTermination) return true;
  if(_isStarted) return false;

  bool result=true;
  for(int i=0; i<tarch::parallel::Node::getInstance().getNumberOfNodes(); i++)
    result &= _currentLoadBuffer[i]==TERMINATE_SIGNAL;

  globalTermination=result;
  return globalTermination;
}

#endif
