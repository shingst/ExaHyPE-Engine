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
    _allreduceBlacklistRequest(MPI_REQUEST_NULL),
    _currentTasks(0),
    _currentTasksSendBuffer(0),
    _remainingTasks(0),
    _tasksPerTimestep(0) {

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  _currentTasksSnapshot = new int[nnodes];
  _currentTasksReceiveBuffer   = new int[nnodes];

  _currentWaitingTimesSnapshot = new int[nnodes*nnodes];
  _currentWaitingTimesReceiveBuffer = new int[nnodes*nnodes];
  _currentWaitingTimesSendBuffer = new int[nnodes];
  _currentWaitingTimes = new int[nnodes];

  _currentBlacklistSnapshot = new double[nnodes*nnodes];
  _currentBlacklistReceiveBuffer = new double[nnodes*nnodes];
  _currentBlacklistSendBuffer = new double[nnodes];
  _currentBlacklist = new double[nnodes];

  std::fill(_currentTasksSnapshot, _currentTasksSnapshot+nnodes, 0);
  std::fill(_currentTasksReceiveBuffer, _currentTasksReceiveBuffer+nnodes, 0);

  std::fill(_currentWaitingTimesSnapshot, _currentWaitingTimesSnapshot+nnodes*nnodes, 0);
  std::fill(_currentWaitingTimesReceiveBuffer, _currentWaitingTimesReceiveBuffer+nnodes*nnodes, 0);
  std::fill(_currentWaitingTimesSendBuffer, _currentWaitingTimesSendBuffer+nnodes, 0);
  std::fill(_currentWaitingTimes, _currentWaitingTimes+nnodes, 0);

  std::fill(_currentBlacklistSnapshot, _currentBlacklistSnapshot+nnodes, 0);
  std::fill(_currentBlacklistReceiveBuffer, _currentBlacklistReceiveBuffer+nnodes, 0);
  std::fill(_currentBlacklistSendBuffer, _currentBlacklistSendBuffer+nnodes, 0);
  std::fill(_currentBlacklist, _currentBlacklist+nnodes, 0);

}

exahype::stealing::PerformanceMonitor::~PerformanceMonitor() {
  delete[] _currentWaitingTimesSendBuffer;
  delete[] _currentWaitingTimesReceiveBuffer;
  delete[] _currentWaitingTimesSnapshot;
  delete[] _currentWaitingTimes;
 
  delete[] _currentTasksSnapshot;
  delete[] _currentTasksReceiveBuffer;
}

void exahype::stealing::PerformanceMonitor::submitWaitingTimeForRank(int waitingTime, int rank) {
//  if(waitingTime>0)
 // logInfo("submitWaitingTimes", "submitting new waiting time "<<waitingTime<< " for rank "<<rank);
  _currentWaitingTimes[rank] =  waitingTime;
}

const int *exahype::stealing::PerformanceMonitor::getWaitingTimesSnapshot() {
  return _currentWaitingTimesSnapshot;
}

void exahype::stealing::PerformanceMonitor::submitBlacklistValueForRank(double bval, int rank) {
  logInfo("submitBlacklistValue", "new value "<<bval<<" for "<<rank);
  _currentBlacklist[rank] = bval;
}

const double *exahype::stealing::PerformanceMonitor::getBlacklistSnapshot() {
  return _currentBlacklistSnapshot;
}

void exahype::stealing::PerformanceMonitor::setTasksPerTimestep(int load) {
   logInfo("setLocalLoadPerTimestep", "setting local load per timestep to "<<load);
  _tasksPerTimestep = load;
  _remainingTasks = _tasksPerTimestep;
}

int exahype::stealing::PerformanceMonitor::getTasksPerTimestep() {
  return _tasksPerTimestep;
}

int exahype::stealing::PerformanceMonitor::getRemainingTasks() {
  return _remainingTasks;
}

const int* exahype::stealing::PerformanceMonitor::getCurrentTasksSnapshot() {
  return _currentTasksSnapshot;
}

exahype::stealing::PerformanceMonitor& exahype::stealing::PerformanceMonitor::getInstance() {
  static PerformanceMonitor perfMon;
  return perfMon;
}

void exahype::stealing::PerformanceMonitor::stop() {
    _isStarted=false;
}

void exahype::stealing::PerformanceMonitor::setCurrentTasks(int num) {
  logInfo("performance monitor", "setting current load to "<<num);
  int myRank = tarch::parallel::Node::getInstance().getRank();
  tarch::multicore::Lock lock(_semaphore);
  _currentTasksSnapshot[myRank] = num;
  _currentTasks = num;
  lock.free();
}

void exahype::stealing::PerformanceMonitor::incCurrentTasks() {
#ifndef StealingStrategyDiffusive
  assertion(_currentTasks>=0);
  _currentTasks++;
#endif
}

void exahype::stealing::PerformanceMonitor::decCurrentTasks() {
#ifndef StealingStrategyDiffusive
  _currentTasks--;
  assertion(_currentTasks>=0);
#endif
}

void exahype::stealing::PerformanceMonitor::decRemainingTasks() {
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

  int completed_tasks = 0;
  int completed_waiting_times = 0;
  int completed_blacklist = 0;
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
    int err = MPI_Test(&_gatherTasksRequest, &completed_tasks, MPI_STATUS_IGNORE); //assert(err==MPI_SUCCESS);
    time += MPI_Wtime();

    //std::string str;
    //for(int i=0;i<nnodes;i++) str=str+" , "+std::to_string(_currentTasksReceiveBuffer[i]);
    //str+="\n";
    //logInfo("performance monitor", "progressing "<<_gatherTasksRequest);
    

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

  if( !isGloballyTerminated() && _gatherWaitingTimesRequest!=MPI_REQUEST_NULL) {
//    double time = - MPI_Wtime();
//    exahype::stealing::StealingProfiler::getInstance().beginCommunication();
    //logInfo("performance monitoŕ()","progressing waiting times");
    int err= MPI_Test(&_gatherWaitingTimesRequest, &completed_waiting_times, MPI_STATUS_IGNORE);
   //assert(err==MPI_SUCCESS);
//    time += MPI_Wtime();_currentTasksSendBuffer
 
  }

  if( !isGloballyTerminated() && _allreduceBlacklistRequest!=MPI_REQUEST_NULL) {
     int err= MPI_Test(&_allreduceBlacklistRequest, &completed_blacklist, MPI_STATUS_IGNORE);
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

  if(completed_tasks) {
  //  logInfo("progressGather","collected new tasks snapshot");
    stealing::StealingProfiler::getInstance().notifyPerformanceUpdate();
    std::copy(&_currentTasksReceiveBuffer[0], &_currentTasksReceiveBuffer[nnodes], &_currentTasksSnapshot[0]);
 
    //std::string str;
    //for(int i=0;i<nnodes;i++) str=str+" , "+std::to_string(_currentTasksSnapshot[i]);
    //str+="\n";
    //logInfo("performance monitor", str);

#if defined(PerformanceAnalysisStealing)
    std::string str="received new update, current load "+std::to_string(_currentTasks.load());
//    if(timeSinceLastGather>0.001) {
//      str=str+ " took too long: "+std::to_string(timeSinceLastGather);
//      stealing::StealingProfiler::getInstance().notifyLatePerformanceUpdate();
//    }

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
    _gatherTasksRequest = MPI_REQUEST_NULL;
  }

  if(completed_waiting_times) {
    //logInfo("progressGather","collected new waiting times snapshot, request"<< _gatherWaitingTimesRequest);
    std::copy(&_currentWaitingTimesReceiveBuffer[0], &_currentWaitingTimesReceiveBuffer[nnodes*nnodes], &_currentWaitingTimesSnapshot[0]);
    _gatherWaitingTimesRequest = MPI_REQUEST_NULL;
      
    //int k = 0;
    //for(int i=0; i<nnodes; i++) {
    //  for(int j=0; j<nnodes; j++) {
    //    logInfo("progressGather()","rank "<<i<<" waiting for "<<_currentWaitingTimesSnapshot[k+j]<<" for rank "<<j);
    //  }
    //  k+= nnodes;
    //}
  }

  if(_gatherWaitingTimesRequest==MPI_REQUEST_NULL && !isGloballyTerminated()) {
    //logInfo("progressGather","post gather waiting times");
    postGatherWaitingTimes();
  }

  if(_gatherTasksRequest==MPI_REQUEST_NULL && !isGloballyTerminated()) {
    //logInfo("progressGather","post gather tasks");
    //postGatherTasks();
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

  if(completed_blacklist) {
    std::copy(&_currentBlacklistReceiveBuffer[0], &_currentBlacklistReceiveBuffer[nnodes], &_currentBlacklistSnapshot[0]);
    _allreduceBlacklistRequest = MPI_REQUEST_NULL;
      
  }

  if(_allreduceBlacklistRequest==MPI_REQUEST_NULL && !isGloballyTerminated()) {
    postAllreduceBlacklist();
  }


  lock.free();
}

void exahype::stealing::PerformanceMonitor::postGatherTasks() {
  int myRank = tarch::parallel::Node::getInstance().getRank();

  _currentTasksSendBuffer = _isStarted ? _currentTasks.load() : TERMINATE_SIGNAL;
  
  //if(!_isStarted)
  //  logInfo("postGatherTasks","posted terminate signal");

  assert(_gatherTasksRequest==MPI_REQUEST_NULL);
  int err = MPI_Iallgather(&_currentTasksSendBuffer, 1, MPI_INTEGER, _currentTasksReceiveBuffer, 1, MPI_INTEGER, exahype::stealing::StealingManager::getInstance().getMPICommunicator(), &_gatherTasksRequest); //assert(err==MPI_SUCCESS);
}

void exahype::stealing::PerformanceMonitor::postGatherWaitingTimes() {
  int nnodes    = tarch::parallel::Node::getInstance().getNumberOfNodes();
  std::copy(&_currentWaitingTimes[0], &_currentWaitingTimes[nnodes], &_currentWaitingTimesSendBuffer[0]);


 // static int nposted = 0;
  //logInfo("postGatherWaitingTimes","posting gather for waiting times");
  //for(int i=0; i<nnodes; i++) {
  //  logInfo("postGatherWaitingTimes","_currentWaitingTimesSendBuffer["<<i<<"]: "<<_currentWaitingTimesSendBuffer[i]);
  //}
  assert(_gatherWaitingTimesRequest==MPI_REQUEST_NULL);
  int err = MPI_Iallgather(&_currentWaitingTimesSendBuffer[0], nnodes, MPI_INTEGER, &_currentWaitingTimesReceiveBuffer[0],
                   nnodes, MPI_INTEGER, exahype::stealing::StealingManager::getInstance().getMPICommunicator(),
                   &_gatherWaitingTimesRequest);// assert(err==MPI_SUCCESS);

  //nposted++;
  //logInfo("postGatherWaitingTimes", "nposted "<<nposted<< " request "<<_gatherWaitingTimesRequest);
}

void exahype::stealing::PerformanceMonitor::postAllreduceBlacklist() {
  int nnodes  = tarch::parallel::Node::getInstance().getNumberOfNodes();
  std::copy(&_currentBlacklist[0], &_currentBlacklist[nnodes], &_currentBlacklistSendBuffer[0]);

  assert(_allreduceBlacklistRequest==MPI_REQUEST_NULL);
  int err = MPI_Iallreduce(&_currentBlacklistSendBuffer[0], &_currentBlacklistReceiveBuffer[0], nnodes, MPI_DOUBLE, MPI_SUM,
             exahype::stealing::StealingManager::getInstance().getMPICommunicator(), &_allreduceBlacklistRequest);
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
