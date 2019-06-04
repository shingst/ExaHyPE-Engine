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

#if defined(SharedTBB) && defined(Parallel)

#include "exahype/offloading/StealingProfiler.h"

#include "tarch/multicore/Core.h"
#include "tarch/parallel/Node.h"

tarch::logging::Log exahype::stealing::StealingProfiler::_log( "exahype::stealing::StealingProfiler" );

exahype::stealing::StealingProfiler::StealingProfiler():
_executedTasksPhase(0),
_executedTasks(0),
_accUsefulCommunicationPhaseTime(0),
_accIdleCommunicationPhaseTime(0),
_accComputationPhaseTime(0),
_accWaitTasksTime(0),
_accWaitTasksPhaseTime(0),
_accUsefulCommunicationTime(0),
_accIdleCommunicationTime(0),
_accComputationTime(0),
_accOffloadPhaseTime(0),
_accOffloadTime(0),
_performanceUpdatesPhase(0),
_performanceUpdates(0),
_latePerformanceUpdatesPhase(0),
_latePerformanceUpdates(0),
_spawnedTasks(0),
_spawnedTasksPhase(0),
_thresholdFailsPhase(0),
_thresholdFails(0),
_stealingDecisionsPhase(0),
_stealingDecisions(0)
{
#if defined(StealingUseProfiler)
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  _offloadedTasksPerRankPhase = new std::atomic<int>[nnodes];
  _offloadedTasksPerRank = new std::atomic<int>[nnodes];
  _targetOffloadedTasksPerRankPhase = new std::atomic<int>[nnodes];
  _targetOffloadedTasksPerRank = new std::atomic<int>[nnodes];
  _receivedTasksPerRankPhase = new std::atomic<int>[nnodes];
  _receivedTasksPerRank = new std::atomic<int>[nnodes];
  std::fill(_offloadedTasksPerRankPhase, _offloadedTasksPerRankPhase+nnodes, 0);
  std::fill(_offloadedTasksPerRank, _offloadedTasksPerRank+nnodes, 0);
  std::fill(_targetOffloadedTasksPerRankPhase, _targetOffloadedTasksPerRankPhase+nnodes, 0);
  std::fill(_targetOffloadedTasksPerRank, _targetOffloadedTasksPerRank+nnodes, 0);
  std::fill(_receivedTasksPerRankPhase, _receivedTasksPerRankPhase+nnodes, 0);
  std::fill(_receivedTasksPerRank, _receivedTasksPerRank+nnodes, 0);
#endif
}

exahype::stealing::StealingProfiler::~StealingProfiler() {
#if defined(StealingUseProfiler)
  delete[] _offloadedTasksPerRank;
  delete[] _offloadedTasksPerRankPhase;
  delete[] _targetOffloadedTasksPerRank;
  delete[] _targetOffloadedTasksPerRankPhase;
  delete[] _receivedTasksPerRank;
  delete[] _receivedTasksPerRankPhase;
#endif
}

exahype::stealing::StealingProfiler& exahype::stealing::StealingProfiler::getInstance() {
  static exahype::stealing::StealingProfiler stealingProfiler;
  return stealingProfiler;
}

void exahype::stealing::StealingProfiler::beginPhase() {
#if defined(StealingUseProfiler)
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  std::fill(_offloadedTasksPerRankPhase, _offloadedTasksPerRankPhase+nnodes, 0);
  std::fill(_receivedTasksPerRankPhase, _receivedTasksPerRankPhase+nnodes, 0);
  std::fill(_targetOffloadedTasksPerRankPhase, _targetOffloadedTasksPerRankPhase+nnodes, 0);
  _executedTasksPhase=0;
  _spawnedTasksPhase=0;
  _thresholdFailsPhase=0;
  _performanceUpdatesPhase=0;
  _latePerformanceUpdatesPhase=0;
  _stealingDecisionsPhase=0;
  _accWaitTasksPhaseTime=0;
  _accUsefulCommunicationPhaseTime=0;
  _accIdleCommunicationPhaseTime=0;
  _accComputationPhaseTime=0;
  _accOffloadPhaseTime=0;
#endif
}

void exahype::stealing::StealingProfiler::endPhase() {
#if defined(StealingUseProfiler)
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  for (int i=0;i<nnodes;i++) {
    _offloadedTasksPerRank[i]+=_offloadedTasksPerRankPhase[i];
    _targetOffloadedTasksPerRank[i]+=_targetOffloadedTasksPerRankPhase[i];
    _receivedTasksPerRank[i]+=_receivedTasksPerRankPhase[i];
  }
  _executedTasks+=_executedTasksPhase;
  _spawnedTasks+= _spawnedTasksPhase;
  _thresholdFails+=_thresholdFailsPhase;
  _performanceUpdates+=_performanceUpdatesPhase;
  _latePerformanceUpdates+=_latePerformanceUpdatesPhase;
  _stealingDecisions+=_stealingDecisionsPhase;
  _accWaitTasksTime+=_accWaitTasksPhaseTime;
  _accUsefulCommunicationTime+=_accUsefulCommunicationPhaseTime;
  _accIdleCommunicationTime+=_accIdleCommunicationPhaseTime;
  _accComputationTime+=_accComputationPhaseTime;
  _accHandlingTime+=_accHandlingPhaseTime;
  _accOffloadTime+=_accOffloadPhaseTime;
#endif
}

void exahype::stealing::StealingProfiler::notifyOffloadedTask(int rank) {
#if defined(StealingUseProfiler)
  _offloadedTasksPerRankPhase[rank]++;
#endif
}

void exahype::stealing::StealingProfiler::notifyTargetOffloadedTask(int ntasks, int rank) {
#if defined(StealingUseProfiler)
  _targetOffloadedTasksPerRankPhase[rank]+=ntasks;
#endif
}

void exahype::stealing::StealingProfiler::notifyReceivedTask(int rank) {
#if defined(StealingUseProfiler)
  _receivedTasksPerRankPhase[rank]++;
#endif
}

void exahype::stealing::StealingProfiler::notifySpawnedTask() {
#if defined(StealingUseProfiler)
  _spawnedTasks++;
#endif
}

void exahype::stealing::StealingProfiler::notifyPerformanceUpdate() {
#if defined(StealingUseProfiler)
  _performanceUpdatesPhase++;
#endif
}

void exahype::stealing::StealingProfiler::notifyLatePerformanceUpdate() {
#if defined(StealingUseProfiler)
  _latePerformanceUpdatesPhase++;
#endif
}

void exahype::stealing::StealingProfiler::notifyStealingDecision() {
#if defined(StealingUseProfiler)
  _stealingDecisionsPhase++;
#endif
}

void exahype::stealing::StealingProfiler::notifyThresholdFail() {
#if defined(StealingUseProfiler)
  _thresholdFailsPhase++;
#endif
}

void exahype::stealing::StealingProfiler::beginComputation() {
#if defined(StealingUseProfiler)
#endif
}

void exahype::stealing::StealingProfiler::endComputation(double elapsed) {
#if defined(StealingUseProfiler)
  const unsigned long long elapsedTime = elapsed*1E6;
  _accComputationPhaseTime+=elapsedTime;
#endif
}

void exahype::stealing::StealingProfiler::beginCommunication() {
#if defined(StealingUseProfiler)
#endif
}

void exahype::stealing::StealingProfiler::endCommunication(bool successful, double elapsed) {
#if defined(StealingUseProfiler)
  const unsigned long long elapsedTime = elapsed *1E6;
  if(successful)
    _accUsefulCommunicationPhaseTime+=elapsedTime;
  else
    _accIdleCommunicationPhaseTime+=elapsedTime;
#endif
}

void exahype::stealing::StealingProfiler::beginHandling() {
#if defined(StealingUseProfiler)
#endif
}

void exahype::stealing::StealingProfiler::endHandling(double elapsed) {
#if defined(StealingUseProfiler)
  const unsigned long long elapsedTime = elapsed*1E6;
#endif
}

//void exahype::stealing::StealingProfiler::beginWaitForBackgroundTasks(exahype::solvers::Solver::JobType type) {
//#if defined(StealingUseProfiler)
//#endif
//}

//void exahype::stealing::StealingProfiler::endWaitForBackgroundTasks(exahype::solvers::Solver::JobType type, double elapsed) {
//#if defined(StealingUseProfiler)
//  unsigned long long elapsedTime;
//  if(type==exahype::solvers::Solver::JobType::EnclaveJob) {
//    elapsedTime = elapsed*1E6;
//    _accWaitEnclaveTasksPhaseTime+=elapsedTime;
//  }
//  else if(type==exahype::solvers::Solver::JobType::SkeletonJob) {
//    elapsedTime = elapsed*1E6;
//    _accWaitSkeletonTasksPhaseTime+=elapsedTime;
//  }

//#endif
//}

void exahype::stealing::StealingProfiler::beginWaitForTasks() {
#if defined(StealingUseProfiler)
#endif
}

void exahype::stealing::StealingProfiler::endWaitForTasks(double elapsed) {
#if defined(StealingUseProfiler)
  const unsigned long long elapsedTime = elapsed*1E6;
  _accWaitTasksPhaseTime+=elapsedTime;
#endif
}

void exahype::stealing::StealingProfiler::beginOffload() {
#if defined(StealingUseProfiler)
#endif
}

void exahype::stealing::StealingProfiler::endOffload(double elapsed) {
#if defined(StealingUseProfiler)
  const unsigned long long elapsedTime = elapsed*1E6;
  _accOffloadPhaseTime+=elapsedTime;
#endif
}

void exahype::stealing::StealingProfiler::printStatistics() {
#if defined(StealingUseProfiler)
  int rank = tarch::parallel::Node::getInstance().getRank();
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  std::string str="Stealing statistics for rank "+std::to_string(rank)+":\n";
  logInfo("exahype::stealing::StealingProfiler", str);
  str="  spawned tasks: "+std::to_string(_spawnedTasks)+"\n";
  logInfo("exahype::stealing::StealingProfiler", str);
  //str+="  executed tasks: "+std::to_string(_executedTasks)+"\n";
  str="  performance updates: "+std::to_string(_performanceUpdates)+"\n";
  logInfo("exahype::stealing::StealingProfiler", str);
  //str+="  late performance updates: "+std::to_string(_latePerformanceUpdates)+"\n";
  str="  stealing decisions: "+std::to_string(_stealingDecisions)+"\n";
  logInfo("exahype::stealing::StealingProfiler", str);
  str="  target offloaded tasks:\n";
  logInfo("exahype::stealing::StealingProfiler", str);
  int totalTargetOffloaded=0;
  for(int i=0;i<nnodes;i++) {
    str="    to rank: "+std::to_string(i)+" : "+std::to_string(_targetOffloadedTasksPerRank[i])+"\n";
    logInfo("exahype::stealing::StealingProfiler", str);
    totalTargetOffloaded+=_targetOffloadedTasksPerRank[i];
  }
  str="  total target offloaded tasks: "+std::to_string(totalTargetOffloaded)+"\n";
  logInfo("exahype::stealing::StealingProfiler", str);
  str="  offloaded tasks:\n";
  logInfo("exahype::stealing::StealingProfiler", str);
  int totalOffloaded=0;
  for(int i=0;i<nnodes;i++) {
    str="    to rank: "+std::to_string(i)+" : "+std::to_string(_offloadedTasksPerRank[i])+"\n";
    logInfo("exahype::stealing::StealingProfiler", str);
    totalOffloaded+=_offloadedTasksPerRank[i];
  }
  str="  total offloaded tasks: "+std::to_string(totalOffloaded)+"\n";
  logInfo("exahype::stealing::StealingProfiler", str);
  str="  offloadable tasks that failed threshold requirement: "+std::to_string(_thresholdFails)+"\n";
  logInfo("exahype::stealing::StealingProfiler", str);
  str="  received tasks:\n";
  logInfo("exahype::stealing::StealingProfiler", str);
  int totalReceived=0;
  for(int i=0;i<nnodes;i++) {
    str="    from rank: "+std::to_string(i)+" : "+std::to_string(_receivedTasksPerRank[i])+"\n";
    logInfo("exahype::stealing::StealingProfiler", str);
    totalReceived+=_receivedTasksPerRank[i];
  }
  str="  total received tasks: "+std::to_string(totalReceived)+"\n";
  logInfo("exahype::stealing::StealingProfiler", str);
  str="  total computation time: "+std::to_string(static_cast<double>(_accComputationTime/1E06))+"\n";
  logInfo("exahype::stealing::StealingProfiler", str);
  str="  total handling time: "+std::to_string(static_cast<double>(_accHandlingTime/1E06))+"\n";
  logInfo("exahype::stealing::StealingProfiler", str);
  str="  total useful communication time: "+std::to_string(static_cast<double>(_accUsefulCommunicationTime/1E06))+"\n";
  logInfo("exahype::stealing::StealingProfiler", str);
  str="  total idle communication time: "+std::to_string(static_cast<double>(_accIdleCommunicationTime/1E06))+"\n";
  logInfo("exahype::stealing::StealingProfiler", str);
  str="  total wait time for tasks: "+std::to_string(static_cast<double>(_accWaitTasksTime/1E06))+"\n";
  logInfo("exahype::stealing::StealingProfiler", str);

#endif
}

#endif
