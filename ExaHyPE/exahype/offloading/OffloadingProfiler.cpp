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

#include "exahype/offloading/OffloadingProfiler.h"

#include "tarch/multicore/Core.h"
#include "tarch/parallel/Node.h"

tarch::logging::Log exahype::offloading::OffloadingProfiler::_log( "exahype::offloading::OffloadingProfiler" );

exahype::offloading::OffloadingProfiler::OffloadingProfiler():
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
_offloadingDecisionsPhase(0),
_offloadingDecisions(0)
{
#if defined(OffloadingUseProfiler)
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

exahype::offloading::OffloadingProfiler::~OffloadingProfiler() {
#if defined(OffloadingUseProfiler)
  delete[] _offloadedTasksPerRank;
  delete[] _offloadedTasksPerRankPhase;
  delete[] _targetOffloadedTasksPerRank;
  delete[] _targetOffloadedTasksPerRankPhase;
  delete[] _receivedTasksPerRank;
  delete[] _receivedTasksPerRankPhase;
#endif
}

exahype::offloading::OffloadingProfiler& exahype::offloading::OffloadingProfiler::getInstance() {
  static exahype::offloading::OffloadingProfiler offloadingProfiler;
  return offloadingProfiler;
}

void exahype::offloading::OffloadingProfiler::beginPhase() {
#if defined(OffloadingUseProfiler)
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  std::fill(_offloadedTasksPerRankPhase, _offloadedTasksPerRankPhase+nnodes, 0);
  std::fill(_receivedTasksPerRankPhase, _receivedTasksPerRankPhase+nnodes, 0);
  std::fill(_targetOffloadedTasksPerRankPhase, _targetOffloadedTasksPerRankPhase+nnodes, 0);
  _executedTasksPhase=0;
  _spawnedTasksPhase=0;
  _thresholdFailsPhase=0;
  _performanceUpdatesPhase=0;
  _latePerformanceUpdatesPhase=0;
  _offloadingDecisionsPhase=0;
  _accWaitTasksPhaseTime=0;
  _accUsefulCommunicationPhaseTime=0;
  _accIdleCommunicationPhaseTime=0;
  _accComputationPhaseTime=0;
  _accOffloadPhaseTime=0;
#endif
}

void exahype::offloading::OffloadingProfiler::endPhase() {
#if defined(OffloadingUseProfiler)
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
  _offloadingDecisions+=_offloadingDecisionsPhase;
  _accWaitTasksTime+=_accWaitTasksPhaseTime;
  _accUsefulCommunicationTime+=_accUsefulCommunicationPhaseTime;
  _accIdleCommunicationTime+=_accIdleCommunicationPhaseTime;
  _accComputationTime+=_accComputationPhaseTime;
  _accHandlingTime+=_accHandlingPhaseTime;
  _accOffloadTime+=_accOffloadPhaseTime;
#endif
}

void exahype::offloading::OffloadingProfiler::notifyOffloadedTask(int rank) {
#if defined(OffloadingUseProfiler)
  _offloadedTasksPerRankPhase[rank]++;
#endif
}

void exahype::offloading::OffloadingProfiler::notifyTargetOffloadedTask(int ntasks, int rank) {
#if defined(OffloadingUseProfiler)
  _targetOffloadedTasksPerRankPhase[rank]+=ntasks;
#endif
}

void exahype::offloading::OffloadingProfiler::notifyReceivedTask(int rank) {
#if defined(OffloadingUseProfiler)
  _receivedTasksPerRankPhase[rank]++;
#endif
}

void exahype::offloading::OffloadingProfiler::notifySpawnedTask() {
#if defined(OffloadingUseProfiler)
  _spawnedTasks++;
#endif
}

void exahype::offloading::OffloadingProfiler::notifyPerformanceUpdate() {
#if defined(OffloadingUseProfiler)
  _performanceUpdatesPhase++;
#endif
}

void exahype::offloading::OffloadingProfiler::notifyLatePerformanceUpdate() {
#if defined(OffloadingUseProfiler)
  _latePerformanceUpdatesPhase++;
#endif
}

void exahype::offloading::OffloadingProfiler::notifyOffloadingDecision() {
#if defined(OffloadingUseProfiler)
  _offloadingDecisionsPhase++;
#endif
}

void exahype::offloading::OffloadingProfiler::notifyThresholdFail() {
#if defined(OffloadingUseProfiler)
  _thresholdFailsPhase++;
#endif
}

void exahype::offloading::OffloadingProfiler::beginComputation() {
#if defined(OffloadingUseProfiler)
#endif
}

void exahype::offloading::OffloadingProfiler::endComputation(double elapsed) {
#if defined(OffloadingUseProfiler)
  const unsigned long long elapsedTime = elapsed*1E6;
  _accComputationPhaseTime+=elapsedTime;
#endif
}

void exahype::offloading::OffloadingProfiler::beginCommunication() {
#if defined(OffloadingUseProfiler)
#endif
}

void exahype::offloading::OffloadingProfiler::endCommunication(bool successful, double elapsed) {
#if defined(OffloadingUseProfiler)
  const unsigned long long elapsedTime = elapsed *1E6;
  if(successful)
    _accUsefulCommunicationPhaseTime+=elapsedTime;
  else
    _accIdleCommunicationPhaseTime+=elapsedTime;
#endif
}

void exahype::offloading::OffloadingProfiler::beginHandling() {
#if defined(OffloadingUseProfiler)
#endif
}

void exahype::offloading::OffloadingProfiler::endHandling(double elapsed) {
#if defined(OffloadingUseProfiler)
  const unsigned long long elapsedTime = elapsed*1E6;
#endif
}

//void exahype::offloading::OffloadingProfiler::beginWaitForBackgroundTasks(exahype::solvers::Solver::JobType type) {
//#if defined(OffloadingUseProfiler)
//#endif
//}

//void exahype::offloading::OffloadingProfiler::endWaitForBackgroundTasks(exahype::solvers::Solver::JobType type, double elapsed) {
//#if defined(OffloadingUseProfiler)
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

void exahype::offloading::OffloadingProfiler::beginWaitForTasks() {
#if defined(OffloadingUseProfiler)
#endif
}

void exahype::offloading::OffloadingProfiler::endWaitForTasks(double elapsed) {
#if defined(OffloadingUseProfiler)
  const unsigned long long elapsedTime = elapsed*1E6;
  _accWaitTasksPhaseTime+=elapsedTime;
#endif
}

void exahype::offloading::OffloadingProfiler::beginOffload() {
#if defined(OffloadingUseProfiler)
#endif
}

void exahype::offloading::OffloadingProfiler::endOffload(double elapsed) {
#if defined(OffloadingUseProfiler)
  const unsigned long long elapsedTime = elapsed*1E6;
  _accOffloadPhaseTime+=elapsedTime;
#endif
}

void exahype::offloading::OffloadingProfiler::printStatistics() {
#if defined(OffloadingUseProfiler)
  int rank = tarch::parallel::Node::getInstance().getRank();
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  std::string str="Offloading statistics for rank "+std::to_string(rank)+":\n";
  logInfo("exahype::offloading::OffloadingProfiler", str);
  str="  spawned tasks: "+std::to_string(_spawnedTasks)+"\n";
  logInfo("exahype::offloading::OffloadingProfiler", str);
  //str+="  executed tasks: "+std::to_string(_executedTasks)+"\n";
  str="  performance updates: "+std::to_string(_performanceUpdates)+"\n";
  logInfo("exahype::offloading::OffloadingProfiler", str);
  //str+="  late performance updates: "+std::to_string(_latePerformanceUpdates)+"\n";
  str="  offloading decisions: "+std::to_string(_offloadingDecisions)+"\n";
  logInfo("exahype::offloading::OffloadingProfiler", str);
  str="  target offloaded tasks:\n";
  logInfo("exahype::offloading::OffloadingProfiler", str);
  int totalTargetOffloaded=0;
  for(int i=0;i<nnodes;i++) {
    str="    to rank: "+std::to_string(i)+" : "+std::to_string(_targetOffloadedTasksPerRank[i])+"\n";
    logInfo("exahype::offloading::OffloadingProfiler", str);
    totalTargetOffloaded+=_targetOffloadedTasksPerRank[i];
  }
  str="  total target offloaded tasks: "+std::to_string(totalTargetOffloaded)+"\n";
  logInfo("exahype::offloading::OffloadingProfiler", str);
  str="  offloaded tasks:\n";
  logInfo("exahype::offloading::OffloadingProfiler", str);
  int totalOffloaded=0;
  for(int i=0;i<nnodes;i++) {
    str="    to rank: "+std::to_string(i)+" : "+std::to_string(_offloadedTasksPerRank[i])+"\n";
    logInfo("exahype::offloading::OffloadingProfiler", str);
    totalOffloaded+=_offloadedTasksPerRank[i];
  }
  str="  total offloaded tasks: "+std::to_string(totalOffloaded)+"\n";
  logInfo("exahype::offloading::OffloadingProfiler", str);
  str="  offloadable tasks that failed threshold requirement: "+std::to_string(_thresholdFails)+"\n";
  logInfo("exahype::offloading::OffloadingProfiler", str);
  str="  received tasks:\n";
  logInfo("exahype::offloading::OffloadingProfiler", str);
  int totalReceived=0;
  for(int i=0;i<nnodes;i++) {
    str="    from rank: "+std::to_string(i)+" : "+std::to_string(_receivedTasksPerRank[i])+"\n";
    logInfo("exahype::offloading::OffloadingProfiler", str);
    totalReceived+=_receivedTasksPerRank[i];
  }
  str="  total received tasks: "+std::to_string(totalReceived)+"\n";
  logInfo("exahype::offloading::OffloadingProfiler", str);
  str="  total computation time: "+std::to_string(static_cast<double>(_accComputationTime/1E06))+"\n";
  logInfo("exahype::offloading::OffloadingProfiler", str);
  str="  total handling time: "+std::to_string(static_cast<double>(_accHandlingTime/1E06))+"\n";
  logInfo("exahype::offloading::OffloadingProfiler", str);
  str="  total useful communication time: "+std::to_string(static_cast<double>(_accUsefulCommunicationTime/1E06))+"\n";
  logInfo("exahype::offloading::OffloadingProfiler", str);
  str="  total idle communication time: "+std::to_string(static_cast<double>(_accIdleCommunicationTime/1E06))+"\n";
  logInfo("exahype::offloading::OffloadingProfiler", str);
  str="  total wait time for tasks: "+std::to_string(static_cast<double>(_accWaitTasksTime/1E06))+"\n";
  logInfo("exahype::offloading::OffloadingProfiler", str);

#endif
}

#endif
