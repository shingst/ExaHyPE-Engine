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

#include "TimeStampAndDubiosityTeamHistory.h"

#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/multicore/Lock.h"
#include "tarch/parallel/Node.h"

#include "tarch/logging/Log.h"

#include <cassert>
#include <string>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <limits>
#include "exahype/reactive/ReactiveContext.h"

tarch::logging::Log exahype::reactive::TimeStampAndDubiosityTeamHistory::_log("exahype::reactive::TimeStampAndLimiterHistory");


exahype::reactive::TimeStampAndDubiosityTeamHistory::TimeStampAndDubiosityTeamHistory() :
  _lastConsistentTimeStepPtr(-1)
{
#if defined(Parallel)
  _timestamps = new std::vector<double> [exahype::reactive::ReactiveContext::getInstance().getTMPINumTeams()];
  _timestepSizes = new std::vector<double> [exahype::reactive::ReactiveContext::getInstance().getTMPINumTeams()];
  _dubiosityStatuses = new std::vector<bool> [exahype::reactive::ReactiveContext::getInstance().getTMPINumTeams()];
#else
  _timestamps = new std::vector<double> [1];
  _timestepSizes = new std::vector<double> [1];
  _dubiosityStatuses = new std::vector<bool> [1];
#endif
}

exahype::reactive::TimeStampAndDubiosityTeamHistory::~TimeStampAndDubiosityTeamHistory() {
  delete[] _dubiosityStatuses;
  delete[] _timestepSizes;
  delete[] _timestamps;
}

exahype::reactive::TimeStampAndDubiosityTeamHistory& exahype::reactive::TimeStampAndDubiosityTeamHistory::getInstance() {
  static TimeStampAndDubiosityTeamHistory history;
  return history;
}

void exahype::reactive::TimeStampAndDubiosityTeamHistory::forwardLastConsistentTimeStepPtr() {
  tarch::multicore::Lock lock(_semaphore, true);
#if defined(Parallel)
  int myTeam = exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber();
  int otherTeam = (myTeam + 1) % 2; //todo: support if more than 2 teams are used

  size_t maxIdx = std::numeric_limits<int>::max();

  for(size_t i=0; i<exahype::reactive::ReactiveContext::getInstance().getTMPINumTeams(); i++) {
    maxIdx = std::min(_timestamps[i].size(), maxIdx);
  }

  unsigned int i = std::max(0, _lastConsistentTimeStepPtr);
  while(i<maxIdx && _timestamps[otherTeam][i]==_timestamps[myTeam][i]
                 && _dubiosityStatuses[otherTeam][i]== _dubiosityStatuses[myTeam][i]  //we may only want to consider a time step consistent if limiter status is both zero (for linear applications)
                 && _timestepSizes[otherTeam][i]==_timestepSizes[myTeam][i]) {
    i++;
  }
  _lastConsistentTimeStepPtr = --i;
  assert(_lastConsistentTimeStepPtr>=-1);
#else
  _lastConsistentTimeStepPtr = _timestamps[0].size();
#endif
  lock.free();
}

void exahype::reactive::TimeStampAndDubiosityTeamHistory::trackTimeStepAndDubiosity(int team, double timeStamp, double timeStepSize, bool triggerActive, double estimatedTimeStepSize) {
  tarch::multicore::Lock lock(_semaphore, true);
  logDebug("trackTimeStepAndLimiterActive", "Tracking time stamp="<<timeStamp<<" for team "<<team);

  int myTeam;
#if defined(Parallel)
  myTeam = exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber();
#else
  myTeam = 0;
#endif
  int size = _timestamps[team].size();
  if(size>0 && _timestamps[team][size-1]==timeStamp) {
    ////only update dubiosity
    _dubiosityStatuses[team][size-1] = _dubiosityStatuses[team][size-1] || triggerActive;
    //assert(_timestepSizes[team][size-1]==std::min(_timestepSizes[team][size-1], timeStepSize));
    _timestepSizes[team][size-1] = std::min(_timestepSizes[team][size-1], timeStepSize);
    if(team==myTeam) {
      _estimatedTimestepSizes[size-1]= std::min(_estimatedTimestepSizes[size-1], estimatedTimeStepSize);
    }
  }
  else if (size>0 && _timestamps[team][size-1]>timeStamp) {
    int idx = size-2;
    while(idx>=0 && _timestamps[team][idx]!=timeStamp) {
      idx--;
    }
    if(idx<0) {
      logError("trackTimeStepAndLimiterActive", "Have received timestamp from another team which seems to be old but has not been observed. This is expected if the other team has done a rollback. Appending new time stamp..");
      _dubiosityStatuses[team].push_back(triggerActive);
      _timestamps[team].push_back(timeStamp);
      _timestepSizes[team].push_back(timeStepSize);
    }
    else {
      _dubiosityStatuses[team][idx] = _dubiosityStatuses[team][idx] || triggerActive;
      _timestepSizes[team][idx] = std::min(_timestepSizes[team][idx], timeStepSize);
      if(team==myTeam) {
        _estimatedTimestepSizes[idx]= std::min(_estimatedTimestepSizes[idx], estimatedTimeStepSize);
      }
    }
  }
  else if (size==0 || _timestamps[team][size-1]<timeStamp){
    _timestamps[team].push_back(timeStamp);
    _timestepSizes[team].push_back(timeStepSize);
    _dubiosityStatuses[team].push_back(triggerActive);
    if(team==myTeam) {
      _estimatedTimestepSizes.push_back(estimatedTimeStepSize);
    }
  }
  //printHistory();
  lock.free();
}

bool exahype::reactive::TimeStampAndDubiosityTeamHistory::checkConsistency() {
#if defined(Parallel)
  int myTeam = exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber();
  int otherTeam = (myTeam + 1) % 2; //todo: support if more than 2 teams are used

  forwardLastConsistentTimeStepPtr();

  tarch::multicore::Lock lock(_semaphore, true);
  size_t maxIdx = std::numeric_limits<int>::max();

  for(unsigned int i=0; i<exahype::reactive::ReactiveContext::getInstance().getTMPINumTeams(); i++) {
    maxIdx = std::min(_timestamps[i].size(), maxIdx);
  }

  logDebug("checkConsistency","checking arrays from " <<_lastConsistentTimeStepPtr<<" until "<<maxIdx);

  bool consistentTimeStamps = true;
  bool consistentTimeStepSizes = true;
  bool consistentDubiosityStatuses = true;

  for(size_t i=_lastConsistentTimeStepPtr; i<maxIdx; i++) {
    std::vector<double> tmp_timestamps;
    std::vector<double> tmp_timestepsizes;
    std::vector<bool> tmp_statuses;
    for(unsigned int t=0; t<exahype::reactive::ReactiveContext::getInstance().getTMPINumTeams(); t++) {
      tmp_timestamps.push_back(_timestamps[t][i]);
      tmp_statuses.push_back(_dubiosityStatuses[t][i]);
      tmp_timestepsizes.push_back(_timestepSizes[t][i]);
    }
    consistentTimeStamps = consistentTimeStamps && std::all_of(tmp_timestamps.begin(), tmp_timestamps.end(), [tmp_timestamps](double x){ return x==tmp_timestamps[0]; });
    consistentDubiosityStatuses = consistentDubiosityStatuses && std::all_of(tmp_statuses.begin(), tmp_statuses.end(), [tmp_statuses](double x){ return x==tmp_statuses[0]; });
    consistentTimeStepSizes = consistentTimeStepSizes && std::all_of(tmp_timestepsizes.begin(), tmp_timestepsizes.end(), [tmp_timestepsizes](double x){ return x==tmp_timestepsizes[0]; });
  }

  if((!consistentTimeStamps || !consistentDubiosityStatuses || !consistentTimeStepSizes)) {
    //logError("checkConsistency"," Time stamps or limiter statuses are diverged between teams! Consistent stamps = "<<consistentTimeStamps
    //    <<" consistent limiter statuses = "<<consistentLimiterStatuses<<" consistent time step sizes="<<consistentTimeStepSizes);
    //logError("checkConsistency","team="<<exahype::reactive::OffloadingContext::getInstance().getTMPITeamNumber()<<": Time stamps or trigger statuses are inconsistent between teams! There must have been a soft error on at least one team.");
    if((long int)_dubiosityStatuses[myTeam].size()> _lastConsistentTimeStepPtr && (long int) _dubiosityStatuses[otherTeam].size()>_lastConsistentTimeStepPtr) {
      if(_dubiosityStatuses[myTeam][_lastConsistentTimeStepPtr+1]==1
         &&_dubiosityStatuses[otherTeam][_lastConsistentTimeStepPtr+1]==0) {
        logError("checkConsistency","team="<<myTeam<<" should be the faulty one, as a dubious result was detected there.");
      }
      else if(_dubiosityStatuses[otherTeam][_lastConsistentTimeStepPtr+1]==1
          &&_dubiosityStatuses[myTeam][_lastConsistentTimeStepPtr+1]==0) {
         logError("checkConsistency","team="<<otherTeam<<" should be the faulty one, as a dubious result was detected there.");
      }
    }

    //printHistory();
  }
  lock.free();

  return consistentTimeStamps && consistentDubiosityStatuses && consistentTimeStepSizes;
#else
  return true;
#endif 
}

bool exahype::reactive::TimeStampAndDubiosityTeamHistory::otherTeamHasTimeStepData(double timeStamp, double timeStep) {
#if defined(Parallel)
  int myTeam = exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber();
  int otherTeam = (myTeam + 1) % 2; //todo: support if more than 2 teams are used

  forwardLastConsistentTimeStepPtr();

  int startIndex = std::max(0, _lastConsistentTimeStepPtr);

  for(size_t i=startIndex; i<_timestamps[otherTeam].size(); i++) {
    if(_timestamps[otherTeam][i]==timeStamp && _timestepSizes[otherTeam][i]==timeStep)
      return true;
  }
  return false;
#else
  return false;
#endif
}

bool exahype::reactive::TimeStampAndDubiosityTeamHistory::otherTeamHasTimeStamp(double timeStamp) {
#if defined(Parallel)
  int myTeam = exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber();
  int otherTeam = (myTeam + 1) % 2; //todo: support if more than 2 teams are used

  forwardLastConsistentTimeStepPtr();

  int startIndex = std::max(0, _lastConsistentTimeStepPtr);

  for(size_t i=startIndex; i<_timestamps[otherTeam].size(); i++) {
    if(_timestamps[otherTeam][i]==timeStamp)
      return true;
  }
  return false;
#else
  return false;
#endif
}

bool exahype::reactive::TimeStampAndDubiosityTeamHistory::otherTeamHasLargerTimeStamp(double timeStamp) {
#if defined(Parallel)
  int myTeam = exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber();
  int otherTeam = (myTeam + 1) % 2; //todo: support if more than 2 teams are used

  forwardLastConsistentTimeStepPtr();

  int startIndex = std::max(0, _lastConsistentTimeStepPtr);

  for(unsigned int i=startIndex; i<_timestamps[otherTeam].size(); i++) {
    if(_timestamps[otherTeam][i]>timeStamp)
      return true;
  }

  return false;
#else
  return false;
#endif
}

bool exahype::reactive::TimeStampAndDubiosityTeamHistory::otherTeamHasLargerTimeStepSizeForStamp(double timeStamp, double timeStep) {
#if defined(Parallel)
  int myTeam = exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber();
  int otherTeam = (myTeam + 1) % 2; //todo: support if more than 2 teams are used

  forwardLastConsistentTimeStepPtr();

  int startIndex = std::max(0, _lastConsistentTimeStepPtr);

  for(unsigned int i=startIndex; i<_timestamps[otherTeam].size(); i++) {
    if(_timestamps[otherTeam][i]==timeStamp && _timestepSizes[otherTeam][i]>timeStep)
      return true;
  }

  return false;
#else
  return false;
#endif
}

void exahype::reactive::TimeStampAndDubiosityTeamHistory::getLastConsistentTimeStepData(double& timestamp, double& timestepSize, double& estimated) {
#if defined(Parallel)
    int myTeam = exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber();
#else
    int myTeam = 0;
#endif

  forwardLastConsistentTimeStepPtr();

  tarch::multicore::Lock lock(_semaphore, true);
  if(_lastConsistentTimeStepPtr>=0) {
    timestamp = _timestamps[myTeam][_lastConsistentTimeStepPtr];
    timestepSize = _timestepSizes[myTeam][_lastConsistentTimeStepPtr];
    estimated = _estimatedTimestepSizes[_lastConsistentTimeStepPtr];
  }
  else {
    timestamp = 0;
    timestepSize = -1;
    estimated = -1;
  }
  lock.free();
}

void exahype::reactive::TimeStampAndDubiosityTeamHistory::resetMyTeamHistoryToLastConsistentTimeStep() {
  tarch::multicore::Lock lock(_semaphore, true);
#if defined(Parallel)
  int myTeam = exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber();
#else
  int myTeam = 0;
#endif

  if(_lastConsistentTimeStepPtr>=0 && _lastConsistentTimeStepPtr<_timestamps[myTeam].size()) {
    _timestamps[myTeam].erase(_timestamps[myTeam].begin()+_lastConsistentTimeStepPtr+1, _timestamps[myTeam].end());
    _timestepSizes[myTeam].erase(_timestepSizes[myTeam].begin()+_lastConsistentTimeStepPtr+1, _timestepSizes[myTeam].end());
    _estimatedTimestepSizes.erase(_estimatedTimestepSizes.begin()+_lastConsistentTimeStepPtr+1, _estimatedTimestepSizes.end());
    _dubiosityStatuses[myTeam].erase(_dubiosityStatuses[myTeam].begin()+_lastConsistentTimeStepPtr+1, _dubiosityStatuses[myTeam].end());
  }

  logInfo("resetMyTeamToLastConsistentTimeStep", "reset (local) team history to last consistent time step. Printing history next...");
  printHistory();
  lock.free();
}


void exahype::reactive::TimeStampAndDubiosityTeamHistory::printHistory() const {
#if defined(Parallel)
  int numTeams = exahype::reactive::ReactiveContext::getInstance().getTMPINumTeams();
  int myTeam = exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber();
#else
  int numTeams = 1;
  int myTeam = 0;
#endif  

  for(unsigned int i=0; i<numTeams; i++) {
    std::ostringstream stream;
    stream <<" Timestamps: "<<std::setprecision(30);
    for(unsigned int j=0; j<_timestamps[i].size(); j++) {
      stream << _timestamps[i][j] << " , ";
    }
    logInfo("printHistory","TimeStampAndLimiterHistory for team "<<i<<" on team "<<myTeam<<":"<<stream.str());

    std::ostringstream stream2;
    stream2 << " TimeStepSizes: "<< std::setprecision(30);
    for(unsigned int j=0; j<_timestepSizes[i].size(); j++) {
      stream2 << _timestepSizes[i][j] << " , ";
    }
    logInfo("printHistory", "TimeStampAndLimiterHistory for team "<<i<<" on team "<<myTeam<<":"<<stream2.str());

    std::ostringstream stream3;
    stream3 << " Trigger statuses ";
    for(unsigned int j=0; j<_dubiosityStatuses[i].size(); j++) {
      stream3 << std::to_string(_dubiosityStatuses[i][j]) + " , ";
    }
    logInfo("printHistory","TimeStampAndLimiterHistory for team "<<i<<" on team "<<myTeam<<":"<<stream3.str());
  }
}
