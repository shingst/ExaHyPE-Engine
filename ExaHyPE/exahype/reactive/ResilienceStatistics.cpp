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

#include "exahype/reactive/ResilienceStatistics.h"

#include "exahype/reactive/ReactiveContext.h"

namespace exahype {
namespace reactive {


tarch::logging::Log ResilienceStatistics::_log( "exahype::reactive::JobTableStatistics" );


ResilienceStatistics::ResilienceStatistics() :
  _spawnedTasks(0),
  _executedTasks(0),
  _savedTasks(0),
  _receivedTasks(0),
  _sentTasks(0),
  _declinedTasks(0),
  _lateTasks(0),
  _recomputedTasks(0),
  _doubleCheckedTasks(0),
  _softErrorsDetected(0),
  _healedTasks(0),
  _softErrorsInjected(0),
  _limitedTasks(0)
{
  // TODO Auto-generated constructor stub

}

ResilienceStatistics::~ResilienceStatistics() {
  // TODO Auto-generated destructor stub
}

ResilienceStatistics& ResilienceStatistics::getInstance() {
  static ResilienceStatistics replicationStats;
  return replicationStats;
}

void ResilienceStatistics::notifyLateTask() {
  _lateTasks++;
}

void ResilienceStatistics::notifyDeclinedTask() {
  _declinedTasks++;
}

void ResilienceStatistics::notifyReceivedTask(){
  _receivedTasks++;
}

void ResilienceStatistics::notifySentTask(){
  _sentTasks++;
}

void ResilienceStatistics::notifySavedTask(){
  _savedTasks++;
}

void ResilienceStatistics::notifySpawnedTask(){
    _spawnedTasks++;
}

void ResilienceStatistics::notifyExecutedTask(){
    _executedTasks++;
}

void ResilienceStatistics::notifyRecomputedTask() {
  _recomputedTasks++;
}

void ResilienceStatistics::notifyDoubleCheckedTask() {
  _doubleCheckedTasks++;
}

void ResilienceStatistics::notifyDetectedError() {
  _softErrorsDetected++;
}

void ResilienceStatistics::notifyHealedTask() {
  _healedTasks++;
}

void ResilienceStatistics::notifyInjectedError() {
  _softErrorsInjected++;
}

void ResilienceStatistics::notifyLimitedTask() {
  _limitedTasks++;
}

void ResilienceStatistics::printStatistics() {
     int team = exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber();
     logInfo("printStatistics", " team "<<team
           <<" spawned tasks = "<<_spawnedTasks
           <<" executed tasks = "<<_executedTasks
           <<" double checked tasks = "<<_doubleCheckedTasks
           <<" soft errors injected = "<<_softErrorsInjected
           <<" detected soft errors = "<<_softErrorsDetected
           <<" limited tasks = "<<_limitedTasks
           <<" healed tasks = "<<_healedTasks
           <<" saved tasks =  "<<_savedTasks
           <<" sent tasks = "<<_sentTasks
           <<" received tasks = "<<_receivedTasks
           <<" declined tasks = "<<_declinedTasks
           <<" late tasks = "<<_lateTasks);
}

} /* namespace offloading */
} /* namespace exahype */
