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

#include "../reactive/JobTableStatistics.h"

#include "OffloadingContext.h"

namespace exahype {
namespace reactive {


tarch::logging::Log JobTableStatistics::_log( "exahype::reactive::JobTableStatistics" );


JobTableStatistics::JobTableStatistics() :
  _spawnedTasks(0),
  _executedTasks(0),
  _savedTasks(0),
  _receivedTasks(0),
  _sentTasks(0),
  _sentKeys(0),
  _receivedKeys(0),
  _declinedTasks(0),
  _lateTasks(0),
  _recomputedTasks(0),
  _doubleCheckedTasks(0)
{
  // TODO Auto-generated constructor stub

}

JobTableStatistics::~JobTableStatistics() {
  // TODO Auto-generated destructor stub
}

JobTableStatistics& JobTableStatistics::getInstance() {
  static JobTableStatistics replicationStats;
  return replicationStats;
}

void JobTableStatistics::notifyLateTask() {
  _lateTasks++;
}

void JobTableStatistics::notifyDeclinedTask() {
  _declinedTasks++;
}

void JobTableStatistics::notifyReceivedTask(){
  _receivedTasks++;
}

void JobTableStatistics::notifySentTask(){
  _sentTasks++;
}

void JobTableStatistics::notifySavedTask(){
  _savedTasks++;
}

void JobTableStatistics::notifySpawnedTask(){
    _spawnedTasks++;
}

void JobTableStatistics::notifyExecutedTask(){
    _executedTasks++;
}

void JobTableStatistics::notifySentKey() {
  _sentKeys++;
}

void JobTableStatistics::notifyReceivedKey() {
  _receivedKeys++;
}

void JobTableStatistics::notifyRecomputedTask() {
  _recomputedTasks++;
}

void JobTableStatistics::notifyDoubleCheckedTask() {
  _doubleCheckedTasks++;
}

void JobTableStatistics::notifyDetectedError() {
  _softErrorsDetected++;
}

void JobTableStatistics::notifyHealedTask() {
  _healedTasks++;
}


void JobTableStatistics::printStatistics() {
     int team = exahype::reactive::OffloadingContext::getInstance().getTMPIInterTeamRank();
     logInfo("printStatistics", " team "<<team
           <<" spawned tasks = "<<_spawnedTasks
           <<" executed tasks = "<<_executedTasks
           <<" double checked tasks = "<<_doubleCheckedTasks
           <<" detected soft errors = "<<_softErrorsDetected
           <<" healed tasks = "<<_healedTasks
           <<" saved tasks =  "<<_savedTasks
           <<" sent tasks = "<<_sentTasks
           <<" received tasks = "<<_receivedTasks
           <<" received keys= "<<_receivedKeys
           <<" sent keys= "<<_sentKeys
           <<" declined tasks = "<<_declinedTasks
           <<" late tasks = "<<_lateTasks
           <<" recomputed tasks = "<<_recomputedTasks);
}

} /* namespace offloading */
} /* namespace exahype */
