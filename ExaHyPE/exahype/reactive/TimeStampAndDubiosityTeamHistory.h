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

#ifndef EXAHYPE_EXAHYPE_REACTIVE_TIMESTAMPANDDUBIOSITYTEAMHISTORY_H_
#define EXAHYPE_EXAHYPE_REACTIVE_TIMESTAMPANDDUBIOSITYTEAMHISTORY_H_

#include "tarch/logging/Log.h"
#include "tarch/multicore/BooleanSemaphore.h"

#include <vector>
#include <limits>

namespace exahype {
 namespace reactive {
   class TimeStampAndDubiosityTeamHistory;
 }
}

/**
 * The time stamp and trigger history keeps track of the history of time stamps and time step sizes
 * and the trigger (indicating whether there was a potential soft error) history.
 * It makes only sense if multiple teams
 * are used.
 * The history of time stamps, time step sizes and trigger values (per time step) can be checked for consistency between the teams in
 * order to find potential soft errors.
 */
class exahype::reactive::TimeStampAndDubiosityTeamHistory {

  private:
    /**
     * Logging device for the trace macros.
     */
    static tarch::logging::Log _log;

    //should probably be a set or similar
    std::vector<double> *_timestamps;
    std::vector<double> *_timestepSizes;
    std::vector<double> _estimatedTimestepSizes;
    std::vector<bool> *_dubiosityStatuses;

    int _lastConsistentTimeStepPtr;

    tarch::multicore::BooleanSemaphore _semaphore;

    void forwardLastConsistentTimeStepPtr();

  public:
    TimeStampAndDubiosityTeamHistory();
    ~TimeStampAndDubiosityTeamHistory();

    static TimeStampAndDubiosityTeamHistory& getInstance();

    void trackTimeStepAndDubiosity(int team, double timeStamp, double timeStepSize, bool limiterActive, double estimated = std::numeric_limits<double>::infinity());
    void trackEstimatedTimeStepSizeLocally(double timeStamp, double estTimeStepSize);

    bool checkConsistency();
    void getLastConsistentTimeStepData(double& timestamp, double& timestepSize, double& estimatedTimeStepSize);
    void resetMyTeamHistoryToLastConsistentTimeStep();

    bool otherTeamHasTimeStepData(double timeStamp, double timeStepSize);
    bool otherTeamHasTimeStamp(double timeStamp);
    bool otherTeamHasLargerTimeStamp(double timeStamp);
    bool otherTeamHasLargerTimeStepSizeForStamp(double timeStamp, double timeStep);

    void printHistory() const;
};

#endif /* EXAHYPE_EXAHYPE_REACTIVE_TIMESTAMPANDDUBIOSITYTEAMHISTORY_H_ */
