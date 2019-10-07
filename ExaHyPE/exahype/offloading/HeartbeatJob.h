/*
 * HeartbeatJob.h
 *
 *  Created on: 05.10.2019
 *      Author: ps659535
 */

#ifndef EXAHYPE_EXAHYPE_OFFLOADING_HEARTBEATJOB_H_
#define EXAHYPE_EXAHYPE_OFFLOADING_HEARTBEATJOB_H_

#include "tarch/multicore/Jobs.h"
#include "tarch/logging/Log.h"
#include <atomic>

namespace exahype {

namespace offloading{

class HeartbeatJob : public tarch::multicore::jobs::Job {
  public:
	static tarch::logging::Log _log;

	static void startHeartbeatJob();
	static void stopHeartbeatJob();
	bool run(bool runOnMasterThread);
  private:
	static HeartbeatJob* _singleton;

	std::atomic <bool> _terminateTrigger;
	std::atomic <bool> _hasTerminated;

	double _lastTimeStampTriggered;

	HeartbeatJob();
	virtual ~HeartbeatJob();
};

}

}
#endif /* EXAHYPE_EXAHYPE_OFFLOADING_HEARTBEATJOB_H_ */
