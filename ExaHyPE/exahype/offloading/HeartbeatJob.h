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

#if !defined(EXAHYPE_OFFLOADING_HEARTBEATJOB_H_) && defined(SharedTBB)  && defined(Parallel)
#define EXAHYPE_OFFLOADING_HEARTBEATJOB_H_

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
#endif /* EXAHYPE_OFFLOADING_HEARTBEATJOB_H_ */
