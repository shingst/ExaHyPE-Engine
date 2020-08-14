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

#ifndef EXAHYPE_EXAHYPE_OFFLOADING_STPSTATSTRACER_H_
#define EXAHYPE_EXAHYPE_OFFLOADING_STPSTATSTRACER_H_

#include <string>

namespace exahype {
namespace offloading {

enum class STPType {ADERDGPrediction, ADERDGOwnMigratable, ADERDGRemoteMigratable, LimitingFusedTimeStep};


class STPStatsTracer {

private:
	std::string _outputDir;
    static tarch::logging::Log     _log;

public:
	STPStatsTracer();
	virtual ~STPStatsTracer();

    void writeTracingEventIteration(unsigned int iterations, STPType type);
    void writeTracingEventRun(unsigned int elapsed, STPType type);
    void writeTracingEventRunIterations(unsigned int elapsed, unsigned int iterations, STPType type);

    void setOutputDir(std::string directory);
    static STPStatsTracer& getInstance();
};

} /* namespace offloading */
} /* namespace exahype */

#endif /* EXAHYPE_EXAHYPE_OFFLOADING_STPSTATSTRACER_H_ */
