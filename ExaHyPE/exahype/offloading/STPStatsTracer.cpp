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

#include <sstream>
#include <fstream>
#include <iostream>

#include "tarch/parallel/Node.h"
#include "tarch/multicore/Core.h"

#include "exahype/offloading/STPStatsTracer.h"
#include "exahype/solvers/Solver.h"

namespace exahype {
namespace offloading {

tarch::logging::Log  exahype::offloading::STPStatsTracer::_log( "exahype::offloading::STPStatsTracer" );

STPStatsTracer::STPStatsTracer() : _outputDir("."), _dumpInterval(1), _dumpCnt(0) {
	// TODO Auto-generated constructor stub

}

STPStatsTracer::~STPStatsTracer() {
	// TODO Auto-generated destructor stub
}

void STPStatsTracer::setOutputDir(std::string outputDir) {
  _outputDir = outputDir;
  std::stringstream stream;
  stream.str(std::string());
  stream<<"mkdir -p "<<outputDir;
  (void) system(stream.str().c_str());
}

void STPStatsTracer::setDumpInterval(int interval) {
  _dumpInterval = interval;
}

void STPStatsTracer::dumpAndResetTraceIfActive() {
  _dumpCnt++;

  int timestep = _dumpCnt/exahype::solvers::Solver::PredictionSweeps;

  if (isActive(timestep)) {
    int rank = tarch::parallel::Node::getInstance().getRank();
#if defined (SharedTBB)
    int threadId = tarch::multicore::Core::getInstance().getThreadNum();
#endif

    for(int type = STPTraceKey::ADERDGPrediction; type<STPTraceKey::LimitingFusedTimeStep; type++) {
      std::stringstream stream;
      stream<<_outputDir<<"/exahype_";

      switch(type) {
        case STPTraceKey::ADERDGPrediction:
          stream<<"solvers_ADERDGSolver_PredictionJob_"; break;
        case STPTraceKey::ADERDGOwnMigratable:
          stream<<"_solvers_ADERDGSolver_OwnMigratableJob_"; break;
        case STPTraceKey::ADERDGRemoteMigratable:
          stream<<"_solvers_ADERDGSolver_AlienMigratableJob_"; break;
        case STPTraceKey::LimitingFusedTimeStep:
          stream<<"_solvers_LimitingADERDGSolver_FusedTimeStepJob_"; break;
      }

      stream<<"_run_rank_"<<rank;
#if defined (SharedTBB)
      stream<<"_"<<threadId;
#endif
      stream<<"_step_"<<timestep;
      stream<<".txt";
      std::string path=stream.str();

      std::ofstream file;
      file.open(path,std::fstream::app);
      file << _elapsed[type][threadId]<< ":" << _iterations[type][threadId]+1 << std::endl;
      file.close();

    }
  }
}


bool STPStatsTracer::isActive(int timestep) {
  return (timestep  % _dumpInterval) == 0;
}

void STPStatsTracer::writeTracingEventIteration(unsigned int iterations, STPTraceKey type) {

#if defined(SharedTBB)
  unsigned int threadId = tarch::multicore::Core::getInstance().getThreadNum();
#endif

  if(threadId>_iterations[type].size()) {
    _iterations[type].resize(threadId);
  }

  _iterations[type][threadId]+= iterations;
}

void STPStatsTracer::writeTracingEventRun(unsigned int elapsed, STPTraceKey type) {

#if defined(SharedTBB)
  unsigned int threadId = tarch::multicore::Core::getInstance().getThreadNum();
#endif

  if(threadId>_elapsed[type].size()) {
    _elapsed[type].resize(threadId);
  }

  _elapsed[type][threadId]+= elapsed;
}

void STPStatsTracer::writeTracingEventRunIterations(unsigned int iterations, unsigned int elapsed, STPTraceKey type) {
#if defined(SharedTBB)
  unsigned int threadId = tarch::multicore::Core::getInstance().getThreadNum();
#endif

  if(threadId>_elapsed[type].size()) {
    _elapsed[type].resize(threadId);
  }

  if(threadId>_iterations[type].size()) {
    _iterations[type].resize(threadId);
  }

  _elapsed[type][threadId]+= elapsed;
  _iterations[type][threadId]+= iterations;
}


/*void STPStatsTracer::writeTracingEventIteration(unsigned int iterations, STPType type) {

  std::stringstream streamIt;
  streamIt.str(std::string());
  streamIt<<_outputDir<<"/exahype_";

  switch(type) {
    case STPType::ADERDGPrediction:
      streamIt<<"_solvers_ADERDGSolver_PredictionJob_"; break;
    case STPType::ADERDGOwnMigratable:
      streamIt<<"_solvers_ADERDGSolver_OwnMigratableJob_"; break;
    case STPType::ADERDGRemoteMigratable:
      streamIt<<"_solvers_ADERDGSolver_AlienMigratableJob_"; break;
    case STPType::LimitingFusedTimeStep:
      streamIt<<"_solvers_LimitingADERDGSolver_FusedTimeStepJob_"; break;
  }

  int rank=tarch::parallel::Node::getInstance().getRank();
  streamIt<<"iterations_rank_"<<rank;
#if defined(SharedTBB)
  int threadId=tarch::multicore::Core::getInstance().getThreadNum();
  streamIt<<"_"<<threadId;
#endif
  streamIt<<".txt";

  std::string pathIt=streamIt.str();

  std::ofstream file;
  file.open(pathIt,std::fstream::app);
  file << (iterations+1) << std::endl;
  file.close();

  //logInfo("writeTracingEventIteration", " wrote to file "<<pathIt);
}

void STPStatsTracer::writeTracingEventRun(unsigned int elapsed, STPType type) {
 //"exahype_solvers_ADERDGSolver_MigratablePredictionJob_run_rank_"
  std::stringstream stream;
  stream<<_outputDir<<"/exahype_";

  switch(type) {
    case STPType::ADERDGPrediction:
      stream<<"solvers_ADERDGSolver_PredictionJob_"; break;
    case STPType::ADERDGOwnMigratable:
      stream<<"_solvers_ADERDGSolver_OwnMigratableJob_"; break;
    case STPType::ADERDGRemoteMigratable:
      stream<<"_solvers_ADERDGSolver_AlienMigratableJob_"; break;
    case STPType::LimitingFusedTimeStep:
      stream<<"_solvers_LimitingADERDGSolver_FusedTimeStepJob_"; break;
  }

  int rank = tarch::parallel::Node::getInstance().getRank();
  stream<<"_run_rank_"<<rank;
#if defined (SharedTBB)
  int threadId=tarch::multicore::Core::getInstance().getThreadNum();
  stream<<"_"<<threadId;
#endif
  stream<<".txt";
  std::string path=stream.str();

  std::ofstream file;
  file.open(path,std::fstream::app);
  file << elapsed << std::endl;
  file.close();
  
  //logInfo("writeTracingEventRun", " wrote to file "<<path);
}

void STPStatsTracer::writeTracingEventRunIterations(unsigned int iterations, unsigned int elapsed, STPType type) {
 //"exahype_solvers_ADERDGSolver_MigratablePredictionJob_run_rank_"
  std::stringstream stream;
  stream<<_outputDir<<"/exahype_";

  switch(type) {
    case STPType::ADERDGPrediction:
      stream<<"solvers_ADERDGSolver_PredictionJob_"; break;
    case STPType::ADERDGOwnMigratable:
      stream<<"_solvers_ADERDGSolver_OwnMigratableJob_"; break;
    case STPType::ADERDGRemoteMigratable:
      stream<<"_solvers_ADERDGSolver_AlienMigratableJob_"; break;
    case STPType::LimitingFusedTimeStep:
      stream<<"_solvers_LimitingADERDGSolver_FusedTimeStepJob_"; break;
  }

  int rank = tarch::parallel::Node::getInstance().getRank();
  stream<<"run_iterations_rank_"<<rank;
#if defined (SharedTBB)
  int threadId=tarch::multicore::Core::getInstance().getThreadNum();
  stream<<"_"<<threadId;
#endif
  stream<<".txt";
  std::string path=stream.str();

  std::ofstream file;
  file.open(path,std::fstream::app);
  file << elapsed<< ":" << iterations+1 << std::endl;
  file.close();

  //logInfo("writeTracingEventRunIterations", " wrote to file "<<path);
}
*/

void STPStatsTracer::writeTracingEventIterationDetailed(unsigned int iterations, STPTraceKey type) {

  std::stringstream streamIt;
  streamIt.str(std::string());
  streamIt<<_outputDir<<"/exahype_";

  switch(type) {
    case STPTraceKey::ADERDGPrediction:
      streamIt<<"_solvers_ADERDGSolver_PredictionJob_"; break;
    case STPTraceKey::ADERDGOwnMigratable:
      streamIt<<"_solvers_ADERDGSolver_OwnMigratableJob_"; break;
    case STPTraceKey::ADERDGRemoteMigratable:
      streamIt<<"_solvers_ADERDGSolver_AlienMigratableJob_"; break;
    case STPTraceKey::LimitingFusedTimeStep:
      streamIt<<"_solvers_LimitingADERDGSolver_FusedTimeStepJob_"; break;
  }

  int rank=tarch::parallel::Node::getInstance().getRank();
  streamIt<<"iterations_rank_"<<rank;
#if defined(SharedTBB)
  int threadId=tarch::multicore::Core::getInstance().getThreadNum();
  streamIt<<"_"<<threadId;
#endif
  streamIt<<".txt";

  std::string pathIt=streamIt.str();

  std::ofstream file;
  file.open(pathIt,std::fstream::app);
  file << (iterations+1) << std::endl;
  file.close();

  //logInfo("writeTracingEventIteration", " wrote to file "<<pathIt);
}

void STPStatsTracer::writeTracingEventRunDetailed(unsigned int elapsed, STPTraceKey type) {
 //"exahype_solvers_ADERDGSolver_MigratablePredictionJob_run_rank_"
  std::stringstream stream;
  stream<<_outputDir<<"/exahype_";

  switch(type) {
    case STPTraceKey::ADERDGPrediction:
      stream<<"solvers_ADERDGSolver_PredictionJob_"; break;
    case STPTraceKey::ADERDGOwnMigratable:
      stream<<"_solvers_ADERDGSolver_OwnMigratableJob_"; break;
    case STPTraceKey::ADERDGRemoteMigratable:
      stream<<"_solvers_ADERDGSolver_AlienMigratableJob_"; break;
    case STPTraceKey::LimitingFusedTimeStep:
      stream<<"_solvers_LimitingADERDGSolver_FusedTimeStepJob_"; break;
  }

  int rank = tarch::parallel::Node::getInstance().getRank();
  stream<<"_run_rank_"<<rank;
#if defined (SharedTBB)
  int threadId=tarch::multicore::Core::getInstance().getThreadNum();
  stream<<"_"<<threadId;
#endif
  stream<<".txt";
  std::string path=stream.str();

  std::ofstream file;
  file.open(path,std::fstream::app);
  file << elapsed << std::endl;
  file.close();

  //logInfo("writeTracingEventRun", " wrote to file "<<path);
}

void STPStatsTracer::writeTracingEventRunIterationsDetailed(unsigned int iterations, unsigned int elapsed, STPTraceKey type) {
 //"exahype_solvers_ADERDGSolver_MigratablePredictionJob_run_rank_"
  std::stringstream stream;
  stream<<_outputDir<<"/exahype_";

  switch(type) {
    case STPTraceKey::ADERDGPrediction:
      stream<<"solvers_ADERDGSolver_PredictionJob_"; break;
    case STPTraceKey::ADERDGOwnMigratable:
      stream<<"_solvers_ADERDGSolver_OwnMigratableJob_"; break;
    case STPTraceKey::ADERDGRemoteMigratable:
      stream<<"_solvers_ADERDGSolver_AlienMigratableJob_"; break;
    case STPTraceKey::LimitingFusedTimeStep:
      stream<<"_solvers_LimitingADERDGSolver_FusedTimeStepJob_"; break;
  }

  int rank = tarch::parallel::Node::getInstance().getRank();
  stream<<"run_iterations_rank_"<<rank;
#if defined (SharedTBB)
  int threadId=tarch::multicore::Core::getInstance().getThreadNum();
  stream<<"_"<<threadId;
#endif
  stream<<".txt";
  std::string path=stream.str();

  std::ofstream file;
  file.open(path,std::fstream::app);
  file << elapsed<< ":" << iterations+1 << std::endl;
  file.close();

  //logInfo("writeTracingEventRunIterations", " wrote to file "<<path);
}

STPStatsTracer& STPStatsTracer::getInstance() {
  static STPStatsTracer tracer;
  return tracer;
}


} /* namespace offloading */
} /* namespace exahype */
