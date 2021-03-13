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

#include "../reactive/STPStatsTracer.h"

#include <sstream>
#include <fstream>
#include <iostream>

#include "tarch/parallel/Node.h"
#include "tarch/multicore/Core.h"

#include "exahype/solvers/Solver.h"

namespace exahype {
namespace reactive {

tarch::logging::Log  exahype::reactive::STPStatsTracer::_log( "exahype::reactive::STPStatsTracer" );

STPStatsTracer::STPStatsTracer() : _outputDir("."), _dumpInterval(1), _dumpCnt(0){
	// TODO Auto-generated constructor stub
  for(int type = STPTraceKey::ADERDGPrediction; type<=STPTraceKey::LimitingFusedTimeStep; type++) {
    _iterations[type].resize(tarch::multicore::Core::getInstance().getNumberOfThreads());
    _elapsed[type].resize(tarch::multicore::Core::getInstance().getNumberOfThreads());
    _semaphores[type] = new tarch::multicore::BooleanSemaphore[tarch::multicore::Core::getInstance().getNumberOfThreads()];

    for(int i=0 ; i<tarch::multicore::Core::getInstance().getNumberOfThreads(); i++) {
      //_semaphores[type].push_back(std::move(tarch::multicore::BooleanSemaphore()));
      _locks[type].push_back(tarch::multicore::Lock(_semaphores[type][i], false));
    }
  }
}

STPStatsTracer::~STPStatsTracer() {
	// TODO Auto-generated destructor stub
	for(int type = STPTraceKey::ADERDGPrediction; type<=STPTraceKey::LimitingFusedTimeStep; type++)
	  delete[] _semaphores[type];
}

void STPStatsTracer::setOutputDir(std::string outputDir) {
  _outputDir = outputDir;
  std::stringstream stream;
  stream.str(std::string());
  stream<<"mkdir -p "<<outputDir;
  int err = system(stream.str().c_str());
  if(err==-1) {
    logError("setOutputDir", "Could not create output directory for STP statistics, exiting...");
    exit(-1);
  }
}

void STPStatsTracer::setDumpInterval(int interval) {
  _dumpInterval = interval;
}

void STPStatsTracer::dumpAndResetTraceIfActive() {
  _dumpCnt++;

  int timestep = _dumpCnt;

  if (isActive(timestep)) {
    int rank = tarch::parallel::Node::getInstance().getRank();

    for(int type = STPTraceKey::ADERDGPrediction; type<=STPTraceKey::LimitingFusedTimeStep; type++) {
    #if defined (SharedTBB)
      //int threadId = tarch::multicore::Core::getInstance().getThreadNum();
      for(size_t threadId = 0; threadId<_elapsed[type].size(); threadId++) {
    #else
      unsigned int threadId = 0;
    #endif
        std::stringstream stream;
        stream<<_outputDir<<"/exahype_";

        switch(type) {
          case STPTraceKey::ADERDGPrediction:
            stream<<"solvers_ADERDGSolver_PredictionJob_"; break;
          case STPTraceKey::ADERDGOwnMigratable:
            stream<<"solvers_ADERDGSolver_OwnMigratableJob_"; break;
          case STPTraceKey::ADERDGRemoteMigratable:
            stream<<"solvers_ADERDGSolver_AlienMigratableJob_"; break;
          case STPTraceKey::LimitingFusedTimeStep:
            stream<<"solvers_LimitingADERDGSolver_FusedTimeStepJob_"; break;
        }

        stream<<"run_rank_"<<rank;
#if defined (SharedTBB)
        stream<<"_"<<threadId;
#endif
        stream<<"_step_"<<timestep;
        stream<<".txt";
        std::string path=stream.str();

        std::ofstream file;
        file.open(path,std::fstream::app);
        _locks[type][threadId].lock();
        file << _elapsed[type][threadId]<< ":" << _iterations[type][threadId]+1 << std::endl;
        file.close();
        
        _elapsed[type][threadId] = 0;
        _iterations[type][threadId] = 0;
        _locks[type][threadId].free();

#if defined (SharedTBB)
      }
#endif
    }
  }
}


bool STPStatsTracer::isActive(int timestep) {
  return (timestep  % _dumpInterval) == 0;
}

void STPStatsTracer::writeTracingEventIteration(unsigned int iterations, STPTraceKey type) {

#if defined(SharedTBB)
  unsigned int threadId = tarch::multicore::Core::getInstance().getThreadNum();
#else
  unsigned int threadId = 0;
#endif

  //if(threadId>_iterations[type].size()) {
  //  _iterations[type].resize(threadId+1);
 // }
  _locks[type][threadId].lock();
  _iterations[type][threadId]+= iterations;
  _locks[type][threadId].free();
}

void STPStatsTracer::writeTracingEventRun(unsigned int elapsed, STPTraceKey type) {

#if defined(SharedTBB)
  unsigned int threadId = tarch::multicore::Core::getInstance().getThreadNum();
#else
  unsigned int threadId = 0;
#endif

//  if(threadId>_elapsed[type].size()) {
//    _elapsed[type].resize(threadId+1);
//  }
  _locks[type][threadId].lock();
  _elapsed[type][threadId]+= elapsed;
  _locks[type][threadId].free();
}

void STPStatsTracer::writeTracingEventRunIterations(unsigned int elapsed, unsigned int iterations, STPTraceKey type) {
#if defined(SharedTBB)
  unsigned int threadId = tarch::multicore::Core::getInstance().getThreadNum();
#else
  unsigned int threadId = 0;
#endif

//  if(threadId>_elapsed[type].size()) {
//    _elapsed[type].resize(threadId+1);
//  }

//  if(threadId>_iterations[type].size()) {
//    _iterations[type].resize(threadId+1);
//  }
  _locks[type][threadId].lock();
  _elapsed[type][threadId]+= elapsed;
  _iterations[type][threadId]+= iterations;
  _locks[type][threadId].free();
}

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
