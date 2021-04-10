#include "ADERDGSolver.h"

#if defined(ScoreP)
#include "scorep/SCOREP_User.h"
#endif

#if defined(FileTrace)
#include <chrono>
#include "exahype/reactive/STPStatsTracer.h"
#endif

#if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
#include <immintrin.h>
#endif

exahype::solvers::ADERDGSolver::PredictionJob::PredictionJob(
    ADERDGSolver&    solver,
    CellDescription& cellDescription,
    const int        cellDescriptionsIndex,
    const int        element,
    const double     predictorTimeStamp,
    const double     predictorTimeStepSize,
    const bool       uncompressBefore,
    const bool       isSkeletonJob,
    const bool       addVolumeIntegralResultToUpdate):
    tarch::multicore::jobs::Job(
        tarch::multicore::jobs::JobType::BackgroundTask,0,
        getTaskPriority(isSkeletonJob)
    ), // ! high priority only if skeleton job
    _solver(solver),
    _cellDescription(cellDescription),
    _cellDescriptionsIndex(cellDescriptionsIndex),
    _element(element),
    _predictorTimeStamp(predictorTimeStamp),
    _predictorTimeStepSize(predictorTimeStepSize),
    _uncompressBefore(uncompressBefore),
    _isSkeletonJob(isSkeletonJob),
    _addVolumeIntegralResultToUpdate(addVolumeIntegralResultToUpdate) {
  if (_isSkeletonJob) {
    NumberOfSkeletonJobs.fetch_add(1);
  } else {
    NumberOfEnclaveJobs.fetch_add(1);
  }
}


bool exahype::solvers::ADERDGSolver::PredictionJob::run(bool runOnMasterThread) {
  #if defined ScoreP
  SCOREP_USER_REGION( "exahype::solvers::ADERDGSolver::PredictionJob::run", SCOREP_USER_REGION_TYPE_FUNCTION ) 
  #endif
   
  #if defined FileTrace
  auto start = std::chrono::high_resolution_clock::now();
  #endif

  int numberIterations =  _solver.predictionAndVolumeIntegralBody(
      _cellDescription,_predictorTimeStamp,_predictorTimeStepSize,
      _uncompressBefore,_isSkeletonJob,_addVolumeIntegralResultToUpdate); // ignore return value

  if (_isSkeletonJob) {
    NumberOfSkeletonJobs.fetch_sub(1);
    assertion( NumberOfSkeletonJobs.load()>=0 );
  } else {
    NumberOfEnclaveJobs.fetch_sub(1);
    assertion( NumberOfEnclaveJobs.load()>=0 );
  }
  #if defined FileTrace
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

  //std::stringstream stream;
  #if defined Picard
  exahype::reactive::STPStatsTracer::getInstance().writeTracingEventRunIterations(duration.count(), numberIterations, exahype::reactive::STPTraceKey::ADERDGPrediction);
  #else
  exahype::reactive::STPStatsTracer::getInstance().writeTracingEventRun(duration.count(), exahype::reactive::STPTraceKey::ADERDGPrediction);
  exahype::reactive::STPStatsTracer::getInstance().writeTracingEventIterations(numberIterations, exahype::reactive::STPTraceKey::ADERDGPrediction);
  #endif 
//  int rank=tarch::parallel::Node::getInstance().getRank();
//  #if defined SharedTBB
//  stream<<rank<<"_";
//  int threadId=tarch::multicore::Core::getInstance().getThreadNum();
//  stream<<threadId;
//  #endif
//  stream<<".txt";
//  std::string path=stream.str();
//
//  std::ofstream file;
//  file.open(path,std::fstream::app);
//  #if defined Picard
//  file << duration.count() <<":"<<(numberIterations+1)<< std::endl;
//  file.close();
//  return false;
//  #else
//  file << duration.count() << std::endl;
//  #endif
//  file.close();
//
//  stream.str(std::string());
//  stream<<"./TraceOutput/exahype_solvers_ADERDGSolver_PredictionJob_iterations_rank_";
//  stream<<rank<<"_";
//  stream<<threadId<<".txt";
//  path=stream.str();
//
//  std::ofstream fileIter;
//  fileIter.open(path,std::fstream::app);
//  fileIter << (numberIterations+1) << std::endl;
//  fileIter.close();
  #endif

  return false;
}

//
// @see PredictionJob
//
void exahype::solvers::ADERDGSolver::PredictionJob::prefetchData() {
#if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
  double* luh  = static_cast<double*>(_cellDescription.getSolution());
  double* lduh = static_cast<double*>(_cellDescription.getUpdate());
  double* lQhbnd = static_cast<double*>(_cellDescription.getExtrapolatedPredictor());
  double* lFhbnd = static_cast<double*>(_cellDescription.getFluctuation());

  _mm_prefetch(luh, _MM_HINT_NTA);
  _mm_prefetch(lduh, _MM_HINT_NTA);
  _mm_prefetch(lQhbnd, _MM_HINT_NTA);
  _mm_prefetch(lFhbnd, _MM_HINT_NTA);
#endif
}

