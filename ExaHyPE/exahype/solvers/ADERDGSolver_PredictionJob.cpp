#include "ADERDGSolver.h"

#if defined(ScoreP)
#include "scorep/SCOREP_User.h"
#endif

#if defined(FileTrace)
#include <iostream>
#include <fstream> 
#include <string>
#include <ctime>
#include <unistd.h>
#include <sys/time.h>
#include <sstream>
#include "tarch/parallel/Node.h"
#include "tarch/multicore/Core.h"
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
  struct timeval start_time, end_time;
  long milli_time, seconds, useconds;
  gettimeofday(&start_time, NULL);
  #endif

  _solver.predictionAndVolumeIntegralBody(
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
  gettimeofday(&end_time, NULL);
  seconds = end_time.tv_sec - start_time.tv_sec; //seconds
  useconds = end_time.tv_usec - start_time.tv_usec; //milliseconds
  milli_time = ((seconds) * 1000 + useconds/1000.0);

  std::stringstream stream;
  stream<<"./TraceOutput/exahype_solvers_ADERDGSolver_PredictionJob_run_rank_";
  int rank=tarch::parallel::Node::getInstance().getRank();
  stream<<rank<<"_";
  //this will only work for 2 cores per Rank
  int threadId=tarch::multicore::Core::getInstance().getThreadNum();
  stream<<threadId<<".txt";
  std::string path=stream.str();

  //char cstr[path.size()];
  //path.copy(cstr,path.size());
  std::ofstream file;
  file.open(path,std::fstream::app);
  file << milli_time << std::endl;
  file.close();
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

