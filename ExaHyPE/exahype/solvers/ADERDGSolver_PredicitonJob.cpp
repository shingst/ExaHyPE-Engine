#include "ADERDGSolver.h"

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
  const bool       isSkeletonJob):
  tarch::multicore::jobs::Job(Solver::getTaskType(isSkeletonJob),0),
  _solver(solver),
  _cellDescription(cellDescription),
  _cellDescriptionsIndex(cellDescriptionsIndex),
  _element(element),
  _predictorTimeStamp(predictorTimeStamp),
  _predictorTimeStepSize(predictorTimeStepSize),
  _uncompressBefore(uncompressBefore),
  _isSkeletonJob(isSkeletonJob) {
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter++;
  }
  lock.free();
}


bool exahype::solvers::ADERDGSolver::PredictionJob::run() {
  _solver.performPredictionAndVolumeIntegralBody(
      _cellDescription,_predictorTimeStamp,_predictorTimeStepSize,
      _uncompressBefore,_isSkeletonJob); // ignore return value
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter--;
    assertion( jobCounter>=0 );
  }
  lock.free();
  return false;
}

void exahype::solvers::ADERDGSolver::PredictionJob::prefetchData() {
  #if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
  const CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,_element);

  double* luh  = static_cast<double*>(cellDescription.getSolution());
  double* lduh = static_cast<double*>(cellDescription.getUpdate());
  double* lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double* lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

  _mm_prefetch(luh, _MM_HINT_T2); // move at least into L3; _MM_HINT_NTA should be tried as well
  _mm_prefetch(lduh, _MM_HINT_T2);
  _mm_prefetch(lQhbnd, _MM_HINT_T2);
  _mm_prefetch(lFhbnd, _MM_HINT_T2);
  #endif
}

