#include "ADERDGSolver.h"

#if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
#include <immintrin.h>
#endif


exahype::solvers::ADERDGSolver::FusedTimeStepJob::FusedTimeStepJob(
  ADERDGSolver&    solver,
  CellDescription& cellDescription,
  const int        cellDescriptionsIndex,
  const int        element,
  const bool       isSkeletonJob):
  tarch::multicore::jobs::Job(Solver::getTaskType(isSkeletonJob),0),
  _solver(solver),
  _cellDescription(cellDescription),
  _cellDescriptionsIndex(cellDescriptionsIndex),
  _element(element),
  _neighbourMergePerformed(cellDescription.getNeighbourMergePerformed()),
  _isSkeletonJob(isSkeletonJob) {
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter++;
  }
  lock.free();
}

bool exahype::solvers::ADERDGSolver::FusedTimeStepJob::run() {
  _solver.fusedTimeStepBody(
      _cellDescription, _cellDescriptionsIndex, _element, false, false, _isSkeletonJob, false, _neighbourMergePerformed );

  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter--;
    assertion( jobCounter>=0 );
  }
  lock.free();
  return false;
}


//
// @see PredictionJob
//
void exahype::solvers::ADERDGSolver::FusedTimeStepJob::prefetchData() {
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
