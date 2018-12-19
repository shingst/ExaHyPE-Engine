#include "ADERDGSolver.h"

#if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
#include <immintrin.h>
#endif


exahype::solvers::ADERDGSolver::FusedTimeStepJob::FusedTimeStepJob(
  ADERDGSolver&    solver,
  CellDescription& cellDescription,
  CellInfo&        cellInfo,
  const bool       isFirstTimeStepOfBatch,
  const bool       isLastTimeStepOfBatch,
  const bool       isSkeletonJob)
  :
  tarch::multicore::jobs::Job(
      isLastTimeStepOfBatch ?
          tarch::multicore::jobs::JobType::RunTaskAsSoonAsPossible :
          Solver::getTaskType(isSkeletonJob),
  0),
  _solver(solver),
  _cellDescription(cellDescription),
  _cellInfo(cellInfo),
  _neighbourMergePerformed(cellDescription.getNeighbourMergePerformed()),
  _isFirstTimeStepOfBatch(isFirstTimeStepOfBatch),
  _isLastTimeStepOfBatch(isLastTimeStepOfBatch),
  _isSkeletonJob(isSkeletonJob) {
  NumberOfReductionJobs++;
  if (_isSkeletonJob) {
    NumberOfSkeletonJobs++;
  } else {
    NumberOfEnclaveJobs++;
  }
}

bool exahype::solvers::ADERDGSolver::FusedTimeStepJob::run() {
  UpdateResult result =
      _solver.fusedTimeStepBody(
          _cellDescription, _cellInfo, _neighbourMergePerformed,
          _isFirstTimeStepOfBatch,_isLastTimeStepOfBatch,
          _isSkeletonJob,false/*mustBeDoneImmediately*/);

    if (_isLastTimeStepOfBatch) {
      tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
      {
        _solver.updateNextMeshUpdateEvent(result._meshUpdateEvent);
        _solver.updateMinNextTimeStepSize(result._timeStepSize);
      }
      lock.free();
    }

  NumberOfReductionJobs--;
  assertion( NumberOfReductionJobs.load()>=0 );
  if (_isSkeletonJob) {
    NumberOfSkeletonJobs--;
    assertion( NumberOfSkeletonJobs.load()>=0 );
  } else {
    NumberOfEnclaveJobs--;
    assertion( NumberOfEnclaveJobs.load()>=0 );
  }
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
