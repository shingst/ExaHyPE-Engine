#include "exahype/solvers/LimitingADERDGSolver.h"

#if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
#include <immintrin.h>
#endif

exahype::solvers::LimitingADERDGSolver::FusedTimeStepJob::FusedTimeStepJob(
  LimitingADERDGSolver& solver,
  SolverPatch&          solverPatch,
  CellInfo&             cellInfo,
  const bool            isFirstTimeStepOfBatch,
  const bool            isLastTimeStepOfBatch,
  const bool            isSkeletonJob):
  tarch::multicore::jobs::Job(
      isLastTimeStepOfBatch ?
          tarch::multicore::jobs::JobType::RunTaskAsSoonAsPossible :
          Solver::getTaskType(isSkeletonJob),
  0),
  _solver(solver),
  _solverPatch(solverPatch),
  _cellInfo(cellInfo),
  _neighbourMergePerformed(solverPatch.getNeighbourMergePerformed()),
  _isFirstTimeStepOfBatch(isFirstTimeStepOfBatch),
  _isLastTimeStepOfBatch(isLastTimeStepOfBatch),
  _isSkeletonJob(isSkeletonJob) {
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    NumberOfReductionJobs++;

    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter++;
  }
  lock.free();
}

bool exahype::solvers::LimitingADERDGSolver::FusedTimeStepJob::run() {
  UpdateResult result =
      _solver.fusedTimeStepBody(
          _solverPatch,_cellInfo,_neighbourMergePerformed,
          _isFirstTimeStepOfBatch,_isLastTimeStepOfBatch,
          _isSkeletonJob,false/*mustBeDoneImmedetially*/);

  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    if (_isLastTimeStepOfBatch) {
      _solver.updateNextMeshUpdateEvent(result._meshUpdateEvent);
      _solver.updateMinNextTimeStepSize(result._timeStepSize);
    }
    NumberOfReductionJobs--;
    assertion( NumberOfReductionJobs>=0 );

    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter--;
    assertion( jobCounter>=0 );
  }
  lock.free();
  return false;
}
