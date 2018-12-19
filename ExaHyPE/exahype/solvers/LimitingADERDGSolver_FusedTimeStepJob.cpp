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
  NumberOfReductionJobs++;
  if (_isSkeletonJob) {
    NumberOfSkeletonJobs++;
  } else {
    NumberOfEnclaveJobs++;
  }
}

bool exahype::solvers::LimitingADERDGSolver::FusedTimeStepJob::run() {
  UpdateResult result =
      _solver.fusedTimeStepBody(
          _solverPatch,_cellInfo,_neighbourMergePerformed,
          _isFirstTimeStepOfBatch,_isLastTimeStepOfBatch,
          _isSkeletonJob,false/*mustBeDoneImmedetially*/);

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
