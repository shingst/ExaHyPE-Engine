#include "exahype/solvers/FiniteVolumesSolver.h"

exahype::solvers::FiniteVolumesSolver::FusedTimeStepJob::FusedTimeStepJob(
  FiniteVolumesSolver& solver,
  CellDescription&     cellDescription,
  CellInfo&            cellInfo,
  const bool           isFirstTimeStepOfBatch,
  const bool           isLastTimeStepOfBatch,
  const bool           isSkeletonJob)
  :
  tarch::multicore::jobs::Job(
      isLastTimeStepOfBatch ?
          tarch::multicore::jobs::JobType::RunTaskAsSoonAsPossible :
          Solver::getTaskType(isSkeletonJob),
  0),
  _solver(solver),
  _cellDescription(cellDescription),
  _cellInfo(cellInfo),
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

bool exahype::solvers::FiniteVolumesSolver::FusedTimeStepJob::run() {
  UpdateResult result =
      _solver.updateBody(
          _cellDescription,_cellInfo,_isFirstTimeStepOfBatch,_isLastTimeStepOfBatch,
          _isSkeletonJob,false/*uncompressBefore*/);

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
