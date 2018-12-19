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

bool exahype::solvers::FiniteVolumesSolver::FusedTimeStepJob::run() {
  UpdateResult result =
      _solver.updateBody(
          _cellDescription,_cellInfo,_neighbourMergePerformed,
          _isFirstTimeStepOfBatch,_isLastTimeStepOfBatch,
          _isSkeletonJob,false/*uncompressBefore*/);

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
