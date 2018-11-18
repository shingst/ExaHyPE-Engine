#include "exahype/solvers/FiniteVolumesSolver.h"

#if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
#include <immintrin.h>
#endif


exahype::solvers::FiniteVolumesSolver::UpdateJob::UpdateJob(
  FiniteVolumesSolver&    solver,
  CellDescription& cellDescription,
  CellInfo&        cellInfo,
  const bool       isAtRemoteBoundary):
  tarch::multicore::jobs::Job(tarch::multicore::jobs::JobType::RunTaskAsSoonAsPossible,0), // !! always high priority
  _solver(solver),
  _cellDescription(cellDescription),
  _cellInfo(cellInfo),
  _neighbourMergePerformed(cellDescription.getNeighbourMergePerformed()),
  _isAtRemoteBoundary(isAtRemoteBoundary) {
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    NumberOfReductionJobs++;
  }
  lock.free();
}

bool exahype::solvers::FiniteVolumesSolver::UpdateJob::run() {
  UpdateResult result =
      _solver.updateBody(
          _cellDescription,_cellInfo,_neighbourMergePerformed,
          true,true,_isAtRemoteBoundary,true);

  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    _solver.updateNextMeshUpdateEvent(result._meshUpdateEvent);
    _solver.updateMinNextTimeStepSize(result._timeStepSize);

    NumberOfReductionJobs--;
    assertion( NumberOfReductionJobs>=0 );
  }
  lock.free();
  return false;
}


//
// @see UpdateJob
//
void exahype::solvers::FiniteVolumesSolver::UpdateJob::prefetchData() {
  #if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
  double* luh    = static_cast<double*>(_cellDescription.getSolution());
  double* luhOld = static_cast<double*>(_cellDescription.getPreviousSolution());

  _mm_prefetch(luh,    _MM_HINT_NTA);
  _mm_prefetch(luhOld, _MM_HINT_NTA);
  #endif
}
