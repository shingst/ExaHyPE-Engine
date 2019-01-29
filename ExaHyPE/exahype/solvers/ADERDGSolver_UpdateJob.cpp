#include "ADERDGSolver.h"

#if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
#include <immintrin.h>
#endif


exahype::solvers::ADERDGSolver::UpdateJob::UpdateJob(
  ADERDGSolver&    solver,
  CellDescription& cellDescription,
  CellInfo&        cellInfo,
  const bool       isAtRemoteBoundary):
  tarch::multicore::jobs::Job(tarch::multicore::jobs::JobType::RunTaskAsSoonAsPossible,0), // !! always high priority
  _solver(solver),
  _cellDescription(cellDescription),
  _cellInfo(cellInfo),
  _isAtRemoteBoundary(isAtRemoteBoundary) {
  NumberOfReductionJobs.fetch_add(1);
}

bool exahype::solvers::ADERDGSolver::UpdateJob::run() {
  UpdateResult result =
      _solver.updateBody(_cellDescription,_cellInfo,_neighbourMergePerformed,_isAtRemoteBoundary);

  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    _solver.updateMeshUpdateEvent(result._meshUpdateEvent);
    _solver.updateMinNextTimeStepSize(result._timeStepSize);
  }
  lock.free();


  NumberOfReductionJobs.fetch_sub(1);
  assertion( NumberOfReductionJobs.load()>=0 );
  return false;
}


//
// @see UpdateJob
//
void exahype::solvers::ADERDGSolver::UpdateJob::prefetchData() {
  #if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
  double* luh    = static_cast<double*>(_cellDescription.getSolution());
  double* luhOld = static_cast<double*>(_cellDescription.getPreviousSolution());
  double* lduh   = static_cast<double*>(_cellDescription.getUpdate());
  double* lQhbnd = static_cast<double*>(_cellDescription.getExtrapolatedPredictor());
  double* lFhbnd = static_cast<double*>(_cellDescription.getFluctuation());

  _mm_prefetch(luh,    _MM_HINT_NTA);
  _mm_prefetch(luhOld, _MM_HINT_NTA);
  _mm_prefetch(lduh,   _MM_HINT_NTA);
  _mm_prefetch(lQhbnd, _MM_HINT_NTA);
  _mm_prefetch(lFhbnd, _MM_HINT_NTA);
  #endif
}
