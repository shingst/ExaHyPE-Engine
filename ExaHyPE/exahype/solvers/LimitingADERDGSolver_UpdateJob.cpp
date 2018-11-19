#include "LimitingADERDGSolver.h"

#if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
#include <immintrin.h>
#endif


exahype::solvers::LimitingADERDGSolver::UpdateJob::UpdateJob(
  LimitingADERDGSolver&  solver,
  SolverPatch&           solverPatch,
  CellInfo&              cellInfo,
  const bool             isAtRemoteBoundary):
  tarch::multicore::jobs::Job(tarch::multicore::jobs::JobType::RunTaskAsSoonAsPossible,0), // !! always high priority
  _solver(solver),
  _solverPatch(solverPatch),
  _cellInfo(cellInfo),
  _isAtRemoteBoundary(isAtRemoteBoundary) {
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    NumberOfReductionJobs++;
  }
  lock.free();
}

bool exahype::solvers::LimitingADERDGSolver::UpdateJob::run() {
  UpdateResult result =
      _solver.updateBody(_solverPatch,_cellInfo,_isAtRemoteBoundary);

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
void exahype::solvers::LimitingADERDGSolver::UpdateJob::prefetchData() {
  #if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
  double* luh    = static_cast<double*>(_solverPatch.getSolution());
  double* luhOld = static_cast<double*>(_solverPatch.getPreviousSolution());
  double* lduh   = static_cast<double*>(_solverPatch.getUpdate());
  double* lQhbnd = static_cast<double*>(_solverPatch.getExtrapolatedPredictor());
  double* lFhbnd = static_cast<double*>(_solverPatch.getFluctuation());

  _mm_prefetch(luh,    _MM_HINT_NTA);
  _mm_prefetch(luhOld, _MM_HINT_NTA);
  _mm_prefetch(lduh,   _MM_HINT_NTA);
  _mm_prefetch(lQhbnd, _MM_HINT_NTA);
  _mm_prefetch(lFhbnd, _MM_HINT_NTA);
  #endif
}
