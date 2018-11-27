#include "exahype/solvers/FiniteVolumesSolver.h"

exahype::solvers::FiniteVolumesSolver::AdjustSolutionDuringMeshRefinementJob::AdjustSolutionDuringMeshRefinementJob(
  FiniteVolumesSolver& solver,
  CellDescription&     cellDescription,
  const bool           isInitialMeshRefinement):
  tarch::multicore::jobs::Job(Solver::getTaskType(false),0),
  _solver(solver),
  _cellDescription(cellDescription),
  _isInitialMeshRefinement(isInitialMeshRefinement)
{
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    NumberOfAMRBackgroundJobs++;
  }
  lock.free();
}

bool exahype::solvers::FiniteVolumesSolver::AdjustSolutionDuringMeshRefinementJob::run() {
  _solver.adjustSolutionDuringMeshRefinementBody(_cellDescription,_isInitialMeshRefinement);

  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    NumberOfAMRBackgroundJobs--;
    assertion( NumberOfAMRBackgroundJobs>=0 );
  }
  lock.free();
  return false;
}
