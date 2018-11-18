#include "ADERDGSolver.h"

exahype::solvers::ADERDGSolver::AdjustSolutionDuringMeshRefinementJob::AdjustSolutionDuringMeshRefinementJob(
  ADERDGSolver&    solver,
  CellDescription& cellDescription,
  const bool       isInitialMeshRefinement):
  tarch::multicore::jobs::Job( Solver::getTaskType(false), 0 ),
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

bool exahype::solvers::ADERDGSolver::AdjustSolutionDuringMeshRefinementJob::run() {
  _solver.ensureNecessaryMemoryIsAllocated(_cellDescription);
  _solver.adjustSolutionDuringMeshRefinementBody(_cellDescription,_isInitialMeshRefinement);

  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    NumberOfAMRBackgroundJobs--;
    assertion( NumberOfAMRBackgroundJobs>=0 );
  }
  lock.free();
  return false;
}
