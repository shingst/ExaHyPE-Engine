#include "exahype/solvers/LimitingADERDGSolver.h"

exahype::solvers::LimitingADERDGSolver::AdjustSolutionDuringMeshRefinementJob::AdjustSolutionDuringMeshRefinementJob(
  LimitingADERDGSolver& solver,
  SolverPatch&          solverPatch,
  CellInfo&             cellInfo,
  const bool            isInitialMeshRefinement):
  tarch::multicore::jobs::Job( Solver::getTaskType(true), 0 ),
  _solver(solver),
  _solverPatch(solverPatch),
  _cellInfo(cellInfo),
  _isInitialMeshRefinement(isInitialMeshRefinement)
{
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    NumberOfAMRBackgroundJobs++;
  }
  lock.free();
}

bool exahype::solvers::LimitingADERDGSolver::AdjustSolutionDuringMeshRefinementJob::run() {
  _solver._solver->ensureNecessaryMemoryIsAllocated(_solverPatch);
  _solver.adjustSolutionDuringMeshRefinementBody(_solverPatch,_cellInfo,_isInitialMeshRefinement);

  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    NumberOfAMRBackgroundJobs--;
    assertion( NumberOfAMRBackgroundJobs>=0 );
  }
  lock.free();
  return false;
}
