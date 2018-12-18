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
  NumberOfAMRBackgroundJobs++;
}

bool exahype::solvers::ADERDGSolver::AdjustSolutionDuringMeshRefinementJob::run() {
  _solver.ensureNecessaryMemoryIsAllocated(_cellDescription);
  _solver.adjustSolutionDuringMeshRefinementBody(_cellDescription,_isInitialMeshRefinement);
  NumberOfAMRBackgroundJobs--;
  assertion( NumberOfAMRBackgroundJobs>=0 );
  return false;
}
