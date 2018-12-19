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
  NumberOfAMRBackgroundJobs++;
}

bool exahype::solvers::FiniteVolumesSolver::AdjustSolutionDuringMeshRefinementJob::run() {
  _solver.adjustSolutionDuringMeshRefinementBody(_cellDescription,_isInitialMeshRefinement);

  NumberOfAMRBackgroundJobs--;
  assertion( NumberOfAMRBackgroundJobs>=0 );
  return false;
}
