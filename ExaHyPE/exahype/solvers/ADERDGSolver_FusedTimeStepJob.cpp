#include "ADERDGSolver.h"
//#include <immintrin.h>


exahype::solvers::ADERDGSolver::FusedTimeStepJob::FusedTimeStepJob(
  ADERDGSolver&    solver,
  CellDescription& cellDescription,
  const int        cellDescriptionsIndex,
  const int        element,
  const bool       isSkeletonJob):
  tarch::multicore::jobs::Job(Solver::getTaskType(isSkeletonJob),0),
  _solver(solver),
  _cellDescription(cellDescription),
  _cellDescriptionsIndex(cellDescriptionsIndex),
  _element(element),
  _neighbourMergePerformed(cellDescription.getNeighbourMergePerformed()),
  _isSkeletonJob(isSkeletonJob) {
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter++;
  }
  lock.free();
}

bool exahype::solvers::ADERDGSolver::FusedTimeStepJob::run() {
  _solver.fusedTimeStepBody(
      _cellDescription, _cellDescriptionsIndex, _element, false, false, _isSkeletonJob, false, _neighbourMergePerformed );

  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter--;
    assertion( jobCounter>=0 );
  }
  lock.free();
  return false;
}


void exahype::solvers::ADERDGSolver::FusedTimeStepJob::prefetchData() {
//  std::cout << "was here (a)" << std::endl;
}
