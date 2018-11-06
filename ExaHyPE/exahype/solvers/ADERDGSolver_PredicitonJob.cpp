#include "ADERDGSolver.h"
//#include <immintrin.h>

exahype::solvers::ADERDGSolver::PredictionJob::PredictionJob(
  ADERDGSolver&    solver,
  CellDescription& cellDescription,
  const int        cellDescriptionsIndex,
  const int        element,
  const double     predictorTimeStamp,
  const double     predictorTimeStepSize,
  const bool       uncompressBefore,
  const bool       isSkeletonJob):
  tarch::multicore::jobs::Job(Solver::getTaskType(isSkeletonJob),0),
  _solver(solver),
  _cellDescription(cellDescription),
  _cellDescriptionsIndex(cellDescriptionsIndex),
  _element(element),
  _predictorTimeStamp(predictorTimeStamp),
  _predictorTimeStepSize(predictorTimeStepSize),
  _uncompressBefore(uncompressBefore),
  _isSkeletonJob(isSkeletonJob) {
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter++;
  }
  lock.free();
}


bool exahype::solvers::ADERDGSolver::PredictionJob::run() {
  _solver.performPredictionAndVolumeIntegralBody(
      _cellDescription,_predictorTimeStamp,_predictorTimeStepSize,
      _uncompressBefore,_isSkeletonJob); // ignore return value
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter--;
    assertion( jobCounter>=0 );
  }
  lock.free();
  return false;
}


/*
void exahype::solvers::ADERDGSolver::PredictionJob::prefetchData() {
  std::cout << "was here (b)" << std::endl;
}
*/
