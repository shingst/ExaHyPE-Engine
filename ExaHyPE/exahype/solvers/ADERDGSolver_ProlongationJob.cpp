#include "ADERDGSolver.h"

#if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
#include <immintrin.h>
#endif

exahype::solvers::ADERDGSolver::ProlongationJob::ProlongationJob(
  ADERDGSolver&     solver,
  CellDescription& cellDescription,
  const CellDescription& parentCellDescription,
  const tarch::la::Vector<DIMENSIONS,int>& subcellIndex):
  tarch::multicore::jobs::Job(Solver::getTaskType(true),0),
  _solver(solver),
  _cellDescription(cellDescription),
  _parentCellDescription(parentCellDescription),
  _subcellIndex(subcellIndex) {
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    NumberOfEnclaveJobs++; // TODO(Dominic): Not sure yet which queue is optimal
  }
  lock.free();
}

bool exahype::solvers::ADERDGSolver::ProlongationJob::run() {
  _solver.prolongateFaceDataToDescendant(
      _cellDescription,_parentCellDescription,_subcellIndex);

  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    NumberOfEnclaveJobs--;
    assertion( NumberOfEnclaveJobs>=0 );
  }
  lock.free();
  return false;
}

void exahype::solvers::ADERDGSolver::ProlongationJob::prefetchData() {
  #if defined(SharedTBB) && !defined(noTBBPrefetchesJobData)
  const int dataPerFace = _solver.getBndFaceSize();
  const int dofsPerFace = _solver.getBndFluxSize();

  bool prefetchedOne = false;
  for (int faceIndex=0; faceIndex<DIMENSIONS_TIMES_TWO; faceIndex++) {
    if ( !prefetchedOne && _cellDescription.getFacewiseCommunicationStatus(faceIndex)>0 ) { // prefetch arrays of first face where we need to prolongate
      double* lQhbnd = static_cast<double*>(_parentCellDescription.getExtrapolatedPredictor()) + faceIndex * dataPerFace; 
      double* lFhbnd = static_cast<double*>(_parentCellDescription.getFluctuation())           + faceIndex * dofsPerFace;
      _mm_prefetch(lQhbnd, _MM_HINT_NTA);
      _mm_prefetch(lFhbnd, _MM_HINT_NTA);
      prefetchedOne = false;
    }
  }
  #endif
}
