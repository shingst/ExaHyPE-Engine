#if defined(Parallel)
#include "exahype/reactive/ResilienceStatistics.h"
#include "exahype/solvers/LimitingADERDGSolver.h"

#include "exahype/solvers/OutcomeDatabase.h"

#include "exahype/reactive/TimeStampAndDubiosityTeamHistory.h"
#include "exahype/reactive/ResilienceTools.h"

#if defined(FileTrace)
#include "exahype/reactive/STPStatsTracer.h"
#include <chrono>
#endif


exahype::solvers::LimitingADERDGSolver::CheckAndCorrectSolutionJob::CheckAndCorrectSolutionJob(
  LimitingADERDGSolver&                                      solver,
  SolverPatch&                                       solverPatch,
  CellInfo&                                          cellInfo,
  const double                                       predictionTimeStamp,
  const double                                       predictionTimeStepSize):
  tarch::multicore::jobs::Job(
      tarch::multicore::jobs::JobType::BackgroundTask,0,
      getTaskPriority(true) //high priority
  ),
  _solver(solver),
  _solverPatch(solverPatch),
  _cellInfo(cellInfo),
  _predictorTimeStamp  (predictionTimeStamp),
  _predictorTimeStepSize(predictionTimeStepSize) {
  NumberOfReductionJobs.fetch_add(1);
}

exahype::solvers::LimitingADERDGSolver::CheckAndCorrectSolutionJob::~CheckAndCorrectSolutionJob() {
  NumberOfReductionJobs.fetch_sub(1);
}

bool exahype::solvers::LimitingADERDGSolver::CheckAndCorrectSolutionJob::run(bool runOnMasterThread) {
#if defined(SharedTBB)
  DeliveryStatus status;
  ADERDGSolver::MigratablePredictionJobData *outcome = nullptr;

  logInfo("run","CheckAndCorrectSolution looking for outcome time stamp= "<< _solverPatch.getTimeStamp()
      <<" previous time step size = "<< _solverPatch.getPreviousTimeStepSize()
      <<" time step size = "<<std::setprecision(30)<< _solverPatch.getTimeStepSize()
      <<" patch "<<_solverPatch.toString() );

  assertion(exahype::reactive::TimeStampAndDubiosityTeamHistory::getInstance().otherTeamHasTimeStepData( _solverPatch.getTimeStamp(),
                                                                                                    _solverPatch.getTimeStepSize()));

  bool found = _solver._solver.get()->tryToFindAndExtractOutcome(_cellInfo._cellDescriptionsIndex, 0 /*is this always correct?*/,
                         _solverPatch.getTimeStamp(), _solverPatch.getTimeStepSize(), status, &outcome);

  if(!found) return true; //reschedule
  else {
    logInfo("run","CheckAndCorrectSolution found outcome time stamp= "<< _solverPatch.getTimeStamp()
        <<" time step size = "<< _solverPatch.getTimeStepSize()
        <<" patch "<<_solverPatch.toString() );

    SDCCheckResult corruption_result = checkAgainstOutcome(outcome);
    UpdateResult update_result;
    switch(corruption_result) {
      case SDCCheckResult::NoCorruption: //run limiter

        update_result._timeStepSize    = _solver.startNewTimeStep(_solverPatch,_cellInfo,true /*isFirstTimeStepOfBatch*/);
        update_result._meshUpdateEvent = _solver.updateRefinementStatusAfterSolutionUpdate(_solverPatch,_cellInfo,true /*troubled*/);
        _solver.reduce(_solverPatch,_cellInfo,update_result);

        _solverPatch.setHasCompletedLastStep(true);
        break;
      case SDCCheckResult::UncorrectableSoftError:  //run limiter
        logError("run()","Soft error detected in cell "<<_solverPatch.toString());
        logError("run()","Cannot correct this error as we don't know which result is the correct one! We run limiter next...");
        exahype::reactive::ResilienceTools::getInstance().setCorruptionDetected(true);

        update_result._timeStepSize    = _solver.startNewTimeStep(_solverPatch,_cellInfo,true /*isFirstTimeStepOfBatch*/);
        update_result._meshUpdateEvent = _solver.updateRefinementStatusAfterSolutionUpdate(_solverPatch,_cellInfo,true /*troubled*/);
        _solver.reduce(_solverPatch,_cellInfo,update_result);

        _solverPatch.setHasCompletedLastStep(true);
        break;
      case SDCCheckResult::OutcomeSaneAsLimiterNotActive:
        logError("run()","Soft error detected in cell "<<_solverPatch.toString());
        // might wanna start healing iterations with solution from sane team
        if(exahype::reactive::ReactiveContext::getInstance().getResilienceStrategy()
            ==exahype::reactive::ReactiveContext::ResilienceStrategy::TaskSharingResilienceCorrection) {
           correctWithOutcomeAndDeleteLimiterStatus(outcome);
           _solver._solver.get()->predictionAndVolumeIntegralBody(
              _solverPatch,
              _predictorTimeStamp,
              _predictorTimeStepSize,
              false/*is uncompressed*/,true,true/*addVolumeIntegralResultToUpdate*/);  //sets has completed!
           //need to insert new prediction outcome for follow-up time step
           _solver._solver.get()->storePendingOutcomeToBeShared(_cellInfo._cellDescriptionsIndex, 0, _predictorTimeStamp, _predictorTimeStepSize);
        }
        else {
          logError("run()","Correction is not activated, so we let the limiter do the job.");
          update_result._timeStepSize    = _solver.startNewTimeStep(_solverPatch,_cellInfo,true /*isFirstTimeStepOfBatch*/);
          update_result._meshUpdateEvent = _solver.updateRefinementStatusAfterSolutionUpdate(_solverPatch,_cellInfo,true /*troubled*/);
          _solver.reduce(_solverPatch,_cellInfo,update_result);
          //skip STP computation
          _solverPatch.setHasCompletedLastStep(true);
        }

        exahype::reactive::ResilienceTools::getInstance().setCorruptionDetected(true);
        break;
    }
    return false;
  }
#else
  assert(false); //todo(Philipp): make it work without TBB, this should not be called
#endif
}

exahype::solvers::LimitingADERDGSolver::CheckAndCorrectSolutionJob::SDCCheckResult exahype::solvers::LimitingADERDGSolver::CheckAndCorrectSolutionJob::checkAgainstOutcome(ADERDGSolver::MigratablePredictionJobData *outcome) {

  double *luh = static_cast<double*>(_solverPatch.getSolution());
  assertion(outcome!=nullptr);

  exahype::reactive::ResilienceStatistics::getInstance().notifyDoubleCheckedTask();

  if(exahype::reactive::ResilienceTools::getInstance().isTrustworthy(outcome->_metadata._errorIndicatorDerivative,
                                                                     outcome->_metadata._errorIndicatorTimeStepSize,
                                                                     outcome->_metadata._errorIndicatorAdmissibility))
  {
    /*logError("checkAgainstOutcome","Limiter was not active previously in other team, we must have a soft error!"
        <<_cellInfo._cellDescriptionsIndex
        <<" stamp "<<_predictorTimeStamp
        <<" step "<<_predictorTimeStepSize
        <<" previous stamp "<<_solverPatch.getPreviousTimeStamp()
        <<" previous step "<<_solverPatch.getPreviousTimeStepSize());*/
     exahype::reactive::ResilienceStatistics::getInstance().notifyDetectedError();
     return SDCCheckResult::OutcomeSaneAsLimiterNotActive;
  }
  else {
    bool equalSolution = exahype::reactive::ResilienceTools::getInstance().isEqual(outcome->_luh.data(), luh, outcome->_luh.size());
    if(!equalSolution) {
      exahype::reactive::ResilienceStatistics::getInstance().notifyDetectedError();
      return SDCCheckResult::UncorrectableSoftError;
    }
    else
      return SDCCheckResult::NoCorruption;
  }
}

void exahype::solvers::LimitingADERDGSolver::CheckAndCorrectSolutionJob::correctWithOutcomeAndDeleteLimiterStatus(ADERDGSolver::MigratablePredictionJobData *outcome){

  //todo: check other data, too?
  double *luh = static_cast<double*>(_solverPatch.getSolution());
  //double *lduh = static_cast<double*>(_solverPatch.getUpdate());
  //double *lQhbnd = static_cast<double*>(_solverPatch.getExtrapolatedPredictor());
  //double *lGradQhbnd = static_cast<double*>(_solverPatch.getExtrapolatedPredictorGradient());
  //double *lFhbnd = static_cast<double*>(_solverPatch.getFluctuation());

  //correct here
  std::memcpy(luh, &outcome->_luh[0], outcome->_luh.size() * sizeof(double));
//  std::memcpy(lduh, &outcome->_lduh[0], outcome->_lduh.size() * sizeof(double));
//  std::memcpy(lQhbnd, &outcome->_lQhbnd[0], outcome->_lQhbnd.size() * sizeof(double));
//  std::memcpy(lFhbnd, &outcome->_lFhbnd[0], outcome->_lFhbnd.size() * sizeof(double));
//#if OffloadingGradQhbnd
//  std::memcpy(lGradQhbnd, &outcome->_lGradQhbnd[0], outcome->_lGradQhbnd.size() * sizeof(double));
//#endif

  logError("correctWithOutcomeAndDeleteLimiterStatus()","We corrected an error and disable limiter for this cell (compute again with ADER-DG).");

  _solverPatch.setRefinementStatus(ADERDGSolver::Keep);

  UpdateResult result;
  result._timeStepSize    = _solver.startNewTimeStep(_solverPatch,_cellInfo,true /*isFirstTimeStepOfBatch*/);
  result._meshUpdateEvent = _solver.updateRefinementStatusAfterSolutionUpdate(_solverPatch,_cellInfo,false /*troubled*/);

  logError("correctWithOutcomeAndDeleteLimiterStatus()","New time step size ="<<result._timeStepSize);

  _solver.reduce(_solverPatch,_cellInfo,result);

  /*if(_solver.getMeshUpdateEvent()==ADERDGSolver::MeshUpdateEvent::IrregularLimiterDomainChangeButMayCorrect) {
    _solver.resetMeshUpdateEvent();
    logInfo("correctWithOutcomeAndDeleteLimiterStatus", "team = "<<exahype::reactive::OffloadingContext::getInstance().getTMPIInterTeamRank()
          <<" corrected patch "<<_solverPatch.toString());
  }*/

  exahype::reactive::ResilienceStatistics::getInstance().notifyHealedTask();
}
#endif
