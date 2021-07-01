#if defined(Parallel) && defined(SharedTBB)
#include "ADERDGSolver.h"

#include "exahype/reactive/ResilienceTools.h"
#include "exahype/reactive/TimeStampAndTriggerTeamHistory.h"


exahype::solvers::ADERDGSolver::CheckAndCorrectSolutionJob::CheckAndCorrectSolutionJob(
    ADERDGSolver & solver,
    CellDescription& solverPatch,
    double predictorTimeStamp,
    double predictorTimeStepSize,
    double confidence)
 : tarch::multicore::jobs::Job(
     tarch::multicore::jobs::JobType::BackgroundTask, 0,
     getTaskPriority(false)),
   _solver(solver),
   _solverPatch(solverPatch),
   _predictorTimeStamp(predictorTimeStamp),
   _predictorTimeStepSize(predictorTimeStepSize),
   _confidence (confidence){

}

bool exahype::solvers::ADERDGSolver::CheckAndCorrectSolutionJob::run(bool isRunOnMaster) {
  //caution, this may delay setting the completed flag and cause slowdowns for the master!
#if !defined(OffloadingUseProgressThread)
  exahype::solvers::ADERDGSolver::progressOffloading(&_solver, false, MAX_PROGRESS_ITS);
#endif
  bool reschedule = true;

  DeliveryStatus status;
  MigratablePredictionJobData *outcome;
  bool found = tryToFindAndExtractEquivalentSharedOutcome(status, &outcome);

  //assert(exahype::reactive::TimeStampAndTriggerTeamHistory::getInstance().otherTeamHasTimeStepData(_predictorTimeStamp, _predictorTimeStepSize));
  if( !found
    && (
        !exahype::reactive::TimeStampAndTriggerTeamHistory::getInstance().otherTeamHasTimeStepData(_predictorTimeStamp, _predictorTimeStepSize)
        &&
        exahype::reactive::TimeStampAndTriggerTeamHistory::getInstance().otherTeamHasLargerTimeStamp(_predictorTimeStamp)))
      // || !exahype::reactive::TimeStampAndTriggerTeamHistory::getInstance().checkConsistency())
  {
    _solverPatch.setHasCompletedLastStep(true);
    reschedule = false;
    exahype::reactive::ResilienceStatistics::getInstance().notifyDetectedError();
  }
  else {

    if(found && status==DeliveryStatus::Received) {
      switch(checkAgainstOutcome(outcome))
      {
      case SDCCheckResult::NoCorruption:
        _solverPatch.setHasCompletedLastStep(true);
        reschedule = false;
        break;
      case SDCCheckResult::UncorrectableSoftError:
        //soft error detected
        reschedule = false;
        _solverPatch.setHasCompletedLastStep(true);
        logError("runCheck"," Detected a soft error and I don't know which outcome is correct. Continuing...");
        break;
      case SDCCheckResult::OutcomeHasHigherConfidence:
         //continue to correction
        logError("runCheck", "I'm correcting with an outcome that has higher confidence!"
                              "My confidence = "<<_confidence<<
			      " other outcome's confidence = "<<outcome->_metadata._confidence);
         //no break!
      case SDCCheckResult::OutcomeSaneAsTriggerNotActive:
        if(exahype::reactive::ReactiveContext::getInstance().getResilienceStrategy()
          ==exahype::reactive::ReactiveContext::ResilienceStrategy::TaskSharingResilienceCorrection) {
          correctWithOutcome(outcome);
          logError("runCheck","Corrected soft error with outcome from other team. Continuing...");
          reschedule = false;
          _solverPatch.setHasCompletedLastStep(true);
        }
        else {
          logError("runCheck","I could probably correct an error, but correction has not been activated. Continuing...");
          reschedule = false;
          _solverPatch.setHasCompletedLastStep(true);
        }
      }
    }
    else {
      logDebug("runCheck"," Not found "<<to_string()<<std::setprecision(30)<<"time stamp ="<<_predictorTimeStamp<<" timestep "<<_predictorTimeStepSize);
      reschedule = true;
    }
  }

  return reschedule;
}

bool exahype::solvers::ADERDGSolver::CheckAndCorrectSolutionJob::tryToFindAndExtractEquivalentSharedOutcome(DeliveryStatus &status, MigratablePredictionJobData **outcome) {

  logDebug("tryToFindAndExtractEquivalentSharedOutcome", "team = "<<exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber()
      <<" looking for "<<_cellDescription.toString())

  return _solver.tryToFindAndExtractOutcome(_solverPatch, _predictorTimeStamp, _predictorTimeStepSize, status, outcome);
}

exahype::solvers::ADERDGSolver::SDCCheckResult exahype::solvers::ADERDGSolver::CheckAndCorrectSolutionJob::checkAgainstOutcome(ADERDGSolver::MigratablePredictionJobData *data) {
  return _solver.checkCellDescriptionAgainstOutcome(_solverPatch, data, _predictorTimeStamp, _predictorTimeStepSize, _confidence);
}

void exahype::solvers::ADERDGSolver::CheckAndCorrectSolutionJob::correctWithOutcome(ADERDGSolver::MigratablePredictionJobData *outcome){

   _solver.correctCellDescriptionWithOutcome(_solverPatch, outcome);
}
#endif

