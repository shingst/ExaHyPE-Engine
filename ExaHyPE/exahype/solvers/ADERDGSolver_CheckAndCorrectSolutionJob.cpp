#if defined(Parallel) && defined(SharedTBB)
#include "ADERDGSolver.h"

#include "exahype/reactive/ResilienceTools.h"
#include "exahype/reactive/TimeStampAndTriggerTeamHistory.h"


exahype::solvers::ADERDGSolver::CheckAndCorrectSolutionJob::CheckAndCorrectSolutionJob(
    ADERDGSolver & solver,
    CellDescription& solverPatch,
    double predictorTimeStamp,
    double predictorTimeStepSize)
 : tarch::multicore::jobs::Job(
     tarch::multicore::jobs::JobType::BackgroundTask, 0,
     getTaskPriority(false)),
   _solver(solver),
   _solverPatch(solverPatch),
   _predictorTimeStamp(predictorTimeStamp),
   _predictorTimeStepSize(predictorTimeStepSize) {

}

bool exahype::solvers::ADERDGSolver::CheckAndCorrectSolutionJob::run(bool isRunOnMaster) {
  //caution, this may delay setting the completed flag and cause slowdowns for the master!
#if !defined(OffloadingUseProgressThread)
  exahype::solvers::ADERDGSolver::progressOffloading(&_solver, false, MAX_PROGRESS_ITS);
#endif
  bool reschedule;

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
        logError("runCheck"," Detected a soft error but I am continuing...");
        break;
        //MPI_Abort(MPI_COMM_WORLD, -1);
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

exahype::solvers::ADERDGSolver::CheckAndCorrectSolutionJob::SDCCheckResult exahype::solvers::ADERDGSolver::CheckAndCorrectSolutionJob::checkAgainstOutcome(ADERDGSolver::MigratablePredictionJobData *data) {

  double *lduh = static_cast<double*>(_solverPatch.getUpdate());
  double *lQhbnd = static_cast<double*>(_solverPatch.getExtrapolatedPredictor());
#if defined(OffloadingGradQhbnd)
  double *lGradQhbnd = static_cast<double*>(_solverPatch.getExtrapolatedPredictorGradient());
#endif
  double *lFhbnd = static_cast<double*>(_solverPatch.getFluctuation());

  tarch::la::Vector<DIMENSIONS, double> center = _solverPatch.getOffset()+0.5*_solverPatch.getSize();

  bool equal = true;
  bool tmp;

  tmp = data->_metadata._predictorTimeStamp == _predictorTimeStamp; equal &= tmp;
  tmp = data->_metadata._predictorTimeStepSize == _predictorTimeStepSize; equal &= tmp;

  tmp = exahype::reactive::ResilienceTools::getInstance().isAdmissibleNumericalError(data->_lQhbnd.data(), lQhbnd, data->_lQhbnd.size()); equal&=tmp;
  if(!tmp) {
    logError("checkAgainstOutcome", "lQhbnd is not (numerically) equal for cell "<<"center[0]="<<_center[0]<<" center[1]="<<_center[1]<<" timestamp "<<_predictorTimeStamp);
  }

#if defined(OffloadingGradQhbnd)
  tmp = exahype::reactive::ResilienceTools::getInstance().isAdmissibleNumericalError(data->_lGradQhbnd.data(), lGradQhbnd, data->_lGradQhbnd.size()); equal&=tmp;
  if(!tmp) {
    logError("checkAgainstOutcome", "lGradQhbnd is not (numerically) equal for cell "<<"center[0]="<<_center[0]<<" center[1]="<<_center[1]<<" timestamp "<<_predictorTimeStamp);
  }
#endif
  tmp = exahype::reactive::ResilienceTools::getInstance().isAdmissibleNumericalError(data->_lFhbnd.data(), lFhbnd, data->_lFhbnd.size()); equal&=tmp;
  if(!tmp) {
    logError("checkAgainstOutcome", "lFhbnd is not  (numerically) equal for cell "<<"center[0]="<<_center[0]<<" center[1]="<<_center[1]<<" timestamp "<<_predictorTimeStamp);
  }
  tmp = exahype::reactive::ResilienceTools::getInstance().isAdmissibleNumericalError(data->_lduh.data(), lduh, data->_lduh.size()); equal&=tmp;
  if(!tmp) {
    logError("checkAgainstOutcome", "lduh is not  (numerically) equal for cell "<<"center[0]="<<_center[0]<<" center[1]="<<_center[1]<<" timestamp "<<_predictorTimeStamp);
  }

  if(!equal) {
    logError("checkAgainstOutcome", "soft error detected: "<<data->_metadata.to_string());
    exahype::reactive::ResilienceStatistics::getInstance().notifyDetectedError();
  }

  exahype::reactive::ResilienceStatistics::getInstance().notifyDoubleCheckedTask();

  if(!equal) {
    return SDCCheckResult::UncorrectableSoftError;
  }
  else
    return SDCCheckResult::NoCorruption;
}


#endif

