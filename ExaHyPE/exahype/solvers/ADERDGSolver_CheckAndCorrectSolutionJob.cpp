#if defined(Parallel) && defined(SharedTBB)
#include "ADERDGSolver.h"

#include "exahype/reactive/ResilienceTools.h"
#include "exahype/reactive/TimeStampAndDubiosityTeamHistory.h"


exahype::solvers::ADERDGSolver::CheckAndCorrectSolutionJob::CheckAndCorrectSolutionJob(
    ADERDGSolver & solver,
    CellDescription& solverPatch,
    double predictorTimeStamp,
    double predictorTimeStepSize,
    double errorIndicatorDerivative,
    double errorIndicatorTimeStepSize,
    double errorIndicatorAdmissibility,
    bool isSkeleton)
 : tarch::multicore::jobs::Job(
     tarch::multicore::jobs::JobType::BackgroundTask, 0,
     getTaskPriorityCheckCorrectJob(false)/2),
   _solver(solver),
   _solverPatch(solverPatch),
   _predictorTimeStamp(predictorTimeStamp),
   _predictorTimeStepSize(predictorTimeStepSize),
   _errorIndicatorDerivative(errorIndicatorDerivative),
   _errorIndicatorTimeStepSize(errorIndicatorTimeStepSize),
   _errorIndicatorAdmissibility(errorIndicatorAdmissibility),
   _startTimeStamp(MPI_Wtime()),
   _printedTimeoutWarning(false){
   NumberOfCheckJobs.fetch_add(1);
}

bool exahype::solvers::ADERDGSolver::CheckAndCorrectSolutionJob::run(bool isRunOnMaster) {
  //caution, this may delay setting the completed flag and cause slowdowns for the master!
#if !defined(OffloadingUseProgressThread)
  if(!isRunOnMaster)
    exahype::solvers::ADERDGSolver::progressOffloading(&_solver, false, MAX_PROGRESS_ITS);
#endif
  bool reschedule = true;

  DeliveryStatus status;
  MigratablePredictionJobData *outcome;
  bool found = tryToFindAndExtractEquivalentSharedOutcome(status, &outcome);
        
  logDebug("runCheck", std::setprecision(30)<<
                      " Checking my outcome: predictorTimeStamp = "<<_predictorTimeStamp<<
                      " predictor time step size = "<<_predictorTimeStepSize<<
                      " My error indicators: derivative = "<<_errorIndicatorDerivative<<
                      " time step size = "<<_errorIndicatorTimeStepSize<<
                      " admissibility = "<<_errorIndicatorAdmissibility);

  bool timeout = (MPI_Wtime()-_startTimeStamp) > exahype::reactive::ReactiveContext::getInstance().getResilienceChecksTimeout();

  if(timeout && !_printedTimeoutWarning) {
    _solverPatch.setHasCompletedLastStep(true);
    reschedule = false;
    //exahype::reactive::ResilienceStatistics::getInstance().notifyDetectedError();
    logWarning("runCheck","Waiting too long for an outcome. You may either have a progression problem or a silent error which could not be corrected so time step sizes between teams may have diverged. Won't check anymore! sent ="<<SentSTPs<< " received "<<ReceivedSTPs<< " CheckJobs "<<NumberOfCheckJobs<< " NumberOfEnclaveJobs "<<NumberOfEnclaveJobs );
    tarch::la::Vector<DIMENSIONS, double> center;
    center = _solverPatch.getOffset() + 0.5 *_solverPatch.getSize();

    logWarning("runCheck", "team = "<<exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber()<<" looking for center[0] = "<<center[0]
                                         <<" center[1] = "<<center[1]
                                         <<" center[2] = "<<center[2]
                                         <<" timestamp = "<<_predictorTimeStamp
                                         <<" time step = "<<_predictorTimeStepSize);
    _printedTimeoutWarning = true;
    return reschedule;
  }

  //assert(exahype::reactive::TimeStampAndTriggerTeamHistory::getInstance().otherTeamHasTimeStepData(_predictorTimeStamp, _predictorTimeStepSize));
  if( !found
    && !exahype::reactive::TimeStampAndDubiosityTeamHistory::getInstance().hasConsistentTeamTimeStampHistories())
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
      case SDCCheckResult::MyOutcomeIsMoreOrEquallyTrustworthy:
        reschedule = false;
        _solverPatch.setHasCompletedLastStep(true);
        logError("runCheck","My result is more trustworthy. Continuing but the error remains which may result in undefined behavior...");
        break;
      case SDCCheckResult::OutcomeIsMoreTrustworthy:
         //continue to correction
        logError("runCheck", "I have an outcome that has lower error indicator!"<<
                              "My error indicators: derivative = "<<_errorIndicatorDerivative<<
                              " time step size = "<<_errorIndicatorTimeStepSize<<
                              " admissibility = "<<_errorIndicatorAdmissibility<<
			                        " other outcome's error indicators: derivative = "<<outcome->_metadata._errorIndicatorDerivative<<
                              " time step size "<<outcome->_metadata._errorIndicatorTimeStepSize<<
                              " admissibility "<<outcome->_metadata._errorIndicatorAdmissibility);
         //no break!
      case SDCCheckResult::OutcomeHasZeroErrorIndicator:
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
          exahype::reactive::ResilienceTools::getInstance().setCorruptionDetected(true);
        }
        break;
      }
      NumberOfCheckJobs.fetch_add(-1);
      delete outcome;
    }
    else {
      logDebug("runCheck"," Not found "<<std::setprecision(30)<<"time stamp ="<<_predictorTimeStamp<<" timestep "<<_predictorTimeStepSize);
      reschedule = true;
    }
  }
  return reschedule;
}

bool exahype::solvers::ADERDGSolver::CheckAndCorrectSolutionJob::tryToFindAndExtractEquivalentSharedOutcome(DeliveryStatus &status, MigratablePredictionJobData **outcome) {
  return _solver.tryToFindAndExtractOutcome(_solverPatch, _predictorTimeStamp, _predictorTimeStepSize, status, outcome);
}

exahype::solvers::ADERDGSolver::SDCCheckResult exahype::solvers::ADERDGSolver::CheckAndCorrectSolutionJob::checkAgainstOutcome(ADERDGSolver::MigratablePredictionJobData *data) {
  return _solver.checkCellDescriptionAgainstOutcome(_solverPatch, data, _predictorTimeStamp, _predictorTimeStepSize,
                                                    _errorIndicatorDerivative,
                                                    _errorIndicatorTimeStepSize,
                                                    _errorIndicatorAdmissibility);
}

void exahype::solvers::ADERDGSolver::CheckAndCorrectSolutionJob::correctWithOutcome(ADERDGSolver::MigratablePredictionJobData *outcome){
   _solver.correctCellDescriptionWithOutcome(_solverPatch, outcome);
}
#endif

