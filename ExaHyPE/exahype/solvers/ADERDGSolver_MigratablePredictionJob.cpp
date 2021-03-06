#if defined(Parallel) && defined(SharedTBB)
#include "ADERDGSolver.h"

#if defined(ScoreP)
#include "scorep/SCOREP_User.h"
#endif

#if defined(FileTrace)
#include "exahype/reactive/STPStatsTracer.h"
#endif

#include "exahype/reactive/PerformanceMonitor.h"
#include "exahype/reactive/OffloadingAnalyser.h"
#include "exahype/reactive/OffloadingProfiler.h"
#include "exahype/reactive/ResilienceStatistics.h"
#include "exahype/reactive/MemoryMonitor.h"
#include "exahype/reactive/NoiseGenerator.h"
#include "exahype/reactive/ResilienceTools.h"
#include "exahype/reactive/TimeStampAndDubiosityTeamHistory.h"
#include "exahype/reactive/RequestManager.h"

#include "tarch/multicore/Core.h"

#if defined(USE_TMPI)
#include "teaMPI.h"
#endif

//#undef assertion
//#define assertion assert

MPI_Datatype exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::_datatype;

exahype::solvers::ADERDGSolver::MigratablePredictionJob::MigratablePredictionJob(
    ADERDGSolver& solver,
    const int cellDescriptionsIndex,
    const int element,
    const double predictorTimeStamp,
    const double predictorTimeStepSize,
    const bool isSkeleton) :
      tarch::multicore::jobs::Job(
          tarch::multicore::jobs::JobType::BackgroundTask, 0,
          getTaskPriorityLocalMigratableJob(cellDescriptionsIndex, element,
              predictorTimeStamp, isSkeleton)),  //this is a locally executed job
      _solver(solver),
      _cellDescriptionsIndex(cellDescriptionsIndex),
      _element(element),
      _predictorTimeStamp(predictorTimeStamp),
      _predictorTimeStepSize(predictorTimeStepSize),
      _originRank(tarch::parallel::Node::getInstance().getRank()),
      _tag(-1),
      _isSkeleton(isSkeleton),
      _luh(nullptr),
      _lduh(nullptr),
      _lQhbnd(nullptr),
      _lFhbnd(nullptr),
      _lGradQhbnd(nullptr),
      _isLocalReplica(false),
      _errorIndicatorTimeStepSize(0),
      _errorIndicatorDerivative(0),
      _errorIndicatorAdmissibility(0),
      _isCorrupted(false)
{
  if (_isSkeleton) {
    NumberOfSkeletonJobs.fetch_add(1);
  } else {
    LocalStealableSTPCounter++;
    NumberOfEnclaveJobs.fetch_add(1);
  }
  exahype::reactive::ResilienceStatistics::getInstance().notifySpawnedTask();
  exahype::reactive::PerformanceMonitor::getInstance().incCurrentNumTasks();

  auto& cellDescription = getCellDescription(cellDescriptionsIndex, element);
  logDebug("MigratablePredictionJob","team "<<exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber()<<" spawning STP for "
        <<cellDescription.toString());

  tarch::la::Vector<DIMENSIONS, double> center = cellDescription.getOffset()+0.5*cellDescription.getSize();
  tarch::la::Vector<DIMENSIONS, double> dx = cellDescription.getSize();

  for (int i = 0; i < DIMENSIONS; i++) {
    _center[i] = center[i];
    _dx[i] = dx[i];
  }

}

exahype::solvers::ADERDGSolver::MigratablePredictionJob::MigratablePredictionJob(
    ADERDGSolver& solver,
    const int cellDescriptionsIndex,
    const int element,
    const double predictorTimeStamp,
    const double predictorTimeStepSize,
    double *luh, double *lduh,
    double *lQhbnd, double *lFhbnd,
    double *lGradQhbnd, double *dx,
    double *center, const int originRank,
    const int tag) :
      tarch::multicore::jobs::Job(
        tarch::multicore::jobs::JobType::BackgroundTask, 0,
        tarch::multicore::DefaultPriority * 4), //this is a remotely executed job -> high priority
      _solver(solver),
      _cellDescriptionsIndex(cellDescriptionsIndex),
      _element(element),
      _predictorTimeStamp(predictorTimeStamp),
      _predictorTimeStepSize(predictorTimeStepSize),
      _originRank(originRank),
      _tag(tag),
      _isSkeleton(false),  //this constructor should never be invoked for skeletons!
      _luh(luh),
      _lduh(lduh),
      _lQhbnd(lQhbnd),
      _lFhbnd(lFhbnd),
      _lGradQhbnd(lGradQhbnd),
      _isLocalReplica(false),
      _errorIndicatorTimeStepSize(0),
      _errorIndicatorDerivative(0),
      _errorIndicatorAdmissibility(0),
      _isCorrupted(false)
{
  assertion(!_isSkeleton);
  assertion(element==0); //todo: Is element still used somewhere? Dominic's code seems to assume it to be zero...

  for (int i = 0; i < DIMENSIONS; i++) {
    _center[i] = center[i];
    _dx[i] = dx[i];
  }

  if (isRemoteJob()) {
    NumberOfStolenJobs++;
  }
  else
    NumberOfEnclaveJobs++;
  exahype::reactive::PerformanceMonitor::getInstance().incCurrentNumTasks();
}

exahype::solvers::ADERDGSolver::MigratablePredictionJob::~MigratablePredictionJob() {
}

//Caution: Compression is not supported yet!
bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::run(
    bool isCalledOnMaster) {

  return runExecution(isCalledOnMaster);
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::runExecution(bool isCalledOnMaster) {
#ifdef USE_ITAC
  VT_begin(event_stp);
#endif

#if defined ScoreP
  SCOREP_USER_REGION( "exahype::solvers::ADERDGSolver::MigratablePredictionJob::run", SCOREP_USER_REGION_TYPE_FUNCTION )
#endif

  bool result = false;
  bool hasComputed = false;
  int curr = std::atomic_fetch_add(&JobCounter, 1);

  if (curr % 1000 == 0) {
    tarch::timing::Watch watch("exahype::MigratablePredictionJob::", "-", false, false);
    watch.startTimer();
    result = handleExecution(isCalledOnMaster, hasComputed);
    watch.stopTimer();

    logDebug("run()","measured time per STP "<<watch.getCalendarTime());
    if (hasComputed) {
      exahype::reactive::OffloadingAnalyser::getInstance().setTimePerSTP(
          watch.getCalendarTime());
    }
  }
  else
    result = handleExecution(isCalledOnMaster, hasComputed);

#ifdef USE_ITAC
  VT_end(event_stp);
#endif
  return result;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::handleExecution(
    bool isRunOnMaster, bool& hasComputed) {

  bool result = false;
#if defined(OffloadingUseProfiler)
  exahype::reactive::OffloadingProfiler::getInstance().beginComputation();
  double time = -MPI_Wtime();
#endif
  //local treatment if this job belongs to the local rank
  if (!isRemoteJob()) {
    result = handleLocalExecution(isRunOnMaster, hasComputed);
    if (!_isLocalReplica) {
      if (_isSkeleton) {
        NumberOfSkeletonJobs.fetch_sub(1);
      } else {
        NumberOfEnclaveJobs.fetch_sub(1);
      }
    }
    assertion( NumberOfEnclaveJobs>=0 );
    assertion( NumberOfSkeletonJobs>=0 );
#if !defined(OffloadingUseProgressThread)
    if (!isRunOnMaster)
      exahype::solvers::ADERDGSolver::progressOffloading(&_solver,
          isRunOnMaster, std::numeric_limits<int>::max());
#endif
  }
  //remote task, need to execute and send it back
  else {
    result = handleRemoteExecution(hasComputed);
    sendBackOutcomeToOrigin();
  }
#if defined(OffloadingUseProfiler)
  time += MPI_Wtime();
  exahype::reactive::OffloadingProfiler::getInstance().endComputation(time);
#endif
  return result;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::handleLocalExecution(bool isRunOnMaster, bool& hasComputed) {
  //fast path
  if(exahype::reactive::ReactiveContext::getInstance().getResilienceStrategy()
      ==exahype::reactive::ReactiveContext::ResilienceStrategy::None) {
    executeLocally();
    hasComputed = true;
    corruptIfActive(false);
    setFinished();
    return false;
  }

  //got some task sharing here
  bool hasFlipped = false;
  //check for outcome and get it + corruption status
  MigratablePredictionJobData *outcome = nullptr;
  DeliveryStatus status;
  if(tryToFindAndExtractEquivalentSharedOutcome(false, status, &outcome)) {
    executeOrCopySTPOutcome(outcome, hasComputed, hasFlipped); //corruption can also happen here if executed
    computeErrorIndicators(hasFlipped);

    if(needToCheckThisSTP(hasComputed)) {
      checkAgainstOutcomeAndCorrectIfActive(outcome);
    }
    if(needToShare(true /*hasOutcome*/,
                   exahype::reactive::ResilienceTools::getInstance().isTrustworthy(
                   _errorIndicatorDerivative,
                   _errorIndicatorTimeStepSize,
                   _errorIndicatorAdmissibility),
                  exahype::reactive::ResilienceTools::getInstance().isTrustworthy(
                   outcome->_metadata._errorIndicatorDerivative,
                   outcome->_metadata._errorIndicatorTimeStepSize,
                   outcome->_metadata._errorIndicatorAdmissibility))) {
      shareSTPImmediatelyOrLater();
    }
    if(needToPutBackOutcome()) {
      MigratablePredictionJobOutcomeKey key(outcome->_metadata.getCenter(),
                                            outcome->_metadata.getPredictorTimeStamp(),
                                            outcome->_metadata.getPredictorTimeStepSize(),
                                            outcome->_metadata.getElement());
      _solver._outcomeDatabase.insertOutcome(key, outcome, DeliveryStatus::Received);
    }
    else {
      delete outcome;
    }
    setFinished(); //have already checked if necessary
    return false;
  }
  //don't have an outcome
  else {
    executeLocally();
    hasComputed = true;
    hasFlipped = corruptIfActive(false);
    computeErrorIndicators(hasFlipped);

    if(exahype::reactive::ReactiveContext::getResilienceStrategy()
      >=exahype::reactive::ReactiveContext::ResilienceStrategy::TaskSharingResilienceChecks) {

        exahype::reactive::TimeStampAndDubiosityTeamHistory::getInstance().trackTimeStepAndDubiosity(
                    exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber(), 
		                                                      _predictorTimeStamp,
							                                            _predictorTimeStepSize,
							      !exahype::reactive::ResilienceTools::getInstance().isTrustworthy(
                                                         _errorIndicatorDerivative,
                                                         _errorIndicatorTimeStepSize,
                                                         _errorIndicatorAdmissibility));
    }

    if(needToShare(false, false, false)) {
      shareSTPImmediatelyOrLater();
    }

    if(needToCheckThisSTP(hasComputed)) {
      if(exahype::reactive::TimeStampAndDubiosityTeamHistory::getInstance().hasConsistentTeamTimeStampHistories()) {
        logDebug("handleLocalExecution","going into check mode "<<to_string());

        CheckAndCorrectSolutionJob *job = new CheckAndCorrectSolutionJob(_solver,
                                           getCellDescription(_cellDescriptionsIndex,_element),
                                           _predictorTimeStamp,
                                           _predictorTimeStepSize,
                                           _errorIndicatorDerivative,
                                           _errorIndicatorTimeStepSize,
                                           _errorIndicatorAdmissibility,
                                           _isSkeleton);
        peano::datatraversal::TaskSet spawnedSet(job);
        return false; //check job will take care next
      }
      else {
        logDebug("handleLocalExecution","Won't be able to find STP anymore as timestamps/timestep sizes have diverged. Timestamp ="
            <<std::setprecision(30)<<_predictorTimeStamp
            <<" time step "<<_predictorTimeStepSize);
        exahype::reactive::ResilienceStatistics::getInstance().notifyDetectedError();
        setFinished();
        return false;
      }
    }

    setFinished();
    return false;
  }
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::handleRemoteExecution(bool& hasComputed) {
#if defined(USE_ITAC)
  VT_begin(event_stp_remote);
#endif
#if defined FileTrace
  auto start = std::chrono::high_resolution_clock::now();
#endif
  bool result = false;

  assertion(_lduh!=nullptr);
  assertion(_lQhbnd!=nullptr);
  assertion(_lGradQhbnd!=nullptr);
  assertion(_lFhbnd!=nullptr);
  assertion(_luh!=nullptr);

  logDebug("handleRemoteExecution",
        " processJob: "<<to_string());

  int iterations=_solver.fusedSpaceTimePredictorVolumeIntegral(
     _lduh,
     _lQhbnd,
     _lGradQhbnd,
     _lFhbnd,
     _luh,
     _center,
     _dx,
     _predictorTimeStamp,
     _predictorTimeStepSize,
     true);

  logDebug("handleRemoteExecution",
        " finished job: "<<to_string());

  hasComputed = true;
#if defined(USE_ITAC)
  VT_end(event_stp_remote);
#endif
#if defined(FileTrace)
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

  exahype::reactive::STPStatsTracer::getInstance().writeTracingEventRunIterations(duration.count(), iterations, exahype::reactive::STPTraceKey::ADERDGRemoteMigratable);
#endif
  return result;
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::setFinished() {
  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
   _element);

  cellDescription.setHasCompletedLastStep(true);
  exahype::reactive::PerformanceMonitor::getInstance().decCurrentNumTasks();
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::corruptIfActive(bool forceCorruption) {
  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
   _element);

  bool shouldCorrupt = false;
  if(!forceCorruption)
    shouldCorrupt = exahype::reactive::ResilienceTools::getInstance().shouldInjectError(_center, _predictorTimeStamp); 
  else 
    shouldCorrupt = true;

  if(shouldCorrupt){
    double *lduh = static_cast<double*>(cellDescription.getUpdate());
    exahype::reactive::ResilienceTools::getInstance().corruptData(
                                                              static_cast<double*>(cellDescription.getSolution()),
                                                              _center,
                                                              DIMENSIONS,
                                                              _predictorTimeStamp,
                                                              lduh,
                                                              _solver.getUpdateSize());
    _isCorrupted = true;
    logError("corruptIfActive","Has corrupted STP job "<<to_string());
  }

  return shouldCorrupt;
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::copyOutcome(MigratablePredictionJobData *outcome) {
  assert(outcome!=nullptr);
  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
   _element);

  double *lduh = static_cast<double*>(cellDescription.getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
#if OffloadingGradQhbnd
  double *lGradQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient());
#endif
  double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

  if(outcome->_metadata._isCorrupted) {
    logError("copyOutcome","Error: we're using a corrupted outcome! Outcome"<<outcome->_metadata.to_string()); //assert(false);
  }

  assertion(outcome!=nullptr);
  std::memcpy(lduh, &outcome->_lduh[0], outcome->_lduh.size() * sizeof(double));
  std::memcpy(lQhbnd, &outcome->_lQhbnd[0], outcome->_lQhbnd.size() * sizeof(double));
  std::memcpy(lFhbnd, &outcome->_lFhbnd[0], outcome->_lFhbnd.size() * sizeof(double));
#if OffloadingGradQhbnd
  std::memcpy(lGradQhbnd, &outcome->_lGradQhbnd[0], outcome->_lGradQhbnd.size() * sizeof(double));
#endif
  logDebug("copyOutcome","team "<<exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber()<<" copied STP outcome for "
        <<to_string());

  exahype::reactive::ResilienceStatistics::getInstance().notifySavedTask();
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::executeOrCopySTPOutcome(MigratablePredictionJobData *outcome, bool& hasComputed, bool& hasFlipped) {
  logDebug("executeOrCopySTPOutcome","team "<<exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber()<<" computing or copying STP for "
        <<to_string()
        <<std::setprecision(30)<<" time step size "<<_predictorTimeStepSize);

  bool shouldInjectError = exahype::reactive::ResilienceTools::getInstance().shouldInjectError(_center, _predictorTimeStamp);
  if( (exahype::reactive::ReactiveContext::getResilienceStrategy()>=exahype::reactive::ReactiveContext::ResilienceStrategy::TaskSharing)
   && exahype::reactive::ReactiveContext::getSaveRedundantComputations()
   && exahype::reactive::ResilienceTools::getInstance().isTrustworthy(
         outcome->_metadata._errorIndicatorDerivative,
         outcome->_metadata._errorIndicatorTimeStepSize,
         outcome->_metadata._errorIndicatorAdmissibility)
   && !shouldInjectError) {
    copyOutcome(outcome);
    hasFlipped = false; //we don't assume any corruption after copying an outcome for now
    hasComputed = false;
  }
  else {
    executeLocally();
    hasFlipped = corruptIfActive(shouldInjectError);
    hasComputed = true;
  }
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::executeLocally() {
  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
   _element);

  double *luh = static_cast<double*>(cellDescription.getSolution());
  double *lduh = static_cast<double*>(cellDescription.getUpdate());
  double *lQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictor());
  double *lGradQhbnd = static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient());
  double *lFhbnd = static_cast<double*>(cellDescription.getFluctuation());

  logDebug("executeLocally","team "<<exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber()<<" computing STP for "
        <<to_string()
        <<std::setprecision(30)<<" time step size "<<_predictorTimeStepSize);

#if defined FileTrace
  auto start = std::chrono::high_resolution_clock::now();
  int iterations = _solver.fusedSpaceTimePredictorVolumeIntegral(lduh,
           lQhbnd,
           lGradQhbnd,
           lFhbnd,
           luh,
           cellDescription.getOffset() + 0.5 * cellDescription.getSize(),
           cellDescription.getSize(),
           _predictorTimeStamp,
           _predictorTimeStepSize,
           true);
#else
  _solver.fusedSpaceTimePredictorVolumeIntegral(lduh,
           lQhbnd,
           lGradQhbnd,
           lFhbnd,
           luh,
           cellDescription.getOffset() + 0.5 * cellDescription.getSize(),
           cellDescription.getSize(),
           _predictorTimeStamp,
           _predictorTimeStepSize,
           true);
#endif

#if defined(FileTrace)
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

  exahype::reactive::STPStatsTracer::getInstance().writeTracingEventRunIterations(duration.count(), iterations, exahype::reactive::STPTraceKey::ADERDGOwnMigratable);
#endif

  exahype::reactive::ResilienceStatistics::getInstance().notifyExecutedTask();
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::computeErrorIndicators(bool flipped) {

  CellDescription& cellDescription = getCellDescription(_cellDescriptionsIndex,
                                                        _element);
    
  //confidence is set later if we are using the limiter checks!
  _errorIndicatorDerivative = 0;
  _errorIndicatorTimeStepSize = 0;
  _errorIndicatorAdmissibility = 0;

  if(exahype::reactive::ReactiveContext::getResilienceStrategy()
    >=exahype::reactive::ReactiveContext::ResilienceStrategy::TaskSharingResilienceChecks) {

    if(exahype::reactive::ResilienceTools::CheckFlipped && flipped) {
      _errorIndicatorDerivative = std::numeric_limits<double>::max();
      _errorIndicatorTimeStepSize = std::numeric_limits<double>::max();
      _errorIndicatorAdmissibility = std::numeric_limits<double>::max();
    }
    else if(exahype::reactive::ResilienceTools::CheckSTPsWithLowConfidence 
    || exahype::reactive::ResilienceTools::CheckAllMigratableSTPs) {
       _solver.computePredictorErrorIndicators(cellDescription,
                              _predictorTimeStepSize,
                              &_errorIndicatorDerivative,
                              &_errorIndicatorTimeStepSize,
                              &_errorIndicatorAdmissibility);
    }
  }

  if(flipped)
    logError("computeErrorIndicators", "Celldesc ="<<_cellDescriptionsIndex<<" errorIndicator derivative "<<std::setprecision(30)<<_errorIndicatorDerivative
                                                                           <<" errorIndicator time step size "<<_errorIndicatorTimeStepSize
                                                                           <<" errorIndicator admissibility "<<_errorIndicatorAdmissibility);
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::shareSTPImmediatelyOrLater() {
  if(exahype::reactive::ResilienceTools::CheckLimitedCellsOnly
     && (exahype::reactive::ReactiveContext::getInstance().getResilienceStrategy()
        >=exahype::reactive::ReactiveContext::ResilienceStrategy::TaskSharingResilienceChecks)) {
    logDebug("handleLocalExecution", "Delaying outcome, as we need to compute corrector first!");
    _solver.storePendingOutcomeToBeShared(this); //delay sharing until we can be sure that trigger has been set (correction must happen first)
  }
  else {
    if(_isCorrupted && exahype::reactive::ResilienceTools::getInstance().isTrustworthy(
                                           _errorIndicatorDerivative,
                                           _errorIndicatorTimeStepSize,
                                           _errorIndicatorAdmissibility)) {
     logWarning("handleLocalExecution","Sharing corrupted outcome but soft error has not been triggered!");
    }
    logDebug("shareSTP"," sharing "<<to_string());
    _solver.sendTaskOutcomeToOtherTeams(this);
  }
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::tryToFindAndExtractEquivalentSharedOutcome(bool previous, DeliveryStatus &status, MigratablePredictionJobData **outcome) {
  return _solver.tryToFindAndExtractOutcome(_cellDescriptionsIndex, _element, _predictorTimeStamp, _predictorTimeStepSize, status, outcome);
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::checkAgainstOutcomeAndCorrectIfActive(MigratablePredictionJobData *outcome) {
  switch(_solver.checkCellDescriptionAgainstOutcome(getCellDescription(_cellDescriptionsIndex,_element), outcome,
                                                            _predictorTimeStamp, _predictorTimeStepSize,
                                                            _errorIndicatorDerivative,
                                                            _errorIndicatorTimeStepSize,
                                                            _errorIndicatorAdmissibility)) {
    case SDCCheckResult::NoCorruption:
      break;
    case SDCCheckResult::OutcomeIsMoreTrustworthy:
      logError("handleLocalExecution", std::setprecision(16)<<"I have an outcome that is more trustworthy!"<<
                            " My error indicators: derivative = "<<_errorIndicatorDerivative<<
                            " time step size = "<<_errorIndicatorTimeStepSize<<
                            " admissibility = "<<_errorIndicatorAdmissibility<<
                            " 0ther outcome's error indicators: derivative = "<<outcome->_metadata._errorIndicatorDerivative<<
                            " time step size = "<<outcome->_metadata._errorIndicatorTimeStepSize<<
                            " admissibility = "<<outcome->_metadata._errorIndicatorAdmissibility);
       //no break!
    case SDCCheckResult::OutcomeHasZeroErrorIndicator:
      logError("handleLocalExecution", "Soft error detected in job execution: "<<to_string());
      if( exahype::reactive::ReactiveContext::getResilienceStrategy()
       == exahype::reactive::ReactiveContext::ResilienceStrategy::TaskSharingResilienceCorrection) {
         _solver.correctCellDescriptionWithOutcome(getCellDescription(_cellDescriptionsIndex,_element), outcome);
      }
      else {
        logError("handleLocalExecution", "Could have corrected soft error but correction is not enabled. Continuing...");
        exahype::reactive::ResilienceTools::getInstance().setCorruptionDetected(true);
      }
      break;
    case SDCCheckResult::MyOutcomeIsMoreOrEquallyTrustworthy:
      logError("handleLocalExecution", "There is disagreement between two outcomes but local result is more or at least equally trustworthy..");
      logError("handleLocalExecution", std::setprecision(16)<<"My error indicators: derivative = "<<_errorIndicatorDerivative<<
                                       " time step size = "<<_errorIndicatorTimeStepSize<<
                                       " admissibility = "<<_errorIndicatorAdmissibility<<
                                       " Other outcome's error indicators: derivative = "<<outcome->_metadata._errorIndicatorDerivative<<
                                       " time step size = "<<outcome->_metadata._errorIndicatorTimeStepSize<<
                                       " admissibility = "<<outcome->_metadata._errorIndicatorAdmissibility);
      break;
  }
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::needToPutBackOutcome() {
  return (exahype::reactive::ReactiveContext::getInstance().getResilienceStrategy()
        >=exahype::reactive::ReactiveContext::ResilienceStrategy::TaskSharingResilienceChecks)
      &&  exahype::reactive::ResilienceTools::CheckLimitedCellsOnly;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::needToCheckThisSTP(bool hasComputed) {
  return (exahype::reactive::ReactiveContext::getInstance().getResilienceStrategy()
        >=exahype::reactive::ReactiveContext::ResilienceStrategy::TaskSharingResilienceChecks)
      &&  hasComputed
      && !exahype::reactive::ResilienceTools::getInstance().isTrustworthy(
         _errorIndicatorDerivative,
         _errorIndicatorTimeStepSize,
         _errorIndicatorAdmissibility);
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJob::needToShare(bool hasOutcome, bool isTrustworthy, bool isOutcomeTrustworthy) {
  if(exahype::reactive::ReactiveContext::getInstance().getResilienceStrategy()
      ==exahype::reactive::ReactiveContext::ResilienceStrategy::None) {
    return false;
  }
  else if(exahype::reactive::ReactiveContext::getInstance().getResilienceStrategy()
          ==exahype::reactive::ReactiveContext::ResilienceStrategy::TaskSharing) {
    return !hasOutcome
        && ((unsigned int) AllocatedSTPsSend
            <= exahype::reactive::PerformanceMonitor::getInstance().getTasksPerTimestep()/((exahype::reactive::ReactiveContext::getInstance().getTMPINumTeams()-1)*tarch::multicore::Core::getInstance().getNumberOfThreads()));
  }
  else {
    return !hasOutcome
         || !isOutcomeTrustworthy
         || exahype::reactive::ResilienceTools::CheckLimitedCellsOnly
         || !isTrustworthy;
  }
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::sendBackOutcomeToOrigin() {
  logDebug("sendBackOutcomeToOrigin",
        " send job outcome: "<<to_string());

  assertion(!_isSkeleton); //skeleton jobs should not be offloaded to a different rank
#if defined(UseSmartMPI)
#if defined(SmartMPINB)
  MPI_Request sendBackRequests[NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME_OFFLOADING];
  _solver.mpiIsendMigratablePredictionJobOutcomeOffload(
       _lduh,
       _lQhbnd,
       _lFhbnd,
       _lGradQhbnd,
       _originRank,
       _tag,
       exahype::reactive::ReactiveContext::getInstance().getMPICommunicatorMapped(),
       sendBackRequests);
  exahype::reactive::ReactiveContext::getInstance().submitRequests(
       sendBackRequests,
       NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME_OFFLOADING,
       _tag,
       _originRank,
       sendBackHandler,
       exahype::reactive::RequestType::SendBack,
       &_solver);
#else
  _solver.mpiSendMigratablePredictionJobOutcomeOffload(
      _lduh,
      _lQhbnd,
      _lFhbnd,
      _lGradQhbnd,
      _originRank,
      _tag,
      exahype::reactive::ReactiveContext::getInstance().getMPICommunicatorMapped()
  );
  MigratablePredictionJob::sendBackHandler(&_solver, _tag, _originRank);
#endif
#else
  MPI_Request sendBackRequests[NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME_OFFLOADING];
  _solver.mpiIsendMigratablePredictionJobOutcome(
       _lduh,
       _lQhbnd,
       _lFhbnd,
       _lGradQhbnd,
       _originRank,
       _tag,
       exahype::reactive::ReactiveContext::getInstance().getMPICommunicatorMapped(),
       sendBackRequests);
  exahype::reactive::RequestManager::getInstance().submitRequests(
       sendBackRequests,
       NUM_REQUESTS_MIGRATABLE_COMM_SEND_OUTCOME_OFFLOADING,
       _tag,
       _originRank,
       sendBackHandler,
       exahype::reactive::RequestType::SendBack,
       &_solver);
#endif
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::packMetaData(MigratablePredictionJobMetaData *metadata) {
#if defined(UseSmartMPI) || defined(OffloadingMetadataPacked)
  char *tmp = metadata->_contiguousBuffer;
  std::memcpy(tmp, _center, DIMENSIONS*sizeof(double)); tmp+= DIMENSIONS*sizeof(double);
  std::memcpy(tmp, _dx, DIMENSIONS*sizeof(double)); tmp+= DIMENSIONS*sizeof(double);
  std::memcpy(tmp, &_predictorTimeStamp, sizeof(double)); tmp+= sizeof(double);
  std::memcpy(tmp, &_predictorTimeStepSize, sizeof(double)); tmp+= sizeof(double);
  std::memcpy(tmp, &_element, sizeof(int)); tmp+= sizeof(int);
  std::memcpy(tmp, &_originRank, sizeof(int)); tmp+= sizeof(int);
  std::memcpy(tmp, &_errorIndicatorDerivative, sizeof(double)); tmp+= sizeof(double);
  std::memcpy(tmp, &_errorIndicatorTimeStepSize, sizeof(double)); tmp+= sizeof(double);
  std::memcpy(tmp, &_errorIndicatorAdmissibility, sizeof(double)); tmp+= sizeof(double);
  std::memcpy(tmp, &_isCorrupted, sizeof(bool)); tmp+= sizeof(bool);
#else
  metadata->_center[0] = _center[0];
  metadata->_center[1] = _center[1];
  #if DIMENSIONS==3
  metadata->_center[2] = _center[2];
  #endif

  metadata->_dx[0] = _dx[0];
  metadata->_dx[1] = _dx[1];
  #if DIMENSIONS==3
  metadata->_dx[2] = _dx[2];
  #endif

  metadata->_predictorTimeStamp = _predictorTimeStamp;
  metadata->_predictorTimeStepSize = _predictorTimeStepSize;
  metadata->_element = _element;
  metadata->_originRank = _originRank;
  metadata->_errorIndicatorDerivative = _errorIndicatorDerivative;
  metadata->_errorIndicatorTimeStepSize = _errorIndicatorTimeStepSize;
  metadata->_errorIndicatorAdmissibility = _errorIndicatorAdmissibility;
  metadata->_isCorrupted = _isCorrupted;
#endif
}

std::string exahype::solvers::ADERDGSolver::MigratablePredictionJob::to_string() const {
  std::string result;

  result =  " cellDescriptionIndex = " + std::to_string(_cellDescriptionsIndex);
  result += " center[0] = " + std::to_string(_center[0]);
  result += " center[1] = " + std::to_string(_center[1]);
  #if DIMENSIONS==3
  result += " center[2] = " + std::to_string(_center[2]);
  #endif

  result += " time stamp = " + std::to_string(_predictorTimeStamp);
  result += " time step = " + std::to_string(_predictorTimeStepSize);
  result += " element = " + std::to_string(_element);
  result += " origin = " + std::to_string(_originRank);
  result += " error indicator derivative = " + std::to_string(_errorIndicatorDerivative);
  result += " error indicator time step size = " + std::to_string(_errorIndicatorTimeStepSize);
  result += " error indicator admissibility  = " + std::to_string(_errorIndicatorAdmissibility);
  result += " isCorrupted = " + std::to_string(_isCorrupted);

  return result;
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::sendHandler(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {
#if defined(PerformanceAnalysisOffloadingDetailed)
  static std::atomic<int> cnt=0;
  cnt++;
#endif
  logDebug("sendHandler","successful send request tag = "<<tag<<" remoteRank = "<<remoteRank);
#if !defined(OffloadingNoEarlyReceiveBacks)
  ADERDGSolver::receiveBackMigratableJob(tag, remoteRank,
      static_cast<exahype::solvers::ADERDGSolver*>(solver));
#endif
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveHandler(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {

  tbb::concurrent_hash_map<std::pair<int, int>, MigratablePredictionJobData*>::accessor a_tagRankToData;
  MigratablePredictionJobData *data;
  bool found = static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToStolenData.find(
      a_tagRankToData, std::make_pair(remoteRank, tag));
  if(!found)
    logError("receiveHandler","Didn't find migrated data, maps are inconsistent!");
  assertion(found);
  data = a_tagRankToData->second;
#if defined(UseSmartMPI) || defined(OffloadingMetadataPacked)
  data->_metadata.unpackContiguousBuffer();
#endif
#if defined(UseSmartMPI)
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToStolenData.erase(a_tagRankToData);
  a_tagRankToData.release();
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToStolenData.insert(std::make_pair(
                                                                                       std::make_pair(data->_metadata.getOrigin(), tag),
                                                                                       data));
#else
  a_tagRankToData.release();
#endif

  exahype::reactive::OffloadingAnalyser::getInstance().notifyReceivedSTPJob();
  MigratablePredictionJob *job =
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->createFromData(data,
          data->_metadata.getOrigin(), tag);
  peano::datatraversal::TaskSet spawnedSet(job);

  logDebug("receiveHandler",
      " received task : "<< data->_metadata.to_string()<<" from "<<data->_metadata.getOrigin()<<" tag = "<<tag);

  exahype::reactive::OffloadingProfiler::getInstance().notifyReceivedTask(data->_metadata.getOrigin());
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::sendBackHandler(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {
  tbb::concurrent_hash_map<std::pair<int, int>, MigratablePredictionJobData*>::accessor a_tagRankToData;

  logDebug("sendBackHandler", "looking for tag = "<<tag<<" remoteRank = "<<remoteRank);

  MigratablePredictionJobData *data;
  bool found =
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToStolenData.find(
          a_tagRankToData, std::make_pair(remoteRank, tag));
  if(!found)
    logError("sendBackHandler","Couldn't find sent back task data. Inconsistent maps! Tag = "<<tag<<" rank = "<<remoteRank);
  assertion(found);
  data = a_tagRankToData->second;
  a_tagRankToData.release();
  delete data;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToStolenData.erase(
      std::make_pair(remoteRank, tag));

  NumberOfStolenJobs--;
  assertion( NumberOfStolenJobs>=0 );
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveBackHandler(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {

  logDebug("receiveBackHandler", "received back STP job tag="<<tag<<" rank="<<remoteRank);
  tbb::concurrent_hash_map<int, CellDescription*>::accessor a_tagToCellDesc;
  bool found =
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToCellDesc.find(
          a_tagToCellDesc, tag);
  if(!found)
    logError("receiveBackHandler", "Received back task but couldn't find cell, inconsistent maps! Tag = "<<tag<<" remote rank = "<<remoteRank);
  assertion(found);
  auto cellDescription = a_tagToCellDesc->second;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToCellDesc.erase(
      a_tagToCellDesc);
  a_tagToCellDesc.release();

  tbb::concurrent_hash_map<int, MigratablePredictionJobMetaData*>::accessor a_tagToMetaData;
  found =
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToMetaData.find(
          a_tagToMetaData, tag);
  assertion(found);
  MigratablePredictionJobMetaData *metadata = a_tagToMetaData->second;
  delete metadata;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToMetaData.erase(
      a_tagToMetaData);
  a_tagToMetaData.release();

  MigratablePredictionJob *job = nullptr;
  tbb::concurrent_hash_map<const CellDescription*, MigratablePredictionJob*>::accessor a_cellDescToJob;
  found = static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapCellDescToJob.find(a_cellDescToJob, cellDescription);
  assertion(found);
  job = a_cellDescToJob->second;
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapCellDescToJob.erase(a_cellDescToJob);
  a_cellDescToJob.release();
  assertion(job!=nullptr);

  tbb::concurrent_hash_map<const CellDescription*, std::pair<int,int>>::accessor a_cellDescToTagRank;
  found = static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapCellDescToTagRank.find(a_cellDescToTagRank, cellDescription);
  assertion(found);
  static_cast<exahype::solvers::ADERDGSolver*> (solver)->_mapCellDescToTagRank.erase(a_cellDescToTagRank);
  a_cellDescToTagRank.release();

  if(exahype::reactive::ReactiveContext::getInstance().getResilienceStrategy()
     >=exahype::reactive::ReactiveContext::ResilienceStrategy::TaskSharingResilienceChecks) {

    MigratablePredictionJobData *outcome = nullptr;
    DeliveryStatus status;
    bool hasOutcome = job->tryToFindAndExtractEquivalentSharedOutcome(false, status, &outcome);

    job->computeErrorIndicators(false);
    bool isTrustworthy = exahype::reactive::ResilienceTools::getInstance().isTrustworthy(
                                                                job->_errorIndicatorDerivative,
                                                                job->_errorIndicatorTimeStepSize,
                                                                job->_errorIndicatorAdmissibility);

    if(hasOutcome) {
      bool outcomeIsTrustworthy =  exahype::reactive::ResilienceTools::getInstance().isTrustworthy(
                                                                outcome->_metadata._errorIndicatorDerivative,
                                                                outcome->_metadata._errorIndicatorTimeStepSize,
                                                                outcome->_metadata._errorIndicatorAdmissibility);
      if(job->needToCheckThisSTP(true)) {
        job->checkAgainstOutcomeAndCorrectIfActive(outcome);
        job->setFinished();
      }

      if(job->needToShare(true /*hasOutcome*/,
                     isTrustworthy,
                     outcomeIsTrustworthy)) {
        job->shareSTPImmediatelyOrLater();
      }
      if(job->needToPutBackOutcome()) {
        MigratablePredictionJobOutcomeKey key(outcome->_metadata.getCenter(), outcome->_metadata.getPredictorTimeStamp(),
                                               outcome->_metadata.getPredictorTimeStepSize(), outcome->_metadata.getElement());
        job->_solver._outcomeDatabase.insertOutcome(key, outcome, DeliveryStatus::Received);
      }
      else {
        delete outcome;
      }
      job->setFinished(); //have already checked if that was necessary, so it's safe to release the job here
    }
    else { //outcome not available
      if(job->needToShare(false, false, false)) {
        job->shareSTPImmediatelyOrLater();
      }

      if(job->needToCheckThisSTP(true)) {
        if(exahype::reactive::TimeStampAndDubiosityTeamHistory::getInstance().hasConsistentTeamTimeStampHistories()) {
          logDebug("handleLocalExecution","going into check mode "<<job->to_string());

          CheckAndCorrectSolutionJob *checkJob = new CheckAndCorrectSolutionJob(job->_solver,
                                             getCellDescription(job->_cellDescriptionsIndex,job->_element),
                                             job->_predictorTimeStamp,
                                             job->_predictorTimeStepSize,
                                             job->_errorIndicatorDerivative,
                                             job->_errorIndicatorTimeStepSize,
                                             job->_errorIndicatorAdmissibility,
                                             job->_isSkeleton);
          peano::datatraversal::TaskSet spawnedSet(checkJob);
        }
        else {
          logDebug("handleLocalExecution","Won't be able to find STP anymore as timestamps/timestep sizes have diverged. Timestamp ="
              <<std::setprecision(30)<<job->_predictorTimeStamp
              <<" time step "<<job->_predictorTimeStepSize);
          exahype::reactive::ResilienceStatistics::getInstance().notifyDetectedError();
          job->setFinished();
        }
      }
    }
  }
  // no resilience checks
  else {
    job->setFinished();
  }

  delete job;

  NumberOfRemoteJobs--;
  NumberOfEnclaveJobs--;

  assertion( NumberOfEnclaveJobs>=0 ); assertion( NumberOfRemoteJobs>=0 );
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::receiveHandlerTaskSharing(
    exahype::solvers::Solver* solver, int tag, int team) {
  logDebug("receiveHandlerTaskSharing","successful receive request");

  tbb::concurrent_hash_map<std::pair<int, int>, MigratablePredictionJobData*>::accessor a_tagRankToData;
  MigratablePredictionJobData *data;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToReplicaData.find(
      a_tagRankToData, std::make_pair(team, tag));
  data = a_tagRankToData->second;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagRankToReplicaData.erase(
      a_tagRankToData);
  a_tagRankToData.release();

#if defined(UseSmartMPI) || defined(OffloadingMetadataPacked)
  data->_metadata.unpackContiguousBuffer();
#endif

  MigratablePredictionJobOutcomeKey key(data->_metadata.getCenter(), data->_metadata.getPredictorTimeStamp(),
                                         data->_metadata.getPredictorTimeStepSize(), data->_metadata.getElement());

  double minTimeStampToKeep = static_cast<exahype::solvers::ADERDGSolver*>(solver)->getPreviousMinTimeStamp();

  if(exahype::reactive::ReactiveContext::getInstance().getResilienceStrategy()>=exahype::reactive::ReactiveContext::ResilienceStrategy::TaskSharingResilienceChecks) {
    exahype::reactive::TimeStampAndDubiosityTeamHistory::getInstance().trackTimeStepAndDubiosity(team, 
		                                              key._timestamp,
							      key._timestepSize,
							      !exahype::reactive::ResilienceTools::getInstance().isTrustworthy(
                                                         data->_metadata.getErrorIndicatorDerivative(),
                                                         data->_metadata.getErrorIndicatorTimeStepSize(),
                                                         data->_metadata.getErrorIndicatorAdmissibility()));
  }

  logDebug("receiveHandlerTaskSharing", "team "
      <<exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber()
      <<" received replica job: "
      <<data->_metadata.to_string()
      <<" min timestamp to keep "<<minTimeStampToKeep);

  if(key._timestamp<minTimeStampToKeep){
    exahype::reactive::ResilienceStatistics::getInstance().notifyLateTask();
    delete data;
    AllocatedSTPsReceive--;
  }
  else {
    MigratablePredictionJobData *data2 = nullptr;
    DeliveryStatus status;
    bool found = static_cast<exahype::solvers::ADERDGSolver*>(solver)->_outcomeDatabase.tryFindAndExtractOutcome(key, &data2, status);

    if(found && status==DeliveryStatus::Transit) {
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_outcomeDatabase.insertOutcome(key, data2, DeliveryStatus::Received);
      logDebug("receiveHandlerTaskSharing", "team "
          <<exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber()
          <<" inserted replica job: "
          <<data2->_metadata.to_string()
          <<std::setprecision(30)<<"timestep "<<data->_metadata.getPredictorTimeStepSize());
    }
    else {
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_outcomeDatabase.insertOutcome(key, data, DeliveryStatus::Received);
      logDebug("receiveHandlerTaskSharing", "team "
          <<exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber()
          <<" inserted replica job: "
          <<data->_metadata.to_string()
          <<std::setprecision(30)<<"timestep "<<data->_metadata.getPredictorTimeStepSize());
    }
    static_cast<exahype::solvers::ADERDGSolver*>(solver)->_allocatedOutcomes.push_back(key);
  }
  exahype::reactive::ResilienceStatistics::getInstance().notifyReceivedTask();

}

void exahype::solvers::ADERDGSolver::MigratablePredictionJob::sendHandlerTaskSharing(
    exahype::solvers::Solver* solver, int tag, int remoteRank) {

  logDebug("sendHandlerReplication","successfully completed send to other teams");
  exahype::reactive::ResilienceStatistics::getInstance().notifySentTask();
  tbb::concurrent_hash_map<int, MigratablePredictionJobData*>::accessor a_tagToData;
  bool found =
      static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToSTPData.find(
          a_tagToData, tag);
  assertion(found);
  MigratablePredictionJobData *data = a_tagToData->second;
  static_cast<exahype::solvers::ADERDGSolver*>(solver)->_mapTagToSTPData.erase(
      a_tagToData);
  delete data;
  AllocatedSTPsSend--;
  CompletedSentSTPs++;
  logDebug("sendHandlerReplication"," allocated stps send "<<AllocatedSTPsSend);
  a_tagToData.release();
}

exahype::solvers::ADERDGSolver::MigratablePredictionJobData::MigratablePredictionJobData(
    ADERDGSolver& solver) :
      _luh(solver.getDataPerCell()),
      _lduh(solver.getUpdateSize()),
      _lQhbnd(solver.getBndTotalSize()),
      _lFhbnd(solver.getBndFluxTotalSize()),
      _lGradQhbnd(solver.getBndGradQTotalSize()){
  AllocatedSTPs++;
}

exahype::solvers::ADERDGSolver::MigratablePredictionJobData::~MigratablePredictionJobData() {
  AllocatedSTPs--;
}

bool exahype::solvers::ADERDGSolver::MigratablePredictionJobData::matches(const MigratablePredictionJobData& a) const {
  bool equal = true;
  bool tmp;
  assertion(a._lQhbnd.size()==_lQhbnd.size());
  assertion(a._lFhbnd.size()==_lFhbnd.size());
  assertion(a._lduh.size()==_lduh.size());
  assertion(a._lGradQhbnd.size()==_lGradQhbnd.size());

  tmp = _metadata == a._metadata; equal&=tmp;

  tmp = exahype::reactive::ResilienceTools::getInstance().isEqual(a._lQhbnd.data(), _lQhbnd.data(), a._lQhbnd.size()); equal&=tmp;
  if(!tmp) {
    logError("matches", "lQhbnd is not (numerically) equal");
  }
#if defined(OffloadingGradQhbnd)
  tmp = exahype::reactive::ResilienceTools::getInstance().isEqual(a._lGradQhbnd.data(), _lGradQhbnd, a._lGradQhbnd.size()); equal&=tmp;
  if(!tmp) {
    logError("matches", "lGradQhbnd is not (numerically) equal");
  }
#endif
  tmp = exahype::reactive::ResilienceTools::getInstance().isEqual(a._lFhbnd.data(), _lFhbnd.data(), a._lFhbnd.size()); equal&=tmp;
  if(!tmp) {
    logError("matches", "lFhbnd is not  (numerically) equal");
  }
  tmp = exahype::reactive::ResilienceTools::getInstance().isEqual(a._lduh.data(), _lduh.data(), a._lduh.size()); equal&=tmp;
  if(!tmp) {
    logError("matches", "lduh is not  (numerically) equal");
  }
  return equal;
}

exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::MigratablePredictionJobMetaData()
: _predictorTimeStamp(0),
  _predictorTimeStepSize(0),
  _element(0),
  _originRank(-1),
  _errorIndicatorTimeStepSize(0),
  _errorIndicatorDerivative(0),
  _errorIndicatorAdmissibility(0),
  _isCorrupted(false),
  _contiguousBuffer(nullptr) {
#if defined(UseSmartMPI) || defined(OffloadingMetadataPacked)
  _contiguousBuffer = (char*) allocate_smartmpi(getMessageLen(),0); //todo: does Alignment help here?
#endif
}

exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::~MigratablePredictionJobMetaData() {
#if defined(UseSmartMPI) || defined(OffloadingMetadataPacked)
  free_smartmpi(_contiguousBuffer);
#endif
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::unpackContiguousBuffer() {
  char *tmp = _contiguousBuffer;
  std::memcpy(_center, tmp, DIMENSIONS*sizeof(double)); tmp+= DIMENSIONS*sizeof(double);
  std::memcpy(_dx, tmp, DIMENSIONS*sizeof(double)); tmp+= DIMENSIONS*sizeof(double);
  std::memcpy(&_predictorTimeStamp, tmp, sizeof(double)); tmp+= sizeof(double);
  std::memcpy(&_predictorTimeStepSize, tmp,  sizeof(double)); tmp+= sizeof(double);
  std::memcpy(&_element, tmp, sizeof(int)); tmp+= sizeof(int);
  std::memcpy(&_originRank, tmp, sizeof(int)); tmp+= sizeof(int);
  std::memcpy(&_errorIndicatorDerivative, tmp, sizeof(double)); tmp+= sizeof(double);
  std::memcpy(&_errorIndicatorTimeStepSize, tmp, sizeof(double)); tmp+= sizeof(double);
  std::memcpy(&_errorIndicatorAdmissibility, tmp, sizeof(double)); tmp+= sizeof(double);
  std::memcpy(&_isCorrupted, tmp, sizeof(bool)); tmp+= sizeof(bool);
}

std::string exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::to_string() const {
  std::string result;

  result = " center[0] = " + std::to_string(_center[0]);
  result += " center[1] = " + std::to_string(_center[1]);
  #if DIMENSIONS==3
  result += " center[2] = " + std::to_string(_center[2]);
  #endif

  result += " time stamp = " + std::to_string(_predictorTimeStamp);
  result += " time step = " + std::to_string(_predictorTimeStepSize);
  result += " element = " + std::to_string(_element);
  result += " origin = " + std::to_string(_originRank);
  result += " error indicator derivative = " + std::to_string(_errorIndicatorDerivative);
  result += " error indicator time step size = " + std::to_string(_errorIndicatorTimeStepSize);
  result += " error indicator admissibility = " + std::to_string(_errorIndicatorAdmissibility);
  result += " isCorrupted = " + std::to_string(_isCorrupted);

  return result;
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::initDatatype() {
#if defined(UseSmartMPI) || defined(OffloadingMetadataPacked)
  _datatype = MPI_BYTE;
#else
  int entries = 2 /*double arrays*/ +8 /*scalars*/;  //be careful here when adapting
  MigratablePredictionJobMetaData dummy;

  int blocklengths[] = {DIMENSIONS, DIMENSIONS, 1, 1, 1, 1, 1, 1, 1, sizeof(bool)};
  MPI_Datatype subtypes[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_BYTE};
  MPI_Aint displs[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  MPI_Aint base;
  MPI_Get_address(&dummy, &base);

  MPI_Get_address(&(dummy._center), &(displs[0]));
  MPI_Get_address(&(dummy._dx), &(displs[1]));
  MPI_Get_address(&(dummy._predictorTimeStamp), &(displs[2]));
  MPI_Get_address(&(dummy._predictorTimeStepSize), &(displs[3]));
  MPI_Get_address(&(dummy._element), &(displs[4]));
  MPI_Get_address(&(dummy._originRank), &(displs[5]));
  MPI_Get_address(&(dummy._errorIndicatorDerivative), &(displs[6]));
  MPI_Get_address(&(dummy._errorIndicatorTimeStepSize), &(displs[7]));
  MPI_Get_address(&(dummy._errorIndicatorAdmissibility), &(displs[8]));
  MPI_Get_address(&(dummy._isCorrupted), &(displs[9]));

  for(int i=0; i<entries; i++) {
    displs[i] = displs[i]-base;
    logDebug("initDatatype" , " displ["<<i<<"]="<<displs[i]);
  }

  int ierr = MPI_Type_create_struct(entries, blocklengths, displs, subtypes, &_datatype); assert(ierr==MPI_SUCCESS);
  ierr = MPI_Type_commit(&_datatype);  assert(ierr==MPI_SUCCESS);
#endif
}

void exahype::solvers::ADERDGSolver::MigratablePredictionJobMetaData::shutdownDatatype() {
#if !(defined(UseSmartMPI) || defined(OffloadingMetadataPacked))
  MPI_Type_free(&_datatype);
#endif
}

#endif
//#undef assertion
//#define assertion(expr) 
