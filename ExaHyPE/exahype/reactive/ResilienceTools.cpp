/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2020  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/


#include "exahype/reactive/ResilienceTools.h"
#include "exahype/reactive/ReactiveContext.h"
#include "exahype/reactive/ResilienceStatistics.h"
//#include "exahype/parser/Parser.h"

#include "tarch/parallel/Node.h"

#ifdef USE_TMPI
#include "teaMPI.h"
#endif

#include <iomanip>

tarch::logging::Log exahype::reactive::ResilienceTools::_log("exahype::reactive::ResilienceTools");

exahype::reactive::ResilienceTools::SoftErrorGenerationStrategy exahype::reactive::ResilienceTools::GenerationStrategy;
bool exahype::reactive::ResilienceTools::CheckAllMigratableSTPs;
bool exahype::reactive::ResilienceTools::CheckLimitedCellsOnly;
bool exahype::reactive::ResilienceTools::CheckFlipped;
bool exahype::reactive::ResilienceTools::CheckSTPsWithLowConfidence;

bool exahype::reactive::ResilienceTools::CheckSTPDerivatives;
bool exahype::reactive::ResilienceTools::CheckSTPTimeSteps;
bool exahype::reactive::ResilienceTools::CheckSTPAdmissibility;
bool exahype::reactive::ResilienceTools::CheckSTPsLazily;

exahype::reactive::ResilienceTools::ResilienceTools()
  : _numFlips(1),
   _cnt(1),
   _numFlipped(0),
   _infNormTol(0.0000001),
   _l1NormTol(0.0000001),
   _l2NormTol(0.0000001),
   _corruptionDetected(false),
   _injectionRank(0),
   _absError(0),
   _relError(0),
   _maxErrorIndicatorDerivatives(0),
   _maxErrorIndicatorTimeStepSizes(0),
   _injectionPosition(),
   _injectionTime(-1) {

  if((_injectionRank>=tarch::parallel::Node::getInstance().getNumberOfNodes() || _injectionRank<0)
     && GenerationStrategy!=SoftErrorGenerationStrategy::None) {
    logWarning("ResilienceTools","The error injection rank is invalid, resetting injection rank to rank 0...");
    _injectionRank = 0;
  }
   	//assertion(_injectionRank<=tarch::parallel::Node::getInstance().getNumberOfNodes());
  //todo: this random shuffling is currently hardcoded
  //we assume 729 cells
  std::random_device r;
  std::default_random_engine generator(r());
  std::uniform_int_distribution<int> un_injection_int(0, 729);

  _injectionInterval = 3000 + un_injection_int(r);
}

exahype::reactive::ResilienceTools::~ResilienceTools() {}

void exahype::reactive::ResilienceTools::configure(double absError,
                    double relError,
                    double injectionTime,
                    tarch::la::Vector<DIMENSIONS, double> injectionPos,
                    int injectionRank,
                    double maxErrorIndicatorDerivatives,
                    double maxErrorIndicatorTimeStepSizes) {
  _absError = absError;
  _relError = relError;
  _maxErrorIndicatorDerivatives = maxErrorIndicatorDerivatives;
  _maxErrorIndicatorTimeStepSizes = maxErrorIndicatorTimeStepSizes;
  
  _injectionTime = injectionTime;
  _injectionPosition = injectionPos;

  if((injectionRank>=tarch::parallel::Node::getInstance().getNumberOfNodes() || injectionRank<0)
     && GenerationStrategy!=SoftErrorGenerationStrategy::None) {
    logWarning("ResilienceTools","The error injection rank is invalid, resetting injection rank to rank 0...");
    _injectionRank = 0;
  }
  else {
    _injectionRank = injectionRank;
  }
}

bool exahype::reactive::ResilienceTools::isTrustworthy(double errorIndicatorDerivatives, double errorIndicatorTimeStepSizes, double errorIndicatorAdmissibility) {
  if(!exahype::reactive::ResilienceTools::CheckSTPsLazily) {
    return !violatesAdmissibility(errorIndicatorAdmissibility)
        && !violatesTimestep(errorIndicatorTimeStepSizes)
        && !violatesDerivatives(errorIndicatorDerivatives);
  }
  else {
    if(!CheckSTPAdmissibility && !CheckSTPAdmissibility && !CheckSTPDerivatives){
          return true;
    }
    bool dubious = true;
    dubious = dubious && (violatesAdmissibility(errorIndicatorAdmissibility) || violatesTimestep(errorIndicatorTimeStepSizes));
    if(dubious && CheckSTPDerivatives) {
      dubious = dubious && violatesDerivatives(errorIndicatorDerivatives);
    }
    return !dubious;
  }
}

bool exahype::reactive::ResilienceTools::violatesCriterion(double value, double threshold) const {
  return value>threshold;
}

bool exahype::reactive::ResilienceTools::violatesAdmissibility(double value) const {
  return CheckSTPAdmissibility && violatesCriterion(value, 0);
}

bool exahype::reactive::ResilienceTools::violatesTimestep(double value) const {
  return CheckSTPTimeSteps && violatesCriterion(value, _maxErrorIndicatorTimeStepSizes);
}

bool exahype::reactive::ResilienceTools::violatesDerivatives(double value) const {
  return CheckSTPDerivatives && violatesCriterion(value, _maxErrorIndicatorDerivatives);
}

bool exahype::reactive::ResilienceTools::shouldInjectError(const double *center, double t) {
  if(GenerationStrategy==SoftErrorGenerationStrategy::OverwriteHardcoded) {
    return tarch::la::equals(center[0], _injectionPosition[0],0.001)
                   && tarch::la::equals(center[1], _injectionPosition[1],0.001)
                   #if DIMENSIONS==3
                   && tarch::la::equals(center[2], _injectionPosition[2],0.001)
                   #endif
                   && tarch::la::equals(t, _injectionTime,0.0001);
  }
  else if(GenerationStrategy==SoftErrorGenerationStrategy::Overwrite
      || GenerationStrategy==SoftErrorGenerationStrategy::Bitflips){

    if(_cnt.fetch_add(1)%_injectionInterval==0
        && _numFlipped<_numFlips
        && tarch::parallel::Node::getInstance().getRank()==_injectionRank) {
      return true;
    }
  }
  return false;
}

exahype::reactive::ResilienceTools& exahype::reactive::ResilienceTools::getInstance() {
  static ResilienceTools singleton;
  return singleton;
}

bool exahype::reactive::ResilienceTools::checkSTPsImmediatelyAfterComputation() {
#if defined(Parallel)
  return (exahype::reactive::ReactiveContext::getInstance().getResilienceStrategy()
      >= exahype::reactive::ReactiveContext::ResilienceStrategy::TaskSharingResilienceChecks)
      && (CheckAllMigratableSTPs || CheckFlipped);
#else
  return false;
#endif
}

void exahype::reactive::ResilienceTools::generateBitflipErrorInDouble(const double *ref, double *center, int dim, double t, double *array, size_t size) {

  std::random_device r;
  std::default_random_engine generator(r());
  std::uniform_int_distribution<int> un_arr(0, size-1);
  std::uniform_int_distribution<int> un_byte(0, sizeof(double)-1);
  std::uniform_int_distribution<int> un_bit(0, 7);

  int idx_array = un_arr(generator);
  int idx_byte = un_byte(generator);
  int idx_bit = un_bit(generator);

  double old_val = array[idx_array];

  double new_val = old_val;
  char * ptr = reinterpret_cast<char*>(&new_val);
  char mask = 1 << idx_bit;
  ptr[idx_byte] = ptr[idx_byte] ^ mask;
  array[idx_array] = new_val;

  logError("generateBitflipErrorInDoubleIfActive()","generating bitflip: pos = "<<idx_array<<" byte = "<<idx_byte<<" bit = "<<idx_bit<< " old ="<<old_val<<" new = "<<new_val);
  _cnt = 0;
  _numFlipped++;
  exahype::reactive::ResilienceStatistics::getInstance().notifyInjectedError();
}

void exahype::reactive::ResilienceTools::overwriteRandomValueInArray(const double *ref, double *center, int dim, double t, double *array, size_t size) {
  std::random_device r;
  std::default_random_engine generator(r());
  std::uniform_int_distribution<int> un_arr(0, size-1);

  int idx_array = un_arr(generator);
  //int idx_array = 0;

  double old_val = array[idx_array];
  double error = (std::abs(_relError)>0) ? (_relError*(ref[idx_array]+old_val)) //introduces relative error into new ref (e.g., new solution), if added to ref
                                : _absError;           //only introduces absolute error


  //overwrite with "random number"
  array[idx_array] += error;
  //std::numeric_limits<double>::max();
  _numFlipped++;

  logError("overwriteDoubleIfActive()", "overwrite double value, pos = "<<idx_array<<std::setprecision(30)
                                     <<" old ="<<old_val
                                     <<" new = "<<array[idx_array]
                                     <<" corresponds to relative error "<<error/(ref[idx_array]+old_val)
                                     <<" corresponds to absolute error "<<error
                                     <<" max error indicator derivatives "<<_maxErrorIndicatorDerivatives
                                     <<" max error indicator timestepsizes "<<_maxErrorIndicatorTimeStepSizes);
  exahype::reactive::ResilienceStatistics::getInstance().notifyInjectedError();
}

void exahype::reactive::ResilienceTools::overwriteHardcoded(const double *ref, double *center, int dim, double t, double *array, size_t size) {

  //logError("overwriteHardcodedIfActive", "center[0]="<<center[0]<<", center[1]="<<center[1]<<", center[2]="<<center[2]<<" t "<<t);

  //without limiter
  /*bool isActive = tarch::la::equals(center[0],3.00,0.001)
                 && tarch::la::equals(center[1],3.00,0.001)
                 && tarch::la::equals(t, 0.0197496,0.0001);*/

  //with limiter
  //exahype::reactive::ResilienceTools::overwriteHardcodedIfActive center[0]=1.8, center[1]=0.12 t 0.138786

  int idx_array = 0;

  double old_val = array[idx_array];
  double error = (std::abs(_relError)>0) ? (_relError*(ref[idx_array]+old_val)) //introduces relative error into new ref (e.g., new solution), if added to ref
                                : _absError;           //only introduces absolute error

  //array[idx_array] = 0.1; //ADER-DG only
  //array[idx_array] = 20; //limiter SWE immediately
  //array[idx_array] = 2; //limiter SWE later
  array[idx_array] += error;
  _numFlipped++;

#if DIMENSIONS==2
  logError("overwriteHardcodedIfActive()", "overwrite double value, pos = "<<idx_array<<" old ="<<old_val<<" new = "<<array[idx_array]
                                                                         <<" in center[0]="<<center[0]
                                                                         <<" center[1]="<<center[1]
                                                                         <<" t ="<<t
                                                                         <<" corresponds to relative error "<<error/(ref[idx_array]+old_val));
#elif DIMENSIONS==3
  logError("overwriteHardcodedIfActive()", "overwrite double value, pos = "<<idx_array<<" old ="<<old_val<<" new = "<<array[idx_array]
                                                                         <<" in center[0]="<<center[0]<<" center[1]="<<center[1]
                                                                         <<" center[2]="<<center[2]
                                                                         <<" t ="<<t
                                                                         <<" corresponds to relative error "<<error/(ref[idx_array]+old_val));
#endif

  exahype::reactive::ResilienceStatistics::getInstance().notifyInjectedError();
}

void exahype::reactive::ResilienceTools::corruptData(const double *ref, double *center, int dim, double t, double *array, size_t size) {
#ifdef USE_TMPI
if(TMPI_IsLeadingRank()) {
#endif
  switch(GenerationStrategy) {
  case SoftErrorGenerationStrategy::None:
    break;
  case SoftErrorGenerationStrategy::Bitflips:
    generateBitflipErrorInDouble(ref, center, dim, t, array, size);
    break;
  case SoftErrorGenerationStrategy::Overwrite:
    overwriteRandomValueInArray(ref, center, dim, t, array, size);
    break;
  case SoftErrorGenerationStrategy::OverwriteHardcoded:
    overwriteHardcoded(ref, center, dim, t, array, size);
    break;
  default:
    break;
  }
#ifdef USE_TMPI
}
#endif
}

void exahype::reactive::ResilienceTools::setSoftErrorGenerationStrategy(SoftErrorGenerationStrategy strat) {
  GenerationStrategy = strat;
}

double exahype::reactive::ResilienceTools::computeInfNormError(const double *a1, const double *a2, size_t length) {
  double result = 0.0;
  for(size_t i = 0; i<length; i++) {
    result = std::max(result, std::abs(a1[i]-a2[i]));
  }
  return result;
}

double exahype::reactive::ResilienceTools::computeL1NormError(const double *a1, const double *a2, size_t length){
  double result = 0.0;
  for(size_t i = 0; i<length; i++) {
    result += std::abs(a1[i]-a2[i]);
  }
  return result;
}

double exahype::reactive::ResilienceTools::computeL2NormError(const double *a1, const double *a2, size_t length){
  double result = 0.0;
  for(size_t i = 0; i<length; i++) {
    result += std::pow((a1[i]-a2[i]),2);
  }
  return std::sqrt(result);
}

double exahype::reactive::ResilienceTools::computeInfNormErrorRel(const double *a1, const double *a2, size_t length) {
  double result = 0.0;
  for(size_t i = 0; i<length; i++) {
    if(a1[i]+a2[i]!=0)
      result = std::max(result, std::abs(a1[i]-a2[i])/(0.5*(a1[i]+a2[i])));
  }
  return result;
}

double exahype::reactive::ResilienceTools::computeL1NormErrorRel(const double *a1, const double *a2, size_t length){
  double result = 0.0;
  for(size_t i = 0; i<length; i++) {
    if(a1[i]+a2[i]!=0)
      result += std::abs(a1[i]-a2[i])/(0.5*(a1[i]+a2[i]));
  }
  return result;
}

double exahype::reactive::ResilienceTools::computeL2NormErrorRel(const double *a1, const double *a2, size_t length){
  double result = 0.0;
  for(size_t i = 0; i<length; i++) {
    if(a1[i]+a2[i]!=0)
     result += std::pow((a1[i]-a2[i])/(0.5*(a1[i]+a2[i])),2);
  }
  return std::sqrt(result);
}

bool exahype::reactive::ResilienceTools::isAdmissibleNumericalError(const double *a1, const double *a2, size_t length) {
  //double infnorm = computeInfNormErrorRel(a1, a2, length);
  //double l1norm = computeL1NormErrorRel(a1, a2, length);
  //double l2norm = computeL2NormErrorRel(a1, a2, length);

  double infnorm = computeInfNormError(a1, a2, length);
  double l1norm = computeL1NormError(a1, a2, length);
  double l2norm = computeL2NormError(a1, a2, length);

  //bool admissible = (infnorm<_infNormTol)
   //         &&  (l1norm<_l1NormTol)
   //         &&  (l2norm<_l2NormTol);
  bool admissible = (infnorm==0)
            &&  (l1norm==0)
            &&  (l2norm==0);

  if(!admissible || infnorm!=0 || l1norm!=0 || l2norm!=0) {
    logError("isAdmissibleNumericalError","We'll likely have a soft error: "
                                         << " inf norm = "<<infnorm
                                         << " l1 norm = "<<l1norm
                                         << " l2 norm = "<<l2norm);
    /*for(size_t i=0; i<length; i++) {
      logError("isAdmissibleNumericalError", "i = "<<i<<" a = "<<a1[i]<<" , b = "<<a2[i]);
    }*/
  }

  return admissible;
}

void exahype::reactive::ResilienceTools::setCorruptionDetected(bool corrupted) {
  _corruptionDetected = true;
}

bool exahype::reactive::ResilienceTools::getCorruptionDetected() {
  return _corruptionDetected;
}

