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
  : _maxNumInjections(1),
   _cntSinceLastInjection(1),
   _numInjected(0),
   _infNormTol(0.0000001),
   _l1NormTol(0.0000001),
   _l2NormTol(0.0000001),
   _corruptionDetected(false),
   _injectionRank(0),
   _injectedErrorVal(0),
   _maxErrorIndicatorDerivatives(0),
   _maxErrorIndicatorTimeStepSizes(0),
   _injectionFrequency(-1),
   _injectionPosition(),
   _injectionTime(-1) {

}

exahype::reactive::ResilienceTools::~ResilienceTools() {}

void exahype::reactive::ResilienceTools::configure(
                    double injectedErrorVal,
                    double injectionTime,
                    tarch::la::Vector<DIMENSIONS, double> injectionPos,
                    int injectionRank,
                    int injectionFrequency,
                    int maxNumInjections,
                    double maxErrorIndicatorDerivatives,
                    double maxErrorIndicatorTimeStepSizes) {
  _injectedErrorVal = injectedErrorVal;

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

  if(injectionFrequency>0) {
    _injectionFrequency = injectionFrequency;
  }
  else {
    //random STP for injection
    std::random_device r;
    std::default_random_engine generator(r());
    std::uniform_int_distribution<int> un_injection_int(0, 729);

    _injectionFrequency = 2000 + un_injection_int(r);
  }

  _maxNumInjections = maxNumInjections;

  _maxErrorIndicatorDerivatives = maxErrorIndicatorDerivatives;
  _maxErrorIndicatorTimeStepSizes = maxErrorIndicatorTimeStepSizes;

}

bool exahype::reactive::ResilienceTools::isTrustworthy(double errorIndicatorDerivatives, double errorIndicatorTimeStepSizes, double errorIndicatorAdmissibility) {
  if(!exahype::reactive::ResilienceTools::CheckSTPsLazily) {
    return !violatesAdmissibility(errorIndicatorAdmissibility)
        && !violatesTimestep(errorIndicatorTimeStepSizes)
        && !violatesDerivatives(errorIndicatorDerivatives);
  }
  else {
    if(!CheckSTPAdmissibility && !CheckSTPTimeSteps && !CheckSTPDerivatives){
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

bool exahype::reactive::ResilienceTools::violatesCriterion(double errorIndicator, double threshold) const {
  return errorIndicator>threshold;
}

bool exahype::reactive::ResilienceTools::violatesAdmissibility(double errorIndicator) const {
  return CheckSTPAdmissibility && violatesCriterion(errorIndicator, 0);
}

bool exahype::reactive::ResilienceTools::violatesTimestep(double errorIndicator) const {
  return CheckSTPTimeSteps && violatesCriterion(errorIndicator, _maxErrorIndicatorTimeStepSizes);
}

bool exahype::reactive::ResilienceTools::violatesDerivatives(double errorIndicator) const {
  return CheckSTPDerivatives && violatesCriterion(errorIndicator, _maxErrorIndicatorDerivatives);
}

bool exahype::reactive::ResilienceTools::shouldInjectError(const double *center, double t) {
  if(GenerationStrategy==SoftErrorGenerationStrategy::OverwriteHardcoded) {
    return tarch::la::equals(center[0], _injectionPosition[0],0.001)
                   && tarch::la::equals(center[1], _injectionPosition[1],0.001)
                   #if DIMENSIONS==3
                   && tarch::la::equals(center[2], _injectionPosition[2],0.001)
                   #endif
                   && (tarch::la::equals(t, _injectionTime,0.0001)
                   ||  _cntSinceLastInjection.fetch_add(1)%_injectionFrequency==0)
                   && _numInjected<_maxNumInjections;
  }
  else if(GenerationStrategy==SoftErrorGenerationStrategy::Overwrite
      || GenerationStrategy==SoftErrorGenerationStrategy::Bitflips
      || GenerationStrategy==SoftErrorGenerationStrategy::OverwriteRandom){

    if(_cntSinceLastInjection.fetch_add(1)%_injectionFrequency==0
        && _numInjected<_maxNumInjections
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

  logError("generateBitflipErrorInDouble()","generating bitflip: pos = "<<idx_array<<" byte = "<<idx_byte<<" bit = "<<idx_bit<< " old ="<<old_val<<" new = "<<new_val);
  _cntSinceLastInjection = 0;
  _numInjected++;
  
  logError("generateBitflipErrorInDouble()", "overwrite double value, pos = "<<idx_array<<std::setprecision(30)
                                     <<" old ="<<old_val
                                     <<" new = "<<array[idx_array]
                                     <<" corresponds to relative error "<<(new_val-old_val)/(ref[idx_array]+old_val)
                                     <<" corresponds to absolute error "<<new_val-old_val
                                     <<" max error indicator derivatives "<<_maxErrorIndicatorDerivatives
                                     <<" max error indicator timestepsizes "<<_maxErrorIndicatorTimeStepSizes);
  exahype::reactive::ResilienceStatistics::getInstance().notifyInjectedError();
}

void exahype::reactive::ResilienceTools::overwriteRandomValueInArrayWithGivenError(const double *ref, double *center, int dim, double t, double *array, size_t size) {
  std::random_device r;
  std::default_random_engine generator(r());
  std::uniform_int_distribution<int> un_arr(0, size-1);

  int idx_array = un_arr(generator);

  double old_val = array[idx_array];
  double error =  _injectedErrorVal;      //introduces error of fixed value

  //overwrite
  array[idx_array] += error;
  _numInjected++;

  logError("overwriteRandomValueInArrayWithGivenError()", "overwrite double value, pos = "<<idx_array<<std::setprecision(30)
                                     <<" old ="<<old_val
                                     <<" new = "<<array[idx_array]
                                     <<" corresponds to relative error "<<error/(ref[idx_array]+old_val)
                                     <<" corresponds to absolute error "<<error
                                     <<" max error indicator derivatives "<<_maxErrorIndicatorDerivatives
                                     <<" max error indicator timestepsizes "<<_maxErrorIndicatorTimeStepSizes);
  exahype::reactive::ResilienceStatistics::getInstance().notifyInjectedError();
}

void exahype::reactive::ResilienceTools::overwriteRandomValueInArrayWithRandomError(const double *ref, double *center, int dim, double t, double *array, size_t size) {
  std::random_device r;
  std::default_random_engine generator(r());
  std::uniform_int_distribution<int> un_arr(0, size-1);
  std::uniform_real_distribution<double> un_err(1e-07, std::numeric_limits<double>::max()-1);
  std::uniform_int_distribution<int> un_sign(0, 1);

  int idx_array = un_arr(generator);
  //int idx_array = 0;

  double old_val = array[idx_array];

  int sign = un_sign(generator);
  sign = (sign==0) ? -1 : 1;
  double error = sign * un_err(generator);

  logError("overwriteRandomValueInArrayWithRandomError()","generated sign = "<<sign<<" error "<<error);

  //overwrite with "random number"
  array[idx_array] += error;
  //std::numeric_limits<double>::max();
  _numInjected++;

  logError("overwriteRandomValueInArrayWithRandomError()", "overwrite double value, pos = "<<idx_array<<std::setprecision(30)
                                     <<" old ="<<old_val
                                     <<" new = "<<array[idx_array]
                                     <<" corresponds to relative error "<<error/(ref[idx_array]+old_val)
                                     <<" corresponds to absolute error "<<error
                                     <<" max error indicator derivatives "<<_maxErrorIndicatorDerivatives
                                     <<" max error indicator timestepsizes "<<_maxErrorIndicatorTimeStepSizes);
  exahype::reactive::ResilienceStatistics::getInstance().notifyInjectedError();
}

void exahype::reactive::ResilienceTools::overwriteHardcoded(const double *ref, double *center, int dim, double t, double *array, size_t size) {

  int idx_array = 0;

  double old_val = array[idx_array];
  double error =  _injectedErrorVal;           //only introduces error of fixed value

  array[idx_array] += error;
  _numInjected++;

#if DIMENSIONS==2
  logError("overwriteHardcoded()", "overwrite double value, pos = "<<idx_array<<" old ="<<old_val<<" new = "<<array[idx_array]
                                                                         <<" in center[0]="<<center[0]
                                                                         <<" center[1]="<<center[1]
                                                                         <<" t ="<<t
                                                                         <<" corresponds to relative error "<<error/(ref[idx_array]+old_val));
#elif DIMENSIONS==3
  logError("overwriteHardcoded()", "overwrite double value, pos = "<<idx_array<<" old ="<<old_val<<" new = "<<array[idx_array]
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
    overwriteRandomValueInArrayWithGivenError(ref, center, dim, t, array, size);
    break;
  case SoftErrorGenerationStrategy::OverwriteRandom:
    overwriteRandomValueInArrayWithRandomError(ref, center, dim, t, array, size);
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

bool exahype::reactive::ResilienceTools::isEqual(const double *a1, const double *a2, size_t length) {
  double infnorm = computeInfNormError(a1, a2, length);
  double l1norm = computeL1NormError(a1, a2, length);
  double l2norm = computeL2NormError(a1, a2, length);

  bool equal = (infnorm==0)
            &&  (l1norm==0)
            &&  (l2norm==0);  //todo: we may need some tolerance but in all experiments,
                              // teams compute exactly (in a bitwise sense) the same results

  return equal;
}

void exahype::reactive::ResilienceTools::setCorruptionDetected(bool corrupted) {
  _corruptionDetected = true;
}

bool exahype::reactive::ResilienceTools::getCorruptionDetected() {
  return _corruptionDetected;
}

