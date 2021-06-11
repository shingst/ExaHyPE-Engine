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

//#if defined(SharedTBB)  && defined(Parallel)
#if  defined(Parallel)

#include "exahype/reactive/ResilienceTools.h"
#include "exahype/reactive/OffloadingContext.h"
//#include "exahype/parser/Parser.h"

#include "tarch/parallel/Node.h"

#ifdef USE_TMPI
#include "teaMPI.h"
#endif


tarch::logging::Log exahype::reactive::ResilienceTools::_log("exahype::reactive::ResilienceTools");

exahype::reactive::ResilienceTools::SoftErrorGenerationStrategy exahype::reactive::ResilienceTools::GenerationStrategy;
bool exahype::reactive::ResilienceTools::CheckAllMigratableSTPs;
bool exahype::reactive::ResilienceTools::CheckLimitedCellsOnly;
bool exahype::reactive::ResilienceTools::CheckFlipped;


exahype::reactive::ResilienceTools::ResilienceTools()
 : _injectionInterval(1000),
   _numFlips(1),
   _cnt(1),
   _numFlipped(0),
   _infNormTol(0.0000001),
   _l1NormTol(0.0000001),
   _l2NormTol(0.0000001),
   _corruptionDetected(false),
   _injectionRank(2)
{
  if(_injectionRank<=tarch::parallel::Node::getInstance().getNumberOfNodes()
     && GenerationStrategy!=SoftErrorGenerationStrategy::None) {
    logWarning("ResilienceTools","The error injection rank is currently hardcoded to a rank bigger than the total number of ranks, no error will be injected...");
  }
   	//assertion(_injectionRank<=tarch::parallel::Node::getInstance().getNumberOfNodes());
}

exahype::reactive::ResilienceTools::~ResilienceTools() {}

exahype::reactive::ResilienceTools& exahype::reactive::ResilienceTools::getInstance() {
  static ResilienceTools singleton;
  return singleton;
}

bool exahype::reactive::ResilienceTools::checkSTPsImmediatelyAfterComputation() {
  return (exahype::reactive::OffloadingContext::getInstance().getResilienceStrategy()
      >= exahype::reactive::OffloadingContext::ResilienceStrategy::TaskSharingResilienceChecks)
      && (CheckAllMigratableSTPs || CheckFlipped);
}

bool exahype::reactive::ResilienceTools::generateBitflipErrorInDoubleIfActive(double *array, size_t size) {
  _cnt++;

  if(_cnt.load()%_injectionInterval==0
      && _numFlipped<_numFlips
      && tarch::parallel::Node::getInstance().getRank()==_injectionRank) {
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

    return true;
  }
  return false;
}

bool exahype::reactive::ResilienceTools::overwriteRandomValueInArrayIfActive(double *array, size_t size) {
  _cnt++;
  if(_cnt.load()%_injectionInterval==0 && _numFlipped<_numFlips
      && tarch::parallel::Node::getInstance().getRank()==_injectionRank) {
    std::random_device r;
    std::default_random_engine generator(r());
    std::uniform_int_distribution<int> un_arr(0, size-1);

    int idx_array = un_arr(generator);

    double old_val = array[idx_array];
    //overwrite with "random number"
    array[idx_array] = 8192; //std::numeric_limits<double>::max();
    _numFlipped++;

    logError("overwriteDoubleIfActive()", "overwrite double value, pos = "<<idx_array<<" old ="<<old_val<<" new = "<<array[idx_array]);
    return true;
  }
  return false;
}

bool exahype::reactive::ResilienceTools::corruptDataIfActive(double *array, size_t size) {
  bool result = false;
#ifdef USE_TMPI
if(TMPI_IsLeadingRank()) {
#endif
  switch(GenerationStrategy) {
  case SoftErrorGenerationStrategy::None:
    result = false;
    break;
  case SoftErrorGenerationStrategy::Bitflips:
    result = generateBitflipErrorInDoubleIfActive(array, size);
    break;
  case SoftErrorGenerationStrategy::Overwrite:
    result = overwriteRandomValueInArrayIfActive(array, size);
    break;
  default:
    result = false;
  }
#ifdef USE_TMPI
}
#endif
  return result;
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
    for(size_t i=0; i<length; i++) {
      logError("isAdmissibleNumericalError", "i = "<<i<<" a = "<<a1[i]<<" , b = "<<a2[i]);
    }
  }

  return admissible;
}

void exahype::reactive::ResilienceTools::setCorruptionDetected(bool corrupted) {
  _corruptionDetected = true;
}

bool exahype::reactive::ResilienceTools::getCorruptionDetected() {
  return _corruptionDetected;
}
#endif
