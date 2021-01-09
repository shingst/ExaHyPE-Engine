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

#include "ResilienceTools.h"

//#include "exahype/parser/Parser.h"

#ifdef TMPI
#include "teaMPI.h"
#endif


tarch::logging::Log exahype::offloading::ResilienceTools::_log("exahype::offloading::ResilienceTools");

bool exahype::offloading::ResilienceTools::GenerateErrors;
bool exahype::offloading::ResilienceTools::TriggerAllMigratableSTPs;
bool exahype::offloading::ResilienceTools::TriggerLimitedCellsOnly;
bool exahype::offloading::ResilienceTools::TriggerFlipped;

exahype::offloading::ResilienceTools::ResilienceTools()
 : _injectionInterval(10000), _numFlips(1), _cnt(1), _numFlipped(0), _infNormTol(0.0000001), _l1NormTol(0.0000001), _l2NormTol(0.0000001)
{}

exahype::offloading::ResilienceTools::~ResilienceTools() {}


exahype::offloading::ResilienceTools& exahype::offloading::ResilienceTools::getInstance() {
  static ResilienceTools singleton;
  return singleton;
}

bool exahype::offloading::ResilienceTools::generateBitflipErrorInDoubleIfActive(double *array, size_t size) {
  if(!GenerateErrors) return false;

  _cnt++;

  if(_cnt.load()%_injectionInterval==0 && _numFlipped<_numFlips) {
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

    logInfo("generateBitflipErrorInDoubleIfActive()","generating bitflip: pos = "<<idx_array<<" byte = "<<idx_byte<<" bit = "<<idx_bit<< " old ="<<old_val<<" new = "<<new_val);
    _cnt = 0;
    _numFlipped++;

    return true;
  }
  return false;
}

double exahype::offloading::ResilienceTools::computeInfNormError(double *a1, double *a2, size_t length) {
  double result = 0.0;
  for(size_t i = 0; i<length; i++) {
    result = std::max(result, std::abs(a1[i]-a2[i]));
  }
  return result;
}

double exahype::offloading::ResilienceTools::computeL1NormError(double *a1, double *a2, size_t length){
  double result = 0.0;
  for(size_t i = 0; i<length; i++) {
    result += std::abs(a1[i]-a2[i]);
  }
  return result;
}

double exahype::offloading::ResilienceTools::computeL2NormError(double *a1, double *a2, size_t length){
  double result = 0.0;
  for(size_t i = 0; i<length; i++) {
    result += std::pow((a1[i]-a2[i]),2);
  }
  return std::sqrt(result);
}

double exahype::offloading::ResilienceTools::computeInfNormErrorRel(double *a1, double *a2, size_t length) {
  double result = 0.0;
  for(size_t i = 0; i<length; i++) {
	if(a1[i]+a2[i]!=0)
      result = std::max(result, std::abs(a1[i]-a2[i])/(0.5*(a1[i]+a2[i])));
  }
  return result;
}

double exahype::offloading::ResilienceTools::computeL1NormErrorRel(double *a1, double *a2, size_t length){
  double result = 0.0;
  for(size_t i = 0; i<length; i++) {
    if(a1[i]+a2[i]!=0)
	 result += std::abs(a1[i]-a2[i])/(0.5*(a1[i]+a2[i]));
  }
  return result;
}

double exahype::offloading::ResilienceTools::computeL2NormErrorRel(double *a1, double *a2, size_t length){
  double result = 0.0;
  for(size_t i = 0; i<length; i++) {
    if(a1[i]+a2[i]!=0)
     result += std::pow((a1[i]-a2[i])/(0.5*(a1[i]+a2[i])),2);
  }
  return std::sqrt(result);
}

bool exahype::offloading::ResilienceTools::isAdmissibleNumericalError(double *a1, double *a2, size_t length) {
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

#endif
