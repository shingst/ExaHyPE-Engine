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

#if defined(SharedTBB)  && defined(Parallel)

#include "exahype/offloading/SoftErrorInjector.h"

//#include "exahype/parser/Parser.h"

#ifdef TMPI
#include "teaMPI.h"
#endif


tarch::logging::Log exahype::offloading::SoftErrorInjector::_log("exahype::offloading::SoftErrorInjector");


exahype::offloading::SoftErrorInjector::SoftErrorInjector()
 : _injectionInterval(10000), _cnt(1)
{}

exahype::offloading::SoftErrorInjector::~SoftErrorInjector() {}


exahype::offloading::SoftErrorInjector& exahype::offloading::SoftErrorInjector::getInstance() {
  static SoftErrorInjector singleton;
  return singleton;
}

void exahype::offloading::SoftErrorInjector::generateBitflipErrorInDoubleIfActive(double *array, size_t size) {
  _cnt++;

  if(_cnt.load()%_injectionInterval==0) {
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
  }
}

#endif
