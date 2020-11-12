
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


#ifndef EXAHYPE_OFFLOADING_SOFTERRORINJECTOR_H_
#define EXAHYPE_OFFLOADING_SOFTERRORINJECTOR_H_

#include <atomic>
#include <random>

#include "tarch/logging/Log.h"

namespace exahype {
  namespace offloading {
    class SoftErrorInjector;
  }
}

class exahype::offloading::SoftErrorInjector {

  private:

  /**
   * The log device of this class.
   */
  static tarch::logging::Log _log;

  int _injectionInterval;
  std::atomic<int> _cnt;

  public:
  SoftErrorInjector();
  static SoftErrorInjector& getInstance();

  void generateBitflipErrorInDoubleIfActive(double *array, size_t length);
  virtual ~SoftErrorInjector();

};


#endif /* EXAHYPE_OFFLOADING_SOFTERRORINJECTOR_H_ */
