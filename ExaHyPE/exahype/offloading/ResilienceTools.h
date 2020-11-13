
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
    class ResilienceTools;
  }
}

class exahype::offloading::ResilienceTools {

  private:

  /**
   * The log device of this class.
   */
  static tarch::logging::Log _log;

  int _injectionInterval;
  std::atomic<int> _cnt;

  double _infNormTol;
  double _l1NormTol;
  double _l2NormTol;

  public:

  static bool GenerateErrors;
  static bool TriggerAllMigratableSTPs;
  static bool TriggerLimitedCellsOnly;
  static bool TriggerFlipped;

  ResilienceTools();
  static ResilienceTools& getInstance();

  bool generateBitflipErrorInDoubleIfActive(double *array, size_t length);

  static double computeInfNormError(double *a1, double *a2, size_t length);
  static double computeL1NormError(double *a1, double *a2, size_t length);
  static double computeL2NormError(double *a1, double *a2, size_t length);

  bool isAdmissibleNumericalError(double *a1, double *a2, size_t length);

  virtual ~ResilienceTools();

};


#endif /* EXAHYPE_OFFLOADING_SOFTERRORINJECTOR_H_ */
