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
#include "tarch/la/Vector.h"

#include "exahype/solvers/Solver.h"

namespace exahype {
  namespace reactive {
    class ResilienceTools;
  }
}

class exahype::reactive::ResilienceTools {
  public:
  enum class SoftErrorGenerationStrategy { None, Bitflips, Overwrite, OverwriteHardcoded };
  private:

  /**
   * The log device of this class.
   */
  static tarch::logging::Log _log;

  /**
   * Specifies how many STPs is a soft-error generated.
   */
  const int _injectionInterval;

  const int _numFlips;

  /**
   * Internal counter of executed STPs since last bitflip.
   */
  std::atomic<int> _cnt;

  /**
   * Counts how many bits have already been flipped.
   */
  std::atomic<int> _numFlipped;

  double _infNormTol;
  double _l1NormTol;
  double _l2NormTol;

  bool _corruptionDetected;

  int _injectionRank;

  double _absError;
  double _relError;
  double _confidenceRequired;

  tarch::la::Vector<DIMENSIONS, double> _injectionPosition;
  double _injectionTime;

  double _minScalingFactorOfDerivative;
  double _maxScalingFactorOfDerivative;

  public:

  static SoftErrorGenerationStrategy GenerationStrategy;
  static bool CheckAllMigratableSTPs;
  static bool CheckSTPsWithLowConfidence;
  static bool CheckLimitedCellsOnly;
  static bool CheckFlipped;

  static bool CheckDerivativesForConfidence;
  static bool CheckTimeStepsForConfidence;
  static bool CheckAdmissibilityForConfidence;


  ResilienceTools();
  static ResilienceTools& getInstance();

  void configure(double absError, double relError, double injectionTime,
                 tarch::la::Vector<DIMENSIONS, double> injectionPos,
                 double confidenceRequired,
                 double minScalingFactorOfDerivative,
                 double maxScalingFactorOfDerivative);

  bool corruptDataIfActive(const double *ref, double *center, int dim, double t, double *array, size_t length);

  bool isTrustworthy(double confidence);

  static void setSoftErrorGenerationStrategy(SoftErrorGenerationStrategy strat);

  static double computeInfNormError(const double *a1, const double *a2, size_t length);
  static double computeL1NormError(const double *a1, const double *a2, size_t length);
  static double computeL2NormError(const double *a1, const double *a2, size_t length);

  static double computeInfNormErrorRel(const double *a1, const double *a2, size_t length);
  static double computeL1NormErrorRel(const double *a1, const double *a2, size_t length);
  static double computeL2NormErrorRel(const double *a1, const double *a2, size_t length);

  static bool checkSTPsImmediatelyAfterComputation();

  bool isAdmissibleNumericalError(const double *a1, const double *a2, size_t length);

  void setCorruptionDetected(bool corrupted);
  bool getCorruptionDetected();

  double getMinDerivativeScalingFactor() const;
  double getMaxDerivativeScalingFactor() const;

  private:
  bool generateBitflipErrorInDoubleIfActive(double *array, size_t length);
  bool overwriteRandomValueInArrayIfActive(const double *ref, double *array, size_t size);
  bool overwriteHardcodedIfActive(const double *ref, double *center, int dim, double t,  double *array, size_t size);

  virtual ~ResilienceTools();

};


#endif /* EXAHYPE_OFFLOADING_SOFTERRORINJECTOR_H_ */
