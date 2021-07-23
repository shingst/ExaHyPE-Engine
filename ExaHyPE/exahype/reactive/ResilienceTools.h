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

/**
 * Utility class for ExaHyPE's reactive resilience mechanisms.
 * Contains functions to inject errors according to specified parameters such as injection position (spatial), frequency and error.
 * It also allows to check double arrays for equality.
 */
class exahype::reactive::ResilienceTools {
  public:
  enum class SoftErrorGenerationStrategy { None, Bitflips, Overwrite, OverwriteHardcoded, OverwriteRandom };
  private:

  /**
   * The log device of this class.
   */
  static tarch::logging::Log _log;

  /**
   * Specifies how many STPs is a soft-error generated.
   */
  int _injectionInterval;

  /**
   * Controls how many bitflips should be injected.
   */
  int _numInjections;

  /**
   * Internal counter of executed STPs since last bitflip.
   */
  std::atomic<int> _cntSinceLastFlip;

  /**
   * Counts how many bits have already been flipped.
   */
  std::atomic<int> _numInjected;

  double _infNormTol;
  double _l1NormTol;
  double _l2NormTol;

  bool _corruptionDetected;

  int _injectionRank;

  double _absError;
  double _relError;
  double _maxErrorIndicatorDerivatives;
  double _maxErrorIndicatorTimeStepSizes;

  tarch::la::Vector<DIMENSIONS, double> _injectionPosition;
  double _injectionTime;

  public:

  static SoftErrorGenerationStrategy GenerationStrategy;
  static bool CheckAllMigratableSTPs;
  static bool CheckSTPsWithLowConfidence;
  static bool CheckLimitedCellsOnly;
  static bool CheckFlipped;

  static bool CheckSTPDerivatives;
  static bool CheckSTPTimeSteps;
  static bool CheckSTPAdmissibility;
  static bool CheckSTPsLazily;

  static void setSoftErrorGenerationStrategy(SoftErrorGenerationStrategy strat);

  static double computeInfNormError(const double *a1, const double *a2, size_t length);
  static double computeL1NormError(const double *a1, const double *a2, size_t length);
  static double computeL2NormError(const double *a1, const double *a2, size_t length);

  static double computeInfNormErrorRel(const double *a1, const double *a2, size_t length);
  static double computeL1NormErrorRel(const double *a1, const double *a2, size_t length);
  static double computeL2NormErrorRel(const double *a1, const double *a2, size_t length);

  static bool checkSTPsImmediatelyAfterComputation();

  ResilienceTools();
  static ResilienceTools& getInstance();

  void configure(double absError, double relError, double injectionTime,
                 tarch::la::Vector<DIMENSIONS, double> injectionPos,
                 int injectionRank,
                 double maxErrorIndicatorDerivatives,
                 double maxErrorIndicatorTimeStepSizes);

  void corruptData(const double *ref, double *center, int dim, double t, double *array, size_t length);

  bool isTrustworthy(double errorIndicatorDerivatives, double errorIndicatorTimeStepSizes, double errorIndicatorAdmissibility);

  /**
   *
   */
  bool violatesCriterion(double errorIndicator, double threshold) const;
  bool violatesAdmissibility(double errorIndicator) const;
  bool violatesDerivatives(double errorIndicator) const;
  bool violatesTimestep(double errorIndicator) const;

  bool shouldInjectError(const double *center, double t);

  bool isEqual(const double *a1, const double *a2, size_t length);

  void setCorruptionDetected(bool corrupted);
  bool getCorruptionDetected();

  private:
  void generateBitflipErrorInDouble(const double *ref, double *center, int dim, double t, double *array, size_t length);
  void overwriteRandomValueInArrayWithGivenError(const double *ref, double *center, int dim, double t, double *array, size_t size);
  void overwriteRandomValueInArrayWithRandomError(const double *ref, double *center, int dim, double t, double *array, size_t size);
  void overwriteHardcoded(const double *ref, double *center, int dim, double t,  double *array, size_t size);

  virtual ~ResilienceTools();

};
#endif /* EXAHYPE_OFFLOADING_SOFTERRORINJECTOR_H_ */
