/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
#ifdef TEST_OPT_KERNEL
 
#ifndef _EXAHYPE_TESTS_OPTIMISED_KERNEL_TEST_H_
#define _EXAHYPE_TESTS_OPTIMISED_KERNEL_TEST_H_

#include "peano/utils/Globals.h"
#include "peano/utils/Loop.h"
#include "tarch/logging/Log.h"
#include "tarch/tests/TestCaseFactory.h"
#include "tarch/tests/TestCase.h"
#include "tarch/compiler/CompilerSpecificSettings.h"

#include "kernels/aderdg/generic/Kernels.h"
#include "kernels/aderdg/optimised/Kernels.h"
#include "kernels/KernelUtils.h"

namespace exahype {
namespace tests {
namespace c {

class OptimisedKernelTest : public tarch::tests::TestCase {
 public:
  OptimisedKernelTest();
  virtual ~OptimisedKernelTest();

  //solver getter
  int getNumberOfVariables();
  int getNodesPerCoordinateAxis(); //_basisSize
  
  //solver methods
  static void adjustedSolutionValues(const double* const x, const double w, const double t, const double dt, double* Q);
  static void flux(const double* const Q, double** F);
  void source(const double* const Q, double* S);
  void ncp(const double* const Q, const double* const gradQ, double* BgradQ);
  void matrixb(const double* const Q, const int normalNonZero, double* Bn);
   
#ifdef Dim3  
  static void fluxSplitted(const double* const Q, double* f, double* g, double* h);
#else 
  static void fluxSplitted(const double* const Q, double* f, double* g);
#endif
  
  //solver method implementation
  static void adjustedSolutionValues_Euler(const double* const x, const double w, const double t, const double dt, double* Q);
  static void flux_Euler(const double* const Q, double** F);
  void source_Euler(const double* const Q, double* S);
  
  //tests
  void run() override;
  void testSolutionAdjustment();
  void testSpaceTimePredictorNonLinear();
  void testSolutionUpdate();

 private:
  static tarch::logging::Log _log;
  static const double eps;  // for quick adaption of the test cases
  static const double eps2;  // for quick adaption of the test cases
  static int _numberOfVariables;
  static int _basisSize;
  static int _order;
  static int _dim;
  static bool _isLinear;
  static const std::string dim; // for log
  
  int _luhSize;
  double* _luh; //goes to the call with generic kernel
  double _dt;  //initialized by testStableTimeStepSize


};

}  // namespace c
}  // namespace tests
}  // namespace exahype


#endif //_EXAHYPE_TESTS_OPTIMISED_KERNEL_TEST_H_

#endif //TEST_OPT_KERNEL