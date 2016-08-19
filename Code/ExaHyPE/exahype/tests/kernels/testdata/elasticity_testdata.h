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

#ifndef _EXAHYPE_TESTS_TESTDATA_ELASTICITY_TESTDATA_H_
#define _EXAHYPE_TESTS_TESTDATA_ELASTICITY_TESTDATA_H_

namespace exahype {
namespace tests {
namespace testdata {
namespace elasticity {

#ifdef Dim2

namespace testSurfaceIntegralLinear {
// extern const double lduh_out_1[80];
// extern const double lduh_out_2[80];
}  // namespace testSurfaceIntegralLinear

namespace testVolumeIntegralLinear {
// extern const double lduh_1[80];
// extern const double lduh_2[80];
}  // namespace testVolumeIntegralLinear

namespace testSpaceTimePredictorLinear {
// extern const double luh[80];
// extern const double lFhi[160];
// extern const double lQhi[80];
// extern const double lQbnd[80];
}  // namespace testSpaceTimePredictorLinear

namespace testRiemannSolverLinear {
extern const double qL_IN[9 * 5];
extern const double qR_IN[9 * 5];
extern const double paramL_IN[3 * 5];
extern const double paramR_IN[3 * 5];
extern const double FL_IN[9 * 5];
extern const double FR_IN[9 * 5];
extern const double FL_OUT[9 * 5];
extern const double FR_OUT[9 * 5];
}  // namespace testRiemannSolverLinear

#endif  // Dim2

}  // namespace elasticity
}  // namespace testdata
}  // namespace tests
}  // namespace exahype

#endif  // _EXAHYPE_TESTS_TESTDATA_ELASTICITY_TESTDATA_H_
