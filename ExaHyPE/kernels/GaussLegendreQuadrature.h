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
#ifndef GAUSSLEGENDRE_H_
#define GAUSSLEGENDRE_H_

#include <set>

namespace kernels {

extern const double gaussLegendreWeights0 [1];
extern const double gaussLegendreWeights1 [2];
extern const double gaussLegendreWeights2 [3];
extern const double gaussLegendreWeights3 [4];
extern const double gaussLegendreWeights4 [5];
extern const double gaussLegendreWeights5 [6];
extern const double gaussLegendreWeights6 [7];
extern const double gaussLegendreWeights7 [8];
extern const double gaussLegendreWeights8 [9];
extern const double gaussLegendreWeights9 [10];
extern const double gaussLegendreWeights10[11];
extern const double gaussLegendreWeights11[12];
extern const double gaussLegendreWeights12[13];
extern const double gaussLegendreWeights13[14];
extern const double gaussLegendreWeights14[15];
extern const double gaussLegendreWeights15[16];
extern const double gaussLegendreWeights16[17];
extern const double gaussLegendreWeights17[18];
extern const double gaussLegendreWeights18[19];
extern const double gaussLegendreWeights19[20];
extern const double gaussLegendreWeights20[21];

extern const double gaussLegendreNodes0 [1];
extern const double gaussLegendreNodes1 [2];
extern const double gaussLegendreNodes2 [3];
extern const double gaussLegendreNodes3 [4];
extern const double gaussLegendreNodes4 [5];
extern const double gaussLegendreNodes5 [6];
extern const double gaussLegendreNodes6 [7];
extern const double gaussLegendreNodes7 [8];
extern const double gaussLegendreNodes8 [9];
extern const double gaussLegendreNodes9 [10];
extern const double gaussLegendreNodes10[11];
extern const double gaussLegendreNodes11[12];
extern const double gaussLegendreNodes12[13];
extern const double gaussLegendreNodes13[14];
extern const double gaussLegendreNodes14[15];
extern const double gaussLegendreNodes15[16];
extern const double gaussLegendreNodes16[17];
extern const double gaussLegendreNodes17[18];
extern const double gaussLegendreNodes18[19];
extern const double gaussLegendreNodes19[20];
extern const double gaussLegendreNodes20[21];

/**
 * The Gauss-Legendre weights mapped onto [0,1]. Array of arrays. The first
 * entry is the order, the second entry the Legendre point.
 */
extern const double* const gaussLegendreWeights[20+1];

/**
 * The Gauss-Legendre nodes mapped onto [0,1]. Array of arrays. The first entry
 * is the order, the second entry the Legendre point.
 */
extern const double* const gaussLegendreNodes[20+1];

}

#endif
