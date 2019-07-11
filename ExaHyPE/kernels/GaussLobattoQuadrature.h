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
#ifndef GAUSSLOBATTO_H_
#define GAUSSLOBATTO_H_

#include <set>

namespace kernels {

extern const double gaussLobattoWeights0 [1];
extern const double gaussLobattoWeights1 [2];
extern const double gaussLobattoWeights2 [3];
extern const double gaussLobattoWeights3 [4];
extern const double gaussLobattoWeights4 [5];
extern const double gaussLobattoWeights5 [6];
extern const double gaussLobattoWeights6 [7];
extern const double gaussLobattoWeights7 [8];
extern const double gaussLobattoWeights8 [9];
extern const double gaussLobattoWeights9 [10];
extern const double gaussLobattoWeights10[11];
extern const double gaussLobattoWeights11[12];
extern const double gaussLobattoWeights12[13];
extern const double gaussLobattoWeights13[14];
extern const double gaussLobattoWeights14[15];
extern const double gaussLobattoWeights15[16];
extern const double gaussLobattoWeights16[17];
extern const double gaussLobattoWeights17[18];
extern const double gaussLobattoWeights18[19];
extern const double gaussLobattoWeights19[20];
extern const double gaussLobattoWeights20[21];

extern const double gaussLobattoNodes0 [1];
extern const double gaussLobattoNodes1 [2];
extern const double gaussLobattoNodes2 [3];
extern const double gaussLobattoNodes3 [4];
extern const double gaussLobattoNodes4 [5];
extern const double gaussLobattoNodes5 [6];
extern const double gaussLobattoNodes6 [7];
extern const double gaussLobattoNodes7 [8];
extern const double gaussLobattoNodes8 [9];
extern const double gaussLobattoNodes9 [10];
extern const double gaussLobattoNodes10[11];
extern const double gaussLobattoNodes11[12];
extern const double gaussLobattoNodes12[13];
extern const double gaussLobattoNodes13[14];
extern const double gaussLobattoNodes14[15];
extern const double gaussLobattoNodes15[16];
extern const double gaussLobattoNodes16[17];
extern const double gaussLobattoNodes17[18];
extern const double gaussLobattoNodes18[19];
extern const double gaussLobattoNodes19[20];
extern const double gaussLobattoNodes20[21];

/**
 * The Gauss-Lobatto weights mapped onto [0,1]. Array of arrays. The first
 * entry is the order, the second entry the Lobatto point.
 */
extern const double* const gaussLobattoWeights[20+1];

/**
 * The Gauss-Lobatto nodes mapped onto [0,1]. Array of arrays. The first entry
 * is the order, the second entry the Lobatto point.
 */
extern const double* const gaussLobattoNodes[20+1];
}

#endif //GAUSSLOBATTO_H_
