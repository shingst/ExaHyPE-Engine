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

#include <algorithm> //copy_n

#include "kernels/SWE_MySWESolver_ADERDG/Kernels.h"
#include "kernels/SWE_MySWESolver_ADERDG/DGMatrices.h"
#include "kernels/SWE_MySWESolver_ADERDG/Quadrature.h"

#include "kernels/SWE_MySWESolver_ADERDG/gemmsCPP.h"

// local help function
inline int powOf3(int exp){
  switch(exp) {
    case 0:
      return 1;
    case 1:
      return 3;
    case 2:
      return 9;
    case 3:
      return 27;
    case 4:
      return 81;
    case 5:
      return 243;
    default:
      int result = 243;
      for (int d=0; d<exp-5; d++) {
        result *= 3;
      }
      return result;
  }
}



void SWE::MySWESolver_ADERDG_kernels::aderdg::faceUnknownsProlongation(
    double* restrict lQhbndFine,
    double* restrict lFhbndFine,
    const double* const restrict lQhbndCoarse,
    const double* const restrict lFhbndCoarse,
    const int coarseGridLevel,
    const int fineGridLevel,
    const int* const subfaceIndex
) {
  const int levelDelta = fineGridLevel - coarseGridLevel;

  // tmp arrays
  double tmpQ[8] __attribute__((aligned(ALIGNMENT)));
  double tmpF[8] __attribute__((aligned(ALIGNMENT)));

  // read only input, start with the function input = coarse
  const double* inputQ = lQhbndCoarse;
  const double* inputF = lFhbndCoarse;

  // output pointer, ensures that the output of the last iteration points to the function output
  double* outputQ;
  double* outputF;
  if (levelDelta % 2 == 0) {
    outputQ = tmpQ;
    outputF = tmpF;
  } else {
    outputQ = lQhbndFine;
    outputF = lFhbndFine;
  }

  int subfaceIndexPrevious_0 = subfaceIndex[0];
  int subfaceIndexCurrent_0;
  int subintervalIndex_0;

  // This loop decodes the elements of subfaceIndex into a tertiary basis
  // starting with the highest significance 3^(levelDelta-1).
  // 
  // Per iteration, the digits corresponding to the current significances then determine
  // the subintervals for the single level prolongation.
  for (int l = 1; l < levelDelta+1; l++) {
    const int significance = powOf3(levelDelta-l);
    subfaceIndexCurrent_0 = subfaceIndexPrevious_0 % significance;
    subintervalIndex_0    = (subfaceIndexPrevious_0 - subfaceIndexCurrent_0)/significance;
    
    // Apply the single level prolongation operator.
    // Use the coarse level unknowns as input in the first iteration.
    // will overwrite outputs, no need to set to 0
    #ifdef USE_IPO
    #pragma forceinline
    #endif
    gemm_4_2_2_face_Q_x(inputQ, fineGridProjector1d[subintervalIndex_0], outputQ);
    #ifdef USE_IPO
    #pragma forceinline
    #endif
    gemm_4_2_2_face_F_x(inputF, fineGridProjector1d[subintervalIndex_0], outputF);
    
    // Prepare next iteration.
    subfaceIndexPrevious_0 = subfaceIndexCurrent_0;

    // Input is previous output
    inputQ = outputQ;
    inputF = outputF;
    
    // Toggle the addresses of the pointers.
    if (outputQ == tmpQ) {
      outputQ = lQhbndFine;
      outputF = lFhbndFine;
    } else {
      outputQ = tmpQ;
      outputF = tmpF;
    }
  }
    
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::faceFluxRestriction(
    double* restrict lFhbndCoarse,
    const double* const restrict lFhbndFine,
    const int* const subfaceIndex,
    const int levelDelta
) {
  
  // tmp array, only allocated if needed (more than one level)
  double tmpF[8] __attribute__((aligned(ALIGNMENT)));

  // read only input, start with the function input = fine
  const double* inputF = lFhbndFine;

  // output pointer, ensures that the output of the last iteration points to the function output
  double* outputF;
  if (levelDelta % 2 == 0) {
    outputF = tmpF;
  } else {
    outputF = lFhbndCoarse;
  }

  int subfaceIndexCurrent_0 = subfaceIndex[0];
  int subintervalIndex_0;
  
  // This loop decodes the indices of subfaceIndex into a tertiary basis
  // starting with the lowest significance 3^0 (in contrast to the prolongation loop).
  //
  // Per iteration, the digits corresponding to the current significances then determine
  // the subintervals for the single level restriction.
  for (int l = 1; l < levelDelta+1; l++) {
    subintervalIndex_0    = subfaceIndexCurrent_0 % 3;  
    subfaceIndexCurrent_0 = (subfaceIndexCurrent_0 - subintervalIndex_0)/3;
    // Apply the single level restriction operator.
    // Use the fine level unknowns as input in the first iteration.
    // will overwrite outputs, no need to set to 0
    #ifdef USE_IPO
    #pragma forceinline
    #endif
    gemm_4_2_2_face_F_x(inputF, fineGridProjector1d_T_weighted[subintervalIndex_0], outputF);

    // Prepare next iteration.
    inputF = outputF;
    // Toggle pointer pairs.
    if (outputF == tmpF) {
      outputF = lFhbndCoarse;
    } else {
      outputF = tmpF;
    }
  }
  
}


void SWE::MySWESolver_ADERDG_kernels::aderdg::volumeUnknownsProlongation(
    double* restrict luhFine,
    const double* const restrict luhCoarse,
    const int coarseGridLevel,
    const int fineGridLevel,
    const int* const subcellIndex
) {
  const int levelDelta = fineGridLevel - coarseGridLevel;

  double tmpLuh[16] __attribute__((aligned(ALIGNMENT)));
  double tmpX[16] __attribute__((aligned(ALIGNMENT)));
  
  // read only input, start with the function input = fine
  const double* inputLuh = luhCoarse;

  // output pointer, ensures that the output of the last iteration points to the function output
  double* outputLuh;
  if (levelDelta % 2 == 0) {
    outputLuh = tmpLuh;
  } else {
    outputLuh = luhFine;
  }

  int subcellIndexPrevious_0 = subcellIndex[0];
  int subcellIndexCurrent_0;
  int subintervalIndex_0;
  int subcellIndexPrevious_1 = subcellIndex[1];
  int subcellIndexCurrent_1;
  int subintervalIndex_1;

  // This loop step by step decodes the elements of subcellIndex into a tertiary basis
  // starting with the highest significance 3^(levelDelta-1).
  // 
  // Per iteration, the digits corresponding to the current significances then determine
  // the subintervals for the single level prolongation.
  for (int l = 1; l < levelDelta+1; l++) {
    const int significance = powOf3(levelDelta-l);
    subcellIndexCurrent_0 = subcellIndexPrevious_0 % significance;
    subintervalIndex_0    = (subcellIndexPrevious_0 - subcellIndexCurrent_0)/significance;
    subcellIndexCurrent_1 = subcellIndexPrevious_1 % significance;
    subintervalIndex_1    = (subcellIndexPrevious_1 - subcellIndexCurrent_1)/significance;

    // will overwrite tmpX, no need to set to 0
    for (int zy = 0; zy < 2; zy++) {
      #ifdef USE_IPO
      #pragma forceinline
      #endif
      gemm_4_2_2_volume_x(inputLuh+zy*8, fineGridProjector1d[subintervalIndex_0], tmpX+zy*8);
    }
    
    for (int x = 0; x < 2; x++) {
      #ifdef USE_IPO
      #pragma forceinline
      #endif
      gemm_4_2_2_volume_y(tmpX+x*4, fineGridProjector1d[subintervalIndex_1], outputLuh+x*4);
    }

    // Prepare next iteration.
    subcellIndexPrevious_0 = subcellIndexCurrent_0;
    subcellIndexPrevious_1 = subcellIndexCurrent_1;

    inputLuh = outputLuh;

    // Toggle pointers.
    if (outputLuh == tmpLuh) {
      outputLuh = luhFine;
    } else {
      outputLuh = tmpLuh;
    }
  }
    
}


void SWE::MySWESolver_ADERDG_kernels::aderdg::volumeUnknownsRestriction(
    double* restrict luhCoarse,
    const double* const restrict luhFine,
    const int coarseGridLevel,
    const int fineGridLevel,
    const int* const subcellIndex
) {
  const int levelDelta = fineGridLevel - coarseGridLevel;
  
  // read only input, start with the function input = fine
  const double* inputLuh = luhFine;
  int subintervalIndex_0 = subcellIndex[0];
  int subintervalIndex_1 = subcellIndex[1];

  double tmpX[16] __attribute__((aligned(ALIGNMENT)));


  if(levelDelta > 1) {
    double tmpLuh[16] __attribute__((aligned(ALIGNMENT)));
    double tmpLuh2[16] __attribute__((aligned(ALIGNMENT)));
    int subcellIndexCurrent_0 = subcellIndex[0];
    int subcellIndexCurrent_1 = subcellIndex[1];
    // output pointer, ensures that the output of the last iteration points to the function output
    double* outputLuh = tmpLuh;
    // This loop step by step decodes the elements of subcellIndex into a tertiary basis
    // starting with the highest significance 3^(levelDelta-1).
    // 
    // Per iteration, the digits corresponding to the current significances then determine
    // the subintervals for the single level prolongation.
    for (int l = 1; l < levelDelta; l++) {
      subintervalIndex_0    = subcellIndexCurrent_0 % 3;
      subcellIndexCurrent_0 = (subcellIndexCurrent_0 - subintervalIndex_0)/3;
      subintervalIndex_1    = subcellIndexCurrent_1 % 3;
      subcellIndexCurrent_1 = (subcellIndexCurrent_1 - subintervalIndex_1)/3;

      // will overwrite tmpX, no need to set to 0
      for (int zy = 0; zy < 2; zy++) {
        #ifdef USE_IPO
        #pragma forceinline
        #endif
        gemm_4_2_2_volume_x(inputLuh+zy*8, fineGridProjector1d_T_weighted[subintervalIndex_0], tmpX+zy*8);
      }
      
      for (int x = 0; x < 2; x++) {
        #ifdef USE_IPO
        #pragma forceinline
        #endif
        gemm_4_2_2_volume_y(tmpX+x*4, fineGridProjector1d_T_weighted[subintervalIndex_1], outputLuh+x*4);
      }

      inputLuh = outputLuh;
      // Toggle pointers.
      if (outputLuh == tmpLuh) {
        outputLuh = tmpLuh2;
      } else {
        outputLuh = tmpLuh;
      }
    }
    
    subintervalIndex_0    = subcellIndexCurrent_0 % 3;
    subintervalIndex_1    = subcellIndexCurrent_1 % 3;
    
  } // if leveldelta>1, 
  
  //now at the case 1 level to do
  
  // will overwrite tmpX, no need to set to 0
  for (int zy = 0; zy < 2; zy++) {
    #ifdef USE_IPO
    #pragma forceinline
    #endif
    gemm_4_2_2_volume_x(inputLuh+zy*8, fineGridProjector1d_T_weighted[subintervalIndex_0], tmpX+zy*8);
  }
  
  // Add to the coarse output
  for (int x = 0; x < 2; x++) {
    #ifdef USE_IPO
    #pragma forceinline
    #endif
    gemm_4_2_2_volume_y_add(tmpX+x*4, fineGridProjector1d_T_weighted[subintervalIndex_1], luhCoarse+x*4);
  }
  

}