#include "totalVariation.h"

#include <vector>

#include "kernels/aderdg/generic/Kernels.h"
#include "kernels/KernelUtils.h"
#include "kernels/GaussLegendreQuadrature.h"

double totalVariation(const double* Q, int order, int numberOfVariables, int numberOfParameters) {
  const auto basisSize = order + 1;
  const auto numberOfData = numberOfVariables + numberOfParameters;
  const auto sizeGradient = basisSize * basisSize * basisSize * numberOfVariables;

  auto idxGradQ = kernels::idx4(basisSize, basisSize, DIMENSIONS, numberOfVariables);
  auto idx_lQi = kernels::idx3(basisSize, basisSize, numberOfData); // idx_lQi(y,x,t,nVar+nPar)

  // TODO(Lukas) Use cell size here!
  auto dx = std::array<double, 2>{1.0, 1.0};
  auto invDx = std::array<double, 2>{1.0, 1.0};

  // First compute gradient.
  // Note that we could avoid this calculation but it's simpler this way.
  // TODO(Lukas): Compute integral of TV directly!
  auto gradQ = std::vector<double>(sizeGradient);
  // x direction (independent from the y derivatives)
  for (int k = 0; k < basisSize; k++) { // k == y
    // Matrix operation
    for (int l = 0; l < basisSize; l++) { // l == x
      for (int m = 0; m < numberOfVariables; m++) {
        for (int n = 0; n < basisSize; n++) { // n == matmul x
          const auto idx = idxGradQ(k,l, /*x*/0,m);
          const auto t = 1.0 * invDx[0] * Q[idx_lQi(k,n,m)] * kernels::dudx[order][l][n];
          gradQ[idx] += t;
        }
      }
    }
  }

  // y direction (independent from the x derivatives)
  for(int k=0; k<basisSize; k++) {
    // Matrix operation
    for (int l = 0; l < basisSize; l++) { // l == y
      for (int m = 0; m < numberOfVariables; m++) {
        for (int n = 0; n < basisSize; n++) { // n = matmul y
          const auto idx = idxGradQ(l, k, /*y*/1, m);
          const auto t = 1.0 * invDx[1] * Q[idx_lQi(n, k, m)] * kernels::dudx[order][l][n]; /* l,n: transpose */
          gradQ[idx] += t;
        }
      }
    }
  }

  // Compute integral of absolute gradient for all dimensions and variables.
  auto tv = 0.0;
  // x:
  const auto& quadratureWeights = kernels::gaussLegendreWeights[order];

  for (int k = 0; k < basisSize; k++) {
    for (int l = 0; l < basisSize; l++) {
      auto w = quadratureWeights[k] * quadratureWeights[l];
        for (int m = 0; m < DIMENSIONS; m++) {
          for (int n = 0; n < numberOfVariables; n++) {
            tv += w * std::abs(gradQ[idxGradQ(k,l,m,n)]);
          }
        }
    }
  }

  // TODO(Lukas) Use correct volume!
  const auto volume = dx[0] * dx[1];
  return volume * tv;
}