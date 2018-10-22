
// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "ErrorWriter.h"
#include "NavierStokesSolverDG_Variables.h"
#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/KernelUtils.h"

#include "peano/utils/Loop.h"

#include "tarch/la/VectorOperations.h"

#include <algorithm>

#include <fstream>
#include <iomanip>
#include <sstream>

// TODO(Lukas) Do not hardcore this, change to proper conv. scenario!
#include "PDE.h"
#include "Scenarios/ConvergenceTest/ConvergenceTest.h"
#include "Scenarios/EntropyWave.h"

NavierStokes::ErrorWriter::ErrorWriter()
    : exahype::plotters::ADERDG2UserDefined::ADERDG2UserDefined(),
      hmin(std::numeric_limits<double>::max()) {}

void NavierStokes::ErrorWriter::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
    double timeStamp) {
  constexpr int numberOfVariables = NavierStokesSolverDG::NumberOfVariables;
  constexpr int numberOfParameters = NavierStokesSolverDG::NumberOfParameters;
  constexpr int numberOfData = numberOfVariables + numberOfParameters;
  constexpr int order = NavierStokesSolverDG::Order;
  constexpr int basisSize = order + 1;
  constexpr int gradSize = DIMENSIONS * numberOfVariables;

  static_assert(DIMENSIONS == 2, "ErrorWriter only supports 2D");
  double x[2] = {0.0, 0.0};

  auto scenario = ConvergenceTest();
  // TODO(Lukas): Change to viscosity that is actually used(?).
  const auto viscosity = 0.1;
  const auto ns = PDE(viscosity, scenario);

  kernels::idx4 idx(basisSize, basisSize, basisSize, numberOfData);
  dfor(i, basisSize) {
    double w_dV = 1.0;
    for (int d = 0; d < DIMENSIONS; d++) {
      x[d] = offsetOfPatch[d] +
             sizeOfPatch[d] * kernels::gaussLegendreNodes[order][i(d)];
      w_dV *= sizeOfPatch[d] * kernels::gaussLegendreWeights[order][i(d)];
      hmin =
          std::min(hmin, sizeOfPatch[d]);  // TODO(Lukas) Is this what we want?
    }

    // Compute analytical solution
    auto uAna = std::array<double, numberOfVariables>{};
    auto uAnaGrad = std::array<double, gradSize>{};
    auto vars = Variables(uAna.data());
    scenario.analyticalSolution(x, timeStamp, ns, vars, uAnaGrad.data());

    const double* uNum;
    if (DIMENSIONS == 3) {
      uNum = u + idx(i(2), i(1), i(0), 0);
    } else {
      uNum = u + idx(0, i(1), i(0), 0);
    }

    for (int v = 0; v < numberOfVariables; v++) {
      const double uDiff = std::abs(uNum[v] - uAna[v]);
      errorL2[v] += uDiff * uDiff * w_dV;
      errorL1[v] += uDiff * w_dV;
      errorLInf[v] = std::max(errorLInf[v], uDiff);

      normL1Ana[v] += std::abs(uAna[v]) * w_dV;
      normL2Ana[v] += uAna[v] * uAna[v] * w_dV;
      normLInfAna[v] = std::max(normLInfAna[v], std::abs(uAna[v]));
    }
  }
}

void NavierStokes::ErrorWriter::startPlotting(double time) {
  timeStamp = time;

  std::fill(errorL1.begin(), errorL1.end(), 0.0);
  std::fill(errorL2.begin(), errorL2.end(), 0.0);
  std::fill(errorLInf.begin(), errorLInf.end(), 0.0);

  std::fill(normL1Ana.begin(), normL1Ana.end(), 0.0);
  std::fill(normL2Ana.begin(), normL2Ana.end(), 0.0);
  std::fill(normLInfAna.begin(), normLInfAna.end(), 0.0);
}

void NavierStokes::ErrorWriter::finishPlotting() {
  constexpr int numberOfVariables =
      AbstractNavierStokesSolverDG::NumberOfVariables;

  // Check whether our file exists already.
  const bool isExisting = [&]() -> bool {
    auto file = std::ifstream(filename);
    return file.good();
  }();

  auto file = std::ofstream(filename, std::ios::app);
  if (!file.is_open()) {
    throw std::runtime_error("ErrorWriter: Error while opening file " +
                             filename + "!");
  }
  assert(file.is_open());

  // TODO(Lukas) Is this enough precision?
  const auto precision = std::numeric_limits<double>::max_digits10 + 2;

  // Write csv header in case of new file.
  if (!isExisting) {
    file << std::setprecision(precision) << std::fixed;
    file << "norm,order,hmin,";
    for (int i = 0; i < numberOfVariables; ++i) {
      file << "var" << i << ",";
    }
    file << "time" << std::endl;
  }

  auto plotRow = [&](const std::string& normType, const Array_t arr) {
    file << std::setprecision(precision) << std::fixed << normType << ","
         << AbstractNavierStokesSolverDG::Order << "," << hmin << ",";

    for (int i = 0; i < numberOfVariables; ++i) {
      file << std::setprecision(precision) << std::fixed << arr[i] << ",";
    }
    file << timeStamp << std::endl;
  };

  // Note that we write the l2-norms without applying the square root.
  // The reason is expected precision loss when using MPI.
  // (need to aggregate results of multiple nodes!)
  plotRow("l1", errorL1);
  plotRow("l2", errorL2);
  plotRow("lInf", errorLInf);

  plotRow("l1Norm", normL1Ana);
  plotRow("l2Norm", normL2Ana);
  plotRow("lInfNorm", normLInfAna);

  // Verify that we actually wrote something.
  file << std::flush;  // Buffering might hide error
  if (!file.good()) {
    throw std::runtime_error("ErrorWriter: Error while writing to file " +
                             filename + "!");
  }

  // Now pretty print for console:
  auto prettyPrintRow = [&](const std::string& normType, const Array_t arr) {
    std::cout << normType << "-Error for order "
              << AbstractNavierStokesSolverDG::Order << " and hmin of " << hmin
              << " at time " << timeStamp << ":" << std::endl;

    for (int i = 0; i < numberOfVariables; ++i) {
      std::cout << std::setprecision(3) << std::scientific << arr[i] << "\t";
    }
    std::cout << std::endl;
  };

  // Compute square roots for pretty printing.
  for (int v = 0; v < numberOfVariables; v++) {
    errorL2[v] = sqrt(errorL2[v]);
    normL2Ana[v] = sqrt(normL2Ana[v]);
  }

  std::cout << "\n\n";
  prettyPrintRow("l1", errorL1);
  prettyPrintRow("l2", errorL2);
  prettyPrintRow("lInf", errorLInf);
  prettyPrintRow("l1Norm", normL1Ana);
  prettyPrintRow("l2Norm", normL2Ana);
  prettyPrintRow("lInfNorm", normLInfAna);
  std::cout << "\n\n";
}

void NavierStokes::ErrorWriter::init(
    const std::string& filename, int orderPlusOne, int solverUnknowns,
    int writtenUnknowns, exahype::parser::ParserView plotterParameters) {
  ADERDG2UserDefined::init(filename, orderPlusOne, solverUnknowns,
                           writtenUnknowns, plotterParameters);

  const auto& node = tarch::parallel::Node::getInstance();
  const auto rank = node.getRank();
  const auto no_nodes = node.getNumberOfNodes();
  isMpi = no_nodes > 1;

  auto fs = std::stringstream();
  fs << filename;
  if (isMpi) {
    fs << "_rank_" << rank;
  }
  fs << ".csv";

  this->filename = fs.str();
}
