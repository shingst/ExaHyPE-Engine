// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "Plotter.h"
#include <kernels/GaussLegendreQuadrature.h>
#include "NavierStokesSolverDG.h"
#include "NavierStokesSolverDG_Variables.h"
#include "PDE.h"

NavierStokes::Plotter::Plotter(NavierStokes::NavierStokesSolverDG& solver) :
        order(solver.Order), solver(&solver) {

}

NavierStokes::Plotter::~Plotter() {
}

void NavierStokes::Plotter::startPlotting( double time) {
  // @TODO Please insert your code here.
}


void NavierStokes::Plotter::finishPlotting() {
  // @TODO Please insert your code here.
}

void NavierStokes::Plotter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  auto vars = Variables(Q);
  constexpr auto writtenUnknowns = vars.Size;

  for (int i = 0; i < writtenUnknowns; ++i){
    outputQuantities[i] = Q[i];
  }

  if (solver->scenarioName == "convergence" ||
          solver->scenarioName == "entropy-wave") {
    // Plot quadrature weights.
    // This is needed to approximate the integral of error norms.
    const auto& weights = kernels::gaussLegendreWeights[order];

    double weight = 1.0;
    for (int i = 0; i < DIMENSIONS; ++i) {
      weight *= weights[pos[i]];
    }
    outputQuantities[vars.Size] = weight;
  } else {
    // For other scenarios, plot the potential temperature.

    const auto& ns = solver->ns;

    const auto pressure = ns.evaluatePressure(vars.E(), vars.rho(), vars.j());
    const auto temperature = ns.evaluateTemperature(vars.rho(), pressure);

    const auto potT = ns.evaluatePotentialTemperature(temperature, pressure);

    // Write potential temperature
    outputQuantities[vars.Size] = potT;
  }
  auto& globalObservables = solver->getGlobalObservables();
  if (timeStamp > 0.0) {
    outputQuantities[vars.Size+1] = globalObservables[0];
    outputQuantities[vars.Size+2] = globalObservables[1];
  } else {
    outputQuantities[vars.Size+1] = 0.0;
    outputQuantities[vars.Size+2] = 0.0;
  }




}