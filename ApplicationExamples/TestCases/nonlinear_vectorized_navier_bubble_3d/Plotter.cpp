// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "Plotter.h"
#include <kernels/GaussLegendreBasis.h>
#include "NavierStokesSolver_ADERDG.h"
#include "NavierStokesSolver_ADERDG_Variables.h"
#include "PDE.h"
//#include "AMR/totalVariation.h"

NavierStokes::Plotter::Plotter(NavierStokes::NavierStokesSolver_ADERDG& solver) :
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
  const auto vars = ReadOnlyVariables{Q};

  const auto rho = NavierStokesSolver_ADERDG_Variables::shortcuts::rho;
  const auto j = NavierStokesSolver_ADERDG_Variables::shortcuts::j;
  const auto E = NavierStokesSolver_ADERDG_Variables::shortcuts::E;
  const auto Z = E + 1; // Only defined if coupling is used!

  // TODO(Lukas) Refactor plotting!
  const bool writePrimitive = true;

  const auto& ns = solver->ns;
  // TODO(Lukas) Fix pott for coupled eq.
  const bool writePotT = !ns.useAdvection && (solver->scenarioName != "lid-driven-cavity");
  const auto pressure = ns.evaluatePressure(vars.E(), vars.rho(), vars.j(),
                                            ns.getZ(Q), ns.getHeight(Q));
  if (writePrimitive) {
    // Primitive variables are density, velocities and pressure.
    outputQuantities[rho] = Q[rho];
    outputQuantities[j+0] = Q[j+0]/Q[rho];
    outputQuantities[j+1] = Q[j+1]/Q[rho];
#if DIMENSIONS == 3
    outputQuantities[j+2] = Q[j+2]/Q[rho];
#endif
    outputQuantities[E] = pressure;
    if (ns.useAdvection) {
      outputQuantities[Z] = Q[Z] / Q[rho];
    }
  } else {
    // Otherwise just plot all conservative variables.
    // (no conversion needed)
    for (int i = 0; i < vars.size(); ++i){
      outputQuantities[i] = Q[i];
    }
  }

  if (writePotT) {
    const auto temperature = ns.evaluateTemperature(vars.rho(), pressure);
    const auto potT = ns.evaluatePotentialTemperature(temperature, pressure);

    // Write potential temperature
    outputQuantities[vars.SizeVariables] = potT;
  } else {
    outputQuantities[vars.SizeVariables] = 0.0;
  }
}