// This file was generated by the ExaHyPE toolkit.
// It will NOT be regenerated or overwritten.
// Please adapt it to your own needs.
// 
// ========================
//   www.exahype.eu
// ========================

#include <cmath>
#include <map>
#include <tuple>

#include "NavierStokesSolver_ADERDG.h"
#include "NavierStokesSolver_ADERDG_Variables.h"

#include "totalVariation.h"
#include "VarianceHelper.h"

#include "stableDiffusiveTimeStepSize.h"
#if DIMENSIONS == 2
#include "diffusiveRiemannSolver2d.h"
#elif DIMENSIONS == 3
#include "diffusiveRiemannSolver3d.h"
#endif

#include "kernels/aderdg/generic/Kernels.h"
#include "kernels/KernelUtils.h"

#include "Scenarios/ScenarioFactory.h"
#include "Scenarios/Atmosphere.h"


tarch::logging::Log NavierStokes::NavierStokesSolver_ADERDG::_log( "NavierStokes::NavierStokesSolver_ADERDG" );

void NavierStokes::NavierStokesSolver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // TODO(Lukas) Refactor init!
  assert(constants.isValueValidString("scenario"));

  double referenceViscosity;
  if (constants.isValueValidString("viscosity") &&
      constants.getValueAsString("viscosity") == "default") {
   throw -1;
  } else {
    assert(constants.isValueValidDouble("viscosity"));
    referenceViscosity = constants.getValueAsDouble("viscosity");
  }

  scenarioName = constants.getValueAsString("scenario");
  scenario = ScenarioFactory::createScenario(scenarioName);

  const auto molecularDiffusionCoeff = scenario->getMolecularDiffusionCoeff();
  auto numberOfNecessaryVariables =
          1 + DIMENSIONS + 1;
  if (scenario->getUseAdvection()) {
    ++numberOfNecessaryVariables;
  }
  if (NumberOfVariables != numberOfNecessaryVariables) {
    throw -1;
  }

  ns = PDE(referenceViscosity, *scenario);
}

void NavierStokes::NavierStokesSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
  if (tarch::la::equals(t, 0.0)) {
    Variables vars(Q);
    scenario->initialValues(x, ns, vars);
    for (int i = 0; i < vars.variables(); ++i) {
      assertion2(std::isfinite(Q[i]), i, Q[i]);
    }
  }
}

void NavierStokes::NavierStokesSolver_ADERDG::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
   scenario->source(x, t, ns, Q, S);
}

void NavierStokes::NavierStokesSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
					const double * const fluxIn,const double* const stateIn, const double* const gradStateIn,
  double *fluxOut,double* stateOut) {
  constexpr auto basisSize = Order + 1;
  constexpr auto gradSize = NumberOfVariables * DIMENSIONS;

  auto gradStateOut = std::array<double, gradSize>{{0.0}};
  kernels::idx2 idxGradQ(DIMENSIONS,NumberOfVariables);

  std::fill_n(fluxOut, NumberOfVariables, 0.0);
  std::fill_n(stateOut, NumberOfVariables, 0.0);

  double _F[DIMENSIONS][NumberOfVariables]={0.0};
#if DIMENSIONS == 2
  double* F[2] = {_F[0], _F[1]};
#elif DIMENSIONS == 3
  double* F[3] = {_F[0], _F[1], _F[2]};
#endif

  for (int i = 0; i < NumberOfVariables; i++) {
    assertion2(std::isfinite(stateIn[i]), stateIn[i], i);
  }
  assertion1(stateIn[0] > 0, stateIn[0]);

  if (scenario->getBoundaryType(faceIndex) == NavierStokes::BoundaryType::analytical) {
    // Integrate over time.
    auto curStateOut = std::array<double, NumberOfVariables>{0.0};
    Variables curVarsOut(curStateOut.data());
    for (int i = 0; i < basisSize; ++i) {
      // TODO(Lukas): Check if we need to reset this data here.
      std::fill(curStateOut.begin(), curStateOut.end(), 0.0);
      std::fill(gradStateOut.begin(), gradStateOut.end(), 0.0);

      const double weight = kernels::gaussLegendreWeights[Order][i];
      const double xi = kernels::gaussLegendreNodes[Order][i];
      const double ti = t + xi * dt;

      scenario->analyticalSolution(x, ti, ns, curVarsOut, gradStateOut.data());

      ns.evaluateFlux(curStateOut.data(), gradStateOut.data(), F);

      for (int j = 0; j < NumberOfVariables; ++j) {
        stateOut[j] += weight * curStateOut[j];
        fluxOut[j] += weight * F[normalNonZero][j];
      }

    }
    return;
  }

  assertion(scenario->getBoundaryType(faceIndex) == BoundaryType::wall ||
                 scenario->getBoundaryType(faceIndex) == BoundaryType::hydrostaticWall);

  // Set no slip wall boundary conditions.
  ReadOnlyVariables varsIn(stateIn);
  Variables varsOut(stateOut);

  // Rho/E extrapolated, velocity mirrored.
  // Leads to zero velocity after Riemann solver.
  std::copy_n(stateIn, NumberOfVariables, stateOut);
  varsOut.j(0) = -varsIn.j(0);
  varsOut.j(1) = -varsIn.j(1);
#if DIMENSIONS == 3
  varsOut.j(2) = -varsIn.j(2);
#endif
  // TODO(Lukas) Refactor these checks.
  if (scenario->getQ0() > 0) {
    varsOut[NumberOfVariables-1] = varsIn[NumberOfVariables-1];
  }

  // Extrapolate gradient.
  std::copy_n(gradStateIn, gradSize, gradStateOut.data());

  // We deal with heat conduction by computing the flux at the boundary without heat conduction,
  // To do this, we reconstruct the incoming flux using the extrapolated/time-averaged state/gradient.
  // Note that this incurs an error.
  // The incoming flux is reconstructed in boundaryConditions.

  // Then compute the outgoing flux.
  if (scenario->getBoundaryType(faceIndex) == BoundaryType::hydrostaticWall) {
    // We need to reconstruct the temperature gradient here.
    // TODO(Lukas) Put those constants into the scenario.
    const auto g = 9.81;
    const auto backgroundPotTemperature = 300;
#if DIMENSIONS == 2
    const auto posZ = x[1];
#else
    const auto posZ = x[2];
#endif

    // In case of flow over a background state that is in hydrostatic equilibrium
    // it becomes necessary to reconstruct the temperature diffusion flux and
    // the energy. Otherwise a small temperature boundary layer forms.
    const double equilibriumTemperatureGradient = computeHydrostaticTemperatureGradient(ns, g, posZ,
            backgroundPotTemperature);

    // TODO(Lukas) Also reconstruct energy?
    const auto pressure = computeHydrostaticPressure(ns, g, posZ, backgroundPotTemperature);
    const auto T = potentialTToT(ns, pressure, backgroundPotTemperature);
    const auto rho = pressure / (ns.gasConstant * T);
    // TODO(Lukas) Is this also an accurate reconstruction for advection-scenarios?
    const auto E = ns.evaluateEnergy(rho, pressure, varsOut.j(), ns.getZ(stateIn));
    varsOut.E() = E;
    varsOut.rho() = rho;

    // Use no viscous effects and use equilibrium temperature gradient.
    ns.evaluateFlux(stateOut, gradStateOut.data(), F, false, true, equilibriumTemperatureGradient);
  } else {
    ns.evaluateFlux(stateOut, gradStateOut.data(), F);
  }

  std::copy_n(F[normalNonZero], NumberOfVariables, fluxOut);

}

bool NavierStokes::NavierStokesSolver_ADERDG::isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,
      const bool wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t, const double dt) const {
    // TODO(Lukas) What about advection-reaction-diffusion?
    return observablesMin[0] > 0.0 && observablesMin[DIMENSIONS + 2] > 0.0;
}

void NavierStokes::NavierStokesSolver_ADERDG::mapDiscreteMaximumPrincipleObservables(double* observables,const int numberOfObservables,const double* const Q) const {
  if (numberOfObservables>0) {
    std::copy_n(Q,numberOfObservables,observables);
  }
}


exahype::solvers::Solver::RefinementControl NavierStokes::NavierStokesSolver_ADERDG::refinementCriterion(
    const double* luh,
    const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx,
    const double t,
    const int level) {
  const bool isAmrScenario =
          scenarioName == "two-bubbles" ||
          scenarioName == "density-current" ||
          scenarioName == "taylor-green" ||
          scenarioName == "coupling-test";
  if (!isAmrScenario || DIMENSIONS != 2 || _globalObservables.size() < 2) {
    return exahype::solvers::Solver::RefinementControl::Keep;
  }

  if (t == 0) {
    // Global observables are reduced after first timestep!
    return exahype::solvers::Solver::RefinementControl::Keep;
  }

  const auto curObservables = mapGlobalObservables(luh, dx);
  const auto curTv = curObservables[0];

  const auto countGlobal = _globalObservables[2];
  const auto meanGlobal = _globalObservables[0];

  // Computed variance applies Bessel's correction.
  // As we use the complete population to compute the statistics
  // this is unnecessary.
  const auto varianceGlobal = (countGlobal - 1)/countGlobal * _globalObservables[1];
  const auto stdGlobal = std::sqrt(varianceGlobal);

  const auto factorRefine = 1.0;
  const auto factorCoarse = 0.5;

  const auto hi = meanGlobal + factorRefine * stdGlobal;
  const auto lo = meanGlobal + factorCoarse * stdGlobal;

  if (curTv > hi) {
    return exahype::solvers::Solver::RefinementControl::Refine;
  }

  if (curTv < lo) {
    return exahype::solvers::Solver::RefinementControl::Erase;
  }

  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void NavierStokes::NavierStokesSolver_ADERDG::eigenvalues(const double* const Q,const int d,double* lambda) {
  ns.evaluateEigenvalues(Q, d, lambda);
}

void NavierStokes::NavierStokesSolver_ADERDG::viscousEigenvalues(const double* const Q,const int d,double* lambda) {
  ns.evaluateDiffusiveEigenvalues(Q, d, lambda);
}

void NavierStokes::NavierStokesSolver_ADERDG::viscousFlux(const double *const Q, const double *const gradQ, double **F) {
  ns.evaluateFlux(Q, gradQ, F, true);
}

double NavierStokes::NavierStokesSolver_ADERDG::stableTimeStepSize(const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx) {
  // TODO(Lukas) Integrate diffusive time step size into standard timestep size!
  return (0.7/0.9) * stableDiffusiveTimeStepSize<NavierStokesSolver_ADERDG>(*static_cast<NavierStokesSolver_ADERDG*>(this),luh,dx);
}

void NavierStokes::NavierStokesSolver_ADERDG::riemannSolver(double* FL,double* FR,const double* const QL,const double* const QR,const double dt,const tarch::la::Vector<DIMENSIONS, double>& lengthScale, const int direction, bool isBoundaryFace, int faceIndex) {
  assertion2(direction>=0,dt,direction);
  assertion2(direction<DIMENSIONS,dt,direction);
  // TODO(Lukas) Integrate Riemann solver changes into standard solver.
  riemannSolverNonlinear<false,NavierStokesSolver_ADERDG>(*static_cast<NavierStokesSolver_ADERDG*>(this),FL,FR,QL,QR,lengthScale, dt,direction);

}

void NavierStokes::NavierStokesSolver_ADERDG::boundaryConditions( double* const fluxIn, const double* const stateIn, const double* const gradStateIn, const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& cellCentre, const tarch::la::Vector<DIMENSIONS,double>&  cellSize, const double t,const double dt, const int direction, const int orientation) {
  constexpr int basisSize     = (Order+1);
  constexpr int sizeStateOut = (NumberOfVariables+NumberOfParameters)*basisSize;
  constexpr int sizeFluxOut  = NumberOfVariables*basisSize;

  constexpr int totalSize = sizeStateOut + sizeFluxOut;
  double* block = new double[totalSize];

  double* memory = block;

  double* stateOut = memory; memory+=sizeStateOut;
  double* fluxOut  = memory; memory+=sizeFluxOut;

  const int faceIndex = 2*direction+orientation;

  kernels::aderdg::generic::c::boundaryConditions<true, NavierStokesSolver_ADERDG>(*static_cast<NavierStokesSolver_ADERDG*>(this),fluxOut,stateOut,fluxIn,stateIn,gradStateIn, cellCentre,cellSize,t,dt,faceIndex,direction);

  if (orientation==0) {
    double* FL = fluxOut; const double* const QL = stateOut;
    double* FR =  fluxIn;  const double* const QR = stateIn;

    riemannSolverNonlinear<false,NavierStokesSolver_ADERDG>(*static_cast<NavierStokesSolver_ADERDG*>(this),FL,FR,QL,QR,cellSize, dt,direction);
  }
  else {
    double* FL =  fluxIn;  const double* const QL = stateIn;
    double* FR = fluxOut; const double* const QR = stateOut;

    riemannSolverNonlinear<false,NavierStokesSolver_ADERDG>(*static_cast<NavierStokesSolver_ADERDG*>(this),FL,FR,QL,QR,cellSize, dt,direction);
  }

  if (scenario->getBoundaryType(faceIndex) == NavierStokes::BoundaryType::wall) {
    static_assert(DIMENSIONS == 2, "BC only implemented for 2D!"); // TODO(Lukas) Implement for 3D
    kernels::idx2 idx_F(Order + 1, NumberOfVariables);
    for (int i = 0; i < (Order + 1); ++i) {
      // Set energy flux to zero!
      fluxIn[idx_F(i, NavierStokesSolver_ADERDG_Variables::shortcuts::E)] = 0.0;
    }
  }

  delete[] block;
}

std::vector<double> NavierStokes::NavierStokesSolver_ADERDG::mapGlobalObservables(const double *const Q,
        const tarch::la::Vector<DIMENSIONS,double>& dx) const {
 if (NumberOfGlobalObservables == 0) return {};

 auto observables = resetGlobalObservables();
 // TODO(Lukas) Implement global observables for 3D!
 const auto idxQ = kernels::idx3(Order+1,Order+1,NumberOfVariables + NumberOfParameters);


 auto tv = 0.0;
 if (scenarioName == "two-bubbles" ||
     scenarioName == "density-current" ||
     scenarioName == "coupling-test") {
   auto computePotT = [this](const double *const Q) {
     const auto vars = ReadOnlyVariables{Q};
     const auto pressure = ns.evaluatePressure(vars.E(), vars.rho(), vars.j(), ns.getZ(Q));
     const auto temperature = ns.evaluateTemperature(vars.rho(), pressure);
     return ns.evaluatePotentialTemperature(temperature, pressure);
   };

   tv = totalVariation(Q, Order, NumberOfVariables, NumberOfParameters, dx, false, computePotT);
 } else {
   auto computeIndicator = [this](const double *const Q) {
     const auto vars = ReadOnlyVariables{Q};
     const auto pressure = ns.evaluatePressure(vars.E(), vars.rho(), vars.j(), ns.getZ(Q));
     return pressure;
   };
   tv = totalVariation(Q, Order, NumberOfVariables, NumberOfParameters, dx, false, computeIndicator);
 }

 return {tv, 0, 1};
}

std::vector<double> NavierStokes::NavierStokesSolver_ADERDG::resetGlobalObservables() const {
  if (NumberOfGlobalObservables == 0) return {};
  return {-1.0, -1.0, 0};
}

void NavierStokes::NavierStokesSolver_ADERDG::reduceGlobalObservables(
        std::vector<double> &reducedGlobalObservables,
        const std::vector<double> &curGlobalObservables) const {
  if (NumberOfGlobalObservables == 0) return;

  assertion2(reducedGlobalObservables.size() == curGlobalObservables.size(),
          reducedGlobalObservables.size(),
          curGlobalObservables.size());

  const auto mean0 = reducedGlobalObservables[0];
  const auto mean1 = curGlobalObservables[0];
  const auto var0 = reducedGlobalObservables[1];
  const auto var1 = curGlobalObservables[1];
  const auto count0 = static_cast<int>(reducedGlobalObservables[2]);
  const auto count1 = static_cast<int>(curGlobalObservables[2]);

  auto mergedMean = 0.0;
  auto mergedVariance = 0.0;

  std::tie(mergedMean, mergedVariance) = mergeVariance(mean0, mean1, var0, var1, count0,
          count1);
  reducedGlobalObservables[0] = mergedMean;
  reducedGlobalObservables[1] = mergedVariance;
  reducedGlobalObservables[2] = count0 + count1;
}
