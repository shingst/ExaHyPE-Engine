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

#include "AMR/Criterion.h"

#include "stableDiffusiveTimeStepSize.h"
#if DIMENSIONS == 2
#include "diffusiveRiemannSolver2d.h"
#elif DIMENSIONS == 3
#include "diffusiveRiemannSolver3d.h"
#endif

#include "kernels/aderdg/generic/Kernels.h"
#include "kernels/KernelUtils.h"

#include "Scenarios/Atmosphere.h"

#include "SetupHelper.h"

tarch::logging::Log NavierStokes::NavierStokesSolver_ADERDG::_log( "NavierStokes::NavierStokesSolver_ADERDG" );

void NavierStokes::NavierStokesSolver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  auto parsedConfig = parseConfig(cmdlineargs, constants, NumberOfVariables, NumberOfParameters, NumberOfGlobalObservables);
  ns = std::move(parsedConfig.ns);
  scenarioName = std::move(parsedConfig.scenarioName);
  scenario = std::move(parsedConfig.scenario);
  amrSettings = std::move(parsedConfig.amrSettings);
}

void NavierStokes::NavierStokesSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
  if (tarch::la::equals(t, 0.0)) {
    // TODO(Lukas) What happens during refinement?
    Q[NumberOfVariables + NumberOfParameters - 1] = 0.0; // TODO(Lukas) Remove!
    ns.setHeight(Q, x[DIMENSIONS-1]);
    ns.setBackgroundState(Q, 0.0, 0.0);

    Variables vars(Q);
    scenario->initialValues(x, ns, vars);
    for (int i = 0; i < vars.variables(); ++i) {
      assertion2(std::isfinite(Q[i]), i, Q[i]);
    }
    //Q[NumberOfVariables+NumberOfParameters-1] = -1;
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
  constexpr auto NumberOfData = NumberOfVariables + NumberOfParameters;

  auto gradStateOut = std::array<double, gradSize>{{0.0}};
  kernels::idx2 idxGradQ(DIMENSIONS,NumberOfVariables);

  std::fill_n(fluxOut, NumberOfVariables, 0.0);
  std::fill_n(stateOut, NumberOfData, 0.0);

  ns.setHeight(stateOut, x[DIMENSIONS-1]);
  ns.setBackgroundState(stateOut, 0.0, 0.0);
  // Need to reconstruct the background state in this case.
  // Extrapolating does not work!
  assert(!ns.useBackgroundState || scenario->getBoundaryType(faceIndex) ==
      BoundaryType::hydrostaticWall);

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
                 scenario->getBoundaryType(faceIndex) == BoundaryType::hydrostaticWall ||
                 scenario->getBoundaryType(faceIndex) == BoundaryType::movingWall ||
                 scenario->getBoundaryType(faceIndex) == BoundaryType::freeSlipWall);

  // Set no slip wall boundary conditions.
  ReadOnlyVariables varsIn(stateIn);
  Variables varsOut(stateOut);

  // Rho/E extrapolated, velocity mirrored.
  std::copy_n(stateIn, NumberOfVariables, stateOut);
  // Extrapolate gradient.
  std::copy_n(gradStateIn, gradSize, gradStateOut.data());

  // TODO(Lukas) Are these gradients here correct?
  if (scenario->getBoundaryType(faceIndex) == BoundaryType::freeSlipWall ||
  scenario->getBoundaryType(faceIndex) == BoundaryType::hydrostaticWall) {
    // Normal velocity zero after Riemann.
    varsOut.j(normalNonZero) = -varsIn.j(normalNonZero);
    for (int i = 0; i < DIMENSIONS; ++i) {
      gradStateOut[idxGradQ(i, 1+normalNonZero)] = -gradStateIn[idxGradQ(i, 1+normalNonZero)];
    }
  } else {
    // No-slip
    // All velocities zero after Riemann.
    for (int i = 0; i < DIMENSIONS; ++i) {
      varsOut.j(i) = -varsIn.j(i);
      for (int j = 0; j < DIMENSIONS; ++j) {
        gradStateOut[idxGradQ(j, 1+i)] = -gradStateIn[idxGradQ(j, 1+i)];
      }
    }
  }

  if (scenario->getBoundaryType(faceIndex) == BoundaryType::movingWall) {
    // Wall speed after Riemann solve
    const auto wallSpeed = 1.0;
    varsOut.j(0) = 2 * wallSpeed - varsIn.j(0);
    // TODO(Lukas) Is this gradient correct?
    for (int i = 0; i < DIMENSIONS; ++i) {
      gradStateOut[idxGradQ(i, 0)] += 2 * wallSpeed;
    }
  }

  // We deal with heat conduction by computing the flux at the boundary without heat conduction,
  // To do this, we reconstruct the incoming flux using the extrapolated/time-averaged state/gradient.
  // Note that this incurs an error.
  // The incoming flux is reconstructed in boundaryConditions.

  // Then compute the outgoing flux.
  if ( scenario->getBoundaryType(faceIndex) == BoundaryType::hydrostaticWall) {
    // TODO(Lukas): Add support for advection-coupling + hydrostatic flows?
    // We need to reconstruct the temperature gradient here.
    const auto posZ = x[DIMENSIONS-1];

    // In case of flow over a background state that is in hydrostatic equilibrium
    // it becomes necessary to reconstruct the temperature diffusion flux and
    // the energy. Otherwise a small temperature boundary layer forms.
    const double equilibriumTemperatureGradient = computeHydrostaticTemperatureGradient(ns, scenario->getGravity(),
            posZ, scenario->getBackgroundPotentialTemperature());

    // We also need to reconstruct the temperature at the border.
    // This corresponds to a heated wall.
    const auto pressure = computeHydrostaticPressure(ns, scenario->getGravity(),
            posZ, scenario->getBackgroundPotentialTemperature());
    const auto T = potentialTToT(ns, pressure, scenario->getBackgroundPotentialTemperature());
    varsOut.rho() = pressure / (ns.gasConstant * T);

    ns.setBackgroundState(stateOut, varsOut.rho(), pressure);
    if (ns.useGravity) {
      varsOut.E() = ns.evaluateEnergy(varsOut.rho(), pressure, varsOut.j(), ns.getZ(stateIn),
              x[DIMENSIONS-1]);
    } else {
      varsOut.E() = ns.evaluateEnergy(varsOut.rho(), pressure, varsOut.j(), ns.getZ(stateIn));
    }

    // Use no viscous effects and use equilibrium temperature gradient.
    ns.evaluateFlux(stateOut, gradStateOut.data(), F, false, true, equilibriumTemperatureGradient);
  } else {
    //ns.evaluateFlux(stateOut, gradStateOut.data(), F);
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
      const double t) const {
    // We now need to do a pointwise check for the primitive variables
    // pressure and Z.
    // TODO(Lukas) At least refactor this. And 3D!
  // return false;
#if DIMENSIONS == 2
    constexpr auto basisSize = Order + 1;
    kernels::idx3 idx(basisSize,basisSize,NumberOfVariables);
    for (int i = 0; i < basisSize; ++i) {
      for (int j = 0; j < basisSize; ++j) {
        const double* const Q = solution + idx(i,j,0);
        auto vars = ReadOnlyVariables{Q};
        const auto Zrho = ns.getZ(Q);
        const auto Z = Zrho / vars.rho();
        const auto height = ns.getHeight(Q);
        const auto pressure = ns.evaluatePressure(vars.E(),
                vars.rho(),
                vars.j(),
                Zrho,
                height
                );
        bool isAdvectionTroubled = ns.useAdvection && (Z < 0.0);
        if (vars.rho() <= 0.0 || pressure < 0.0 || isAdvectionTroubled) {
          return false;
        }

        // Surprisingly, this is necessary.
        for (int v = 0; v < NumberOfVariables; v++) {
          if (!std::isfinite(solution[v])) {
            return false;
          }
        }
      }
    }
#else
    // TODO(Lukas) Limiting in 3D!
    std::abort();
#endif
  return true;
}

void NavierStokes::NavierStokesSolver_ADERDG::mapDiscreteMaximumPrincipleObservables(double* observables,const double* const Q) const {
  if (NumberOfDMPObservables > 0) {
    // TODO(Lukas) Remove this.
    std::fill_n(observables, NumberOfDMPObservables, 0.0);
    assert(NumberOfDMPObservables >= 2);
    observables[0] = Q[0];
    const auto vars = ReadOnlyVariables{Q};
    observables[1] = ns.evaluatePressure(vars.E(),
                                         vars.rho(),
                                         vars.j(),
                                         ns.getZ(Q),
                                         ns.getHeight(Q));
    if (ns.useAdvection) {
      observables[2] = ns.getZ(Q) / Q[0];
    }
  }
}


exahype::solvers::Solver::RefinementControl NavierStokes::NavierStokesSolver_ADERDG::refinementCriterion(
    const double* luh,
    const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx,
    const double t,
    const int level) {
  if (!amrSettings.useAMR) {
    // Default: Delete cells.
    // This is useful when one wants to use limiting-guided refinement
    // without another source of AMR.
    return exahype::solvers::Solver::RefinementControl::Erase;
  }

  if (t == 0) {
    // Global observables are reduced after first timestep!
    return exahype::solvers::Solver::RefinementControl::Keep;
  }

  const auto curObservables = mapGlobalObservables(luh, dx);
  const auto curTv = curObservables[0];

  const auto countGlobal = _globalObservables[2];
  const auto meanGlobal = _globalObservables[0];
  // Merging computes sample variance (Bessel's correction), we need population variance.
  const auto varianceGlobal = ((countGlobal - 1)/countGlobal) * _globalObservables[1];
  const auto stdGlobal = std::sqrt(varianceGlobal);

  const auto factorRefine = amrSettings.factorRefine;
  const auto factorCoarse = amrSettings.factorErase;

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
  return (0.7/0.9) * kernels::aderdg::generic::c::stableTimeStepSize<NavierStokesSolver_ADERDG,true>(*static_cast<NavierStokesSolver_ADERDG*>(this),luh,dx);
}

void NavierStokes::NavierStokesSolver_ADERDG::riemannSolver(double* FL,double* FR,const double* const QL,const double* const QR,const double t, const double dt,const tarch::la::Vector<DIMENSIONS, double>& dx, const int direction, bool isBoundaryFace, int faceIndex) {
  assertion2(direction>=0,dt,direction);
  assertion2(direction<DIMENSIONS,dt,direction);
  kernels::aderdg::generic::c::riemannSolverNonlinear<false,true, NavierStokesSolver_ADERDG>(*static_cast<NavierStokesSolver_ADERDG*>(this),FL,FR,QL,QR,t,dt,dx,direction);
}

void NavierStokes::NavierStokesSolver_ADERDG::boundaryConditions( double* const fluxIn, const double* const stateIn, const double* const gradStateIn, const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& cellCentre, const tarch::la::Vector<DIMENSIONS,double>&  cellSize, const double t,const double dt, const int direction, const int orientation) {
#if DIMENSIONS == 2
  constexpr int basisSize     = (Order+1);
#else
  constexpr int basisSize     = (Order+1) * (Order+1);
#endif
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

    riemannSolver(FL,FR,QL,QR,t,dt,cellSize,direction,true,faceIndex);
  }
  else {
    double* FL =  fluxIn;  const double* const QL = stateIn;
    double* FR = fluxOut; const double* const QR = stateOut;

    riemannSolver(FL,FR,QL,QR,t,dt,cellSize,direction,true,faceIndex);
  }

  if (scenario->getBoundaryType(faceIndex) == NavierStokes::BoundaryType::wall ||
      scenario->getBoundaryType(faceIndex) == NavierStokes::BoundaryType::freeSlipWall ||
      scenario->getBoundaryType(faceIndex) == NavierStokes::BoundaryType::movingWall) {
#if DIMENSIONS == 2
    kernels::idx2 idx_F(Order + 1, NumberOfVariables);
    for (int i = 0; i < (Order + 1); ++i) {
      // Set energy flux to zero!
      fluxIn[idx_F(i, NavierStokesSolver_ADERDG_Variables::shortcuts::E)] = 0.0;
      //ns.setZ(fluxIn + idx_F(i, 0), 0.0);
    }
#else
   // TODO(Lukas) Is this correct for 3D? Untested!
    kernels::idx3 idx_F(Order + 1, Order + 1, NumberOfVariables);
    for (int i = 0; i < (Order + 1); ++i) {
      for (int j = 0; j < (Order + 1); ++j) {
        // Set energy flux to zero!
        fluxIn[idx_F(i, j, NavierStokesSolver_ADERDG_Variables::shortcuts::E)] = 0.0;
      }
    }
#endif
  }

  delete[] block;
}

std::vector<double> NavierStokes::NavierStokesSolver_ADERDG::mapGlobalObservables(const double *const Q,
        const tarch::la::Vector<DIMENSIONS,double>& dx) const {
  return ::NavierStokes::mapGlobalObservables(Q, dx, scenarioName, ns, amrSettings,
          Order, NumberOfVariables, NumberOfParameters, NumberOfGlobalObservables);
}

std::vector<double> NavierStokes::NavierStokesSolver_ADERDG::resetGlobalObservables() const {
    return ::NavierStokes::resetGlobalObservables(NumberOfGlobalObservables);
}

void NavierStokes::NavierStokesSolver_ADERDG::reduceGlobalObservables(
        std::vector<double> &reducedGlobalObservables,
        const std::vector<double> &curGlobalObservables) const {
    ::NavierStokes::reduceGlobalObservables(reducedGlobalObservables,
            curGlobalObservables, NumberOfGlobalObservables);
}
