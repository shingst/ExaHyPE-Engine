#include "PotentialEulerSolver.h"
#include "PotentialEulerSolver_Variables.h"
#include "Primitives.h"
#include "tarch/la/Vector.h"

tarch::logging::Log PotentialHydro::PotentialEulerSolver::_log( "PotentialHydro::PotentialEulerSolver" );

constexpr int numberOfVariables = PotentialHydro::PotentialEulerSolver::NumberOfVariables;

void PotentialHydro::PotentialEulerSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
}

bool PotentialHydro::PotentialEulerSolver::useAdjustSolution(const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, const double t, const double dt) const {
  return tarch::la::equals(t,0.0);
}

void PotentialHydro::PotentialEulerSolver::adjustSolution(const double* const x,const double w,const double t,const double dt, double* Q) {
  // initial data: stable matter clumb  
  using namespace tarch::la;
  typedef Vector<DIMENSIONS,double> vec;

  vec xvec(x[0], x[1]);
  vec v0(0.0, 0);
  vec x0(0.0, 0.0);

  double width = 0.20;
  double V[numberOfVariables];
  V[0] = 0.2 + 0.4 * exp(-norm2(xvec - x0 - v0 * t) /  pow(width, DIMENSIONS));
  V[1] = v0[0];
  V[2] = v0[1];
  V[3] = 0.; // vz
  V[4] = 0.; // energy
  V[5] = 0.; // Potential

  prim2con(Q, V);
}

exahype::solvers::Solver::RefinementControl PotentialHydro::PotentialEulerSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void PotentialHydro::PotentialEulerSolver::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
  /*** Ordinary EulerSolver eigenvalues ***/
  // For the time being, set all eigenvalues to one. Is pessmistic but okay.
  for(int i=0; i<numberOfVariables; i++)
    lambda[i] = 1.0;
}

void PotentialHydro::PotentialEulerSolver::flux(const double* const Q, double** F) {
  /* Ordinary EulerSolver flux */
  ReadOnlyVariables vars(Q);
  Fluxes f(F);

  tarch::la::Matrix<3,3,double> I;
  I = 1, 0, 0,
      0, 1, 0,
      0, 0, 1;

  const double irho = 1./vars.rho();
  const double p = (eos_gamma-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  f.rho ( vars.j()                                 );
  f.j   ( irho * outerDot(vars.j(),vars.j()) + p*I );
  f.E   ( irho * (vars.E() + p) * vars.j()         );
  
  // the potential has no flux
  f.Pot (0, 0);
}



void PotentialHydro::PotentialEulerSolver::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* stateOutside) {
  // Dimensions             = 2
  // Number of variables    = 5 + #parameters

  // Copy Bc for FV solver
  for(int i=0; i<numberOfVariables; i++)
	stateOutside[i] = stateInside[i];
}


//void PotentialHydro::algebraicSource(const double* const Q, double* S) {}

void PotentialHydro::PotentialEulerSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  // for testing, first set it to zero
  for(int i=0; i<numberOfVariables; i++)
    BgradQ[i] = 0;
}
