#include "LimitingADERDG_ADERDG.h"

#include "InitialData.h"
#include "LimitingADERDG_ADERDG_Variables.h"


tarch::logging::Log Euler::LimitingADERDG_ADERDG::_log( "Euler::LimitingADERDG_ADERDG" );


void Euler::LimitingADERDG_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

void Euler::LimitingADERDG_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
  if (tarch::la::equals(t,0.0)) {
    Euler::initialData(x,Q);
  }
}

void Euler::LimitingADERDG_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut)
  fluxOut[0] = fluxIn[0];
  fluxOut[1] = fluxIn[1];
  fluxOut[2] = fluxIn[2];
  fluxOut[3] = fluxIn[3];
  fluxOut[4] = fluxIn[4];

  stateOut[0] = stateIn[0];
  stateOut[1] = stateIn[1];
  stateOut[2] = stateIn[2];
  stateOut[3] = stateIn[3];
  stateOut[4] = stateIn[4];
}

exahype::solvers::Solver::RefinementControl Euler::LimitingADERDG_ADERDG::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  if (level>getCoarsestMeshLevel()) {
    return exahype::solvers::Solver::RefinementControl::Erase;
  }

  return exahype::solvers::Solver::RefinementControl::Keep;
}


void Euler::LimitingADERDG_ADERDG::mapDiscreteMaximumPrincipleObservables(double* observables, const double* const Q) const {
  assertion(numberOfObservables==2);
  ReadOnlyVariables vars(Q);

  observables[0]=vars.rho(); //extract density

  const double GAMMA = 1.4;
  const double irho = 1./vars.rho();
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );
  observables[1]=p; //extract pressure
}


bool Euler::LimitingADERDG_ADERDG::isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,
      const bool wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t) const {
  
//  if ((center[0]-0.5)*(center[0]-0.5)+(center[1]-0.5)*(center[1]-0.5)<0.25*dx[0]*dx[0]) return false;

  if (observablesMin[0] <= 0.0) return false;
  if (observablesMin[1] < 0.0) return false;
  return true;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void Euler::LimitingADERDG_ADERDG::eigenvalues(const double* const Q,const int d,double* lambda) {
  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  constexpr double GAMMA = 1.4;
  const     double irho  = 1./vars.rho();
  const     double p     = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  double u_n = vars.j(d) * irho;
  double c   = std::sqrt(GAMMA * p * irho);

  eigs.rho()=u_n - c;
  eigs.E()  =u_n + c;
  eigs.j(u_n,u_n,u_n);
}


void Euler::LimitingADERDG_ADERDG::flux(const double* const Q,double** F) {
  ReadOnlyVariables vars(Q);
  Fluxes f(F);

  tarch::la::Matrix<3,3,double> I;
  I = 1, 0, 0,
      0, 1, 0,
      0, 0, 1;

  const double GAMMA = 1.4;
  const double irho = 1./vars.rho();
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  f.rho ( vars.j()                                 );
  f.j   ( irho * outerDot(vars.j(),vars.j()) + p*I );
  f.E   ( irho * (vars.E() + p) * vars.j()         );
  
}





