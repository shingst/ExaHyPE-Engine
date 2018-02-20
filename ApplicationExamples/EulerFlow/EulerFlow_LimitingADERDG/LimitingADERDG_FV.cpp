#include "LimitingADERDG_FV.h"

#include "InitialData.h"
#include "LimitingADERDG_FV_Variables.h"


tarch::logging::Log Euler::LimitingADERDG_FV::_log( "Euler::LimitingADERDG_FV" );


void Euler::LimitingADERDG_FV::init(const std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

void Euler::LimitingADERDG_FV::adjustSolution(const double* const x,const double t,const double dt, double* Q) {
  initialData(x,Q);
}

void Euler::LimitingADERDG_FV::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const double GAMMA = 1.4;
  const double irho = 1./vars.rho();
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  double u_n = vars.j(dIndex) * irho;
  double c   = std::sqrt(GAMMA * p * irho);

  eigs.rho()=u_n - c;
  eigs.E()  =u_n + c;
  eigs.j(u_n,u_n,u_n);
}

void Euler::LimitingADERDG_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* stateOutside) {
  for (int i=0; i<NumberOfVariables+NumberOfParameters; i++) {
    stateOutside[i] = stateInside[i];
  }
}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit


void Euler::LimitingADERDG_FV::flux(const double* const Q,double** F) {
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




