#include "MyEulerSolver.h"
#include "MyEulerSolver_Variables.h"

#include "picopng.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>


tarch::logging::Log Euler::MyEulerSolver::_log( "Euler::MyEulerSolver" );

picopng::Image Image;

void Euler::MyEulerSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::Parser::ParserView& constants) {
	// Reading in Initial Data at startup
	// You might pass your own PNG file on the command line next to the specfile.
	const char* imagename = cmdlineargs.size() > 1 ? cmdlineargs.back().c_str() : "logo.png";
	int ret = Image.loadPNG(imagename);
	if(ret) {
		fprintf(stderr, "Please ensure image '%s' is readable and a proper PNG file (err %d)\n", imagename, ret);
		exit(-1);
	}
}

bool Euler::MyEulerSolver::useAdjustSolution(const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, const double t, const double dt) const {
  return tarch::la::equals(t,0.0);
}

void Euler::MyEulerSolver::adjustSolution(const double* const x,const double w,const double t,const double dt, double* Q) {
  Variables vars(Q);
  
  tarch::la::Vector<DIMENSIONS,double> myX( x[0], 1.0-x[1] );
  myX *= static_cast<double>(Image.width);
  tarch::la::Vector<DIMENSIONS,int>    myIntX( myX(0), myX(1) );

  double Energy = 0.1;

  if (
    myIntX(0) < static_cast<int>(Image.width)
    &&
    myIntX(1) < static_cast<int>(Image.height)
  ) {
    Energy += (
        Image.pixel_data[myIntX(1)*Image.width*4+myIntX(0)*4+0]
      + Image.pixel_data[myIntX(1)*Image.width*4+myIntX(0)*4+1]
      + Image.pixel_data[myIntX(1)*Image.width*4+myIntX(0)*4+2]) / 3.0 / 256.0;
  }
  else {
    Energy += (
        Image.pixel_data[0]
      + Image.pixel_data[1]
      + Image.pixel_data[2]) / 3.0 / 256.0;
  }

  vars.rho() = 1.0;
  vars.E()   = Energy;
  vars.j(0,0,0);
}

exahype::solvers::Solver::RefinementControl Euler::MyEulerSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) {
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void Euler::MyEulerSolver::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const double GAMMA = 1.4;
  const double irho = 1./vars.rho();
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  double u_n = Q[dIndex + 1] * irho;
  double c  = std::sqrt(GAMMA * p * irho);

  eigs.rho()=u_n - c;
  eigs.E()  =u_n + c;
  eigs.j(u_n,u_n,u_n);
}

void Euler::MyEulerSolver::flux(const double* const Q, double** F) {
  ReadOnlyVariables vars(Q);
  Fluxes fluxes(F);

  tarch::la::Matrix<3,3,double> I;
  I = 1, 0, 0,
      0, 1, 0,
      0, 0, 1;

  const double GAMMA = 1.4;
  const double irho = 1./vars.rho();
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  fluxes.rho ( vars.j()                                 );
  fluxes.j   ( irho * outerDot(vars.j(),vars.j()) + p*I );
  fluxes.E   ( irho * (vars.E() + p) * vars.j()         );
}



void Euler::MyEulerSolver::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* stateOutside) {
  ReadOnlyVariables varsInside(stateInside);
  Variables         varsOutside(stateOutside);

  varsOutside = varsInside;
}
