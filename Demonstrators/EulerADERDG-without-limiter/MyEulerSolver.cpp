#include "MyEulerSolver.h"

#include "MyEulerSolver_Variables.h"

#include "LogoExaHyPE.h"

#include "peano/utils/Loop.h"
#include "kernels/KernelUtils.h"

#include <limits>


tarch::logging::Log EulerADERDG::MyEulerSolver::_log( "EulerADERDG::MyEulerSolver" );


void EulerADERDG::MyEulerSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}



void EulerADERDG::MyEulerSolver::flux(const double* const Q,double** const F) {
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


void EulerADERDG::MyEulerSolver::eigenvalues(const double* const Q,const int d,double* const lambda) {
  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const double GAMMA = 1.4;
  const double irho = 1./vars.rho();
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  double u_n = Q[d + 1] * irho;
  double c  = std::sqrt(GAMMA * p * irho);

  eigs.rho()=u_n - c;
  eigs.E()  =u_n + c;
  eigs.j(u_n,u_n,u_n);
}

double EulerADERDG::MyEulerSolver::getInitialProfile(const double* const x) {
#ifdef Dim2
  tarch::la::Vector<DIMENSIONS,double> myX( x[0] - 0.06, 1.0-x[1] - 0.25 ); // translate
  myX *= static_cast<double>(LogoExaHyPE.width);
  tarch::la::Vector<DIMENSIONS,int>    myIntX( 1.2*myX(0) , 1.2*myX(1) );  // scale
#else
  tarch::la::Vector<DIMENSIONS,double> myX( x[0] - 0.06, 1.0-x[1] - 0.25, 1.0-x[2] - 0.25 ); // translate
  myX *= static_cast<double>(LogoExaHyPE.width);
  tarch::la::Vector<DIMENSIONS,int>    myIntX( 1.2*myX(0) , 1.2*myX(1) , 1.2*myX(2) );  // scale
#endif

  if (
    myIntX(0) > 0 && myIntX(0) < static_cast<int>(LogoExaHyPE.width)
    &&
    myIntX(1) > 0 && myIntX(1) < static_cast<int>(LogoExaHyPE.height)
  ) {
    double value =  LogoExaHyPE.pixel_data[myIntX(1)*LogoExaHyPE.width+myIntX(0)];
    assertion(value>=0.0);
    assertion(value<=1.0);
    return (2.0-value);
  }
  else return 1.0;
}

void EulerADERDG::MyEulerSolver::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
  Variables vars(Q);
  if ( tarch::la::equals( t,0.0 ) ) {
    vars.rho() = 1.0;
    vars.j(0,0,0);
    vars.E() = getInitialProfile(x);
  }
}


void EulerADERDG::MyEulerSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  std::copy_n(stateIn, NumberOfVariables, stateOut);
  stateOut[1+normalNonZero] =  -stateOut[1+normalNonZero];
  double _F[3][NumberOfVariables]={0.0};
  double* F[3] = {_F[0], _F[1], _F[2]};
  flux(stateOut,F);
  std::copy_n(F[normalNonZero], NumberOfVariables, fluxOut);
}


exahype::solvers::Solver::RefinementControl EulerADERDG::MyEulerSolver::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  double largestE   = -std::numeric_limits<double>::max();
  double smallestE  = +std::numeric_limits<double>::max();

  kernels::idx3 idx_luh(Order+1,Order+1,NumberOfVariables);
  dfor(i,Order+1) {
    ReadOnlyVariables vars(luh + idx_luh(i(1),i(0),0));

    largestE  = std::max (largestE,   vars.E());
    smallestE = std::min (smallestE,  vars.E());
  }

  assertion(largestE>=smallestE);

  if (largestE-smallestE > 0.8) {
    return exahype::solvers::Solver::RefinementControl::Refine;
  }
  else if (largestE-smallestE < 0.2 && level > getCoarsestMeshLevel() && t>1e-8) {
    return exahype::solvers::Solver::RefinementControl::Erase;
  }
  else return exahype::solvers::Solver::RefinementControl::Keep;
}
