#include "MySWESolver.h"
#include "../InitialData.h"
#include "MySWESolver_Variables.h"

const double grav = 9.81;


tarch::logging::Log SWE::MySWESolver::_log( "SWE::MySWESolver" );

void SWE::MySWESolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

void SWE::MySWESolver::adjustSolution(const double* const x,const double t,const double dt, double* Q) {
    // Dimensions             = 2
    // Number of variables    = 4 + #parameters

    if (tarch::la::equals(t, 0.0)) {
//        MySWESolver::Variables vars(Q);
//
//        if (x[0] <= 5){
//            vars.h() = x[0];
//        }
//        else {
//            vars.h() = 10 - x[0];
//        }
//        vars.hu() = 0.0;
//        vars.hv() = 0.0;
//        vars.b() = 0.0;

        MySWESolver::Variables vars(Q);

        if((x[0] -5) *(x[0] -5) + (x[1] -5) *(x[1] -5) < 2) {
            vars.h() = 4.0;
            vars.hu()= 0.0;
            vars.hv()= 0.0;
            vars.b() = 0;
        } else {
            vars.h() = 3.0;
            vars.hu()= 0.0;
            vars.hv()= 0.0;
            vars.b() = 0.0;
        }
    }
}

void SWE::MySWESolver::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
  // Dimensions             = 2
  // Number of variables    = 4 + #parameters

  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const double c= std::sqrt(grav*vars.h());
  const double ih = 1./vars.h();
  double u_n = Q[dIndex + 1] *ih;

  eigs.h() = u_n + c ;
  eigs.hu()= u_n -c;
  eigs.hv()= u_n ;
}

void SWE::MySWESolver::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* stateOutside) {
  // Dimensions             = 2
  // Number of variables    = 4 + #parameters

  stateOutside[0] = stateInside[0];
  stateOutside[1] = stateInside[1];
  stateOutside[2] = stateInside[2];
  stateOutside[3] = stateInside[3];

    //for WALL BCs
    stateOutside[d+1]=-stateInside[d+1];
}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit


void SWE::MySWESolver::flux(const double* const Q,double** F) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 3 + 1

  ReadOnlyVariables vars(Q);

  const double ih = 1./vars.h();

  double* f = F[0];
  double* g = F[1];

  f[0] = vars.hu();
  f[1] = vars.hu()*vars.hu()*ih + 0.5*grav*vars.h()*vars.h();
  f[2] = vars.hu()*vars.hv()*ih;

  g[0] = vars.hv();
  g[1] = vars.hu()*vars.hv()*ih;
  g[2] = vars.hv()*vars.hv()*ih + 0.5*grav*vars.h()*vars.h();
}




