#include "MySWESolver.h"
#include "InitialData.h"
#include "MySWESolver_Variables.h"

#include "kernels/KernelUtils.h"

using namespace kernels;

double grav;
int scenario;

tarch::logging::Log SWE::MySWESolver::_log( "SWE::MySWESolver" );

void SWE::MySWESolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  if (constants.isValueValidDouble( "grav" )) {
    grav = constants.getValueAsDouble("grav");
  }
  if (constants.isValueValidInt( "scenario" )) {
    scenario = constants.getValueAsInt( "scenario" );
  }
}

void SWE::MySWESolver::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {
  // Dimensions             = 2
  // Number of variables    = 4 + #parameters

  if (tarch::la::equals(t,0.0)) {
    initialData(x, Q);
  }
}

void SWE::MySWESolver::eigenvalues(const double* const Q, const int dIndex, double* const lambda) {
  // Dimensions             = 2
  // Number of variables    = 4 + #parameters

  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const double c = std::sqrt(grav*vars.h());
  const double ih = 1./vars.h();
  double u_n = Q[dIndex + 1] * ih;

  eigs.h() = u_n + c;
  eigs.hu() = u_n - c;
  eigs.hv() = u_n;
  eigs.b() = 0.0;
}

void SWE::MySWESolver::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* const stateOutside) {
  // Dimensions             = 2
  // Number of variables    = 4 + #parameters

  //Outflow
  stateOutside[0] = stateInside[0];
  stateOutside[1] = stateInside[1];
  stateOutside[2] = stateInside[2];
  stateOutside[3] = stateInside[3];

  //Wall
  stateOutside[d+1]=-stateInside[d+1];


}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit


void SWE::MySWESolver::flux(const double* const Q,double** const F) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 4 + 0

  ReadOnlyVariables vars(Q);

  const double ih = 1./vars.h();

  double* f = F[0];
  double* g = F[1];

  f[0] = vars.hu();
  f[1] = vars.hu()*vars.hu()*ih + 0.5*grav*vars.h()*vars.h();
  f[2] = vars.hu()*vars.hv()*ih;
  f[3] = 0.0;

  g[0] = vars.hv();
  g[1] = vars.hu()*vars.hv()*ih;
  g[2] = vars.hv()*vars.hv()*ih + 0.5*grav*vars.h()*vars.h();
  g[3] = 0.0;
  
}



void  SWE::MySWESolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  idx2 idx_gradQ(DIMENSIONS,NumberOfVariables);

  BgradQ[0] = 0.0;
  BgradQ[1] = grav*Q[0]*gradQ[idx_gradQ(0,3)];
  BgradQ[2] = grav*Q[0]*gradQ[idx_gradQ(1,3)];
  BgradQ[3] = 0.0;
}

