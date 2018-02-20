#include "MySWESolver.h"
#include "InitialData.h"
#include "MySWESolver_Variables.h"

#include "kernels/KernelUtils.h"

using namespace kernels;

const double grav= 9.81;


//tarch::logging::Log SWE::MySWESolver::_log( "SWE::MySWESolver" );



void SWE::MySWESolver::init(const std::vector<std::string>& cmdlineargs) {
  //logInfo( "init(...)", "SWE is called with these parameters:" );
  for(size_t i=0; i<cmdlineargs.size(); i++) {
  //  logInfo( "init(...)", "- argument " << i << ": " << cmdlineargs[i] );
  }
}


bool SWE::MySWESolver::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
  return tarch::la::equals(t,0.0);
}


void SWE::MySWESolver::adjustSolution(const double* const x,const double w,const double t,const double dt,double* Q) {
  assertion(tarch::la::equals(t, 0.0));

  initialData(x,Q);
}


void SWE::MySWESolver::eigenvalues(const double* const Q,const int normalNonZeroIndex,double* lambda) {
  // Dimensions             = 2
  // Number of variables    = 3 (#unknowns + #parameters)
  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);  


  const double c= std::sqrt(grav*vars.h());
  const double ih = 1./vars.h();

  double u_n = Q[normalNonZeroIndex + 1] * ih;

  eigs.h() = u_n + c ;
  eigs.hu()= u_n -c;
  eigs.hv()= u_n ;
}


void SWE::MySWESolver::flux(const double* const Q,double** F) {
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


void SWE::MySWESolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const stateIn,double* stateOut)
 {
  // Dimensions             = 2
  // Number of variables    = 3 + 1

  //for OUTFLOW and WALL
  stateOut[0] = stateIn[0];
  stateOut[1] = stateIn[1];
  stateOut[2] = stateIn[2];
  stateOut[3] = stateIn[3];

  //for WALL BCs
  stateOut[normalNonZero+1]=-stateIn[normalNonZero+1];
}

  return true;
}

void SWE::MySWESolver::coefficientMatrix(const double* const Q, const int d, double* Bn)
{
  idx2 idx_Bn(NumberOfVariables+NumberOfParameters,NumberOfVariables+NumberOfParameters);

  Bn[0] = 0.0;
  Bn[1] = 0.0;
  Bn[2] = 0.0;
  Bn[3] = 0.0;
  Bn[4] = 0.0;
  Bn[5] = 0.0;
  Bn[6] = 0.0;
  Bn[7] = 0.0;
  Bn[8] = 0.0;
  Bn[9] = 0.0;
  Bn[10]= 0.0;
  Bn[11]= 0.0;
  Bn[12]= 0.0;
  Bn[13]= 0.0;
  Bn[14]= 0.0;
  Bn[15]= 0.0;    

  Bn[idx_Bn(3,d+1)]=grav*Q[0]; //g*h
}

exahype::solvers::Solver::RefinementControl SWE::MySWESolver::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}
