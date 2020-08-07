#include "LimitingADERDG_FV.h"

#include "LimitingADERDG_FV_Variables.h"
#include "LimitingADERDG_ADERDG.h"

#include "PDE.h"

tarch::logging::Log Euler::LimitingADERDG_FV::_log( "Euler::LimitingADERDG_FV" );

void Euler::LimitingADERDG_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {}

void Euler::LimitingADERDG_FV::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {
  if (tarch::la::equals(t,0.0)) {
      initialdata(x,t,Q);
  }
  for(int i=0; i< 9; i++){
      if(!std::isfinite(Q[i])){
          std::cout << "Q[" << i << "]=" << Q[i] << std::endl;
          throw(2);
      }
  }
}

void Euler::LimitingADERDG_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int normalNonZero,
    const double* const stateIn,
    double* const stateOut) {

  const int nVar = NumberOfVariables;
  std::copy_n(stateIn,nVar,stateOut);
  stateOut[1+normalNonZero] = -stateIn[1+normalNonZero];
    if(faceIndex == 0){ // inflow
        stateOut[1+normalNonZero] = stateIn[1+normalNonZero];
    }
    if(faceIndex == 1) { //outflow
        stateOut[1+normalNonZero] = stateIn[1+normalNonZero];
        stateOut[0] = 1.0;
        stateOut[4] = 2.5;
    }
}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

void Euler::LimitingADERDG_FV::eigenvalues(const double* const Q, const int direction, double* const lambda) {
    PDEEigenvalues(Q,direction, lambda);
}

void Euler::LimitingADERDG_FV::flux(const double* const Q,double** const F) {
    PDEflux(Q,F);
}


void Euler::LimitingADERDG_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ){
   PDEncp(Q,gradQ,BgradQ);
}

void Euler::LimitingADERDG_FV::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S){
    for(int i = 0; i < 9; i++) S[i]=0.0;
}

double Euler::LimitingADERDG_FV::riemannSolver(
        double* fnL, double *fnR, const double* qL, const double* qR,
        const double* gradQL, const double* gradQR,
        const double* const cellSize,
        int normalNonZero) {
  constexpr int numberOfVariables  = Euler::AbstractLimitingADERDG_FV::NumberOfVariables;
  constexpr int numberOfParameters = Euler::AbstractLimitingADERDG_FV::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;

  // Compute penalty contribution of convective eigenvalues.
  double sL[numberOfVariables];
  double sR[numberOfVariables];

  eigenvalues(qL, normalNonZero, sL);
  eigenvalues(qR, normalNonZero, sR);

  double maxEigenvalue = -1.0;
  for (int i = 0; i < numberOfVariables; i++) {
    maxEigenvalue = std::max( std::abs(sL[i]), maxEigenvalue);
  }
  for (int i = 0; i < numberOfVariables; i++) {
    maxEigenvalue = std::max(std::abs(sR[i]), maxEigenvalue);
  }

  double s_max = maxEigenvalue;
  double s_max_dt = maxEigenvalue;

  // determine BgradQ from ncp
  double ncp[numberOfData] = {0.0};
  double gradQ[DIMENSIONS][numberOfData] = {0.0};
    double Qavg[numberOfData];
    for(int k=0; k < numberOfData; k++) {
       Qavg[k] = (qR[k] + qL[k]) / 2;
       gradQ[normalNonZero][k] = qR[k] - qL[k];
    }
    nonConservativeProduct(Qavg, gradQ[0], ncp);

  double FL2[DIMENSIONS][numberOfVariables] = {0.0}; // Q: Can we skip this memset?
  double FR2[DIMENSIONS][numberOfVariables] = {0.0};
  double* FL[DIMENSIONS]={FL2[0], FL2[1]};
  double* FR[DIMENSIONS]={FR2[0], FR2[1]};
    flux(qL, FL);
    flux(qR, FR);

  for (int i = 0; i < numberOfVariables; i++) {
    fnL[i] = 0.5 * s_max * (qL[i] - qR[i]);
      fnL[i] += 0.5 * (FL2[normalNonZero][i] + FR2[normalNonZero][i]);

      fnR[i] = fnL[i] - 0.5 * ncp[i];
      fnL[i] = fnL[i] + 0.5 * ncp[i];
  }
  //fnL[5] = 0.0;
  //fnR[5] = 0.0;

  return s_max_dt;  
}
