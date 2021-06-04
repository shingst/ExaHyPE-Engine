// This file was generated by the ExaHyPE toolkit.
// It will NOT be regenerated or overwritten.
// Please adapt it to your own needs.
// 
// ========================
//   www.exahype.eu
// ========================

#include "MySWESolver_ADERDG.h"
#include "InitialData.h"
#include "MySWESolver_ADERDG_Variables.h"

#include "peano/utils/Loop.h"
#include "kernels/KernelUtils.h"
#include "../../../ExaHyPE/kernels/GaussLegendreBasis.h"

using namespace kernels;

double grav_DG;
double epsilon_DG;
int scenario_DG;

tarch::logging::Log SWE::MySWESolver_ADERDG::_log( "SWE::MySWESolver_ADERDG" );

#if defined USE_TMPI && defined(GenerateError)
#include "teaMPI.h"
#endif

void SWE::MySWESolver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  if (constants.isValueValidDouble( "grav" )) {
    grav_DG = constants.getValueAsDouble("grav")*1.0e-3;
  }
  if (constants.isValueValidDouble( "epsilon" )) {
    epsilon_DG = constants.getValueAsDouble( "epsilon" )*1.0e-3;
  }
  if (constants.isValueValidInt( "scenario" )) {
    initialData= new InitialData(constants.getValueAsInt( "scenario" ));
  }

}

void SWE::MySWESolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 4 + 0
  if (tarch::la::equals(t,0.0)) {
    initialData->getInitialData(x, Q); 
  }

#if defined(GenerateError)
  if (tarch::la::equals(t,0.927,0.01)
     && tarch::la::equals(x[0],2.6319,0.01) 
     && tarch::la::equals(x[1],5.428,0.01)
     //2.6319 x[1]=5.42863
#if defined USE_TMPI
    && TMPI_IsLeadingRank()
#endif	     
    ) {
    logError("adjustPointSolution", "Introducing error into solution x[0]="<<x[0]<<" x[1]="<<x[1]<<" t= "<<t);
    Q[0]=Q[0]*10;
  }
#endif
}

void SWE::MySWESolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 4 + 0

  //Outflow (not working)
//    for(int i=0; i < NumberOfVariables; i++) {
//        fluxOut[i] = fluxIn[i];
//        stateOut[i] = stateIn[i];
//    }

  //Wall
  std::copy_n(stateIn, NumberOfVariables, stateOut);
  stateOut[1+normalNonZero] =  -stateOut[1+normalNonZero];
  double _F[2][NumberOfVariables]={0.0};
  double* F[2] = {_F[0], _F[1]};
  flux(stateOut,F);
  std::copy_n(F[normalNonZero], NumberOfVariables, fluxOut);
}

exahype::solvers::Solver::RefinementControl SWE::MySWESolver_ADERDG::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
    double largestH = -std::numeric_limits<double>::max();
    double smallestH = std::numeric_limits<double>::max();

    kernels::idx3 idx_luh(Order+1,Order+1,NumberOfVariables);
    dfor(i,Order+1) {
        ReadOnlyVariables vars(luh + idx_luh(i(1),i(0),0));
        largestH = std::max (largestH, vars.h());
        smallestH = std::min(smallestH, vars.h());
    }

    //gradient
//    if (largestH - smallestH > 5e-2){
//        return exahype::solvers::Solver::RefinementControl::Refine;
//    }

    //height
//    if (smallestH < 3.5 && level > getCoarsestMeshLevel() + 1) {
//        return exahype::solvers::Solver::RefinementControl::Refine;
//    }
//    if (smallestH < 3.7 && level > getCoarsestMeshLevel()) {
//        return exahype::solvers::Solver::RefinementControl::Refine;
//    }
//
//    if (smallestH < 3.9 && level == getCoarsestMeshLevel()) {
//        return exahype::solvers::Solver::RefinementControl::Refine;
//    }
//
//    if (level > getCoarsestMeshLevel())
//        return exahype::solvers::Solver::RefinementControl::Erase;
    return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void SWE::MySWESolver_ADERDG::eigenvalues(const double* const Q,const int d,double* const lambda) {
  /// Dimensions                        = 2
  // Number of variables + parameters  = 4 + 0
  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const double c = std::sqrt(grav_DG*vars.h());
  const double ih = 1./vars.h();
  double u_n = Q[d + 1] * ih;


  if (vars.h() < epsilon_DG){
    eigs.h() = epsilon_DG;
    eigs.hu() = 0.0;
    eigs.hv() = 0.0;
    eigs.b() = 0.0;
    //    //std::cout << 0.0 << std::endl;
  }
  else {
    eigs.h() = u_n + c;
    eigs.hu() = u_n - c;
    eigs.hv() = u_n;
    eigs.b() = 0.0;
    //    //std::cout << eigs.h() + std::abs(c) << std::endl;
  }
}


void SWE::MySWESolver_ADERDG::flux(const double* const Q,double** const F) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 4 + 0

  ReadOnlyVariables vars(Q);

  const double ih = 1./vars.h();

  double* f = F[0];
  double* g = F[1];

  if (Q[0] < epsilon_DG){
    f[0] = 0.0;
    f[1] = 0.0;
    f[2] = 0.0;
    f[3] = 0.0;

    g[0] = 0.0;
    g[1] = 0.0;
    g[2] = 0.0;
    g[3] = 0.0;
  }
  else {
    f[0] = vars.hu();
    // Moved hydrostatic pressure to ncp for well balancedness
    //f[1] = vars.hu() * vars.hu() * ih + 0.5 * grav_DG * vars.h() * vars.h();
    f[1] = vars.hu() * vars.hu() * ih;
    f[2] = vars.hu() * vars.hv() * ih;
    f[3] = 0.0;

    g[0] = vars.hv();
    g[1] = vars.hu() * vars.hv() * ih;
    // Moved hydrostatic pressure to ncp for well balancedness
    //g[2] = vars.hv() * vars.hv() * ih + 0.5 * grav_DG * vars.h() * vars.h();
    g[2] = vars.hv() * vars.hv() * ih;
    g[3] = 0.0;
  }
}



void  SWE::MySWESolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  idx2 idx_gradQ(DIMENSIONS,NumberOfVariables);

  BgradQ[0] = 0.0;
  BgradQ[1] = grav_DG*Q[0]*gradQ[idx_gradQ(0,3)] + grav_DG*Q[0]*gradQ[idx_gradQ(0,0)]; 
  BgradQ[2] = grav_DG*Q[0]*gradQ[idx_gradQ(1,3)] + grav_DG*Q[0]*gradQ[idx_gradQ(1,0)];
  BgradQ[3] = 0.0;
}


bool SWE::MySWESolver_ADERDG::isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,
      const bool wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t) const {

  //limiter along coast
  double bMin;
  double hMin;
  idx3 id(Order+1,Order+1,NumberOfVariables);
  bMin=std::abs(solution[id(0,0,3)]);
  hMin=         solution[id(0,0,0)];
 
  for(int i = 0 ; i < Order+1 ; i++){
    for(int j = 0 ; j < Order+1 ; j++){
      bMin=std::min(bMin, std::abs(solution[id(i,j,3)]));
      hMin=std::min(hMin, solution[id(i,j,0)]);
    }
  }

/*  if(bMin < 0.100){
    //    //std::cout <<center[0] << "," << center[1] <<": false" << std::endl;
    return false;
  }*/

  if(hMin < epsilon_DG * 10.0){
    return false;
  }
  
  /*  for(int i = 0 ; i < Order+1 ; i++){
    for(int j = 0 ; j < Order+1 ; j++){
      if(solution[id(i,j,0)] < epsilon_DG * 10.0) {
	return false;
      }
    }
    }*/
  //  //std::cout <<center[0] << "," << center[1] <<": true" << std::endl;

  return true;

}

void SWE::MySWESolver_ADERDG::riemannSolver(double* const FL,double* const FR,const double* const QL,const double* const QR,const double* gradQL, const double* gradQR, const double dt,const int direction,bool isBoundaryFace, int faceIndex) {
  constexpr int numberOfVariables  = NumberOfVariables;
  constexpr int numberOfData       = numberOfVariables;
  constexpr int order              = Order;
  constexpr int basisSize          = order+1;
  bool useNCP = true;

  // Compute the average variables and parameters from the left and the right
  double QavL[numberOfData] = {0.0}; // ~(numberOfVariables+numberOfParameters)
  double QavR[numberOfData] = {0.0}; // ~(numberOfVariables+numberOfParameters)
  {
    idx2 idx_QLR(basisSize, numberOfData);
    for (int j = 0; j < basisSize; j++) {
      const double weight = kernels::legendre::weights[order][j];

      for (int k = 0; k < numberOfData; k++) {
        QavL[k] += weight * QL[idx_QLR(j, k)];
        QavR[k] += weight * QR[idx_QLR(j, k)];
      }
    }
  }

  double LL[numberOfVariables] = {0.0}; // do not need to store material parameters
  double LR[numberOfVariables] = {0.0};
  eigenvalues(QavL, direction, LL);
  eigenvalues(QavR, direction, LR);

  //std::cout<<"compute riemann "<<std::endl;

  // skip parameters
  std::transform(LL, LL + numberOfVariables, LL, std::abs<double>);
  std::transform(LR, LR + numberOfVariables, LR, std::abs<double>);
  const double* smax_L = std::max_element(LL, LL + numberOfVariables);
  const double* smax_R = std::max_element(LR, LR + numberOfVariables);
  const double smax    = std::max(*smax_L, *smax_R);

  // compute fluxes (and fluctuations for non-conservative PDEs)
  double Qavg[numberOfData];
  idx2 idx_gradQ(DIMENSIONS, numberOfVariables);
  double gradQ[DIMENSIONS][numberOfVariables] = {0.0};
  double ncp[numberOfVariables]               = {0.0};
    idx2 idx_FLR(basisSize, numberOfVariables);
    idx2 idx_QLR(basisSize, numberOfData);

    for (int j = 0; j < basisSize; j++) {

      if(useNCP) { // we don't use matrixB but the NCP call here.
        for(int l=0; l < numberOfVariables; l++) {
          gradQ[direction][l] = QR[idx_QLR(j, l)] - QL[idx_QLR(j, l)];
          Qavg[l] = 0.5 * (QR[idx_QLR(j, l)] + QL[idx_QLR(j, l)]);
        }

        nonConservativeProduct(Qavg, gradQ[0], ncp);
      }

      // skip parameters
      for (int k = 0; k < numberOfVariables; k++) {

        FL[idx_FLR(j, k)] =
            0.5 * (FR[idx_FLR(j, k)] + FL[idx_FLR(j, k)]) -
            0.5 * smax * (QR[idx_QLR(j, k)] - QL[idx_QLR(j, k)]);

	//to consider the bathmetry
	if(k == 0){
	  FL[idx_FLR(j, k)] -= 0.5 * smax * (QR[idx_QLR(j, 3)] - QL[idx_QLR(j, 3)]);
	}

        if(useNCP) {
          FR[idx_FLR(j, k)] = FL[idx_FLR(j, k)] - 0.5 * ncp[k];
          FL[idx_FLR(j, k)] = FL[idx_FLR(j, k)] + 0.5 * ncp[k];
        } else {
          FR[idx_FLR(j, k)] = FL[idx_FLR(j, k)];
        }

	if(k == 3){
	  //bathymetry doesn't change
          FR[idx_FLR(j, k)] = 0.0;
          FL[idx_FLR(j, k)] = 0.0;
	}
      }
    }

//   std::cout << "QL" << std::endl;
//   for (int j = 0; j < basisSize; j++) {
//       for (int k = 0; k < numberOfVariables; k++) {
// 	std::cout <<QL[idx_FLR(j,k)] << " ";
//       } std::cout << std::endl;
//   }
//   std::cout << "QR" << std::endl;
//   for (int j = 0; j < basisSize; j++) {
//       for (int k = 0; k < numberOfVariables; k++) {
// 	std::cout <<QR[idx_FLR(j,k)] << " ";
//       } std::cout << std::endl;
//   }

//   std::cout << "FL" << std::endl;
//   for (int j = 0; j < basisSize; j++) {
//       for (int k = 0; k < numberOfVariables; k++) {
// 	std::cout <<FL[idx_FLR(j,k)] << " ";
//       } std::cout << std::endl;
//   }

//   std::cout << "FR" << std::endl;
//   for (int j = 0; j < basisSize; j++) {
//       for (int k = 0; k < numberOfVariables; k++) {
// 	std::cout <<FR[idx_FLR(j,k)] << " ";
//       } std::cout << std::endl;
//   }



}

// void SWE::MySWESolver_ADERDG::mapDiscreteMaximumPrincipleObservables(double* const observables, const double* const Q){

//}


