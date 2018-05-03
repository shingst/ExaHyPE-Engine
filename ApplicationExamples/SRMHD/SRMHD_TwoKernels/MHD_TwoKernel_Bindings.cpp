/**
 * This is the fancy FORTRAN and C calling kernel generated stuff.
 *
 **/
#include "MHDSolver.h"

// includes fortran and c kernels
#include "kernels/aderdg/generic/Kernels.h"
#include "kernels/KernelUtils.h"

// std::abs
#include <climits>
#include <cstdlib>
#include <stdio.h>


struct check {
	FILE* out;
	int counter;
	int outputAll;
	const char* name;
	
	check(const char* name, const char* filename) {
		counter = 0;
		outputAll = 100;
		this->name = name;
		out = fopen(filename, "w");
		if(out == NULL) {
			fprintf(stderr, "Checker: Could not open output file '%s'\n", filename);//, strerror(errno));
			exit(-1);
		}
		printf("Check on %s\n", name);
	}
	
	void count() {
		counter++;
	}

	void operator()(
	   const char* d1Name, const double* const d1,
	   const char* d2Name, const double* const d2,
	   kernels::index idx, double accepterror,
	   bool printPointErrors=true,
	   bool exitOnError=false
	) {
		int errors = 0;
		double minError = 1234, maxError = 0;
		for(int n=0; n < idx.size; n++) {
			double error = std::abs(d1[n] - d2[n]);
			if(error > accepterror) {
				if(printPointErrors) {
					fprintf(out, "%d. idx%s %s[%i] = %+.20e %s[%i] = %+.20e error=%+e\n",
						counter, idx.revStr(n).c_str(),
						d1Name ? d1Name : "d1", n, d1[n],
						d2Name ? d2Name : "d2", n, d2[n],
						error);
				}
				errors++;
			}
			// log for statistics
			minError = std::min(minError, error);
			maxError = std::max(maxError, error);
		}

		if(errors > 0 && counter % outputAll == 0) {
			fprintf(out, "%10d. %s: %d/%d errors > %e (minError=%e, maxError=%e)\n", counter, name, errors, idx.size, accepterror, minError, maxError);
			if(exitOnError) exit(-1);
		}
	}
};


MHDSolver::MHDSolver::MHDSolver(double maximumMeshSize,exahype::solvers::Solver::TimeStepping timeStepping,std::vector<std::string>& cmdlineargs, exahype::parser::ParserView constants):
  exahype::solvers::ADERDGSolver("MHDSolver",nVar /* numberOfUnknowns */,nParams /* numberOfParameters */,order + 1  /* nodesPerCoordinateAxis */,maximumMeshSize,timeStepping) {
  init(cmdlineargs, constants);
}

check spaceTimePredictorCheck("spaceTimePredictor", "checks/spaceTimePredictor.txt");
void MHDSolver::MHDSolver::spaceTimePredictor(double* lQhbnd,double* lFhbnd,double** tempSpaceTimeUnknowns,double** tempSpaceTimeFluxUnknowns,double* tempUnknowns,double* tempFluxUnknowns,double* tempStateSizedVectors,const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& dx,const double dt) {
  kernels::darray fQ(2 * DIMENSIONS, basisSize, basisSize, nVar);
  kernels::darray fF(2 * DIMENSIONS, basisSize, basisSize, nVar);
  std::memcpy(fQ.data, lQhbnd, fQ.idx.size * sizeof(double));
  std::memcpy(fF.data, lFhbnd, fF.idx.size * sizeof(double));

  kernels::aderdg::generic::fortran::spaceTimePredictorNonlinear<MHDSolver>(*this,fQ.data,fF.data,tempSpaceTimeUnknowns,tempSpaceTimeFluxUnknowns,tempUnknowns,tempFluxUnknowns,tempStateSizedVectors,luh,dx,dt);
  kernels::aderdg::generic::c::spaceTimePredictorNonlinear<MHDSolver>(*this,lQhbnd,lFhbnd,tempSpaceTimeUnknowns,tempSpaceTimeFluxUnknowns,tempUnknowns,tempFluxUnknowns,tempStateSizedVectors,luh,dx,dt);

  spaceTimePredictorCheck.count();
  spaceTimePredictorCheck("lQhFORT", fQ.data, "lQhCPP", lQhbnd, fQ.idx, 1e-10, true);
  spaceTimePredictorCheck("lFhFORT", fF.data, "lFCPP", lFhbnd, fF.idx, 1e-10, true);
}


check solutionUpdateCheck("solutionUpdate", "checks/solutionUpdate.txt");
void MHDSolver::MHDSolver::solutionUpdate(double* luh,const double* const lduh,const double dt) {
  kernels::darray fort_luh(basisSize, basisSize, basisSize, nVar);
  std::memcpy(fort_luh.data, lduh, fort_luh.idx.size * sizeof(double));

  kernels::aderdg::generic::fortran::solutionUpdate(fort_luh.data,lduh,dt,getNumberOfVariables(),getNumberOfParameters(),getNodesPerCoordinateAxis());
  kernels::aderdg::generic::c::solutionUpdate(luh,lduh,dt,getNumberOfVariables(),getNumberOfParameters(),getNodesPerCoordinateAxis());

  solutionUpdateCheck.count();
  solutionUpdateCheck("luhFORT", fort_luh.data, "luhCPP", lduh, fort_luh.idx, 1e-10, true);
}


check volumeIntegralCheck("volumeIntegral", "checks/volumeIntegral.txt");
void MHDSolver::MHDSolver::volumeIntegral(double* lduh,const double* const lFhi,const tarch::la::Vector<DIMENSIONS,double>& dx) {
  kernels::darray fort_lduh(basisSize, basisSize, basisSize, nVar);
  std::memcpy(fort_lduh.data, lduh, fort_lduh.idx.size * sizeof(double));
  
  kernels::aderdg::generic::fortran::volumeIntegralNonlinear(fort_lduh.data,lFhi,dx,getNumberOfVariables(),getNumberOfParameters(),getNodesPerCoordinateAxis());
  kernels::aderdg::generic::c::volumeIntegralNonlinear(lduh,lFhi,dx,getNumberOfVariables(),getNumberOfParameters(),getNodesPerCoordinateAxis());
  
  volumeIntegralCheck.count();
  volumeIntegralCheck("lduhFORT", fort_lduh.data, "lduhCPP", lduh, fort_lduh.idx, 1e-10, true);
}

void MHDSolver::MHDSolver::surfaceIntegral(double* lduh,const double* const lFhbnd,const tarch::la::Vector<DIMENSIONS,double>& dx) {
  kernels::aderdg::generic::c::surfaceIntegralNonlinear(lduh,lFhbnd,dx,getNumberOfVariables(),getNodesPerCoordinateAxis());
}

check riemannSolverCheck("riemannSolver", "checks/riemannSolver.txt");
void MHDSolver::MHDSolver::riemannSolver(double* FL,double* FR,const double* const QL,const double* const QR,double* tempFaceUnknownsArray,double** tempStateSizedVectors,double** tempStateSizedSquareMatrices,const double dt,const int normalNonZeroIndex) {
  assertion2(normalNonZeroIndex>=0,dt,normalNonZeroIndex);
  assertion2(normalNonZeroIndex<DIMENSIONS,dt,normalNonZeroIndex);

  // First call Fortran as we know it works and doesn't corrupt data.
  // output quantities
  kernels::darray fortFL(basisSize, basisSize, nVar);
  kernels::darray fortFR(basisSize, basisSize, nVar);
  std::memcpy(fortFL.data, FL, fortFL.idx.size * sizeof(double));
  std::memcpy(fortFL.data, FR, fortFL.idx.size * sizeof(double));
  /*
  // input quantities. May be copied for security. Not used do far.
  double fort_QL[basissize*basissize*nVar];
  double fort_QR[basissize*basissize*nVar];
  // temporary storage for computation
  double *fort_tempFaceUnknownsArray = NULL; // This is never used.
  const int numTempStateSizedVectors = 5; // cf. ADERDGSolver.h, line 738ff.
  double *fort_tempStateSizedVectors[numTempStateSizedVectors];
  double v1[numTempStateSizedVectors][nVar];
  for(int i=0; i<numTempStateSizedVectors; i++) {
    fort_tempStateSizedVectors[i] = v1[i];
  }
  double **fort_tempStateSizedSquareMatrices = NULL; // Never used
  */

  kernels::aderdg::generic::fortran::riemannSolverNonlinear<MHDSolver>(*this,fortFL.data,fortFL.data,QL,QR,tempFaceUnknownsArray,tempStateSizedVectors,tempStateSizedSquareMatrices,dt,normalNonZeroIndex);
  kernels::aderdg::generic::c::riemannSolverNonlinear<MHDSolver>(*this,FL,FR,QL,QR,tempFaceUnknownsArray,tempStateSizedVectors,tempStateSizedSquareMatrices,dt,normalNonZeroIndex);

  riemannSolverCheck.count();
  riemannSolverCheck("fortFL", fortFL.data, "FL", FL, fortFL.idx, 1e-10, true);
  riemannSolverCheck("fortFR", fortFR.data, "FR", FR, fortFR.idx, 1e-10, true);
}


//check boundaryConditionsCheck("boundaryConditions", "checks/boundaryConditions.txt");
void MHDSolver::MHDSolver::boundaryConditions(double* fluxOut,double* stateOut,const double* const fluxIn,const double* const stateIn,const tarch::la::Vector<DIMENSIONS,double>& cellCentre,const tarch::la::Vector<DIMENSIONS,double>& cellSize,const double t,const double dt,const int faceIndex,const int normalNonZero) {
  //kernels::darray F_fluxOut(nVar), F_stateOut(nVar);

  // There is no FORTRAN version of the boundaryConditions!
  kernels::aderdg::generic::c::boundaryConditions<MHDSolver>(*this,fluxOut,stateOut,fluxIn,stateIn,cellCentre,cellSize,t,dt,faceIndex,normalNonZero);
  
  //boundaryConditionsCheck.count();
  //boundaryConditionsCheck("FORT_fluxOut", F_fluxOut.data, "fluxOut", fluxOut, F_fluxOut.idx, 1e-10, true);
  //boundaryConditionsCheck("FORT_stateOut", F_stateOut.data, "stateOut", stateOut, F_stateOut.idx, 1e-10, true);
}



double MHDSolver::MHDSolver::stableTimeStepSize(const double* const luh,double* tempEigenvalues,const tarch::la::Vector<DIMENSIONS,double>& dx) {
  double d = kernels::aderdg::generic::c::stableTimeStepSize<MHDSolver>(*this,luh,tempEigenvalues,dx);
  return d;
}



void MHDSolver::MHDSolver::solutionAdjustment(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {
  kernels::aderdg::generic::c::solutionAdjustment<MHDSolver>(*this,luh,center,dx,t,dt);
}



void MHDSolver::MHDSolver::faceUnknownsProlongation(double* lQhbndFine,double* lFhbndFine,const double* lQhbndCoarse,const double* lFhbndCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1,int>& subfaceIndex) {
  kernels::aderdg::generic::c::faceUnknownsProlongation(lQhbndFine,lFhbndFine,lQhbndCoarse,lFhbndCoarse,coarseGridLevel,fineGridLevel,subfaceIndex,getNumberOfVariables(),getNodesPerCoordinateAxis());
}



void MHDSolver::MHDSolver::faceUnknownsRestriction(double* lQhbndCoarse,double* lFhbndCoarse,const double* lQhbndFine,const double* lFhbndFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1,int>& subfaceIndex) {
  kernels::aderdg::generic::c::faceUnknownsRestriction(lQhbndCoarse,lFhbndCoarse,lQhbndFine,lFhbndFine,coarseGridLevel,fineGridLevel,subfaceIndex,getNumberOfVariables(),getNodesPerCoordinateAxis());
}



void MHDSolver::MHDSolver::volumeUnknownsProlongation(double* luhFine,const double* luhCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) {
  kernels::aderdg::generic::c::volumeUnknownsProlongation(luhFine,luhCoarse,coarseGridLevel,fineGridLevel,subcellIndex,getNumberOfVariables(),getNodesPerCoordinateAxis());
}



void MHDSolver::MHDSolver::volumeUnknownsRestriction(double* luhCoarse,const double* luhFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) {
  kernels::aderdg::generic::c::volumeUnknownsRestriction(luhCoarse,luhFine,coarseGridLevel,fineGridLevel,subcellIndex,getNumberOfVariables(),getNodesPerCoordinateAxis());
}