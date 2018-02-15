// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "ErrorWriter.h"
#include "PDE.h"
#include "kernels/GaussLegendreQuadrature.h"
#include <cmath>

#include "kernels/aderdg/generic/c/sizes.cpph"
#include "kernels/KernelUtils.h" // matrix indexing


//#include <iostream.h>

GPR::ErrorWriter::ErrorWriter(GPR::GPRSolver_ADERDG& solver) : ErrorWriter(){
	plotForADERSolver = true;
}

GPR::ErrorWriter::ErrorWriter(GPR::GPRSolver_FV& solver) : ErrorWriter(){
	plotForADERSolver = false;
}

GPR::ErrorWriter::ErrorWriter() : exahype::plotters::ADERDG2UserDefined::ADERDG2UserDefined(),
	errors("output/error-"){
  // @TODO Please insert your code here.
	char name[10];
	for(int m=0; m < nVar; m++) {
		pdevarname_(name,&m);
		//printf("***********************************************************");
		//printf(name);
		//printf("\n");
		errors.add(m, name);
    }

}


void GPR::ErrorWriter::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
    double timeStamp) {

  // @TODO Please insert your code here.

	kernels::idx3 id_xy_dof(basisSize,basisSize,nVar);
	double dV=sizeOfPatch[0]*sizeOfPatch[1];
  	//double dV=sizeOfPatch[0]*sizeOfPatch[1]*sizeOfPatch[2];
	/*printf("=====================================================\n");
	for(int i = 0; i < basisSize; i++){
		for(int j = 0; j <basisSize ; j++){
			printf("%d %d : ", i,j);
			for(int k = 0; k <nVar ; k++){
				printf("%e ", u[id_xy_dof(i,j,:)]);
				//printf(" . ");
			}
		printf("\n");
		}
	}
	printf("=====================================================\n");*/
	
	double localError[nVar]={0.};
	double w_x, w_y;
	double pos[DIMENSIONS] = {0.};
	
	for(int i = 0; i < basisSize; i++){
		w_x= kernels::gaussLegendreWeights[order][i];
		pos[0] = kernels::gaussLegendreNodes[order][i]*sizeOfPatch[0]+offsetOfPatch[0];
		for(int j = 0; j <basisSize ; j++){
			w_y = kernels::gaussLegendreWeights[order][j];
			pos[1] = kernels::gaussLegendreNodes[order][j]*sizeOfPatch[1]+offsetOfPatch[1];
			
			double numerical[nVar];
			//getNumericalSolution(numerical,&u[id_xy_dof(i,j,0)]);
			getnumericalsolution_(numerical,&u[id_xy_dof(i,j,0)]);
			double exact[nVar];
			//getExactSolution(exact,pos,timeStamp);
			getexactsolution_(exact,pos,&timeStamp);
						
			for(int k = 0; k <nVar ; k++){
				localError[k] += std::abs(numerical[k]-exact[k]) * w_x * w_y;
			}	
		}
	}
	
	/*if(plotForADERSolver) {
		const int order = GPR::AbstractGPRSolver_ADERDG::Order;
		dV = kernels::ADERDGVolume(order, sizeOfPatch, pos);
	} else {
		const int patchSize = GPR::AbstractGPRSolver_FV::PatchSize;
		dV = tarch::la::volume(sizeOfPatch)/patchSize; // correct is probably (patchSize+1)
	}*/
	// printf("dV=%e\n", dV);
	// printf("localError=%e\n", localError[0]);
  	errors.addValue(localError, dV);
}


void GPR::ErrorWriter::startPlotting( double time) {
  // @TODO Please insert your code here.
	printf("***********************************************************\n");
	printf("I am starting the error writer:");
    errors.startRow(time);
   printf("Done \n");
   printf("***********************************************************\n");
}


void GPR::ErrorWriter::finishPlotting() {
  // @TODO Please insert your code here.
   printf("***********************************************************\n");
   printf("I am finishing the error writer:");
   errors.finishRow();
   printf("Done \n");
   printf("***********************************************************\n");
}

