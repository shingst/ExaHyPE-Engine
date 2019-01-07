// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "ErrorWriter.h"
#include "Tools.h"
#include "PDE.h"
#include "GRMHDbSolver_ADERDG.h"
#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/KernelUtils.h"
#include <cmath>

#include <fstream>
#include <iostream>
using namespace std;

#include "peano/utils/Loop.h"
#include "kernels/aderdg/generic/c/sizes.cpph"
#include "kernels/KernelUtils.h" // matrix indexing

#include "tarch/la/VectorOperations.h"

#include <algorithm>

#include <iomanip>

// see ch.18 Post-processing of the Guidebook.
//#include <mpi.h> 
#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"

GRMHDb::ErrorWriter::ErrorWriter() : exahype::plotters::ADERDG2UserDefined::ADERDG2UserDefined(){
  // @TODO Please insert your code here.
}


void GRMHDb::ErrorWriter::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
    double timeStamp) {
  // @TODO Please insert your code here.
  constexpr int numberOfVariables  = AbstractGRMHDbSolver_ADERDG::NumberOfVariables;
  constexpr int numberOfParameters = AbstractGRMHDbSolver_ADERDG::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
  constexpr int basisSize          = AbstractGRMHDbSolver_ADERDG::Order+1;
  constexpr int order              = basisSize-1;

  double x[DIMENSIONS];

  kernels::idx4 idx(basisSize,basisSize,basisSize,numberOfData);
  dfor(i,basisSize) {
     double w_dV = 1.0;
     for (int d=0; d<DIMENSIONS; d++) {
       x[d]  = offsetOfPatch[d] + sizeOfPatch[d] * kernels::gaussLegendreNodes[order][i(d)];
       w_dV *= sizeOfPatch[d] * kernels::gaussLegendreWeights[order][i(d)];
     }

     double uAna[numberOfVariables];
     GRMHDbSolver_ADERDG::referenceSolution(x,timeStamp,uAna);

     const double* uNum = u + idx ( (DIMENSIONS==3) ? i(2) : 0, i(1), i(0), 0);

     for (int v=0; v<numberOfVariables; v++) {
        const double uDiff = std::abs(uNum[v]-uAna[v]);
        errorL2[v]   += uDiff*uDiff * w_dV;
        errorL1[v]   += uDiff * w_dV;
        errorLInf[v]  = std::max( errorLInf[v], uDiff );

        normL1Ana[v]  += std::abs(uAna[v]) * w_dV;
        normL2Ana[v]  += uAna[v] * uAna[v] * w_dV;
        normLInfAna[v] = std::max( normLInfAna[v], std::abs(uAna[v]) );
     }
  }
}

void GRMHDb::ErrorWriter::startPlotting( double time) {
  // @TODO Please insert your code here.
  _timeStamp = time;

  std::fill_n(errorL1,  AbstractGRMHDbSolver_ADERDG::NumberOfVariables, 0.0);
  std::fill_n(errorL2,  AbstractGRMHDbSolver_ADERDG::NumberOfVariables, 0.0);
  std::fill_n(errorLInf,AbstractGRMHDbSolver_ADERDG::NumberOfVariables, 0.0);
  
  std::fill_n(normL1Ana,  AbstractGRMHDbSolver_ADERDG::NumberOfVariables, 0.0);
  std::fill_n(normL2Ana,  AbstractGRMHDbSolver_ADERDG::NumberOfVariables, 0.0);
  std::fill_n(normLInfAna,AbstractGRMHDbSolver_ADERDG::NumberOfVariables, 0.0);
}

void GRMHDb::ErrorWriter::finishPlotting() {
  constexpr int numberOfVariables = AbstractGRMHDbSolver_ADERDG::NumberOfVariables;
  ofstream myfile;
  int mpirank = tarch::parallel::Node::getInstance().getRank();
  const int myMessageTagUseForTheReduction = tarch::parallel::Node::getInstance().reserveFreeTag("GRMHDb::ErrorWriter::finishPlotting()");


//#define Parallel
#ifdef Parallel 
  //int NumberOfNodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  const int numberOfAvailableRanks =
      tarch::parallel::Node::getInstance().getNumberOfNodes();
  const int numberOfData = 3 * numberOfVariables;
  const int numberOfDataTot = 3 * numberOfVariables * numberOfAvailableRanks;
  int tag;
  //int l1tag;
  //int l2tag;
  //int linftag;
  MPI_Request rq_recv[numberOfAvailableRanks];
  MPI_Request rq_send[numberOfAvailableRanks];
  //MPI_Request rq_L2  [numberOfAvailableRanks];
  //MPI_Request rq_Linf[numberOfAvailableRanks];
  double	  receivedValues[numberOfAvailableRanks][numberOfData];
  //double      receivedValueL1  [numberOfAvailableRanks][numberOfData];
  //double      receivedValueL2  [numberOfAvailableRanks][numberOfData];
  //double      receivedValueLinf[numberOfAvailableRanks][numberOfData]; 
  double SentValues[numberOfData];
  //std::fill_n(SentValues, numberOfData, 0.0);
  //std::fill_n(receivedValues, numberOfDataTot, 0.0);
  tag = 333;
  //l1tag = 111;
  //l2tag = 222;
  //linftag = 333;
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
	  // then this is the master rank
    for (int rank = 1; rank < numberOfAvailableRanks; rank++) {
	  //
      if (!tarch::parallel::NodePool::getInstance().isIdleNode(rank)) {
		//
        MPI_Irecv(&receivedValues[rank-1][0]  , numberOfData, MPI_DOUBLE, rank, tag+rank  , tarch::parallel::Node::getInstance().getCommunicator(), &rq_recv[rank-1]);
        //MPI_Irecv(&receivedValueL2[rank-1][0]  , numberOfVariables, MPI_DOUBLE, rank, l2tag  , tarch::parallel::Node::getInstance().getCommunicator(), &rq_L2[rank-1]);
        //MPI_Irecv(&receivedValueLinf[rank-1][0], numberOfVariables, MPI_DOUBLE, rank, linftag, tarch::parallel::Node::getInstance().getCommunicator(), &rq_Linf[rank-1]);
      }
    }
    for (int rank = 1; rank < numberOfAvailableRanks; rank++) {
      //
      if (!tarch::parallel::NodePool::getInstance().isIdleNode(rank)) {
        //
        MPI_Wait(&rq_recv[rank-1], MPI_STATUS_IGNORE);
        //MPI_Wait(&rq_L2[rank-1], MPI_STATUS_IGNORE);
        //MPI_Wait(&rq_Linf[rank-1], MPI_STATUS_IGNORE);
        for (int i = 0; i < numberOfVariables; i++) {
          errorL1[i] += receivedValues[rank - 1][3*i];
          errorL2[i] += receivedValues[rank - 1][3 * i + 1];
          errorLInf[i] = std::max(receivedValues[rank - 1][3 * i + 2],errorLInf[i]);
			//errorL2  [i] += receivedValueL2  [rank-1][i];
            //errorLInf[i] =std::max(receivedValueLinf[rank-1][i],errorLInf[i]); 
		}
      }
    }
  } else {
		// then this is a PDE-worker rank.
		for (int i = 0; i < numberOfVariables; i++) {
		              SentValues[3*i] = errorL1[i];
		              SentValues[3*i+1] = errorL2[i];
		              SentValues[3*i+2] = errorLInf[i];
		} 
		MPI_Isend(&SentValues, numberOfData, MPI_DOUBLE,tarch::parallel::Node::getGlobalMasterRank(), tag+mpirank  , tarch::parallel::Node::getInstance().getCommunicator(), &rq_send[mpirank-1]);
        MPI_Wait(&rq_send[mpirank-1], MPI_STATUS_IGNORE);
		//MPI_Isend(&errorL2  , numberOfVariables, MPI_DOUBLE,tarch::parallel::Node::getGlobalMasterRank(), l2tag  , tarch::parallel::Node::getInstance().getCommunicator(), &rq_L2[mpirank-1]);
		//MPI_Isend(&errorLInf, numberOfVariables, MPI_DOUBLE,tarch::parallel::Node::getGlobalMasterRank(), linftag, tarch::parallel::Node::getInstance().getCommunicator(), &rq_Linf[mpirank-1]);
  }
  tarch::parallel::Node::getInstance().releaseTag(myMessageTagUseForTheReduction);
#endif

  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
	    //
		for (int v = 0; v < numberOfVariables; v++) {
		  errorL2[v] = sqrt(errorL2[v]);
		  normL2Ana[v] = sqrt(normL2Ana[v]);
		}
		//
		if (tarch::la::equals(_timeStamp, 0.0)) {
			myfile.open("./output/ErrorNorms.dat", ios::trunc);
			myfile << "*********************************************" << std::endl;
			myfile << "**Errors for ADER-DG solver with order=" << AbstractGRMHDbSolver_ADERDG::Order << "**" << std::endl;
			myfile << "*********************************************" << std::endl;
			myfile << "---------------------------------------------" << std::endl;
                        myfile << "variable:\t";
			for (int v = 0; v < numberOfVariables; v++) {
				myfile << v << " \t ";
			}
			myfile << std::endl;
			myfile.close();
			//
		}

		myfile.open("./output/ErrorNorms.dat", ios::app);
		myfile << "*********************************************" << std::endl;
		myfile << "t_eval : " << _timeStamp << std::endl;
		/*
		std::cout << "**Errors for ADER-DG solver with order="<<AbstractGRMHDbSolver_ADERDG::Order<<"**" << std::endl;
		std::cout << "t_eval : "<<_timeStamp << std::endl;
		std::cout << "variable     : ";
		for (int v=0; v<numberOfVariables; v++) {
		  std::cout << v << " \t ";
		}
		std::cout << std::endl;

		std::cout << "absErrorL1   : ";
		for (int v=0; v<numberOfVariables; v++) {
		  std::cout << std::setprecision(2) << errorL1[v] << " \t ";
		}
		std::cout << std::endl;
		*/
		/*
		  this is for the output to file
		*/

		myfile << "ErrorL1:\t";
		for (int v = 0; v < numberOfVariables; v++) {
		  myfile << std::scientific << errorL1[v] << " \t ";
		}
		myfile << std::endl;

		myfile << "ErrorL2:\t";
		for (int v = 0; v < numberOfVariables; v++) {
		  myfile << std::scientific << errorL2[v] << " \t ";
		}
		myfile << std::endl;

		myfile << "errorLInf:\t";
		for (int v = 0; v < numberOfVariables; v++) {
		  myfile << std::scientific << errorLInf[v] << " \t ";
		}
		myfile << std::endl;
		 
		myfile.close();  
	}
}

