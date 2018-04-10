/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#include "exahype/solvers/TemporaryVariables.h"

#include "exahype/solvers/Solver.h"
#include "exahype/solvers/LimitingADERDGSolver.h"

#include <algorithm>
#include <mm_malloc.h> //g++
#include <cstring> //memset

#include "tarch/Assertions.h"
#include "tarch/multicore/MulticoreDefinitions.h"

exahype::solvers::PredictionTemporaryVariables::PredictionTemporaryVariables() {}


double* exahype::solvers::allocateArray( std::vector<int>& heapIndices, const int size ) {
  //don't need to allocate anything
  if(size < 1) {
    return nullptr;
  }

  tarch::multicore::Lock lock(exahype::HeapSemaphore);
  const int heapIndex = exahype::DataHeap::getInstance().createData(size,size,
    exahype::DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired
  );
  auto& vector = exahype::DataHeap::getInstance().getData(heapIndex);
  lock.free();

  assertionEquals(static_cast<int>(vector.size()),size);
  std::fill(vector.begin(),vector.end(),0.0);

  heapIndices.push_back(heapIndex);
  return vector.data();
}


void exahype::solvers::freeArrays( std::vector<int>& heapIndices ) {
  for (int i : heapIndices) {
    tarch::multicore::Lock lock(exahype::HeapSemaphore);
    assertion(exahype::DataHeap::getInstance().isValidIndex(i));
    const bool doRecycle = true;
    exahype::DataHeap::getInstance().deleteData(i,doRecycle);
    lock.free();
  }

  heapIndices.clear();
}

void exahype::solvers::initialiseTemporaryVariables(exahype::solvers::PredictionTemporaryVariables& temporaryVariables) {
  assertion(temporaryVariables._dataHeapIndices.size()==0);
  assertion(temporaryVariables._tempSpaceTimeUnknowns    ==nullptr);
  assertion(temporaryVariables._tempSpaceTimeFluxUnknowns==nullptr);
  assertion(temporaryVariables._tempUnknowns             ==nullptr);
  assertion(temporaryVariables._tempFluxUnknowns         ==nullptr);
  assertion(temporaryVariables._tempPointForceSources    ==nullptr);

  const int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
  temporaryVariables._dataHeapIndices.reserve(numberOfSolvers*7); // count the arrays

  temporaryVariables._tempSpaceTimeUnknowns     = new double**[numberOfSolvers]; // == lQi, rhs
  temporaryVariables._tempSpaceTimeFluxUnknowns = new double**[numberOfSolvers]; // == lFi, gradQ
  temporaryVariables._tempUnknowns              = new double* [numberOfSolvers]; // == lQhi
  temporaryVariables._tempFluxUnknowns          = new double* [numberOfSolvers]; // == lFhi
  temporaryVariables._tempPointForceSources     = new double**[numberOfSolvers]; // == PSi, forceVectorSourceN 

  exahype::solvers::ADERDGSolver* aderdgSolver = nullptr;

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    switch( solver->getType() ) {
    case exahype::solvers::Solver::Type::ADERDG:
      aderdgSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);
      break;
    case exahype::solvers::Solver::Type::LimitingADERDG:
      aderdgSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver().get();
      break;
    default:
      aderdgSolver = nullptr;
      break;
    }

    if (aderdgSolver!=nullptr) {
      temporaryVariables._tempSpaceTimeUnknowns[solverNumber] = new double*[2];
      for (int i=0; i<2; ++i) { // max; see spaceTimePredictorNonlinear
        temporaryVariables._tempSpaceTimeUnknowns[solverNumber][i] = allocateArray( temporaryVariables._dataHeapIndices,
                                                                                    aderdgSolver->getTempSpaceTimeUnknownsSize() );
      }
      //
      temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber] = new double*[2];
      //lFi
      temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber][0] = allocateArray( temporaryVariables._dataHeapIndices,
                                                                                        aderdgSolver->getTempSpaceTimeFluxUnknowns0Size() );
      //gradQ, might be not required if ncp is not used
      temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber][1] = allocateArray( temporaryVariables._dataHeapIndices,
                                                                                        aderdgSolver->getTempSpaceTimeFluxUnknowns1Size() );
      //
      temporaryVariables._tempUnknowns[solverNumber] = allocateArray( temporaryVariables._dataHeapIndices,
                                                                      aderdgSolver->getTempUnknownsSize() );
      //
      temporaryVariables._tempFluxUnknowns[solverNumber] = allocateArray( temporaryVariables._dataHeapIndices,
                                                                          aderdgSolver->getTempFluxUnknownsSize() );
      //
      temporaryVariables._tempPointForceSources[solverNumber] = new double*[2];
      temporaryVariables._tempPointForceSources[solverNumber][0] = allocateArray( temporaryVariables._dataHeapIndices,
                                                                               aderdgSolver->getTempPointForceSourcesSize() ); //PSi     
      temporaryVariables._tempPointForceSources[solverNumber][1] = allocateArray( temporaryVariables._dataHeapIndices,
                                                                               aderdgSolver->getTempPointForceSourcesSize() ); //forceVectorSourceN
    } else {
      temporaryVariables._tempSpaceTimeUnknowns    [solverNumber] = nullptr;
      temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber] = nullptr;
      temporaryVariables._tempUnknowns             [solverNumber] = nullptr;
      temporaryVariables._tempFluxUnknowns         [solverNumber] = nullptr;
      temporaryVariables._tempPointForceSources    [solverNumber] = nullptr;
    }
  }
}

void exahype::solvers::deleteTemporaryVariables(exahype::solvers::PredictionTemporaryVariables& temporaryVariables) {
  if (temporaryVariables._tempSpaceTimeUnknowns!=nullptr) {
    // assertion(!temporaryVariables._dataHeapIndices.empty()); This does not hold if all solvers are FV
    assertion(temporaryVariables._tempSpaceTimeFluxUnknowns!=nullptr);
    assertion(temporaryVariables._tempPointForceSources    !=nullptr);

    exahype::solvers::ADERDGSolver* aderdgSolver = nullptr;

    // release memory
    freeArrays( temporaryVariables._dataHeapIndices );

    // set the pointers to null
    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

      switch( solver->getType() ) {
      case exahype::solvers::Solver::Type::ADERDG:
        aderdgSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);
        break;
      case exahype::solvers::Solver::Type::LimitingADERDG:
        aderdgSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver().get();
        break;
      default:
        aderdgSolver = nullptr;
        break;
      }

      if (aderdgSolver!=nullptr) {
          //
          for (int i=0; i<2; ++i) {
            temporaryVariables._tempSpaceTimeUnknowns[solverNumber][i] = nullptr;
          }
          delete[] temporaryVariables._tempSpaceTimeUnknowns[solverNumber];
          temporaryVariables._tempSpaceTimeUnknowns[solverNumber] = nullptr;
          //
          for (int i=0; i<2; ++i) {
            temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber][i] = nullptr;
          }
          delete[] temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber];
          temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber] = nullptr;
          //
          temporaryVariables._tempUnknowns[solverNumber] = nullptr;
          //
          temporaryVariables._tempFluxUnknowns[solverNumber] = nullptr;
          //
           for (int i=0; i<2; ++i) {
            temporaryVariables._tempPointForceSources[solverNumber][i] = nullptr;
          }
          delete[] temporaryVariables._tempPointForceSources[solverNumber];
          temporaryVariables._tempPointForceSources[solverNumber] = nullptr;
      }
    }

    delete[] temporaryVariables._tempSpaceTimeUnknowns;
    delete[] temporaryVariables._tempSpaceTimeFluxUnknowns;
    delete[] temporaryVariables._tempUnknowns;
    delete[] temporaryVariables._tempFluxUnknowns;
    delete[] temporaryVariables._tempPointForceSources;
    temporaryVariables._tempSpaceTimeUnknowns     = nullptr;
    temporaryVariables._tempSpaceTimeFluxUnknowns = nullptr;
    temporaryVariables._tempUnknowns              = nullptr;
    temporaryVariables._tempFluxUnknowns          = nullptr;
    temporaryVariables._tempPointForceSources     = nullptr;
  }
}

exahype::solvers::MergingTemporaryVariables::MergingTemporaryVariables() {
}

void exahype::solvers::initialiseTemporaryVariables(exahype::solvers::MergingTemporaryVariables& temporaryVariables) {
  assertion(temporaryVariables._dataHeapIndices.size()==0);
  assertion(temporaryVariables._tempFaceUnknowns==nullptr);

  const int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
  temporaryVariables._dataHeapIndices.reserve(numberOfSolvers*3);

  temporaryVariables._tempFaceUnknowns = new double**[numberOfSolvers];

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    int numberOfFaceUnknowns = 0;
    int lengthOfFaceUnknowns = 0;

    switch (solver->getType()) {
      case exahype::solvers::Solver::Type::ADERDG:
        numberOfFaceUnknowns       = 3; // See exahype::solvers::ADERDGSolver::applyBoundaryConditions
        lengthOfFaceUnknowns       =
            static_cast<exahype::solvers::ADERDGSolver*>(solver)->getBndFaceSize(); // == getDataPerFace() + eventual padding
        break;
      case exahype::solvers::Solver::Type::LimitingADERDG:
        // Needs the same temporary variables as the normal ADER-DG scheme.
        numberOfFaceUnknowns       = 3;
        lengthOfFaceUnknowns       = std::max(
            static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver()->getBndFaceSize(), // == getDataPerFace() + eventual padding
            static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiter()->getDataPerPatchFace() );
        break;
      case exahype::solvers::Solver::Type::FiniteVolumes:
        numberOfFaceUnknowns = 2; // See exahype::solvers::FiniteVolumesSolver::mergeWithBoundaryData
        lengthOfFaceUnknowns =
            static_cast<exahype::solvers::FiniteVolumesSolver*>(solver)->getDataPerPatchFace();
        break;
    }
    //
    temporaryVariables._tempFaceUnknowns[solverNumber] = nullptr;
    if (numberOfFaceUnknowns>0) {
      temporaryVariables._tempFaceUnknowns[solverNumber] = new double*[numberOfFaceUnknowns];
      for (int i=0; i<numberOfFaceUnknowns; ++i) { // see ADERDGSolver::applyBoundaryConditions(...)
        temporaryVariables._tempFaceUnknowns[solverNumber][i] = allocateArray( temporaryVariables._dataHeapIndices, lengthOfFaceUnknowns );
      }
    }
  }
}

void exahype::solvers::deleteTemporaryVariables(exahype::solvers::MergingTemporaryVariables& temporaryVariables) {
  if (temporaryVariables._tempFaceUnknowns!=nullptr) {
    assertion(!temporaryVariables._dataHeapIndices.empty());
    assertion(temporaryVariables._tempFaceUnknowns!=nullptr);

    // release memory
    freeArrays(temporaryVariables._dataHeapIndices);

    // set the pointers to null
    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

      int numberOfFaceUnknowns = 0;
      switch (solver->getType()) {
        case exahype::solvers::Solver::Type::ADERDG:
          numberOfFaceUnknowns = 3; // See exahype::solvers::ADERDGSolver::applyBoundaryConditions
          break;
        case exahype::solvers::Solver::Type::LimitingADERDG:
          numberOfFaceUnknowns = 3; // See exahype::solvers::ADERDGSolver::applyBoundaryConditions
          break;
        case exahype::solvers::Solver::Type::FiniteVolumes:
          numberOfFaceUnknowns = 2; // See exahype::solvers::FiniteVolumesSolver::mergeWithBoundaryData
          break;
      }

      if (numberOfFaceUnknowns>0) {
        for (int i=0; i<numberOfFaceUnknowns; ++i) { // see ADERDGSolver::applyBoundaryConditions(...)
          temporaryVariables._tempFaceUnknowns[solverNumber][i] = nullptr;
        }
        delete[] temporaryVariables._tempFaceUnknowns[solverNumber];
        temporaryVariables._tempFaceUnknowns[solverNumber] = nullptr;
      }
    }

    delete[] temporaryVariables._tempFaceUnknowns;
    temporaryVariables._tempFaceUnknowns = nullptr;
  }
}
