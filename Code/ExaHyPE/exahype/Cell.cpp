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
 
#include "exahype/Cell.h"
#include "exahype/State.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "tarch/la/ScalarOperations.h"

#include "kernels/KernelCalls.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

#include "exahype/records/ADERDGCellDescription.h"

tarch::logging::Log exahype::Cell::_log("exahype::Cell");

exahype::Cell::Cell() : Base() {
  _cellData.setCellDescriptionsIndex(
      multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
}

exahype::Cell::Cell(const Base::DoNotCallStandardConstructor& value)
    : Base(value) {
  _cellData.setCellDescriptionsIndex(
      multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
}

exahype::Cell::Cell(const Base::PersistentCell& argument) : Base(argument) {
  // do nothing
}

int exahype::Cell::getCellDescriptionsIndex() const {
  return _cellData.getCellDescriptionsIndex();
}

void exahype::Cell::setupMetaData() {
  assertion1(!ADERDGCellDescriptionHeap::getInstance().isValidIndex(_cellData.getCellDescriptionsIndex()),toString());

  const int CellDescriptionIndex = ADERDGCellDescriptionHeap::getInstance().createData(0, 0);
  FiniteVolumesCellDescriptionHeap::getInstance().createDataForIndex(CellDescriptionIndex,0,0);

  _cellData.setCellDescriptionsIndex(CellDescriptionIndex);
}

void exahype::Cell::shutdownMetaData() {
  assertion1(
    ADERDGCellDescriptionHeap::getInstance().isValidIndex(
      _cellData.getCellDescriptionsIndex()
    ),
    toString());

  ADERDGCellDescriptionHeap::getInstance().deleteData(_cellData.getCellDescriptionsIndex());
  FiniteVolumesCellDescriptionHeap::getInstance().deleteData(_cellData.getCellDescriptionsIndex());

  _cellData.setCellDescriptionsIndex(multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
}


bool exahype::Cell::isInitialised() const {
  if (_cellData.getCellDescriptionsIndex() != multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex) {
    assertion( ADERDGCellDescriptionHeap::getInstance().isValidIndex(_cellData.getCellDescriptionsIndex()) );
    assertion( FiniteVolumesCellDescriptionHeap::getInstance().isValidIndex(_cellData.getCellDescriptionsIndex()) );
  }  // dead code elimination will get rid of this loop if Asserts flag is not set

  return _cellData.getCellDescriptionsIndex() != multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex;
}


void exahype::Cell::addNewCellDescription(
    const int solverNumber,
    const exahype::records::FiniteVolumesCellDescription::Type cellType,
/*
    const exahype::records::FiniteVolumesCellDescription::RefinementEvent
        refinementEvent,
*/
    const int level, const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, double>& size,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre) {
  if (_cellData.getCellDescriptionsIndex() == multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex) {
    setupMetaData();
  }

  assertion1(
    FiniteVolumesCellDescriptionHeap::getInstance().isValidIndex(
      _cellData.getCellDescriptionsIndex()
    ),
    toString());

  assertion2(static_cast<unsigned int>(solverNumber) <
                 solvers::RegisteredSolvers.size(),
             solverNumber, exahype::solvers::RegisteredSolvers.size());

  //const solvers::Solver* solver = solvers::RegisteredSolvers[solverNumber];
  assertion(solver->getType()==exahype::solvers::Solver::Type::FiniteVolumes);

  exahype::records::FiniteVolumesCellDescription newCellDescription;
  newCellDescription.setSolverNumber(solverNumber);

  // Default AMR settings
  newCellDescription.setType(cellType);
  newCellDescription.setLevel(level);
  //newCellDescription.setRefinementEvent(refinementEvent);

  // Pass geometry information to the cellDescription description
  newCellDescription.setSize(size);
  newCellDescription.setOffset(cellCentre);

  // Default field data indices
  newCellDescription.setSolution(-1);

  FiniteVolumesCellDescriptionHeap::getInstance()
      .getData(_cellData.getCellDescriptionsIndex())
      .push_back(newCellDescription);

}


void exahype::Cell::addNewCellDescription(
  const int                                     solverNumber,
  const exahype::records::ADERDGCellDescription::Type cellType,
  const exahype::records::ADERDGCellDescription::RefinementEvent refinementEvent,
  const int                                     level,
  const int                                     parentIndex,
  const tarch::la::Vector<DIMENSIONS, double>&  size,
  const tarch::la::Vector<DIMENSIONS, double>&  cellCentre) {
  if (_cellData.getCellDescriptionsIndex() == multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex) {
    setupMetaData();
  }

  assertion1(
    ADERDGCellDescriptionHeap::getInstance().isValidIndex(
      _cellData.getCellDescriptionsIndex()
    ),
    toString());

  assertion2(parentIndex == -1 ||
             parentIndex != _cellData.getCellDescriptionsIndex(),
             parentIndex, _cellData.getCellDescriptionsIndex());

  assertion2(parentIndex != _cellData.getCellDescriptionsIndex(),
             parentIndex, _cellData.getCellDescriptionsIndex());

  assertion2(static_cast<unsigned int>(solverNumber) <
                 solvers::RegisteredSolvers.size(),
             solverNumber, exahype::solvers::RegisteredSolvers.size());

  //const solvers::Solver* solver = solvers::RegisteredSolvers[solverNumber];
  assertion(solver->getType()==exahype::solvers::Solver::Type::ADER_DG);

  exahype::records::ADERDGCellDescription newCellDescription;
  newCellDescription.setSolverNumber(solverNumber);

  // Default AMR settings
  newCellDescription.setType(cellType);
  newCellDescription.setParentIndex(parentIndex);
  newCellDescription.setLevel(level);
  newCellDescription.setRefinementEvent(refinementEvent);

  std::bitset<DIMENSIONS_TIMES_TWO>
      riemannSolvePerformed;  // default construction: no bit set
  newCellDescription.setRiemannSolvePerformed(riemannSolvePerformed);

  // Pass geometry information to the cellDescription description
  newCellDescription.setSize(size);
  newCellDescription.setOffset(cellCentre);

  // Default field data indices
  newCellDescription.setSpaceTimePredictor(-1);
  newCellDescription.setSpaceTimeVolumeFlux(-1);
  newCellDescription.setPredictor(-1);
  newCellDescription.setVolumeFlux(-1);
  newCellDescription.setSolution(-1);
  newCellDescription.setUpdate(-1);
  newCellDescription.setExtrapolatedPredictor(-1);
  newCellDescription.setFluctuation(-1);

  // Limiter meta data (oscillations identificator)
  newCellDescription.setSolutionMin(std::numeric_limits<double>::min());
  newCellDescription.setSolutionMax(std::numeric_limits<double>::max());

  ADERDGCellDescriptionHeap::getInstance()
      .getData(_cellData.getCellDescriptionsIndex())
      .push_back(newCellDescription);
}


void exahype::Cell::ensureNecessaryMemoryIsAllocated(const int solverNumber) {
  // Ensure that -1 value is invalid index (cf. addNewCellDescription)
  assertion(!DataHeap::getInstance().isValidIndex(-1));

  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                 _cellData.getCellDescriptionsIndex()),
             toString());

  assertion3(static_cast<unsigned int>(solverNumber) <
                 solvers::RegisteredSolvers.size(),
             solverNumber, solvers::RegisteredSolvers.size(), toString());

  const solvers::Solver* solver = solvers::RegisteredSolvers[solverNumber];

  switch (solver->getType()) {
    case exahype::solvers::Solver::Type::ADER_DG:
      for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(
               _cellData.getCellDescriptionsIndex())) {
        if (solverNumber == p.getSolverNumber()) {
          switch (p.getType()) {
            case exahype::records::ADERDGCellDescription::Cell:
              if (!DataHeap::getInstance().isValidIndex(p.getSolution())) {
                assertion(!DataHeap::getInstance().isValidIndex(
                    p.getSpaceTimePredictor()));
                assertion(!DataHeap::getInstance().isValidIndex(
                    p.getSpaceTimeVolumeFlux()));
                assertion(
                    !DataHeap::getInstance().isValidIndex(p.getPredictor()));
                assertion(!DataHeap::getInstance().isValidIndex(p.getUpdate()));
                assertion(
                    !DataHeap::getInstance().isValidIndex(p.getVolumeFlux()));

                const int spaceTimeUnknownsPerCell =
                    static_cast<const exahype::solvers::ADERDGSolver*>(solver)->getSpaceTimeUnknownsPerCell();
                const int spaceTimeFluxUnknownsPerCell =
                    static_cast<const exahype::solvers::ADERDGSolver*>(solver)->getSpaceTimeFluxUnknownsPerCell();
                const int unknownsPerCell =
                    static_cast<const exahype::solvers::ADERDGSolver*>(solver)->getUnknownsPerCell();
                const int fluxUnknownsPerCell =
                    static_cast<const exahype::solvers::ADERDGSolver*>(solver)->getFluxUnknownsPerCell();

                // Allocate space-time DoF
                p.setSpaceTimePredictor(DataHeap::getInstance().createData(
                    spaceTimeUnknownsPerCell, spaceTimeUnknownsPerCell));
                p.setSpaceTimeVolumeFlux(DataHeap::getInstance().createData(
                    spaceTimeFluxUnknownsPerCell,
                    spaceTimeFluxUnknownsPerCell));

                // Allocate volume DoF
                p.setPredictor(DataHeap::getInstance().createData(
                    unknownsPerCell, unknownsPerCell));
                p.setVolumeFlux(DataHeap::getInstance().createData(
                    fluxUnknownsPerCell, fluxUnknownsPerCell));
                p.setUpdate(DataHeap::getInstance().createData(
                    unknownsPerCell, unknownsPerCell));
                p.setSolution(DataHeap::getInstance().createData(
                    unknownsPerCell, unknownsPerCell));

                // TODO(Dominic): Material parameters what are the correct array sizes?
              }
              break;
            default:
              break;
          }
          // Allocate face DoF
          switch (p.getType()) {
            case exahype::records::ADERDGCellDescription::Cell:
            case exahype::records::ADERDGCellDescription::Ancestor:
            case exahype::records::ADERDGCellDescription::Descendant:
              if (!DataHeap::getInstance().isValidIndex(
                      p.getExtrapolatedPredictor())) {
                assertion(
                    !DataHeap::getInstance().isValidIndex(p.getFluctuation()));

                const int unknownsPerCellBoundary =
                    static_cast<const exahype::solvers::ADERDGSolver*>(solver)->getUnknownsPerCellBoundary();

                p.setExtrapolatedPredictor(DataHeap::getInstance().createData(
                    unknownsPerCellBoundary, unknownsPerCellBoundary));
                p.setFluctuation(DataHeap::getInstance().createData(
                    unknownsPerCellBoundary, unknownsPerCellBoundary));
              }
              break;
            default:
              break;
          }
        }
      }
      break;
    case exahype::solvers::Solver::Type::FiniteVolumes:
      for (auto& p : FiniteVolumesCellDescriptionHeap::getInstance().getData(_cellData.getCellDescriptionsIndex())) {
        if (solverNumber == p.getSolverNumber()) {
          switch (p.getType()) {
            case exahype::records::FiniteVolumesCellDescription::Cell:
              const int unknownsPerCell =
                  static_cast<const exahype::solvers::FiniteVolumesSolver*>(solver)->getUnknownsPerCell();
              assertion(unknownsPerCell>0);
              p.setSolution(DataHeap::getInstance().createData(unknownsPerCell, unknownsPerCell));
              logDebug( "initialiseCellDescription(...)", "allocated " << unknownsPerCell << " records for a cell" );
              break;
          }
        }
      }
      break;
    default:
      logDebug("initialiseCellDescription(...)",
               "solver is not associated with any cell descriptions of this "
               "cell. cell="
                   << toString());
      break;
  }
}



void exahype::Cell::ensureNoUnnecessaryMemoryIsAllocated(const int solverNumber) {
  // Ensure that -1 value is invalid index (cf. addNewCellDescription)
  assertion(!DataHeap::getInstance().isValidIndex(-1));

  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                 _cellData.getCellDescriptionsIndex()),
             toString());

  assertion3(static_cast<unsigned int>(solverNumber) <
                 solvers::RegisteredSolvers.size(),
             solverNumber, solvers::RegisteredSolvers.size(), toString());

  const solvers::Solver* solver = solvers::RegisteredSolvers[solverNumber];

  switch (solver->getType()) {
    case exahype::solvers::Solver::Type::ADER_DG:
      for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(
               _cellData.getCellDescriptionsIndex())) {
        if (solverNumber == p.getSolverNumber()) {
          if (DataHeap::getInstance().isValidIndex(p.getSolution())) {
            switch (p.getType()) {
              case exahype::records::ADERDGCellDescription::Erased:
              case exahype::records::ADERDGCellDescription::EmptyAncestor:
              case exahype::records::ADERDGCellDescription::EmptyDescendant:
              case exahype::records::ADERDGCellDescription::Ancestor:
              case exahype::records::ADERDGCellDescription::Descendant:
                assertion1(
                    DataHeap::getInstance().isValidIndex(p.getSolution()),
                    p.getSolution());
                assertion1(DataHeap::getInstance().isValidIndex(
                               p.getSpaceTimePredictor()),
                           p.getSpaceTimePredictor());
                assertion(DataHeap::getInstance().isValidIndex(
                    p.getSpaceTimeVolumeFlux()));
                assertion(
                    DataHeap::getInstance().isValidIndex(p.getPredictor()));
                assertion(
                    DataHeap::getInstance().isValidIndex(p.getVolumeFlux()));
                assertion(DataHeap::getInstance().isValidIndex(p.getUpdate()));

                DataHeap::getInstance().deleteData(p.getSpaceTimePredictor());
                DataHeap::getInstance().deleteData(p.getSpaceTimeVolumeFlux());
                DataHeap::getInstance().deleteData(p.getPredictor());
                DataHeap::getInstance().deleteData(p.getVolumeFlux());
                DataHeap::getInstance().deleteData(p.getUpdate());
                DataHeap::getInstance().deleteData(p.getSolution());

                p.setSpaceTimePredictor(-1);
                p.setSpaceTimeVolumeFlux(-1);
                p.setPredictor(-1);
                p.setVolumeFlux(-1);
                p.setSolution(-1);
                p.setUpdate(-1);
                break;
              default:
                break;
            }
          }

          if (DataHeap::getInstance().isValidIndex(
                  p.getExtrapolatedPredictor())) {
            switch (p.getType()) {
              case exahype::records::ADERDGCellDescription::Erased:
              case exahype::records::ADERDGCellDescription::EmptyAncestor:
              case exahype::records::ADERDGCellDescription::EmptyDescendant:
                assertion(
                    DataHeap::getInstance().isValidIndex(p.getFluctuation()));

                DataHeap::getInstance().deleteData(
                    p.getExtrapolatedPredictor());
                DataHeap::getInstance().deleteData(p.getFluctuation());

                p.setExtrapolatedPredictor(-1);
                p.setFluctuation(-1);
                break;
              default:
                break;
            }
          }
        }
      }
      break;
    default:
      logDebug("cleanCellDescription(...)",
               "solver is not associated with any cell descriptions of this "
               "cell. cell="
                   << toString());
      break;
  }
}

exahype::Cell::SubcellPosition
exahype::Cell::computeSubcellPositionOfCellOrAncestor(
    const exahype::records::ADERDGCellDescription& pChild) const {
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                 _cellData.getCellDescriptionsIndex()),
             toString());
  assertion1(solvers::RegisteredSolvers[pChild.getSolverNumber()]->getType() ==
                 exahype::solvers::Solver::Type::ADER_DG,
             toString());
  assertion1(
      pChild.getType() == exahype::records::ADERDGCellDescription::Cell ||
          pChild.getType() == exahype::records::ADERDGCellDescription::Ancestor,
      toString());

  exahype::Cell::SubcellPosition subcellPosition;
  // Initialisation.
  subcellPosition.parentIndex = pChild.getParentIndex();
  exahype::records::ADERDGCellDescription* pParent = 0;
  for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(
           subcellPosition.parentIndex)) {  // Loop over cell descriptions
    if (p.getSolverNumber() == pChild.getSolverNumber()) {
      pParent = &p;
      assertion(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
          subcellPosition.parentIndex));
    }
  }

  if (pParent != 0) {
    // Iterative determining of the top most parent that might hold data.
    while (pParent->getType() ==
               exahype::records::ADERDGCellDescription::EmptyAncestor &&
           ADERDGCellDescriptionHeap::getInstance().isValidIndex(
               pParent->getParentIndex())) {
      const int currentParentIndex =
          pParent->getParentIndex();  // Value must be fixed. We update pParent
                                      // within the loop.

      for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(
               currentParentIndex)) {  // Loop over cell descriptions
        if (p.getSolverNumber() == pChild.getSolverNumber()) {
          subcellPosition.parentIndex = pParent->getParentIndex();
          pParent = &p;
        }
      }
    }
    assertion(pParent->getType() ==
                  exahype::records::ADERDGCellDescription::Ancestor ||
              exahype::records::ADERDGCellDescription::EmptyAncestor);

    // compute subcell index
    double scaling = tarch::la::aPowI(pChild.getLevel() - 1, 3);
    for (int xi = 0; xi < DIMENSIONS; xi++) {
      assertion((pChild.getOffset(xi) - pParent->getOffset(xi)) >= 0);
      subcellPosition.subcellIndex[xi] = tarch::la::round(
          scaling * (pChild.getOffset(xi) - pParent->getOffset(xi)));
    }
  }

  return subcellPosition;
}


void exahype::Cell::validateNoNansInADERDGSolver(
  int                                  number,
  exahype::Cell&                       fineGridCell,
  const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
  const std::string&                   methodTraceOfCaller
) {
  auto& p = fineGridCell.getADERDGCellDescription(number);

  assertion1(DataHeap::getInstance().isValidIndex(p.getSpaceTimePredictor()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getSpaceTimeVolumeFlux()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getSolution()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getUpdate()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getPredictor()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getVolumeFlux()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getExtrapolatedPredictor()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getFluctuation()),p.toString());

  assertionEquals4(p.getPredictorTimeStepSize(),p.getPredictorTimeStepSize(),
                   fineGridVerticesEnumerator.toString(),
                   p.toString(),fineGridCell.toString(),methodTraceOfCaller);

  exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
      exahype::solvers::RegisteredSolvers[p.getSolverNumber()]);
  #ifdef Dim2
  const int numberOfVolumeEntries = solver->getNodesPerCoordinateAxis() * solver->getNodesPerCoordinateAxis()
  #else
  const int numberOfVolumeEntries = solver->getNodesPerCoordinateAxis() * solver->getNodesPerCoordinateAxis() * solver->getNodesPerCoordinateAxis()
  #endif
  * (solver->getNumberOfParameters()+solver->getNumberOfVariables());

  //double* luh = DataHeap::getInstance().getData(p.getSolution()).data();
  for (int i=0; i<numberOfVolumeEntries; i++)
    assertionEquals5(luh[i], luh[i],
                   fineGridVerticesEnumerator.toString(),
                   p.toString(),fineGridCell.toString(),methodTraceOfCaller,i);

  //double* lQi = DataHeap::getInstance().getData(p.getSpaceTimePredictor()).data();
  //double* lFi = DataHeap::getInstance().getData(p.getSpaceTimeVolumeFlux()).data();
  //double* lQhi = DataHeap::getInstance().getData(p.getPredictor()).data();
  //double* lFhi = DataHeap::getInstance().getData(p.getVolumeFlux()).data();
  //double* lQhbnd = DataHeap::getInstance().getData(p.getExtrapolatedPredictor()).data();
  //double* lFhbnd = DataHeap::getInstance().getData(p.getFluctuation()).data();


  for (int i=0; i<numberOfVolumeEntries; i++) {
    assertionEquals11(lQi[i], lQi[i],
     fineGridVerticesEnumerator.toString(),
     p.toString(),fineGridCell.toString(),methodTraceOfCaller,i,
     (long int)(lQi), (long int)(lFi), (long int)(lQhi), (long int)(lFhi), (long int)(lQhbnd), (long int)(lFhbnd) );
  }

  assertionEquals4(lFi[0], lFi[0],
                   fineGridVerticesEnumerator.toString(),
                   p.toString(),fineGridCell.toString(),methodTraceOfCaller);

  assertionEquals4(lQhi[0], lQhi[0],
                   fineGridVerticesEnumerator.toString(),
                   p.toString(),fineGridCell.toString(),methodTraceOfCaller);

  assertionEquals4(lFhi[0], lFhi[0],
                   fineGridVerticesEnumerator.toString(),
                   p.toString(),fineGridCell.toString(),methodTraceOfCaller);

  assertionEquals4(lQhbnd[0], lQhbnd[0],
                   fineGridVerticesEnumerator.toString(),
                   p.toString(),fineGridCell.toString(),methodTraceOfCaller);

  assertionEquals4(lFhbnd[0], lFhbnd[0],
                   fineGridVerticesEnumerator.toString(),
                   p.toString(),fineGridCell.toString(),methodTraceOfCaller);
}


exahype::Cell::SubcellPosition
exahype::Cell::computeSubcellPositionOfDescendant(
    const exahype::records::ADERDGCellDescription& pChild) const {
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                 _cellData.getCellDescriptionsIndex()),
             toString());
  assertion1(solvers::RegisteredSolvers[pChild.getSolverNumber()]->getType() ==
                 exahype::solvers::Solver::Type::ADER_DG,
             toString());
  assertion1(
      pChild.getType() == exahype::records::ADERDGCellDescription::Descendant,
      toString());

  exahype::Cell::SubcellPosition subcellPosition;

  // Initialisation.
  assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                 pChild.getParentIndex()),
             pChild.getParentIndex());
  subcellPosition.parentIndex = pChild.getParentIndex();
  exahype::records::ADERDGCellDescription* pParent = 0;

  for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(
           pChild.getParentIndex())) {
    if (p.getSolverNumber() == pChild.getSolverNumber()) {
      pParent = &p;
    }
  }

  if (pParent != 0) {
    // recursion
    while (pParent->getType() ==
           exahype::records::ADERDGCellDescription::EmptyDescendant) {
      const int currentParentIndex =
          pParent->getParentIndex();  // Value must be fixed. We update pParent
                                      // within the loop.
      assertion1(ADERDGCellDescriptionHeap::getInstance().isValidIndex(
                     currentParentIndex),
                 currentParentIndex);
      for (auto& p : ADERDGCellDescriptionHeap::getInstance().getData(
               currentParentIndex)) {  // Loop over cell descriptions
        if (p.getSolverNumber() == pChild.getSolverNumber()) {
          subcellPosition.parentIndex = pParent->getParentIndex();
          pParent = &p;
        }
      }
    }

    assertion(pParent->getType() ==
                  exahype::records::ADERDGCellDescription::Descendant ||
              pParent->getType() ==
                  exahype::records::ADERDGCellDescription::Cell);

    // compute subcell index
    double scaling = tarch::la::aPowI(pChild.getLevel() - 1, 3);

    for (int xi = 0; xi < DIMENSIONS; xi++) {
      assertion((pChild.getOffset(xi) - pParent->getOffset(xi)) >= 0);
      subcellPosition.subcellIndex[xi] = tarch::la::round(
          scaling * (pChild.getOffset(xi) - pParent->getOffset(xi)));
    }
  } else {
    std::cerr << "exahype::Cell::computeSubcellPositionOfDescendant: parent of "
                 "descendant could not be found!"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  return subcellPosition;
}


int exahype::Cell::getNumberOfADERDGCellDescriptions() const {
  return ADERDGCellDescriptionHeap::getInstance().getData(
      getCellDescriptionsIndex()).size();
}


int exahype::Cell::getNumberOfFiniteVolumeCellDescriptions() const {
  return FiniteVolumesCellDescriptionHeap::getInstance().getData(
      getCellDescriptionsIndex()).size();
}


#ifdef Parallel
void exahype::Cell::clearLoadBalancingWorkloads() {
  if (isRefined()) {
    _cellData.setLocalWorkload(0.0);
    _cellData.setGlobalWorkload(0.0);
  }
  else {
    // @todo really insert here the number of real solvers and weight them accordingly
    _cellData.setLocalWorkload(1.0);
    _cellData.setGlobalWorkload(1.0);
  }
}


void exahype::Cell::restrictLoadBalancingWorkloads(const Cell& childCell, bool isRemote) {
  if (isRemote) {
    // _cellData.setLocalWorkload(  _cellData.getLocalWorkload()  + childCell._cellData.getLocalWorkload() );
    _cellData.setGlobalWorkload(
      std::max(_cellData.getLocalWorkload(), childCell._cellData.getGlobalWorkload())
    );
  }
  else {
    _cellData.setLocalWorkload(  _cellData.getLocalWorkload()  + childCell._cellData.getLocalWorkload() );
    _cellData.setGlobalWorkload( _cellData.getGlobalWorkload() + childCell._cellData.getGlobalWorkload() );
  }

  if ( ADERDGCellDescriptionHeap::getInstance().isValidIndex(getCellDescriptionsIndex()) ) {
    const double embeddedADERDGCells = getNumberOfADERDGCellDescriptions();
    const double embeddedFVPatches   = getNumberOfFiniteVolumeCellDescriptions();

    // @todo this will require further tuning and it might become necessary to
    //       take the order or patch size into account.
    double additionalWeight = 4.0 * embeddedADERDGCells + 1.0 * embeddedFVPatches;

    assertion(additionalWeight>=0.0);

    _cellData.setLocalWorkload(  _cellData.getLocalWorkload()  + additionalWeight );
    _cellData.setGlobalWorkload( _cellData.getGlobalWorkload() + additionalWeight );
  }
}


double exahype::Cell::getLocalWorkload() const {
  return _cellData.getLocalWorkload();
}


double exahype::Cell::getGlobalWorkload() const {
  return _cellData.getGlobalWorkload();
}
#endif


bool exahype::Cell::setSolutionMinMaxAndAnalyseValidity( double min, double max, int solverIndex ) {
  assertion( ADERDGCellDescriptionHeap::getInstance().isValidIndex( getCellDescriptionsIndex() ) ) ;
  assertion( ADERDGCellDescriptionHeap::getInstance().getData( getCellDescriptionsIndex() ).size()>solverIndex ) ;
  assertion( max>=min );

  exahype::records::ADERDGCellDescription& cellDescription = ADERDGCellDescriptionHeap::getInstance().getData(getCellDescriptionsIndex())[solverIndex];

  double cellMinimum = std::numeric_limits<double>::max();
  double cellMaximum = std::numeric_limits<double>::min();

  for (int faceNumber = 0; faceNumber<DIMENSIONS_TIMES_TWO; faceNumber++ ) {
    cellMinimum = std::min(cellMinimum,cellDescription.getSolutionMin(faceNumber));
    cellMaximum = std::max(cellMaximum,cellDescription.getSolutionMax(faceNumber));
  }

  bool isValidNewCombination = tarch::la::greater(min,cellMinimum)
                             & tarch::la::greater(cellMaximum,max);

  for (int faceNumber = 0; faceNumber<DIMENSIONS_TIMES_TWO; faceNumber++ ) {
    cellDescription.setSolutionMin(min);
    cellDescription.setSolutionMax(max);
  }

  return isValidNewCombination;
}


void exahype::Cell::mergeSolutionMinMaxOnFace(
  const int cellDescriptionsIndexOfLeftCell,
  const int cellDescriptionsIndexOfRightCell,
  const int faceIndexForLeftCell,
  const int faceIndexForRightCell
) {
  assertion( ADERDGCellDescriptionHeap::getInstance().isValidIndex( cellDescriptionsIndexOfLeftCell ) ) ;
  assertion( ADERDGCellDescriptionHeap::getInstance().isValidIndex( cellDescriptionsIndexOfRightCell ) ) ;

  for (auto& leftCellDescription: ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionsIndexOfLeftCell))
  for (auto& rightCellDescription: ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionsIndexOfRightCell)) {
    if (
      leftCellDescription.getType() == exahype::records::ADERDGCellDescription::Cell
      &&
      rightCellDescription.getType() == exahype::records::ADERDGCellDescription::Cell
    ) {
      double min = std::min( leftCellDescription.getSolutionMin(faceIndexForLeftCell), rightCellDescription.getSolutionMin(faceIndexForRightCell) );
      double max = std::max( leftCellDescription.getSolutionMax(faceIndexForLeftCell), rightCellDescription.getSolutionMax(faceIndexForRightCell) );

      leftCellDescription.setSolutionMin(faceIndexForLeftCell,min);
      rightCellDescription.setSolutionMin(faceIndexForRightCell,min);
      leftCellDescription.setSolutionMax(faceIndexForLeftCell,max);
      rightCellDescription.setSolutionMax(faceIndexForRightCell,max);
    }
    else {
      assertionMsg( false, "Dominic, please implement" );
    }
  }
}


void exahype::Cell::mergeSolutionMinMaxOnFace(
  const int cellDescriptionsIndex,
  double min, double max,
  int faceNumber,
  int ADERDGSolverNumber
) {
  assertion( ADERDGCellDescriptionHeap::getInstance().isValidIndex( cellDescriptionsIndex ) ) ;
  assertion( ADERDGCellDescriptionHeap::getInstance().getData( cellDescriptionsIndex ).size()>ADERDGSolverNumber ) ;
  assertion( max>=min );

  records::ADERDGCellDescription& cellDescription = ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionsIndex)[ADERDGSolverNumber];

  if (cellDescription.getType() == exahype::records::ADERDGCellDescription::Cell) {
    cellDescription.setSolutionMin( faceNumber, std::min( cellDescription.getSolutionMin(faceNumber),min ) );
    cellDescription.setSolutionMax( faceNumber, std::max( cellDescription.getSolutionMax(faceNumber),max ) );
  }
}
