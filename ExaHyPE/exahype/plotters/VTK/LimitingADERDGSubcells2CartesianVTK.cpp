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

#include "LimitingADERDGSubcells2CartesianVTK.h"
#include "tarch/parallel/Node.h"

// @todo 16/05/03:Dominic Etienne Charreir Plotter depends now on kernels.
// Should thus be placed in kernel module or the solver
// should provide a function that computes solution values
// at equidistant grid points
#include "kernels/DGMatrices.h"
#include "peano/utils/Loop.h"


#include "tarch/plotter/griddata/unstructured/vtk/VTKTextFileWriter.h"
#include "tarch/plotter/griddata/unstructured/vtk/VTKBinaryFileWriter.h"


#include "kernels/DGBasisFunctions.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

std::string exahype::plotters::LimitingADERDGSubcells2CartesianCellsVTKAscii::getIdentifier() {
  return "vtk::Cartesian::subcells::limited::ascii";
}


exahype::plotters::LimitingADERDGSubcells2CartesianCellsVTKAscii::LimitingADERDGSubcells2CartesianCellsVTKAscii(
    exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
    const int ghostLayerWidth):
    LimitingADERDGSubcells2CartesianVTK(postProcessing,ghostLayerWidth,false) {
}


std::string exahype::plotters::LimitingADERDGSubcells2CartesianCellsVTKBinary::getIdentifier() {
 return "vtk::Cartesian::subcells::limited::binary";
}


exahype::plotters::LimitingADERDGSubcells2CartesianCellsVTKBinary::LimitingADERDGSubcells2CartesianCellsVTKBinary(
    exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
    const int ghostLayerWidth):
    LimitingADERDGSubcells2CartesianVTK(postProcessing,ghostLayerWidth,true) {
}

tarch::logging::Log exahype::plotters::LimitingADERDGSubcells2CartesianVTK::_log("exahype::plotters::LimitingADERDGSubcells2CartesianVTK");

exahype::plotters::LimitingADERDGSubcells2CartesianVTK::LimitingADERDGSubcells2CartesianVTK(
    exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
    const int ghostLayerWidth,
    const bool isBinary)
  :
  Device(postProcessing),
  _fileCounter(-1),
  _isBinary(isBinary),
  _order(-1),
  _solverUnknowns(-1),
  _writtenUnknowns(-1),
  _ghostLayerWidth(ghostLayerWidth),
  _gridWriter(nullptr),
  _patchWriter(nullptr),
  _vertexDataWriter(nullptr),
  _cellDataWriter(nullptr),
  _timeStampCellDataWriter(nullptr),
  _cellRefinementStatusWriter(nullptr),
  _cellPreviousRefinementStatusWriter(nullptr)
{}

void exahype::plotters::LimitingADERDGSubcells2CartesianVTK::init(
  const std::string& filename,
  int                orderPlusOne,
  int                unknowns,
  int                writtenUnknowns,
  exahype::parser::ParserView plotterParameters
) {
  _filename          = filename;
  _order             = orderPlusOne-1;
  _solverUnknowns    = unknowns;
  _plotterParameters            = plotterParameters;
  _patchWriter       = nullptr;
  _writtenUnknowns   = writtenUnknowns;

  _slicer = Slicer::bestFromSelectionQuery(plotterParameters);

  if(_slicer) {
    logInfo("init", "Plotting selection "<<_slicer->toString()<<" to Files "<<filename);
  }
}


void exahype::plotters::LimitingADERDGSubcells2CartesianVTK::startPlotting( double time ) {
  _fileCounter++;

  assertion( _patchWriter==nullptr );

  if (_writtenUnknowns>0) {
    if (_isBinary) {
      _patchWriter =
        new tarch::plotter::griddata::blockstructured::PatchWriterUnstructured(
          new tarch::plotter::griddata::unstructured::vtk::VTKBinaryFileWriter());
    }
    else {
      _patchWriter =
        new tarch::plotter::griddata::blockstructured::PatchWriterUnstructured(
          new tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter());
    }

    _gridWriter                = _patchWriter->createSinglePatchWriter();

    _cellDataWriter            = _patchWriter->createCellDataWriter("Q", _writtenUnknowns);
    _vertexDataWriter          = nullptr;

    _cellRefinementStatusWriter   = _patchWriter->createCellDataWriter("RefinementStatus", 1);
    _cellPreviousRefinementStatusWriter   = _patchWriter->createCellDataWriter("PreviousRefinementStatus", 1);

    _timeStampCellDataWriter   = _patchWriter->createCellDataWriter("time", 1);

    assertion( _patchWriter!=nullptr );
    assertion( _gridWriter!=nullptr );
    assertion( _timeStampCellDataWriter!=nullptr );
  }

  _postProcessing->startPlotting( time );
}


void exahype::plotters::LimitingADERDGSubcells2CartesianVTK::finishPlotting() {
  _postProcessing->finishPlotting();

  if (_writtenUnknowns>0) {
    assertion( _patchWriter!=nullptr );
    assertion( _gridWriter!=nullptr );
    assertion( _timeStampCellDataWriter!=nullptr );

    _gridWriter->close();
//    if (_timeStampCellDataWriter!=nullptr) _timeStampCellDataWriter->close();
    if (_vertexDataWriter!=nullptr)                  _vertexDataWriter->close();
    if (_cellDataWriter!=nullptr)                    _cellDataWriter->close();
    if (_cellRefinementStatusWriter!=nullptr)           _cellRefinementStatusWriter->close();
    if (_cellPreviousRefinementStatusWriter!=nullptr)   _cellPreviousRefinementStatusWriter->close();
    _timeStampCellDataWriter->close();

    std::ostringstream snapshotFileName;
    snapshotFileName << _filename
                     << "-" << _fileCounter;

    const bool hasBeenSuccessful =
      _patchWriter->writeToFile(snapshotFileName.str());
    if (!hasBeenSuccessful) {
      exit(-1);
    }
  }

  if (_vertexDataWriter!=nullptr)          delete _vertexDataWriter;
  if (_cellDataWriter!=nullptr)            delete _cellDataWriter;
  if (_timeStampCellDataWriter!=nullptr)   delete _timeStampCellDataWriter;
  if (_cellRefinementStatusWriter!=nullptr)   delete _cellRefinementStatusWriter;
  if (_cellPreviousRefinementStatusWriter!=nullptr)   delete _cellPreviousRefinementStatusWriter;
  if (_gridWriter!=nullptr)                delete _gridWriter;
  if (_patchWriter!=nullptr)               delete _patchWriter;

  _vertexDataWriter                  = nullptr;
  _cellDataWriter                    = nullptr;
  _patchWriter                       = nullptr;
  _timeStampCellDataWriter           = nullptr;
  _cellRefinementStatusWriter           = nullptr;
  _cellPreviousRefinementStatusWriter   = nullptr;
  _gridWriter                = nullptr;
}



exahype::plotters::LimitingADERDGSubcells2CartesianVTK::~LimitingADERDGSubcells2CartesianVTK() {
}

void exahype::plotters::LimitingADERDGSubcells2CartesianVTK::plotPatch(const int solverNumber,const solvers::Solver::CellInfo& cellInfo) {
  const int element = solvers::Solver::indexOfCellDescription(cellInfo._ADERDGCellDescriptions,solverNumber);
  auto& solverPatch  = cellInfo._ADERDGCellDescriptions[element];

  if (solverPatch.getType()==exahype::solvers::ADERDGSolver::CellDescription::Type::Cell) {
    assertion(exahype::solvers::RegisteredSolvers[solverPatch.getSolverNumber()]->getType()==
        exahype::solvers::Solver::Type::LimitingADERDG);
    auto* limitingADERDG =
        static_cast<exahype::solvers::LimitingADERDGSolver*>(
            exahype::solvers::RegisteredSolvers[solverPatch.getSolverNumber()]);

    // ignore limiter status on coarser mesh levels
    int refinementStatus         = solverPatch.getRefinementStatus();
    int previousRefinementStatus = solverPatch.getPreviousRefinementStatus();
    assertion(static_cast<unsigned int>(solverPatch.getSolverNumber())
        <exahype::solvers::RegisteredSolvers.size());
    if (solverPatch.getLevel()<limitingADERDG->getMaximumAdaptiveMeshLevel()) {
      refinementStatus         = 0;
      previousRefinementStatus = 0;
    }

    if (refinementStatus>=limitingADERDG->getSolver()->getMinimumRefinementStatusForActiveFVPatch()) {
      auto& limiterPatch = limitingADERDG->
              getLimiterPatchForSolverPatch(solverPatch,cellDescriptionsIndex);

      double* limiterSolution = DataHeap::getInstance().getData(limiterPatch.getSolution()).data();
      plotFiniteVolumesPatch(
          limiterPatch.getOffset(),
          limiterPatch.getSize(), limiterSolution,
          limiterPatch.getTimeStamp(),
          refinementStatus,
          previousRefinementStatus);
    }
  }
}

void exahype::plotters::LimitingADERDGSubcells2CartesianVTK::plotFiniteVolumesPatch(
  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
  const double* const u,
  const double timeStamp,
  const int RefinementStatus, const int previousRefinementStatus) {
  if (!_slicer || _slicer->isPatchActive(offsetOfPatch, sizeOfPatch)) {
    logDebug("plotPatch(...)","offset of patch: "<<offsetOfPatch
    <<", size of patch: "<<sizeOfPatch
    <<", time stamp: "<<timeStamp);

    assertion( _writtenUnknowns==0 || _patchWriter!=nullptr );
    assertion( _writtenUnknowns==0 || _gridWriter!=nullptr );
    assertion( _writtenUnknowns==0 || _timeStampCellDataWriter!=nullptr );

    const int numberOfCellsPerAxis = 2*_order+1;

    int cellIndex = _writtenUnknowns==0 ? -1 : _gridWriter->plotPatch(offsetOfPatch, sizeOfPatch, numberOfCellsPerAxis).second;

    double* sourceValue = new double[_solverUnknowns];
    double* value       = _writtenUnknowns==0 ? nullptr : new double[_writtenUnknowns];

    dfor(i,numberOfCellsPerAxis+_ghostLayerWidth) {
      if (tarch::la::allSmaller(i,numberOfCellsPerAxis+_ghostLayerWidth)
          && tarch::la::allGreater(i,_ghostLayerWidth-1)) {
        if (_writtenUnknowns>0) {
          _timeStampCellDataWriter->plotCell(cellIndex, timeStamp);
          _cellRefinementStatusWriter->plotCell(cellIndex, RefinementStatus);
          _cellPreviousRefinementStatusWriter->plotCell(cellIndex, previousRefinementStatus);
        }

        for (int unknown=0; unknown < _solverUnknowns; unknown++) {
          sourceValue[unknown] =
            u[peano::utils::dLinearisedWithoutLookup(i,numberOfCellsPerAxis+2*_ghostLayerWidth)*_solverUnknowns+unknown];
        } // !!! Be aware of the "2*_ghostLayerWidth" !!!

        assertion(sizeOfPatch(0)==sizeOfPatch(1));

        _postProcessing->mapQuantities(
          offsetOfPatch,
          sizeOfPatch,
          offsetOfPatch + (i-_ghostLayerWidth).convertScalar<double>()* (sizeOfPatch(0)/(numberOfCellsPerAxis)),
          i-_ghostLayerWidth,
          sourceValue,
          value,
          timeStamp
        );

        if (_writtenUnknowns>0) {
          _cellDataWriter->plotCell(cellIndex, value, _writtenUnknowns);
        }
        cellIndex++;
      }
    }

    delete[] sourceValue;
    delete[] value;
  }
}
