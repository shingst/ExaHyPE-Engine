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
 *
 * @authors: Sven Koeppel
 **/

#include "exahype/plotters/CarpetHDF5/FiniteVolume2CarpetHDF5.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

#include <cstdlib>
#include <stdio.h>
#include <sstream>
#include <memory>

typedef tarch::la::Vector<DIMENSIONS,int> ivec;


std::string exahype::plotters::FiniteVolume2CarpetHDF5::getIdentifier() {
	return std::string("Carpet::Cartesian::Vertices::HDF5");
}


// my small C++11 to_string-independent workaround.
template <typename T> std::string toString( T Number ) {
	std::ostringstream ss; ss << Number; return ss.str();
}

typedef tarch::la::Vector<DIMENSIONS, double> dvec;

#ifndef HDF5
/*************************************************************************************************
 * FiniteVolume2CarpetHDF5 Dummy implementation in case HDF5 support is skipped.
 *************************************************************************************************/

exahype::plotters::FiniteVolume2CarpetHDF5::FiniteVolume2CarpetHDF5(
  exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
  const int _ghostLayerWidth) : Device(postProcessing), ghostLayerWidth(_ghostLayerWidth) {
	printf("ERROR: Compile with HDF5, otherwise you cannot use the HDF5 plotter.\n");
	abort();
}

// all other methods are stubs
exahype::plotters::FiniteVolume2CarpetHDF5::~FiniteVolume2CarpetHDF5() {}
void exahype::plotters::FiniteVolume2CarpetHDF5::init(const std::string& filename, int orderPlusOne, int solverUnknowns, int writtenUnknowns, exahype::parser::ParserView plotterParameters) {}
void exahype::plotters::FiniteVolume2CarpetHDF5::plotPatch(const int solverNumber,solvers::Solver::CellInfo& cellInfo) {}
void exahype::plotters::FiniteVolume2CarpetHDF5::plotPatch(const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,double timeStamp) {}
void exahype::plotters::FiniteVolume2CarpetHDF5::startPlotting(double time) {}
void exahype::plotters::FiniteVolume2CarpetHDF5::finishPlotting() {}

#else
#include "exahype/plotters/CarpetHDF5/CarpetHDF5Writer.h"
#include "kernels/KernelUtils.h" // indexing
#include "peano/utils/Loop.h" // dfor
#include "kernels/DGMatrices.h"
#include "exahype/solvers/ADERDGSolver.h"
#include "kernels/DGBasisFunctions.h"
#include "tarch/logging/Log.h"
#include <sstream>

/*************************************************************************************************
 * FiniteVolume2CarpetHDF5 non-dummy implementation
 *************************************************************************************************/

exahype::plotters::FiniteVolume2CarpetHDF5::FiniteVolume2CarpetHDF5(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
  const int _ghostLayerWidth) :  Device(postProcessing), ghostLayerWidth(_ghostLayerWidth) { writer = nullptr; }

exahype::plotters::FiniteVolume2CarpetHDF5::~FiniteVolume2CarpetHDF5() {
	if(writer) delete writer;
}

void exahype::plotters::FiniteVolume2CarpetHDF5::init(const std::string& filename, int _numberOfCellsPerAxis, int _solverUnknowns, int writtenUnknowns, exahype::parser::ParserView  plotterParameters) {
	// Determine names of output fields
	char **writtenQuantitiesNames = new char*[writtenUnknowns];
	std::fill_n(writtenQuantitiesNames, writtenUnknowns, nullptr);
	_postProcessing->writtenQuantitiesNames(writtenQuantitiesNames);
	
	numberOfCellsPerAxis = _numberOfCellsPerAxis;
	numberOfVerticesPerAxis = _numberOfCellsPerAxis + 1;
	solverUnknowns = _solverUnknowns;

	const int basisSize = numberOfVerticesPerAxis;
	writer = new exahype::plotters::CarpetHDF5Writer(filename, basisSize, solverUnknowns, writtenUnknowns, plotterParameters,
		writtenQuantitiesNames);	

	if(writer->slicer && writer->slicer->getIdentifier() == "CartesianSlicer") {
		throw std::domain_error("FiniteVolume2CarpetHDF5 currently does not support the CartesianSlicer.");
	}
}

void exahype::plotters::FiniteVolume2CarpetHDF5::plotPatch(const int solverNumber,solvers::Solver::CellInfo& cellInfo) {
  const int element = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);
  auto& cellDescription  = cellInfo._FiniteVolumesCellDescriptions[element];
	if (cellDescription.getType()==exahype::solvers::FiniteVolumesSolver::CellDescription::Type::Cell) {
		double* solution = static_cast<double*>(cellDescription.getSolution());

		plotPatch(
			cellDescription.getOffset(),
			cellDescription.getSize(), solution,
			cellDescription.getTimeStamp());
	}
}

void exahype::plotters::FiniteVolume2CarpetHDF5::plotPatch(
  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
  double* u, /* unknown */
  double timeStamp) {

	if(writer->slicer && !writer->slicer->isPatchActive(offsetOfPatch, sizeOfPatch)) return;
	
	double* mappedCell;

	if(writer->slicer && writer->slicer->getIdentifier() == "CartesianSlicer") {
		mappedCell = new double[writer->writtenFieldsSize];
		
		interpolateCartesianSlicedVertexPatch(offsetOfPatch, sizeOfPatch, u, mappedCell, timeStamp,
			static_cast<exahype::plotters::CartesianSlicer&>(*writer->slicer));
	} else {
		mappedCell  = new double[writer->patchFieldsSize];
		interpolateVertexPatch(offsetOfPatch, sizeOfPatch, u, mappedCell, timeStamp);
	}
	delete[] mappedCell;
}

void exahype::plotters::FiniteVolume2CarpetHDF5::interpolateVertexPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    double* u, /* ingoing unknowns in cell, size numberOfCellsPerAxis^DIMENSIONS */
    double *mappedCell, /* outgoing mapped cell, size writtenUnknowns^DIMENSIONS */
    double timeStamp) {

	double* vertexValue = new double[solverUnknowns];

	// the following assumes quadratic cells.
	assertion(sizeOfPatch(0)==sizeOfPatch(1));
	kernels::dindex patchPos(numberOfCellsPerAxis + 2*ghostLayerWidth); // including ghost zones

	dfor(ivertex, numberOfVerticesPerAxis) {
		// We do no smearing, so we only take into account the 2 nearest neighbours.
		constexpr int neighbourCellsPerAxis = 2;
		//constexpr int neighbourCellsMax = std::pow(neighbourCellsPerAxis, DIMENSIONS); // maximum possible cells (ie. 4 in 2D)
		std::fill_n(vertexValue, solverUnknowns, 0.0);
		int neighbourCells = 0; // actual neighbour cells taken into account
		dfor(icells, neighbourCellsPerAxis) {
			ivec icell = ghostLayerWidth + ivertex + (icells - neighbourCellsPerAxis / 2);
			
			// if the target cell position in the patch is *not* in the ghost layers:
			if (tarch::la::allSmaller(icell,numberOfCellsPerAxis+ghostLayerWidth)
			&& tarch::la::allGreater(icell,ghostLayerWidth-1)) {
				double *cell = u + patchPos.rowMajor(icell)*solverUnknowns;
				for (int unknown=0; unknown < solverUnknowns; unknown++) {
					vertexValue[unknown] += cell[unknown];
				}
				neighbourCells++;
			}
		}

		// normalize value
		for (int unknown=0; unknown < solverUnknowns; unknown++) {
			vertexValue[unknown] = vertexValue[unknown] / neighbourCells;
		}
		
		/*
		// The following code could be used instead of the neighbour contributions as
		// above and was used for the start. Just one  cell.
		// This works and shows how badly it is if we rely on ghost zones.
		double *cell = u + patchPos.rowMajor(ghostLayerWidth + ivertex)*solverUnknowns;
		for (int unknown=0; unknown < solverUnknowns; unknown++) {
			vertexValue[unknown] = 42; // cell[unknown];
		}
		*/
		

		double *outputValue = mappedCell + (DIMENSIONS == 3 ?
			writer->writtenCellIdx->get(ivertex(2),ivertex(1),ivertex(0),0) :
			writer->writtenCellIdx->get(ivertex(1),ivertex(0),0));
		_postProcessing->mapQuantities(
		offsetOfPatch,
		sizeOfPatch,
		offsetOfPatch + ivertex.convertScalar<double>()* (sizeOfPatch(0)/(numberOfVerticesPerAxis)), // coordinate of vertex
		ivertex,
		vertexValue,
		outputValue,
		timeStamp
		);
	}

	delete[] vertexValue;
	
	dvec dx = 1./numberOfCellsPerAxis * sizeOfPatch;
	writer->plotPatch(offsetOfPatch, sizeOfPatch, dx, mappedCell, timeStamp);
}

void exahype::plotters::FiniteVolume2CarpetHDF5::interpolateCartesianSlicedVertexPatch(
  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
  double *u,
  double *mappedCell,
  double timeStamp,
  const exahype::plotters::CartesianSlicer& slicer) {
	
	
	
}

void exahype::plotters::FiniteVolume2CarpetHDF5::startPlotting(double time) {
	_postProcessing->startPlotting(time);
	writer->startPlotting(time);
}

void exahype::plotters::FiniteVolume2CarpetHDF5::finishPlotting() {
	_postProcessing->finishPlotting();
	writer->finishPlotting();
}

#endif /* HDF5 */
