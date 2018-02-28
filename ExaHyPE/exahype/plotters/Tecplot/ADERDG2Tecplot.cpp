#include "exahype/plotters/Tecplot/ADERDG2Tecplot.h"


#include "tarch/logging/Log.h"


#include <cstdlib>
#include <stdio.h>
#include <sstream>
#include <memory>
#include <limits> // signaling_NaN


std::string exahype::plotters::ADERDG2Tecplot::getIdentifier() {
	return std::string("Tecplot::Binary");
}


tarch::logging::Log exahype::plotters::ADERDG2Tecplot::_log("exahype::plotters::ADERDG2Tecplot");


#ifndef TECPLOT
/*************************************************************************************************
 * ADERDG2Tecplot Dummy implementation in case TECPLOT support is skipped.
 * Probably such a section is (except the constructor) not neccessary as the methods are never
 * referenced/called.
 *************************************************************************************************/

exahype::plotters::ADERDG2Tecplot::ADERDG2Tecplot(
  exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing, exahype::solvers::Solver::Type type) : Device(postProcessing) {
	logError("ADERDG2Tecplot()", "ERROR: Compile with TECPLOT, otherwise you cannot use the Tecplot plotter.");
	// if(std::getenv("EXAHYPE_STRICT"))
	abort();
}

// all other methods are stubs
exahype::plotters::ADERDG2Tecplot::~ADERDG2Tecplot() {}

void exahype::plotters::ADERDG2Tecplot::init(const std::string& filename, int orderPlusOne, int solverUnknowns, int writtenUnknowns, const std::string& select) {
	logError("init()","ERROR: Compile with TECPLOT, otherwise you cannot use the Tecplot plotter. There will be no output going to " << filename << " today.");
	logError("init()", "Will fail gracefully. If you want to stop the program in such a case, please set the environment variable EXAHYPE_STRICT=\"Yes\".");
}
void exahype::plotters::ADERDG2Tecplot::plotPatch(const int cellDescriptionsIndex, const int element) {}
void exahype::plotters::ADERDG2Tecplot::startPlotting(double time) {
	logError("startPlotting()", "Skipping HDF5 output due to missing support.");
}
void exahype::plotters::ADERDG2Tecplot::finishPlotting() {}

#else

// Define the Fortran routines here.




#include "exahype/solvers/ADERDGSolver.h"

/*************************************************************************************************
 * ADERDG2Tecplot non-dummy implementation
 *************************************************************************************************/

exahype::plotters::ADERDG2Tecplot::ADERDG2Tecplot(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing, exahype::solvers::Solver::Type type) :
    Device(postProcessing), _solverType(type) {}

exahype::plotters::ADERDG2Tecplot::~ADERDG2Tecplot() {}

void exahype::plotters::ADERDG2Tecplot::init(const std::string& filename, int basisSize, int solverUnknowns, int writtenUnknowns, const std::string& select) {
	
	// Here you can register the filename, DG order, etc.
	
}

void exahype::plotters::ADERDG2Tecplot::plotPatch(
        const int cellDescriptionsIndex,
        const int element) {
  auto& aderdgCellDescription = exahype::solvers::ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

  // you can also use _solverType to see whether you have an ordinary ADERDG solver
  // or a limiting solver
  
  // this if allows you to understand whether you have an ADERDG cell or
  // a limiting cell or whatever.
  if (aderdgCellDescription.getType()==exahype::solvers::ADERDGSolver::CellDescription::Type::Cell) {
    // The internal cell structure is (order,order,order,nVar) in C and 3D.
    double* solverSolution = DataHeap::getInstance().getData(aderdgCellDescription.getSolution()).data();

    // vectors of length DIMENSIONS:
    double* cellOffset = aderdgCellDescription.getOffset().data();
    double* cellSize = aderdgCellDescription.getSize().data();
    
    double time = aderdgCellDescription.getCorrectorTimeStamp();
    
    int order     = _orderPlusOne - 1;
    // compute the dx as vector for instance with:
    // vector  dx    = 1./order * sizeOfPatch;

  } // if cellldescription == ADERDG cell
}

void exahype::plotters::ADERDG2Tecplot::startPlotting(double time) {
	_postProcessing->startPlotting(time);
	
	// This is called when a new plotting grid swipe happens
}

void exahype::plotters::ADERDG2Tecplot::finishPlotting() {
	_postProcessing->finishPlotting();
	
	// This is called when a plotting grid swipe ends.

}

#endif /* TECPLOT */
