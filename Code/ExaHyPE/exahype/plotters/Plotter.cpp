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
 
#include "exahype/plotters/Plotter.h"

#include "exahype/plotters/ADERDG2CartesianVTK.h"
#include "exahype/plotters/ADERDG2LegendreVTK.h"
#include "exahype/plotters/ADERDG2LegendreCSV.h"
#include "exahype/plotters/ADERDG2ProbeAscii.h"
#include "exahype/plotters/FiniteVolumes2VTKAscii.h"
#include "exahype/plotters/LimitingADERDG2CartesianVTK.h"
#include "exahype/solvers/LimitingADERDGSolver.h"


std::vector<exahype::plotters::Plotter*> exahype::plotters::RegisteredPlotters;

tarch::logging::Log exahype::plotters::Plotter::_log( "exahype::solvers::Plotter" );


exahype::plotters::Plotter::Plotter(
  int solver, int plotterCount,
  const exahype::Parser& parser, UserOnTheFlyPostProcessing* postProcessing)
    : _solver(solver),
      _identifier(parser.getIdentifierForPlotter(solver, plotterCount)),
      _writtenUnknowns(parser.getUnknownsForPlotter(solver, plotterCount)),
      _time(parser.getFirstSnapshotTimeForPlotter(solver, plotterCount)),
      _repeat(parser.getRepeatTimeForPlotter(solver, plotterCount)),
      _filename(parser.getFilenameForPlotter(solver, plotterCount)),
      _select(parser.getSelectorForPlotter(solver, plotterCount)),
      _isActive(false),
      _device(nullptr) {
  if (_time < 0.0) {
    logError("Plotter(...)",
      "plotter's first snapshot time is set to negative value "
      << _time << ". Plotter configuration=" << toString() );
  }
  if (_repeat < 0.0) {
    logError("Plotter(...)", "plotter's repeat time is set to negative value "
      << _repeat << ". Plotter configuration=" << toString() );
  }
  logInfo("Plotter(...)", "write snapshot to file "
    << _filename << " every " << _repeat
    << " time units with first snapshot at " << _time
    << ". plotter type is " << _identifier << ". Plotter configuration=" << toString() );

  assertion(_solver < static_cast<int>(solvers::RegisteredSolvers.size()));

  switch (solvers::RegisteredSolvers[_solver]->getType()) {
    case exahype::solvers::Solver::Type::ADER_DG:
      /**
       * This is actually some kind of switch expression though switches do
       * not work for strings, so we map it onto an if-then-else cascade.
       */
      if (_identifier.compare( ADERDG2CartesianVerticesVTKAscii::getIdentifier() ) == 0) {
        _device = new ADERDG2CartesianVerticesVTKAscii(postProcessing);
      }
      if (_identifier.compare( ADERDG2CartesianVerticesVTKBinary::getIdentifier() ) == 0) {
        _device = new ADERDG2CartesianVerticesVTKBinary(postProcessing);
      }
      if (_identifier.compare( ADERDG2CartesianCellsVTKAscii::getIdentifier() ) == 0) {
        _device = new ADERDG2CartesianCellsVTKAscii(postProcessing);
      }
      if (_identifier.compare( ADERDG2CartesianCellsVTKBinary::getIdentifier() ) == 0) {
        _device = new ADERDG2CartesianCellsVTKBinary(postProcessing);
      }
      if (_identifier.compare( ADERDG2LegendreVerticesVTKAscii::getIdentifier() ) == 0) {
        _device = new ADERDG2LegendreVerticesVTKAscii(postProcessing);
      }
      if (_identifier.compare( ADERDG2LegendreVerticesVTKBinary::getIdentifier() ) == 0) {
        _device = new ADERDG2LegendreVerticesVTKBinary(postProcessing);
      }
      if (_identifier.compare( ADERDG2LegendreCellsVTKAscii::getIdentifier() ) == 0) {
        _device = new ADERDG2LegendreCellsVTKAscii(postProcessing);
      }
      if (_identifier.compare( ADERDG2LegendreCellsVTKBinary::getIdentifier() ) == 0) {
        _device = new ADERDG2LegendreCellsVTKBinary(postProcessing);
      }
      if (_identifier.compare( ADERDG2ProbeAscii::getIdentifier() ) == 0) {
        _device = new ADERDG2ProbeAscii(postProcessing);
      }
      if (_identifier.compare( ADERDG2LegendreCSV::getIdentifier() ) == 0) {
        _device = new ADERDG2LegendreCSV(postProcessing);
      }
    break;
    case exahype::solvers::Solver::Type::FiniteVolumes:
      /**
       * This is actually some kind of switch expression though switches do
       * not work for strings, so we map it onto an if-then-else cascade.
       */
      if (_identifier.compare( FiniteVolumes2VTKAscii::getIdentifier() ) == 0) {
        _device = new FiniteVolumes2VTKAscii(
            postProcessing,static_cast<exahype::solvers::FiniteVolumesSolver*>(
                solvers::RegisteredSolvers[_solver])->getGhostLayerWidth());
      }
/*
      else if (_identifier.compare( ADERDG2VTKBinary::getIdentifier() ) == 0) {
        _device = new ADERDG2VTKBinary();
      }
      else if (_identifier.compare( ADERDG2ProbeAscii::getIdentifier() ) == 0) {
        _device = new ADERDG2ProbeAscii();
      }
*/
    break;
    case exahype::solvers::Solver::Type::LimitingADERDG:
      /**
       * Plotters for amy ADER-DG scheme.
       *
       * This is actually some kind of switch expression though switches do
       * not work for strings, so we map it onto an if-then-else cascade.
       */
      if (_identifier.compare( ADERDG2CartesianVerticesVTKAscii::getIdentifier() ) == 0) {
        _device = new ADERDG2CartesianVerticesVTKAscii(postProcessing);
      }
      if (_identifier.compare( ADERDG2CartesianVerticesVTKBinary::getIdentifier() ) == 0) {
        _device = new ADERDG2CartesianVerticesVTKBinary(postProcessing);
      }
      if (_identifier.compare( ADERDG2CartesianCellsVTKAscii::getIdentifier() ) == 0) {
        _device = new ADERDG2CartesianCellsVTKAscii(postProcessing);
      }
      if (_identifier.compare( ADERDG2CartesianCellsVTKBinary::getIdentifier() ) == 0) {
        _device = new ADERDG2CartesianCellsVTKBinary(postProcessing);
      }
      if (_identifier.compare( ADERDG2LegendreVerticesVTKAscii::getIdentifier() ) == 0) {
        _device = new ADERDG2LegendreVerticesVTKAscii(postProcessing);
      }
      if (_identifier.compare( ADERDG2LegendreVerticesVTKBinary::getIdentifier() ) == 0) {
        _device = new ADERDG2LegendreVerticesVTKBinary(postProcessing);
      }
      if (_identifier.compare( ADERDG2LegendreCellsVTKAscii::getIdentifier() ) == 0) {
        _device = new ADERDG2LegendreCellsVTKAscii(postProcessing);
      }
      if (_identifier.compare( ADERDG2LegendreCellsVTKBinary::getIdentifier() ) == 0) {
        _device = new ADERDG2LegendreCellsVTKBinary(postProcessing);
      }
      if (_identifier.compare( ADERDG2ProbeAscii::getIdentifier() ) == 0) {
        _device = new ADERDG2ProbeAscii(postProcessing);
      }
      if (_identifier.compare( ADERDG2LegendreCSV::getIdentifier() ) == 0) {
        _device = new ADERDG2LegendreCSV(postProcessing);
      }

      /**
       * Plotters specifically for the limiting ADER-DG scheme.
       *
       * This is actually some kind of switch expression though switches do
       * not work for strings, so we map it onto an if-then-else cascade.
       */
      if (_identifier.compare( LimitingADERDG2CartesianVerticesVTKAscii::getIdentifier() ) == 0) {
//        _device = new LimitingADERDG2CartesianVerticesVTKAscii(
//            postProcessing,static_cast<exahype::solvers::LimitingADERDGSolver*>(
//                solvers::RegisteredSolvers[_solver])->_limiter->getGhostLayerWidth());
      }
      if (_identifier.compare( LimitingADERDG2CartesianVerticesVTKBinary::getIdentifier() ) == 0) {
//        _device = new LimitingADERDG2CartesianVerticesVTKBinary(
//            postProcessing,static_cast<exahype::solvers::LimitingADERDGSolver*>(
//            solvers::RegisteredSolvers[_solver])->_limiter->getGhostLayerWidth());
      }
      if (_identifier.compare( LimitingADERDG2CartesianCellsVTKAscii::getIdentifier() ) == 0) {
//        _device = new LimitingADERDG2CartesianCellsVTKAscii(
//            postProcessing,static_cast<exahype::solvers::LimitingADERDGSolver*>(
//                solvers::RegisteredSolvers[_solver])->_limiter->getGhostLayerWidth());
      }
      if (_identifier.compare( LimitingADERDG2CartesianCellsVTKBinary::getIdentifier() ) == 0) {
//        _device = new LimitingADERDG2CartesianCellsVTKBinary(
//            postProcessing,static_cast<exahype::solvers::LimitingADERDGSolver*>(
//                solvers::RegisteredSolvers[_solver])->_limiter->getGhostLayerWidth());
      }
//      if (_identifier.compare( ADERDG2LegendreVerticesVTKAscii::getIdentifier() ) == 0) {
//        _device = new LimitingADERDG2LegendreVerticesVTKAscii(
//            postProcessing,static_cast<exahype::solvers::LimitingADERDGSolver*>(
//                solvers::RegisteredSolvers[_solver])->_limiter->getGhostLayerWidth());
//      }
//      if (_identifier.compare( ADERDG2LegendreVerticesVTKBinary::getIdentifier() ) == 0) {
//        _device = new LimitingADERDG2LegendreVerticesVTKBinary(
//            postProcessing,static_cast<exahype::solvers::LimitingADERDGSolver*>(
//                solvers::RegisteredSolvers[_solver])->_limiter->getGhostLayerWidth());
//      }
//      if (_identifier.compare( ADERDG2LegendreCellsVTKAscii::getIdentifier() ) == 0) {
//        _device = new LimitingADERDG2LegendreCellsVTKAscii(
//            postProcessing,static_cast<exahype::solvers::LimitingADERDGSolver*>(
//                solvers::RegisteredSolvers[_solver])->_limiter->getGhostLayerWidth());
//      }
//      if (_identifier.compare( ADERDG2LegendreCellsVTKBinary::getIdentifier() ) == 0) {
//        _device = new LimitingADERDG2LegendreCellsVTKBinary(
//            postProcessing,static_cast<exahype::solvers::LimitingADERDGSolver*>(
//                solvers::RegisteredSolvers[_solver])->_limiter->getGhostLayerWidth());
//      if (_identifier.compare( ADERDG2ProbeAscii::getIdentifier() ) == 0) {
//        _device = new LimitingADERDG2ProbeAscii(
//            postProcessing,static_cast<exahype::solvers::LimitingADERDGSolver*>(
//            solvers::RegisteredSolvers[_solver])->_limiter->getGhostLayerWidth());
//      }
//      if (_identifier.compare( ADERDG2LegendreCSV::getIdentifier() ) == 0) {
//        _device = new LimitingADERDG2LegendreCSV(postProcessing);
//      }
    break;
  }


  if (_device!=nullptr) {
    _device->init(
        _filename,
        solvers::RegisteredSolvers[_solver]->getNodesPerCoordinateAxis(),
        solvers::RegisteredSolvers[_solver]->getNumberOfVariables(),
        _writtenUnknowns,
        _select
    );
  }
  else if (_identifier=="notoken") {
    logError(
      "Plotter(...)",
      "unable to set up " << (plotterCount+1) << "th plotter for the "
      << (_solver+1) << "th solver. Ensure number of plot sections "
      << "equals number of plotters originally passed to toolkit and "
      << "validate that plot syntax is correct"
    );
  }
  else {
    logError(
      "Plotter(...)",
      "unknown plotter type "
          << _identifier << " for "
          << solvers::RegisteredSolvers[_solver]->getIdentifier()
	  << ". Potential reasons: you have not specified a valid identifier following the plot keyword or you have specified a plotter in the ExaHyPE toolkit and later removed this plotter from the config"
    );
  }
}


std::string exahype::plotters::Plotter::toString() const {
  std::ostringstream msg;

  msg << "(solver no=" << _solver
      << ",plotter identifier (type)=" << _identifier
      << ",written unknowns=" << _writtenUnknowns
      << ",time=" << _time
      << ",repeat=" << _repeat
      << ",file name=" << _filename
      << ",select statement=" << _select
      << ",device configured=" << (_device!=nullptr)
      << ")";

  return msg.str();
}


exahype::plotters::Plotter::~Plotter() {
  if (_device!=nullptr) {
    delete _device;
    _device = nullptr;
  }
}


double exahype::plotters::Plotter::getNextPlotTime() const {
  return _time;
}


bool exahype::plotters::Plotter::checkWetherSolverBecomesActive(double currentTimeStamp) {
  if ((_time >= 0.0) && tarch::la::greaterEquals(currentTimeStamp, _time)) {
    
    if (_device==nullptr){
      logError(
        "checkWetherSolverBecomesActive(double)",
        "unknown plotter type " << _identifier << " piping into file " << _filename
      );
    }
    else {
      assertion(_device!=nullptr);
      _isActive = true;
      _device->startPlotting(currentTimeStamp);
    }
  }

  return isActive();
}


bool exahype::plotters::Plotter::isActive() const {
  return _isActive;
}


bool exahype::plotters::Plotter::plotDataFromSolver(int solver) const {
  return isActive() && _solver == solver;
}


void exahype::plotters::Plotter::plotPatch(
  const int cellDescriptionsIndex,
  const int element) {
  assertion(_device != nullptr);
  if (_device!=nullptr) {
    _device->plotPatch(cellDescriptionsIndex,element);
  }
}


void exahype::plotters::Plotter::finishedPlotting() {
  assertion(isActive());
  if (_repeat > 0.0) {
    _time += _repeat;
  } else {
    _time = -1.0;
  }
  if (_device!=nullptr) {
    _device->finishPlotting();
  }
  _isActive = false;
}


bool exahype::plotters::isAPlotterActive(double currentTimeStep) {
  bool result = false;
  for (const auto& p : RegisteredPlotters) {
    result |= p->checkWetherSolverBecomesActive(currentTimeStep);
  }
  return result;
}


double exahype::plotters::getTimeOfNextPlot() {
  double result = std::numeric_limits<double>::max();
  for (const auto& p : RegisteredPlotters) {
    result = std::min(result,p->getNextPlotTime());
  }
  return result;
}


void exahype::plotters::finishedPlotting() {
  for (auto& p : RegisteredPlotters) {
    if (p->isActive()) {
      p->finishedPlotting();
    }
  }
}


std::string exahype::plotters::Plotter::getFileName() const {
  return _filename;
}
