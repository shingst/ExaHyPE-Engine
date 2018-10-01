// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
// ========================
//   www.exahype.eu
// ========================
#ifndef POSTPROCESSING_Plotter_CLASS_HEADER_
#define POSTPROCESSING_Plotter_CLASS_HEADER_

#include "exahype/plotters/Plotter.h"

namespace NavierStokes {
  class NavierStokesSolverDG;
  class Plotter;
}

class NavierStokes::Plotter : public exahype::plotters::Plotter::UserOnTheFlyPostProcessing {
public:
  Plotter(NavierStokes::NavierStokesSolverDG& solver);
  virtual ~Plotter();

  void startPlotting(double time) override;
  void finishPlotting() override;
  void mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp) override;
private:
    int order;
    NavierStokesSolverDG* solver;
};

#endif /* POSTPROCESSING_Plotter_CLASS_HEADER_ */