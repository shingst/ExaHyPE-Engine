// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
// ========================
//   www.exahype.eu
// ========================
#ifndef POSTPROCESSING_ConservedWriter_CLASS_HEADER_
#define POSTPROCESSING_ConservedWriter_CLASS_HEADER_

#include "exahype/plotters/Plotter.h"

namespace MHD {
  class MHDSolver;
  class ConservedWriter;
}

class MHD::ConservedWriter : public exahype::plotters::Plotter::UserOnTheFlyPostProcessing {
public:
  ConservedWriter(MHD::MHDSolver& solver);
  virtual ~ConservedWriter();

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
};

#endif /* POSTPROCESSING_ConservedWriter_CLASS_HEADER_ */