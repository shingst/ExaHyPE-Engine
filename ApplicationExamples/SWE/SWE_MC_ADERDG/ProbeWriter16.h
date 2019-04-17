// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
// ========================
//   www.exahype.eu
// ========================
#ifndef POSTPROCESSING_ProbeWriter16_CLASS_HEADER_
#define POSTPROCESSING_ProbeWriter16_CLASS_HEADER_

#include "exahype/plotters/Plotter.h"

namespace SWE {
  class MySWESolver;
  class ProbeWriter16;
}

class SWE::ProbeWriter16 : public exahype::plotters::Plotter::UserOnTheFlyPostProcessing {
public:
  ProbeWriter16(SWE::MySWESolver& solver);
  virtual ~ProbeWriter16();

  void startPlotting(double time) override;
  void finishPlotting() override;
  void mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp) override;
};

#endif /* POSTPROCESSING_ProbeWriter16_CLASS_HEADER_ */