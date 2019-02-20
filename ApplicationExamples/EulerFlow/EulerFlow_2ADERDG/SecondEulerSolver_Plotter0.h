// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include "exahype/plotters/Plotter.h"
namespace Euler{
  class SecondEulerSolver_Plotter0;

  /**
   * Forward declaration
   */
  class SecondEulerSolver;
}




class Euler::SecondEulerSolver_Plotter0: public exahype::plotters::Plotter::UserOnTheFlyPostProcessing{
  public:
  SecondEulerSolver_Plotter0(SecondEulerSolver&  solver);
  virtual ~SecondEulerSolver_Plotter0();
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
