// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include "exahype/plotters/Plotter.h"
namespace GRMHD{
  class IntegralsWriter;

  /**
   * Forward declaration
   */
  class GRMHDSolver_ADERDG;
}

#include "exahype/plotters/ascii/MultipleReductionsWriter.h"
#include "GRMHDSolver_ADERDG.h"

class GRMHD::IntegralsWriter: public exahype::plotters::Plotter::UserOnTheFlyPostProcessing{
  public:
  bool plotForADERSolver;
  static const int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;
  exahype::plotters::ascii::MultipleReductionsWriter conserved;
  exahype::plotters::ascii::MultipleReductionsWriter primitives;
  exahype::plotters::ascii::MultipleReductionsWriter errors;
  exahype::plotters::ascii::ReductionsWriter statistics;
  exahype::plotters::ascii::ReductionsWriter masschange;
  
  IntegralsWriter();
  IntegralsWriter(GRMHDSolver_ADERDG&  solver);

  virtual ~IntegralsWriter();
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