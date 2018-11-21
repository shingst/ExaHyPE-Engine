// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
// ========================
//   www.exahype.eu
// ========================
#ifndef ErrorWriter_CLASS_HEADER_
#define ErrorWriter_CLASS_HEADER_

#include "NavierStokesSolverDG.h"
#include "exahype/plotters/ADERDG2UserDefined.h"
#include <array>

namespace NavierStokes {
class ErrorWriter;
}

class NavierStokes::ErrorWriter : public exahype::plotters::ADERDG2UserDefined {
 private:
  double timeStamp;
  std::string filename;
  bool isMpi;
  NavierStokesSolverDG* solver;

  using Array_t = std::array<double, NavierStokesSolverDG::NumberOfVariables>;
  double hmin;
  Array_t errorL1, errorL2, errorLInf;
  Array_t normL1Ana, normL2Ana, normLInfAna;

 public:
  void init(const std::string& filename, int orderPlusOne, int solverUnknowns,
            int writtenUnknowns,
            exahype::parser::ParserView plotterParameters) override;

  /**
   * Constructor.
   *
   * \note ExaHyPE does not increment file counters for
   * you if you use user defined plotting. You have
   * to declare and manage such member variables yourself.
   */
  ErrorWriter(NavierStokesSolverDG& solver);

  /**
   * This method is invoked every time a cell
   * is touched by the plotting device.
   *
   * \note Use the protected variables _order, _variables to
   * determine the size of u.
   * The array u has the size _variables * (_order+1)^DIMENSIONS.
   *
   * \param[in] offsetOfPatch the offset of the cell/patch.
   * \param[in] sizeOfPatch the offset of the cell/patch.
   * \param[in] u the degrees of freedom "living" inside of the patch.
   */
  void plotPatch(const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
                 const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
                 double* u, double timeStamp) override;

  /**
   * This method is called at the beginning of the plotting.
   * You can use it to reset member variables, e.g., those
   * used for calculations, or to increment file counters.
   *
   * \param[in] time a characteristic solver time stamp.
   *            Usually the global minimum.
   */
  void startPlotting(double time) override;

  /**
   * This method is called at the end of the plotting.
   * You can use it to reset member variables, finalise calculations (compute
   * square roots etc.), or to increment file counters
   */
  void finishPlotting() override;
};

#endif /* ErrorWriter_CLASS_HEADER_ */
