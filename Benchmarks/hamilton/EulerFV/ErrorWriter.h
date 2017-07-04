// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
// ========================
//   www.exahype.eu
// ========================
#ifndef ErrorWriter_CLASS_HEADER_
#define ErrorWriter_CLASS_HEADER_

#include "exahype/plotters/FiniteVolumes2UserDefined.h"

#include "EulerSolver.h"

namespace EulerFV {
  class ErrorWriter;
}

class EulerFV::ErrorWriter : public exahype::plotters::FiniteVolumes2UserDefined {
 private:
   double _timeStamp;

   double errorL1  [EulerSolver::NumberOfVariables]; 
   double errorL2  [EulerSolver::NumberOfVariables]; 
   double errorLInf[EulerSolver::NumberOfVariables];
   
   double normL1Ana  [EulerSolver::NumberOfVariables];  
   double normL2Ana  [EulerSolver::NumberOfVariables];  
   double normLInfAna[EulerSolver::NumberOfVariables]; 

 public:
  /**
   * Constructor.
   * 
   * \note ExaHyPE does not increment file counters for
   * you if you use user defined plotting. You have
   * to declare and manage such member variables yourself. 
   */
  ErrorWriter();

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
  void plotPatch(
      const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
      double timeStamp) override;

  /** 
   * This method is called at the beginning of the plotting.
   * You can use it to reset member variables, e.g., those
   * used for calculations, or to increment file counters.
   *
   * \param[in] time a characteristic solver time stamp.
   *            Usually the global minimum.
   */
  void startPlotting( double time) override;
  
  /** 
   * This method is called at the end of the plotting.
   * You can use it to reset member variables, finalise calculations (compute square roots etc.),
   * or to increment file counters
   */
  void finishPlotting() override;
};

#endif /* ErrorWriter_CLASS_HEADER_ */
