#ifndef _EXAHYPE_KERNELS_ADERDG_GENERIC_PDEFLUXES_H_
#define _EXAHYPE_KERNELS_ADERDG_GENERIC_PDEFLUXES_H_

#include "tarch/la/Vector.h"
#include "peano/utils/Globals.h"

namespace kernels {
  namespace aderdg {
    namespace generic {

      // @todo Dominic Etienne Charrier
      // Inconsistent ordering of inout and in arguments for
      // template argument functions and non-template argument function.
      template <void PDEFlux(const double* const Q,double* f,double* g)>
      void spaceTimePredictor(
          double* lQi,
          double* lFi,
          double* lQhi,
          double* lFhi,
          double* lQhbnd,
          double* lFhbnd,
          const double* const luh,
          const tarch::la::Vector<DIMENSIONS,double>&  dx,
          const double predictorTimeStepSize,
          const int numberOfVariables,
          const int basisSize
      );

      // @todo Dominic Etienne Charrier
      // Inconsistent ordering of inout and in arguments for
      // template argument functions and non-template argument function.
      template <void PDEFlux(const double* const Q,double* f,double* g)>
      void spaceTimePredictor(
          double* lQi,
          double* lFi,
          const double* const luh,
          const tarch::la::Vector<DIMENSIONS,double>&  dx,
          const double predictorTimeStepSize,
          const int numberOfVariables,
          const int basisSize
      );

      /**
       * (At the moment, we always evaluate the time averaged space-time
       * predictor unknowns.)
       * todo docu
       */
      void predictor(
          double* lQhi,
          double* lFhi,
          const double* const lQi,
          const double* const lFi,
          const double predictorTimeStepSize,
          const int numberOfVariables,
          const int basisSize
      );

      /**
       * @todo Dominic Etienne Charrier
       * This is just a "parent" function that
       * invokes the function going by the same
       * name 2*dim times.
       */
      void extrapolatedPredictor(
          double* lQhbnd,
          double* lFhbnd,
          const double* const lQhi,
          const double* const lFhi,
          const double predictorTimeStepSize,
          const int numberOfVariables,
          const int basisSize
      );

      /**
       * @todo Dominic Etienne Charrier
       * docu
       * Note that we need to replace lQhi and LFhi by
       * the space-time predictor unknowns if we want to employ
       * local/anarchic time stepping. Since we will have
       * to perform evaluations of the extrapolated boundary fluxes
       * at various appropriate times in this case.
       * The evaluation of the extrapolated predictor requires a time integration
       * of the space-time predictor unknowns. Clearly,
       * \p evaluationTimeStepSize must be smaller than or equal to
       * \p predictorTimeStepSize.
       * At the moment, we always evaluate the time averaged space-time
       * predictor unknowns. Thus it is not necessary to pass these values.
       */
      void extrapolatedPredictorXDirection(
          double* lQhbnd,
          double* lFhbnd,
          const double* const lQhi,
          const double* const lFhi,
          const int facePosition, // 0 for "left" face, 1 far "right" face
          const double evaluationTimeStepSize,
          const double predictorTimeStepSize,
          const int numberOfVariables,
          const int basisSize
      );

      void extrapolatedPredictorYDirection(
          double* lQhbnd,
          double* lFhbnd,
          const double* const lQhi,
          const double* const lFhi,
          const int facePosition, // 0 for "left" face, 1 far "right" face
          const double evaluationTimeStepSize,
          const double predictorTimeStepSize,
          const int numberOfVariables,
          const int basisSize
      );
#if DIMENSIONS == 3

      void extrapolatedPredictorZDirection(
          double* lQhbnd,
          double* lFhbnd,
          const double* const lQhi,
          const double* const lFhi,
          const int facePosition, // 0 for "left" face, 1 far "right" face
          const double evaluationTimeStepSize,
          const double predictorTimeStepSize,
          const int numberOfVariables,
          const int basisSize
      );
#endif

      // todo Dominic Etienne Charrier:
      // The DIMENSIONS depending mesh size vector enables overloading at the moment.
      // If we replace it by scalar mesh size, we have to add a template argument "int dim".

      void solutionUpdate(
          double* luh,
          const double* const lduh,
          const double dt,
          const int numberOfVariables,
          const int basisSize
      );

      void volumeIntegral(
          double* lduh,
          const double* const lFhi,
          const tarch::la::Vector<DIMENSIONS,double>&  dx,
          const int numberOfVariables,
          const int basisSize

      );


      // todo 10/02/16: Dominic
      // Keep only one surfaceIntegral.
      void surfaceIntegral(
          double* lduh,
          const double* const lFhbnd,
          const tarch::la::Vector<DIMENSIONS,double>&  dx,
          const int numberOfVariables,
          const int basisSize
      );

      /*void surfaceIntegral2(
          double* lduh,
          const double* const lFhbnd,
          const tarch::la::Vector<DIMENSIONS,double>&  dx,
          const int numberOfVariables,
          const int basisSize
      );*/

      void surfaceIntegralXDirection(
          double * lduh,
          const double * const lFhbnd,
          const double area,
          const int facePosition,  // 0 for "left" face, 1 for "right" face.
          const double updateSign, // -1 for "left" face, 1 for "right" face.
          const int numberOfVariables,
          const int basisSize
      );

      void surfaceIntegralYDirection(
          double * lduh,
          const double * const lFhbnd,
          const double area,
          const int facePosition,  // 0 for "left" face, 1 for "right" face.
          const double updateSign, // -1 for "left" face, 1 for "right" face.
          const int numberOfVariables,
          const int basisSize
      );

#if DIMENSIONS == 3

      void surfaceIntegralZDirection(
          double * lduh,
          const double * const lFhbnd,
          const double area,
          const int facePosition,  // 0 for "left" face, 1 for "right" face.
          const double updateSign, // -1 for "left" face, 1 for "right" face.
          const int numberOfVariables,
          const int basisSize
      );
#endif

      // @todo Dominic Etienne Charrier
      // Inconsistent ordering of inout and in arguments for
      // template argument functions and non-template argument function.
      template <void PDEInitialValues(const double* const x,double* Q)>
      void initialCondition(
          double* luh,
          const tarch::la::Vector<DIMENSIONS,double>& center,
          const tarch::la::Vector<DIMENSIONS,double>& dx,
          const int numberOfVariables,
          const int basisSize
      );

      // @todo Dominic Etienne Charrier
      // Inconsistent ordering of inout and in arguments
      // template argument functions and non-template argument function.
      template <void PDEEigenvalues(const double* const Q,const int normalNonZero,double* lambda)>
      void riemannSolver(
          double* FL,
          double* FR,
          const double* const QL,
          const double* const QR,
          const double dt,
          const int normalNonZero,
          const int numberOfVariables,
          const int basisSize
      );

      // @todo Dominic Etienne Charrier
      // Inconsistent ordering of inout and in arguments for
      // template argument functions and non-template argument function.
      template <void PDEEigenvalues(const double* const Q,const int normalNonZero,double* lambda)>
      double stableTimeStepSize(
          const double* const luh,
          const tarch::la::Vector<DIMENSIONS,double>& dx,
          const int numberOfVariables,
          const int basisSize
      );
    }
  }
}


#include "kernels/aderdg/generic/Kernels.cpph"

#endif
