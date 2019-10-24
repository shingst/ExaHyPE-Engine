#ifndef HEADER_MAINVARIABLES
#define HEADER_MAINVARIABLES

! DIM Parameters
  
  MODULE MainVariables 
    IMPLICIT NONE  
    PUBLIC  

    ! This typesDef.f90 is a relict still from the old Fortran interface in Exahype.
    ! However, the following two variables are needed in the Fortran code. They
    ! should provided by the glue code generator in an appropriate way.

    ! If you modify the SRHD.exahype, please do the following mapping by hand:
    !
    ! solver ADER-DG SRHDSolver / unknowns   ->  goes to ->  nVar
    ! computational-domain / dimension       ->  goes to ->  nDim
    !

    ! Here, we obtain DIMENSIONS from the C++ preprocessor
#if defined(Dim3)
    INTEGER, PARAMETER             	:: nDim = 3                   ! The number of space dimensions
#elif defined(Dim2)
    INTEGER, PARAMETER             	:: nDim = 2                   ! The number of space dimensions
#endif
	INTEGER, PARAMETER             	:: nAux = 0
    INTEGER, PARAMETER             	:: nVar = 96                           ! The number of variables of the PDE system 
    INTEGER, PARAMETER 				:: nLin = 0
	INTEGER, PARAMETER				:: nParam=0
	INTEGER, PARAMETER				:: d=3
    !CHARACTER(LEN=20), PARAMETER	:: ICType='NLOPRUPTURE'
	CHARACTER(LEN=20)				:: ICType
  TYPE tEquations 
      REAL(8)    :: gamma, Pi, c0, g = 9.81, friction = 1.0     
      REAL(8)    :: CCZ4k1, CCZ4k2, CCZ4k3, CCZ4eta, CCZ4itau, CCZ4f, CCZ4g, CCZ4xi, CCZ4e, CCZ4c, CCZ4mu, CCZ4ds, CCZ4sk, CCZ4bs  
      REAL(8)    :: CCZ4GLMc0 = 0.5, CCZ4GLMc = 0.75, CCZ4GLMd = 0.75, CCZ4GLMepsD = 1e-2, CCZ4GLMepsA = 1e-2, CCZ4GLMepsP = 1e-2, cs, alpha, beta, lambda, cv, rho0, p0, tau1, tau2, mu, kappa ,tau 
      INTEGER :: CCZ4LapseType, EinsteinAutoAux = 0, ReferenceDepth = 1.0    
      REAL(8)    :: DivCleaning_a = 1.0 
  END TYPE tEquations 
	
	TYPE(tEquations) :: EQN
	
  ! 3-point Gaussian quadrature 
  REAL, PARAMETER     :: sGP3(3) = (/ 0.5-sqrt(15.)/10.,0.5,0.5+sqrt(15.)/10. /) 
  REAL, PARAMETER     :: wGP3(3) = (/ 5./18., 8./18., 5./18. /) 
  INTEGER, PARAMETER  :: nGP3 = 3 

  END MODULE MainVariables  
#endif
