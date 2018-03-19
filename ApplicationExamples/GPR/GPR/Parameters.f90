! DIM Parameters
  
  MODULE Parameters 
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
    INTEGER, PARAMETER             	:: nVar = 17                           ! The number of variables of the PDE system 
	REAL, PARAMETER 				:: rho0 = 1.0
	REAL, PARAMETER 				:: cs = 5.0
	REAL, PARAMETER 				:: cv = 1.0
	REAL, PARAMETER 				:: gamma = 1.4
	REAL, PARAMETER 				:: alpha = 5.0
	REAL, PARAMETER 				:: p0 = 0.0
	REAL, PARAMETER 				:: tau1 =1.e-2! 1.e+20
	REAL, PARAMETER 				:: tau2 = 0.3111111111e-2!1.e+20
  END MODULE Parameters  
