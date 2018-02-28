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
    INTEGER, PARAMETER             	:: nVar = 30                           ! The number of variables of the PDE system
	REAL, PARAMETER 				:: gamma = 5./3.
	REAL, PARAMETER 				:: cs = 0.6
	REAL, PARAMETER 				:: tau =0.1! 1.e+20
	REAL, PARAMETER 				:: mu = 1./6.*cs**2*tau


	!REAL, PARAMETER 				:: tau2 = 0.3111111111e-2!1.e+20
	REAL, PARAMETER    ::  p_floor = 1.0e-11, rho_floor = 1.0e-10
	REAL, PARAMETER    :: tol = 1e-12
  END MODULE Parameters  
