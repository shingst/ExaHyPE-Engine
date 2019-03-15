! MHD Parameters: Riemann problem, Rotor problem, Blast wave
  
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
    INTEGER, PARAMETER :: d = 3
#if defined(Dim3)
    INTEGER, PARAMETER             :: nDim = 3                   ! The number of space dimensions
#elif defined(Dim2)
    INTEGER, PARAMETER             :: nDim = 2                   ! The number of space dimensions
#endif
#ifdef VECTOR    
#ifdef AVX512 
  INTEGER, PARAMETER :: nVar = 24                           ! The number of variables of the PDE system 
#else   
  INTEGER, PARAMETER :: nVar = 20                           ! The number of variables of the PDE system 
#endif 
#else
  INTEGER, PARAMETER :: nVar = 19                           ! The number of variables of the PDE system 
#endif
    INTEGER, PARAMETER :: nparam = 0
    ! Ideal EOS:
    ! 4/3 used in ADERDG3D-SRMHD-Z4 by MD, 01.Nov2016
    ! 2.0 used for TOV stars
    !REAL, PARAMETER                :: gamma = 4.0/3.0
    REAL, PARAMETER                :: gamma = 2.0
    REAL, PARAMETER                :: igamma = 1.0/gamma ! 4.0/3.0
    !
    ! Divergence cleaning:
    REAL, PARAMETER :: DivCleaning_a = 1.0
   !Variables for NSTOV module 
#ifdef SPHERICAL  
    INTEGER, PARAMETER :: NSTOV_nODE = 3
#else
    INTEGER, PARAMETER :: NSTOV_nODE = 4
#endif
    INTEGER, PARAMETER :: NSTOV_nODE_p = 3, NSTOV_ATMO = 0 ! 1 => Landau atmo
    REAL, PARAMETER :: NSTOV_rho_c = 1.28e-3
    REAL, PARAMETER :: NSTOV_kappa = 100
    REAL, PARAMETER :: NSTOV_t_atm=1.0
    !
    TYPE tNSTOVVar
        INTEGER :: Computed
        REAL :: Int
        !INTEGER, PARAMETER :: nODE = 3
        REAL :: Mass, radius,r1,dr1,r2,dr2,rW,drW,rUni,C,p_R ,rho_R,lapse_C
        INTEGER :: n_B,iradius,ir1,ir2,irW 
        REAL, DIMENSION (:,:),     ALLOCATABLE :: q,dq
        REAL, DIMENSION (:),     ALLOCATABLE :: r,dr
        REAL, DIMENSION (:),     ALLOCATABLE :: qloc
    END TYPE tNSTOVVar
    TYPE(tNSTOVVar) :: NSTOVVar,NSTOVVar_bar,NSTOVVar_barNew
    !
    REAL, PARAMETER :: CoordTol =1e-11
    ! 
	REAL, PARAMETER :: p_floor = 1.0e-16
	REAL, PARAMETER :: rho_floor = 1.0e-10 
	REAL, PARAMETER :: NSTOV_p_atmo = 1e-15
    REAL, PARAMETER :: NSTOV_rho_atmo = (NSTOV_p_atmo/NSTOV_kappa)**igamma  

  REAL, PARAMETER :: ExcisionRadius=1.5
	!
	REAL, PARAMETER :: aom = 0.0, Mbh = 1.0
  !CHARACTER(LEN=200), PARAMETER  :: ICType = TRIM('GRMHD-SphericalAccretion')
  !CHARACTER(LEN=200), PARAMETER  :: ICType2 = TRIM('GRMHD-SphericalAccretion')
!  CHARACTER(LEN=200), PARAMETER  :: ICType2 = TRIM('GRMHD-BondiAccretion')
  CHARACTER(LEN=200), PARAMETER  :: ICType = TRIM('GRMHDTOV')
  CHARACTER(LEN=200), PARAMETER  :: ICType2 = TRIM('GRMHD-SphericalAccretion')
	!ICType = 'GRMHDTOV' 
    !
  END MODULE Parameters  
