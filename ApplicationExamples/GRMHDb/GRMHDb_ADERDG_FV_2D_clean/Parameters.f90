! DIM Parameters
  
  MODULE Parameters 
    IMPLICIT NONE  
    PUBLIC
    ! These are the PDE-dependent parameters, used in the PDE.f90 and related subroutines.
	!
    ! If you modify the <name>.exahype, please do the following mapping by hand:
    !
    ! solver ADER-DG SRHDSolver / unknowns   ->  goes to ->  nVar
    ! computational-domain / dimension       ->  goes to ->  nDim
    !
	INTEGER, PARAMETER :: N = 3   ! used in VSMM.f90                  
	INTEGER :: myrank_f90
	INTEGER, PARAMETER :: StrandID = 1   ! any postitive integer number: time strand for tecplot output; 0: inactive
    INTEGER, PARAMETER :: nElem_max=1000000     ! max element available for tecplot output.
	!
#if defined(AVX512)    
	INTEGER, PARAMETER :: VECTORLENGTH = 8                    ! Define length of vector registers for AVX 512 
#else 
	INTEGER, PARAMETER :: VECTORLENGTH = 4                    ! Define length of vector registers for AVX 256 
#endif   
	!
#ifdef VECTOR
#ifdef AVX512 
	INTEGER, PARAMETER :: nVar = 24                           ! The number of variables of the PDE system 
#else   
	INTEGER, PARAMETER :: nVar = 20                           ! The number of variables of the PDE system 
#endif 
#else
	INTEGER, PARAMETER :: nVar = 19                           ! The number of variables of the PDE system 
#endif 
  INTEGER, PARAMETER :: nLin = nVar - 2 
	INTEGER, PARAMETER :: nParam = 0                          ! The number of material parameters for the PDE system 
    	! Here, we obtain DIMENSIONS from the C++ preprocessor
	INTEGER, PARAMETER :: d=3
#if defined(Dim3)
    INTEGER, PARAMETER             	:: nDim = 3                   ! The number of space dimensions
#elif defined(Dim2)
    INTEGER, PARAMETER             	:: nDim = 2                   ! The number of space dimensions
#endif
	!
  ! 3-point Gaussian quadrature 
  REAL, PARAMETER     :: sGP3(3) = (/ 0.5-sqrt(15.)/10.,0.5,0.5+sqrt(15.)/10. /) 
  REAL, PARAMETER     :: wGP3(3) = (/ 5./18., 8./18., 5./18. /) 
  INTEGER, PARAMETER  :: nGP3 = 3 
    !CHARACTER(LEN=200), PARAMETER :: ICType      = 'GRMHDAccretion'
    !CHARACTER(LEN=200), PARAMETER :: ICType2     = 'GRMHD-MichelAccretion' 
    !CHARACTER(LEN=200), PARAMETER :: ICType      = 'GRMHDTOV'
    !CHARACTER(LEN=200), PARAMETER :: ICType2     = 'GRMHDTOV' 
    CHARACTER(LEN=200), PARAMETER :: ICType      = 'Riemann0'
    CHARACTER(LEN=200), PARAMETER :: ICType2     = 'Riemann0' 
    !CHARACTER(LEN=200), PARAMETER :: ICType      = 'Sod'
    !CHARACTER(LEN=200), PARAMETER :: ICType2     = 'Sod' 
	
	REAL, PARAMETER 				:: Exc_radius = -1.0  
  
#if defined(SPHERICAL)
    INTEGER, PARAMETER :: nAux = 9
#elif defined(CYLINDRICAL)
    INTEGER, PARAMETER :: nAux = 9
#else
    INTEGER, PARAMETER :: nAux = 0
#endif
	!
	INTEGER ::  nAdd = 0, nMul = 0
	! Important info and parameters concerning the governing PDE system 
	!
	TYPE tEquations 
		REAL(8)    :: gamma, Pi, c0   
		REAL(8)    :: CCZ4k1, CCZ4k2, CCZ4k3, CCZ4eta, CCZ4itau, CCZ4f, CCZ4g, CCZ4xi, CCZ4e, CCZ4c, CCZ4mu, CCZ4ds, CCZ4sk, CCZ4bs  
		REAL(8)    :: cs, alpha, cv, rho0, p0, tau1, tau2, mu, kappa ,tau 
		INTEGER :: CCZ4LapseType, EinsteinAutoAux = 0   
		REAL(8)    :: DivCleaning_a = 1.0 
	END TYPE tEquations 
	! 
	TYPE(tEquations)   :: EQN 
	!  
	REAL(8), DIMENSION(VECTORLENGTH) :: p_floor_VECTOR, rho_floor_VECTOR
	REAL(8), PARAMETER    :: aom = 0.0, Mbh = 1.0
	REAL(8), PARAMETER    :: P_eps = 1e-4  
	!REAL(8), PARAMETER    :: b0_torus_par = 0.0, al_torus_par = 3.8   ! magnetic field parameter and angular momentum parameter for the torus  
	!   
	REAL(8)    :: ExcisionRadius 
	! 
    !Variables for NSTOV module 
#ifdef SPHERICAL  
    INTEGER, PARAMETER :: NSTOV_nODE = 3
#else
    INTEGER, PARAMETER :: NSTOV_nODE = 4
#endif
    INTEGER, PARAMETER :: NSTOV_nODE_p = 3, NSTOV_ATMO = 0
    REAL(8), PARAMETER :: AV00 = 1.0, AV01 = 1e-3, rho01 = 5.0e-5, rho02 = 1.0e-4, irho00 = 1.0e-3, irho01 = 1.0e-6
    REAL(8), PARAMETER :: NSTOV_rho_c = 1.28e-3
    REAL(8), PARAMETER :: NSTOV_kappa = 100
    REAL(8), PARAMETER :: p_floor = 1.0e-16, rho_floor = 1.0e-10
    REAL(8) ::NSTOV_p_atmo, NSTOV_rho_atmo != 1e-10
    REAL(8), PARAMETER :: NSTOV_t_atm=1.0, NSTOV_p_zero = 1e-14
    !REAL(8), PARAMETER :: NSTOV_p_atmo = 1e-13, NSTOV_t_atm=1.0
    !
    TYPE tNSTOVVar
        INTEGER :: Computed
        REAL(8) :: Int
        !INTEGER, PARAMETER :: nODE = 3
        REAL(8) :: Mass, radius,r1,dr1,r2,dr2,rW,drW,rUni,C,p_R ,rho_R,lapse_C
        INTEGER :: n_B,iradius,ir1,ir2,irW 
        REAL(8), DIMENSION (:,:),     ALLOCATABLE :: q,dq
        REAL(8), DIMENSION (:),     ALLOCATABLE :: r,dr
        REAL(8), DIMENSION (:),     ALLOCATABLE :: qloc
    END TYPE tNSTOVVar
    TYPE(tNSTOVVar) :: NSTOVVar,NSTOVVar_bar,NSTOVVar_barNew
    !
    REAL(8), PARAMETER :: CoordTol =1e-11
  END MODULE Parameters  
