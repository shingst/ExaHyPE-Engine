!
! The typesDef module is a global storage for runtime variables and compile time parameters.
! Here in the CCZ4 ExaHyPE application, it is a remnant of the much larger typesDef from the
! Trento code.
! By using ISO C bindings, we can access these global variables with a well-defined name from C.
! They are subsequently set at initialization.
!
MODULE typesDef
  USE ISO_C_BINDING
  IMPLICIT NONE 
  PUBLIC 

  ! The number of space dimensions that we actually want to simulate 
  ! The CCZ4 application is currently written always in 3D.
  ! Always use 3D when compiling it.
#if defined(Dim3)
    INTEGER, PARAMETER             :: nDim = 3
#elif defined(Dim2)
    INTEGER, PARAMETER             :: nDim = 2
#endif

  INTEGER, PARAMETER :: nVar = 59     ! The number of variables of the PDE system 
  INTEGER, PARAMETER :: d = 3         ! Maximal number of dimensions supported by the code. Always 3.

  ! A derived type in Fortran which translates to a structure in C.
  ! It holds the parameters governing the CCZ4 PDE system
  TYPE, BIND(C) :: tEquations 
      REAL(C_DOUBLE) :: CCZ4k1        ! CCZ4 damping parameter k1 
      REAL(C_DOUBLE) :: CCZ4k2        ! CCZ4 damping parameter k2  
      REAL(C_DOUBLE) :: CCZ4k3        ! CCZ4 damping parameter k3 
      REAL(C_DOUBLE) :: CCZ4eta       ! CCZ4 damping parameter for the PDE of b^i in the gamma driver 
      REAL(C_DOUBLE) :: CCZ4itau      ! inverse relaxation time for the relaxed approach in order to enforce the unit determinant of the conformal metric and the trace-free matrix A_ij
      REAL(C_DOUBLE) :: CCZ4f         ! set f=0.75 or f=1.0 for the gamma driver. Typical BSSNOK value: f=0.75. Set f=0 to avoid shift evolution, i.e. to have d/dt beta^i=0.  
      REAL(C_DOUBLE) :: CCZ4g         ! not used at the moment: reserved for the generalized harmonic shift   
      REAL(C_DOUBLE) :: CCZ4xi        ! set to zero to switch off the gamma driver and to avoid shift evolution (i.e. to have d/dt b^i = 0) 
      REAL(C_DOUBLE) :: CCZ4e         ! cleaning speed e>=1 for the Hamiltonian constraint. Typical value for CCZ4 is e=1. However, e>1 gives better constraints and better hyperbolicity!  
      REAL(C_DOUBLE) :: CCZ4c         ! set c=0 to remove the algebraic source terms of the type -2*Theta 
      REAL(C_DOUBLE) :: CCZ4mu        ! mu=1 adds second order ordering constraints. Has shown to be important for the gamma driver, but too large values can also become a problem... 
      REAL(C_DOUBLE) :: CCZ4ds        ! set this value always to ds=1, unless you know what you are doing. It allows to increase the cleaning speed for the momentum constraints, but in CCZ4 this does not seem to do what it should do...  
      REAL(C_DOUBLE) :: CCZ4sk        ! setting sk=0 removes the contribution of the shift also in the auxiliary variables. If you want to evolve the shift, set sk=1. 
      REAL(C_DOUBLE) :: CCZ4bs        ! set bs=1 if you want to activate the shift convection for beta, b and B (standard CCZ4 formulation). set it to bs=0 to switch off shift convection for those quantities 
      INTEGER(C_INT) :: CCZ4LapseType ! LapseType = 0 is harmonic lapse, LapseType = 1 is 1+log slicing.   
  END TYPE tEquations

  ! A global "singleton" instance of the tEquations structure.
  ! Access from Fortran is like EQN%CCZ4k1 after "USE typesDef"
  ! Access from C is like "typesDef_eqn.k1" when defined as external struct.
  TYPE(tEquations), bind(C, name="typesDef_eqn") :: EQN 

  ! To distinguish the initial data, before we had a
  !!  CHARACTER(LEN=200) :: ICType  
  ! However, while Fortran<->C strings work in practice, compilers print plenty of warnings.
  ! Therefore, we use a safer ENUM idiom here. Fortran 2003 does not yet have enums, so we
  ! use this kind of "run-time enum" structure which also serves as initial data registry.
  TYPE, BIND(C) :: tInitialDataEnum
      INTEGER(C_INT) :: CCZ4GaugeWave = 10
      INTEGER(C_INT) :: CCZ4Puncture = 20
      INTEGER(C_INT) :: CCZ4TwoPunctures = 30
      INTEGER(C_INT) :: CCZ4GRMHDAccretion = 40
      INTEGER(C_INT) :: CCZ4Kerr2D = 50
  END TYPE tInitialDataEnum

  TYPE(tInitialDataEnum), bind(C, name="typesDef_ICType") :: ICType

  ! The storage is not a tInitialDataEnum, but legit values are members of the tInitialDataEnum.
  ! Access from Fortran is like IC = ICType%CCZ4GaugeWave
  ! Access from C is like typesDef_ic = typesDef_ICType.CCZ4GaugeWave after having defined the struct.
  INTEGER(C_INT), bind(C, name="typesDef_ic") :: IC

END MODULE typesDef 
