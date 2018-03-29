
! was ADERDGInit
RECURSIVE SUBROUTINE EQNInit
    USE typesDef
    IMPLICIT NONE
    ! Local variables 
    ! 
    !
    ! Some info about the PDE system 
    !
    
    IC = ICType%CCZ4GaugeWave
    
    ! no shift for black hole simulation 
    ! Damit klappt Kerr-Wuerfel
    ! Gamma driver for the black hole simulation! Gamma driver for the black hole simulation 
    IF(IC.EQ.ICType%CCZ4TwoPunctures.OR.IC.EQ.ICType%CCZ4Puncture) THEN
    EQN%CCZ4k1   = 0.1      ! CCZ4 damping parameter k1 
    EQN%CCZ4k2   = 0.1      ! CCZ4 damping parameter k2  
    EQN%CCZ4k3   = 0.1      ! CCZ4 damping parameter k3 
    EQN%CCZ4eta  = 0.0      ! CCZ4 damping parameter for the PDE of b^i in the gamma driver 
    EQN%CCZ4itau = 0.0      ! inverse relaxation time for the relaxed approach in order to enforce the unit determinant of the conformal metric and the trace-free matrix A_ij
    EQN%CCZ4f    = 1.0      ! set f=0.75 or f=1.0 for the gamma driver. Typical BSSNOK value: f=0.75. Set f=0 to avoid shift evolution, i.e. to have d/dt beta^i=0.  
    EQN%CCZ4g    = 0.0      ! not used at the moment: reserved for the generalized harmonic shift   
    EQN%CCZ4xi   = 1.0      ! set to zero to switch off the gamma driver and to avoid shift evolution (i.e. to have d/dt b^i = 0) 
    EQN%CCZ4e    = 1.0      ! cleaning speed e>=1 for the Hamiltonian constraint. Typical value for CCZ4 is e=1. However, e>1 gives better constraints and better hyperbolicity!  
    EQN%CCZ4c    = 1.0      ! set c=0 to remove the algebraic source terms of the type -2*Theta 
    EQN%CCZ4mu   = 0.2      ! mu=1 adds second order ordering constraints. Has shown to be important for the gamma driver, but too large values can also become a problem... 
    EQN%CCZ4ds   = 1.0      ! set this value always to ds=1, unless you know what you are doing. It allows to increase the cleaning speed for the momentum constraints, but in CCZ4 this does not seem to do what it should do...  
    EQN%CCZ4sk   = 1.0      ! setting sk=0 removes the contribution of the shift also in the auxiliary variables. If you want to evolve the shift, set sk=1. 
    EQN%CCZ4bs   = 1.0      ! set bs=1 if you want to activate the shift convection for beta, b and B (standard CCZ4 formulation). set it to bs=0 to switch off shift convection for those quantities 
    EQN%CCZ4LapseType   = 1 ! LapseType = 0 is harmonic lapse, LapseType = 1 is 1+log slicing.   
    
    ELSE IF(IC.EQ.ICType%CCZ4GRMHDAccretion) THEN
    
    EQN%CCZ4k1   = 0.1      ! CCZ4 damping parameter k1 
    EQN%CCZ4k2   = 0.1      ! CCZ4 damping parameter k2  
    EQN%CCZ4k3   = 0.1      ! CCZ4 damping parameter k3 
    EQN%CCZ4eta  = 0.0      ! CCZ4 damping parameter for the PDE of b^i in the gamma driver 
    EQN%CCZ4itau = 0.0      ! inverse relaxation time for the relaxed approach in order to enforce the unit determinant of the conformal metric and the trace-free matrix A_ij
    EQN%CCZ4f    = 1.0      ! set f=0.75 or f=1.0 for the gamma driver. Typical BSSNOK value: f=0.75. Set f=0 to avoid shift evolution, i.e. to have d/dt beta^i=0.  
    EQN%CCZ4g    = 0.0      ! not used at the moment: reserved for the generalized harmonic shift   
    EQN%CCZ4xi   = 1.0      ! set to zero to switch off the gamma driver and to avoid shift evolution (i.e. to have d/dt b^i = 0) 
    EQN%CCZ4e    = 1.0      ! cleaning speed e>=1 for the Hamiltonian constraint. Typical value for CCZ4 is e=1. However, e>1 gives better constraints and better hyperbolicity!  
    EQN%CCZ4c    = 1.0      ! set c=0 to remove the algebraic source terms of the type -2*Theta 
    EQN%CCZ4mu   = 0.2      ! mu=1 adds second order ordering constraints. Has shown to be important for the gamma driver, but too large values can also become a problem... 
    EQN%CCZ4ds   = 1.0      ! set this value always to ds=1, unless you know what you are doing. It allows to increase the cleaning speed for the momentum constraints, but in CCZ4 this does not seem to do what it should do...  
    EQN%CCZ4sk   = 1.0      ! setting sk=0 removes the contribution of the shift also in the auxiliary variables. If you want to evolve the shift, set sk=1. 
    EQN%CCZ4bs   = 0.0      ! set bs=1 if you want to activate the shift convection for beta, b and B (standard CCZ4 formulation). set it to bs=0 to switch off shift convection for those quantities 
    EQN%CCZ4LapseType   = 1 ! LapseType = 0 is harmonic lapse, LapseType = 1 is 1+log slicing.   
    
    ELSE IF(IC.EQ.ICType%CCZ4GaugeWave) THEN
    ! Parameters for the Gauge wave 
    EQN%CCZ4k1  = 0.0
    EQN%CCZ4k2  = 0.0 
    EQN%CCZ4k3  = 0.0 
    EQN%CCZ4eta = 0.0 
    EQN%CCZ4f   = 0.0 
    EQN%CCZ4g   = 0.0 
    EQN%CCZ4xi  = 0.0 
    EQN%CCZ4e   = 2.0 
    EQN%CCZ4c   = 0.0 
    EQN%CCZ4mu  = 0.0 
    EQN%CCZ4ds  = 1.0 
    EQN%CCZ4sk  = 0.0
    EQN%CCZ4bs  = 0.0
    EQN%CCZ4LapseType   = 0 ! harmonic lapse 
    
    ELSE
        PRINT *, "Unkown ICtype in EQNInit"
        ! CALL ABORT
    END IF
END SUBROUTINE EQNInit

