! MHD Initial Data
 
 SUBROUTINE MinimumTreeDepth(depth)
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE 
    ! Argument list 
    INTEGER, INTENT(OUT)              :: depth        ! maximal depth of tree recursion
    
    depth = 4
    
END SUBROUTINE MinimumTreeDepth 

SUBROUTINE HasToAdjustSolution(time, refine)
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE 
    ! Argument list 
    REAL   , INTENT(IN)               :: time        ! 

    LOGICAL, INTENT(OUT)              :: refine      ! 
    
    IF(time<0.000000001) THEN
      refine = .TRUE.
    ELSE
      refine = .FALSE.
    ENDiF
    
END SUBROUTINE HasToAdjustSolution


SUBROUTINE AdjustedSolutionValues(x, w, t, dt, Q)
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar, nDim
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: x(nDim)        ! 
	REAL, INTENT(IN)               :: w           ! 
	REAL, INTENT(IN)               :: t           ! 
	REAL, INTENT(IN)               :: dt          ! 

	REAL, INTENT(OUT)              :: Q(nVar)        ! 
	
	IF ( t < 1e-15 ) THEN
		CALL InitialData(x, Q)
	ENDIF
END SUBROUTINE AdjustedSolutionValues

SUBROUTINE InitialData(x, Q)
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar, nDim
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: x(nDim)        ! 
	REAL, INTENT(OUT)              :: Q(nVar)        ! 

	! We call a C++ function which helps us to get access to the
	! exahype specification file constants
	INTERFACE
		SUBROUTINE InitialDataByExaHyPESpecFile(x,Q) BIND(C)
			USE, INTRINSIC :: ISO_C_BINDING
			USE Parameters, ONLY : nVar, nDim
			IMPLICIT NONE
			REAL, INTENT(IN)               :: x(nDim)
			REAL, INTENT(OUT)              :: Q(nVar)
		END SUBROUTINE InitialDataByExaHyPESpecFile
	END INTERFACE
	
	! Call here one of
	! CALL InitialBlast(x, Q)
  ! Call InitialAlfenWave(x, Q)
	Call InitialRotor(x,Q)
	! Call InitialBlast(x, Q)
	! Call InitialOrsagTang(x, Q)

	! CALL InitialDataByExaHyPESpecFile(x,Q)
END SUBROUTINE InitialData

SUBROUTINE InitialAlfenWave(x, Q)
    ! Get the AlfenWave for t=0
    ! This subroutine is neccessary for a common argument interface InitialFooBar(x,Q)
    ! for all initial data.
    USE, INTRINSIC :: ISO_C_BINDING
    Use Parameters, ONLY : nDim, nVar
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 

    CALL AlfenWave(x, Q, 0.0)
END SUBROUTINE InitialAlfenWave

SUBROUTINE AlfenWave(x, Q, t)
    ! Computes the AlfenWave at a given time t.
    ! Use it ie. with t=0 for initial data
    ! Use it for any other time ie. for comparison

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim, gamma
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 

    REAL :: rho0, p0, eta, B0, hh, tempaa, tempab, tempac, va2, vax
    REAL :: V(nVar), BV(3), VV(3), Pi = ACOS(-1.0)

    rho0 = 1.
    p0   = 1.
    eta  = 1.
    B0   = 1.0 
    !
    hh = 1.0 + gamma / ( gamma - 1.0) * p0 / rho0
    tempaa = rho0 * hh + B0**2 * ( 1.0 + eta**2)
    tempab = 2.0 * eta * B0**2 / tempaa
    tempac = 0.5 * ( 1.0 + sqrt ( 1.0 - tempab**2))
    va2 = b0**2 / ( tempaa * tempac)
    vax = sqrt ( va2)
    !
    BV(1) = B0
    BV(2) = eta * B0 * COS(2*Pi*( x(1) - vax*t))
    BV(3) = eta * B0 * SIN(2*Pi*( x(1) - vax*t))
    !
    VV(1)   = 0.0
    VV(2:3) = - vax * BV(2:3) / B0
    !
    ! Now convert to conservative variables
    !
    V = (/ rho0, VV(1:3), p0, BV(1:3), 0.0 /)
    CALL PDEPrim2Cons(Q,V)
END SUBROUTINE AlfenWave

SUBROUTINE InitialBlast(x, Q)
    ! Blast wave:
    ! Simulation domain:  -6 .. +6

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 
    
    REAL :: rho0, p0, p1, rho1, r, r0, r1, taper
    REAL :: V(nVAR)


    p0   = 5.0e-4
    p1   = 1.0
    rho0 = 1.0e-4 
    rho1 = 1.0e-2
    ! Here, we set the primitive variables 
    V(1) = 0.0   ! rho
    V(2) = 0.0   ! vx
    V(3) = 0.0   ! vy
    V(4) = 0.0   ! vz
    V(5) = 0.0   ! p
    V(6) = 0.1   ! bx
    V(7) = 0.0   
    V(8) = 0.0   
    V(9) = 0.0   ! preserving
    r = SQRT(x(1)**2 + x(2)**2)
    r0 = 0.8
    r1 = 1.0
    !
    IF(r .LE. 0.8) THEN
        V(1) = rho1
        V(5) = p1  
    ELSEIF ( r.GT.0.8 .AND. r.LE.1.0) THEN
        taper = (1.-(r-r0)/(r1-r0))
        V(1) = rho0 + taper*(rho1-rho0)
        taper = (1.-(r-r0)/(r1-r0))
        V(5) = p0 + taper*(p1-p0)
    ELSE  
        V(1) = rho0
        V(5) = p0 
    ENDIF

    CALL PDEPrim2Cons(Q,V)
END SUBROUTINE

SUBROUTINE InitialOrsagTang(x, Q)
    ! Simulation Domain: 0 .. 2*pi = 6.283185307179586

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 

    REAL :: rho0, p0, vel0, B0
    REAL :: V(nVAR), VV(3), BV(3)

    ! RMHDOrszagTang
    rho0 = 1.0
    p0   = 1.0
    vel0 = 0.75
    B0   = 1.0
    !
    VV(1) = -1./SQRT(2.0)*vel0*SIN(x(2))
    VV(2) =  1./SQRT(2.0)*vel0*SIN(x(1))
    VV(3) = 0. 
    !
    BV(1) = -B0*SIN(x(2))
    BV(2) =  B0*SIN(2.0*x(1))
    BV(3) =  0.
    !
    V = (/ rho0, VV(1:3), p0, BV(1:3), 0.0 /)

    Call PDEPrim2Cons(Q,V)
END SUBROUTINE InitialOrsagTang

SUBROUTINE InitialRotor(x,Q)
    ! see http://plutocode.ph.unito.it/Doxygen/Test_Problems/_m_h_d_2_rotor_2init_8c.html
    ! Domain: -0.6 .. 0.6  (square domain)

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 
    
    REAL ::  rho, p, r
    REAL :: rho0, p0, rvel0, rho1, p1, B0, v0
    REAL :: V(nVAR), VV(3), BV(3)
    
    REAL, PARAMETER :: MHDRotomega = 10.0
    REAL :: EPCenter(2), EPRadius, TaperRadius, f
    
    EPCenter = (/ 0.0, 0.0 /)
    EPRadius    = 0.1
    TaperRadius = 0.105
    
    rho0 = 1.0
    p0   = 1.0
    rho1 = 10.0
    p1   = 1.0
    B0   = 1.0
    v0   = MHDRotomega  ! this is omega
    r = sqrt ( (x(1)-EPCenter(1))**2 + (x(2)-EPCenter(2))**2)   ! cylindrical

    IF  (r < EPRadius) THEN
        rho   =  rho1                   
        VV(1) = -v0 * x(2) ! angular velocity v_phi = w * r * e_phi , with e_phi = (-y/r, x/r )^T
        VV(2) =  v0 * x(1)
        VV(3) =  0.
        p     =  p1
    ELSE IF (r.GE.EPRadius.AND.r.LT.TaperRadius) THEN
        f = (TaperRadius - r) / (EPRadius-TaperRadius) ! linear interpolant;

        rho   = rho0 + (rho1-rho0) * f
        VV(1) = -v0 * f * x(2)
        VV(2) =  v0 * f * x(1)
        VV(3) =  0.
    ELSE
        rho   = rho0
        VV    = 0.
        p     = p0
    ENDIF
    !
    BV(1) = B0
    BV(2) = 0.
    BV(3) = 0.
    !
    V = (/ rho, VV(1:3), p, BV(1:3), 0.0 /)
    CALL PDEPrim2Cons(Q, V)
END SUBROUTINE InitialRotor



