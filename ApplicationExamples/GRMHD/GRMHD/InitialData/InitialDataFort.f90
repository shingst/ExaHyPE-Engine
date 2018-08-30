! GRMHD Initial Data


RECURSIVE SUBROUTINE InitialData(x, t, Q)
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar, nDim
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: x(nDim), t        ! 
	REAL, INTENT(OUT)              :: Q(nVar)        ! 

	! Call here one of
	! CALL InitialBlast(x, 0.0,  Q)
	Call AlfenWave(x, t, Q)
	! Call InitialRotor(x, 0.0, Q)
	! Call InitialBlast(x, 0.0, Q)
	! Call InitialOrsagTang(x, 0.0 , Q)
	
	! CALL InitialAccretionDisc(x, 0.0,  Q)
!	CALL InitialAccretionDisc3D(x, 0.0, Q)
END SUBROUTINE InitialData


RECURSIVE SUBROUTINE AlfenWave(x, t, Q)
    ! Computes the AlfenWave conserved variables (Q) at a given time t.
    ! Use it ie. with t=0 for initial data
    ! Use it for any other time ie. for comparison
    
    ! GRID FOR ALFENWAVE:
    !     dimension const                = 2
    !     width                          = 1.0, 0.3
    !     offset                         = 0.0, 0.0
    !     end-time                       = 2.1
    !
    !  maximum-mesh-size              = 0.04

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim, gamma
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 
    REAL, PARAMETER                :: t_offset = 1.0 ! time offset

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
!    BV(2) = eta * B0 * COS(2*Pi*( x(1) - vax*(t-t_offset)))
!    BV(3) = eta * B0 * SIN(2*Pi*( x(1) - vax*(t-t_offset)))
    BV(2) = eta * B0 * COS(x(1) - vax*t)
    BV(3) = eta * B0 * SIN(x(1) - vax*t)
    print *, eta, B0, x(1), vax, t , BV(2)   
    print *, eta, B0, x(1), vax, t , BV(3)   
    !
    VV(1)   = 0.0
    VV(2:3) = - vax * BV(2:3) / B0
    !
    ! Now convert to conservative variables
    !
    V(1:9) = (/ rho0, VV(1:3), p0, BV(1:3), 0.0 /) ! psi
    !            lapse, shift(3),    gamma11, ...  gamma22, -, gamma33
    V(10:19) = (/ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0 /)
    
    ! create primitive ID
    !Q = V
    CALL PDEPrim2Cons(Q,V)
END SUBROUTINE AlfenWave

RECURSIVE SUBROUTINE InitialBlast(x, t, Q)
    ! Blast wave initial data in conserved variables (Q):
    ! Simulation domain:  -6 .. +6

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
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
    V(10:19) = (/ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0 /)

    CALL PDEPrim2Cons(Q,V)
END SUBROUTINE

RECURSIVE SUBROUTINE InitialOrsagTang(x, t, Q)
    ! Orsang-Tang initial data in conserved variables (Q)
    ! Simulation Domain: 0 .. 2*pi = 6.283185307179586

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
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
    V(1:9) = (/ rho0, VV(1:3), p0, BV(1:3), 0.0 /)
    V(10:19) = (/ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0 /)

    Call PDEPrim2Cons(Q,V)
END SUBROUTINE InitialOrsagTang

RECURSIVE SUBROUTINE InitialRotor(x,t,Q)
    ! GRMHD Rotor initial data in conserved variables (Q)
    ! Needs Limiting, otherwise crash.
    ! Domain: 0.0 .. 1.0 square

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 
    
    REAL ::  rho, p, r
    REAL :: rho0, p0, rvel0, rho1, p1, B0, v0
    REAL :: V(nVAR), VV(3), BV(3)
    
    REAL, PARAMETER :: MHDRotomega = 0.95
    REAL :: EPCenter(2), EPRadius
    
    EPCenter = (/ 0.5, 0.5 /)
    EPRadius = 0.1
    
    rho0 = 1.0
    p0   = 1.0
    rho1 = 10.0
    p1   = 1.0
    B0   = 1.0
    v0   = MHDRotomega  ! this is omega
    r = sqrt ( (x(1)-EPCenter(1))**2 + (x(2)-EPCenter(2))**2)   ! cylindrical
    IF  (r < EPRadius) THEN
        rho   =  rho1                   
        VV(1) = -v0 * x(2)
        VV(2) =  v0 * x(1) 
        VV(3) =  0.
        p     =  p1
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
    V(1:9) = (/ rho, VV(1:3), p, BV(1:3), 0.0 /)
    V(10:19) = (/ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0 /)
    CALL PDEPrim2Cons(Q, V)
END SUBROUTINE InitialRotor



RECURSIVE SUBROUTINE InitialAccretionDisc(x,t,Q)
    ! 2D Accretion disk, with simulation domain:
    !    dimension const                = 2
    !    width                          = 8.5, 2.0
    !    offset                         = 1.5, 0.5
    !    end-time                       = 2.1
    ! mesh:     maximum-mesh-size              = 0.1
    ! Timestep size is around repeat = 0.000785674, for plotter.

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim, gamma
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 
    ! Local variables
    
    REAL :: rho0, p0, eta, B0, hh, tempaa, tempab, tempac, va2, vax
    REAL :: V(nVar), BV(3), VV(3), Pi = ACOS(-1.0)
    REAL :: r, zz, urc, vc2, tc, pc,tt, c1, c2, urr, f
    REAL :: df, dt, ut, LF, vr, vtheta, vphi, rho, p, VV_cov(3), g_cov(3,3), g_contr(3,3)
    REAL :: phi_c
    REAL :: gp, gm, shift(3), lapse
    
    ! PARAMETERS:
    REAL :: rhoc = 0.0625  ! Critical radius
    REAL :: rc = 8.0
    INTEGER :: MAXNEWTON = 50, iNewton
    REAL :: ng = 1.0 / (gamma-1.0)

    CALL METRIC ( x, lapse, gp, gm, shift, g_cov, g_contr,phi_c)

    
     ! The Following is for Kerr-Schild spherical coordinates
       r      = x(1)
       !
       zz     = 2.0/r               ! we are computing the solution at theta=pi/2
       !
       urc = sqrt(1.0 / (2.0*rc))
       vc2 = urc**2 / (1.0 - 3.0*urc**2)
       tc  = ng*vc2 / ((1.0 + ng)*(1.0 - ng*vc2))
       pc  = rhoc*tc
      
       c1 = urc*tc**ng*rc**2
       c2 = (1.0 + ( 1.0 + ng)*tc)**2*(1.0 - 2.0/rc+urc**2)
       !
       tt = tc
       DO iNewton = 1, MAXNEWTON  
          urr = c1 / (r**2*tt**ng)
          f   = (1.0 + (1.0 + ng)*tt)**2*(1.0 - 2.0/r + urr**2) - c2
          df  = 2.0 * (1.0 + ng)*(1.0 + (1.0 + ng)*tt)*(1.0 - 2.0/r + urr**2) - 2.0*ng*urr**2/tt*(1.0 + (1.0 + ng)*tt)**2
          dt  = -f/df
          IF (abs(dt) < 1.e-10) EXIT
          tt = tt + dt
       ENDDO
       ut     = (-zz*urr + sqrt(urr**2 - zz + 1.0))/(zz - 1.0)
       LF     = lapse*ut
       vr     = ( urr / LF + shift(1) / lapse)
       vtheta = 0.0
       vphi   = 0.0
       !
       VV(1:3) = (/ vr, vtheta, vphi /)
       ! Convert to covariant velocities
       VV_cov = MATMUL(g_cov, VV)
       !
       rho = rhoc*(tt/tc)**ng
       p   = rho*tt      

       V(1:9) = (/ rho, VV_cov(1:3), p, 0., 0., 0., 0. /)
       V(10:19) = (/ 1., 0., 0., 0., 1., 0., 0., 1., 0., 1. /)
       CALL PDEPrim2Cons(Q,V)
END SUBROUTINE InitialAccretionDisc


RECURSIVE SUBROUTINE InitialAccretionDisc3D(x,t,Q)
    ! 3D Accretion disk, with simulation domain
    !     width                          = 2.0, 2.0, 2.0
    !     offset                         = 0.0, 0.0, 0.0
    ! 

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim, gamma
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 
    ! Local variables
    
    REAL :: rho0, p0, eta, B0, hh, tempaa, tempab, tempac, va2, vax
    REAL :: V(nVar), BV(3), BV_contr(3), VV(3), Pi = ACOS(-1.0)
    REAL :: r, zz, urc, vc2, tc, pc,tt, c1, c2, urr, f
    REAL :: df, dt, ut, LF, vr, vtheta, vphi, rho, p, VV_cov(3), g_cov(3,3), g_contr(3,3)
    REAL :: gp, gm, shift(3), lapse, gammaij(6), betaru, g_tt, phi, theta, vx, vy, vz
    REAL :: phi_c
    REAL :: detgamma    
    ! PARAMETERS:
    REAL :: rhoc = 0.625  ! Critical radius
   ! REAL :: rhoc = 0.0625  ! Critical radius
    
    REAL :: rc = 8.0
    INTEGER :: MAXNEWTON = 50, iNewton
    REAL :: ng = 1.0 / (gamma-1.0)

    
     CALL METRIC ( x, lapse, gp, gm, shift, g_cov, g_contr,phi_c)

     ng     = 1.0/(gamma - 1.0)

       ! The Following is for Kerr-Schild Cartesian coordinates
       r      = SQRT( x(1)**2 + x(2)**2 + x(3)**2)
       !
       !IF ( x(1) .LT. 0.1) RETURN
       !IF ( x(2) .LT. 0.1) RETURN
       !IF ( x(3) .LT. 0.1) RETURN
       IF ( r .LT. 0.8) THEN
		! To avoid division by zero, never used for evolution or BC
		rho = 1.0 !1.0
		VV_cov(1:3) = 0.0
		p = 1.0
		BV = 0.0 
		V(1:9) = (/ rho, VV_cov(1:3), p, BV(1:3), 0. /)
		V(10:19) = (/ 1.0, 0.0,0.0,0.0, 1.0,0.0,0.0,1.0,0.0,1.0 /)
		CALL PDEPrim2Cons(Q,V)       
		RETURN
       ENDIF
       !
       theta  = ACOS( x(3)/r)
       phi    = ATAN2( x(2), x(1))
       !phi    = ACOS( x(1) / (r*SIN(theta)))
       zz     = 2.0/r   ! we are computing the solution at theta=pi/2
       betaru = zz/(1.0 + zz)
       g_tt   = zz - 1.0
       !
       urc = sqrt(1.0 / (2.0*rc))
       vc2 = urc**2 / (1.0 - 3.0*urc**2)
       tc  = ng*vc2 / ((1.0 + ng)*(1.0 - ng*vc2))
       pc  = rhoc*tc
      
       c1 = urc*tc**ng*rc**2
       c2 = (1.0 + ( 1.0 + ng)*tc)**2*(1.0 - 2.0/rc+urc**2)
       !
       tt = tc
       DO iNewton = 1, MAXNEWTON  
          urr = c1 / (r**2*tt**ng)
          f   = (1.0 + (1.0 + ng)*tt)**2*(1.0 - 2.0/r + urr**2) - c2
          df  = 2.0 * (1.0 + ng)*(1.0 + (1.0 + ng)*tt)*(1.0 - 2.0/r + urr**2) - 2.0*ng*urr**2/tt*(1.0 + (1.0 + ng)*tt)**2
          dt  = -f/df
          IF (abs(dt) < 1.e-10) EXIT
          tt = tt + dt
       ENDDO
       ut     = (-zz*urr + sqrt(urr**2 - zz + 1.0))/(zz - 1.0)
       LF     = lapse*ut
       vr     = ( urr / LF + betaru / lapse)
       vtheta = 0.0
       vphi   = 0.0
       !
       vx  = SIN(theta)*COS(phi)*vr
       vy  = SIN(theta)*SIN(phi)*vr
       vz  = COS(theta)*         vr
       !
       VV(1:3) = (/ vx, vy, vz /)
       ! Convert to covariant velocities
       VV_cov = MATMUL(g_cov, VV)
       !
       rho = rhoc*(tt/tc)**ng
       p   = rho*tt
       !
       lapse = lapse
       !
       !shift_contr = MATMUL(g_contr,shift)  !shift is controvariant. See fluxes....
       !
       gammaij(1) = g_cov(1,1)
       gammaij(2) = g_cov(1,2)
       gammaij(3) = g_cov(1,3)
       gammaij(4) = g_cov(2,2)
       gammaij(5) = g_cov(2,3)
       gammaij(6) = g_cov(3,3)

       !If hydro (Clasical) Michel accretion                                                    

       !BV(1:3) = 0.                                                                            

       ! MHD  Michel accretion                                                                  

!       detgamma= gammaij(1)*( gammaij(4)*gammaij(6)-gammaij(5)*gammaij(5)) &
!               - gammaij(2)*( gammaij(2)*gammaij(6)-gammaij(5)*gammaij(3)) &
!               + gammaij(3)*( gammaij(2)*gammaij(5)-gammaij(3)*gammaij(4))
!       detgamma = sqrt(detgamma)




       ! B0 = (2.2688/M)*(sqrt(b^2/rho0)), here sqrt(b^2/rho0)=4.0                              

       B0 = 2.2688*2.0

       !---------------------------------------------------------------------- 
       BV_contr = B0*x(1:3)/(sqrt(gp) * r*r*r)

!       BV = BV_contr
       
       BV = MATMUL(g_cov,BV_contr)
       

       V(1:9) = (/ rho, VV_cov(1:3), p, BV(1:3), 0. /)
       V(10:19) = (/ lapse, shift(1:3), gammaij(1:6) /)
       CALL PDEPrim2Cons(Q,V)

 END SUBROUTINE InitialAccretionDisc3D
