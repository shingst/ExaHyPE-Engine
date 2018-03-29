! Initial Data for benchmarking the Z4 formulation


!!!!! THE GAUGE WAVE HAS NOT BEEN PORTED YET TOT CCZ4
RECURSIVE SUBROUTINE Z4_GaugeWaveTest(xGp, u0, tGP)
	USE, INTRINSIC :: ISO_C_BINDING
	USE typesDef, ONLY : nVar, nDim
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: xGP(nDim), tGP      ! 
	REAL, INTENT(OUT)              :: u0(nVar)        ! 
	! Local variables
	REAL :: ICA, Pi
	REAL :: HH, dXh, Kxx, V0(nVar), dxphi, phi, traceK
	Pi = ACOS(-1.0)
	ICA = 0.1
	
	V0 = 0.0
        !
        ICA = 0.1 
        ! Gauge wave propagating along x 
        HH     = 1.0 - ICA*SIN( 2.0*Pi*( xGP(1) - tGP)) 
        dxH    = -2.0*Pi*ICA*COS( 2.0 * Pi*(xGP(1) - tGP)) 
        dxphi  = - HH**(-7.0/6.0)*dxH/6.0
        phi    = ( 1.0 / HH)**(1.0/6.0)
        Kxx    = - Pi*ICA*COS( 2.0 * Pi*(xGP(1) - tGP))/SQRT( 1.0 - ICA*SIN( 2.0*Pi*( xGP(1) - tGP))  )  ! extrinsic curvature  
        traceK = Kxx/HH
        !
!        HH = 1.0 
!        phi = 1.0 
!        traceK = 0.0
!        Kxx = 0.0 
!        traceK = 0.0 
!        dxH = 0.0
!        dxphi = 0.0 
        !
        V0(:)  = 0.0
        V0(1)  = phi**2*HH                          ! \tilde\gamma_xx
        V0(4)  = phi**2                             ! \tilde\gamma_yy
        V0(6)  = phi**2                             ! \tilde\gamma_zz
        V0(7)  = phi**2*(Kxx - 1.0/3.0*traceK*HH )  ! \tilde A_{xx}
        V0(10) = phi**2*(0.0 - 1.0/3.0*traceK*1.0)  ! \tilde A_{yy}
        V0(12) = phi**2*(0.0 - 1.0/3.0*traceK*1.0)  ! \tilde A_[zz}
        V0(17) = SQRT(HH)
        V0(14) = 2.0/(3.0*HH**(5.0/3.0))*dxH        ! Gtilde 
        !
        ! Auxiliary variables
        V0(24) = 1.0/(2.0*HH)*dxH                ! A_x  
        V0(36) = HH**(-1.0/3.0)*dxH/3.0          ! D_xxx
        V0(39) = phi*dxphi                       ! D_xyy
        V0(41) = phi*dxphi                       ! D_xzz
        ! 
        V0(54) = traceK
        V0(55) = phi
        V0(56) = dxphi/phi                       ! P_x
 
	CALL PDEPrim2Cons(u0,V0) 

END SUBROUTINE Z4_GaugeWaveTest

RECURSIVE SUBROUTINE CCZ4_MinkowskiRobustScability(x, u0, t)
	USE, INTRINSIC :: ISO_C_BINDING
	USE typesDef, ONLY : nVar, nDim
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: x(nDim), t      ! 
	REAL, INTENT(OUT)              :: u0(nVar)        ! 
	! local variables
	REAL :: V0(nVar)

	! Alle Quantities random stoeren
	CALL RANDOM_NUMBER(V0)  ! Fortran intrinsic
	V0 = V0*1e-4   
	! 
	V0(1)  = V0(1)  + 1.0 
	V0(4)  = V0(4)  + 1.0 
	V0(6)  = V0(6)  + 1.0 
	V0(17) = V0(17) + 1.0 
	V0(55) = V0(55) + 1.0 
	!
	! V0(13) = 0.1*EXP(-0.5*( xGP(1)**2+xGP(2)**2)/0.1**2 ) 
	
	CALL PDEPrim2Cons(u0,V0) 

END SUBROUTINE CCZ4_MinkowskiRobustScability

