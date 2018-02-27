! DIM Initial Data


RECURSIVE SUBROUTINE InitialData(x, t, Q)
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar, nDim
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: x(nDim), t        ! 
	REAL, INTENT(OUT)              :: Q(nVar)        ! 
		
	call CurvedGaussian(x, t, Q);
	!call SmoothShock(x, t, Q);
END SUBROUTINE InitialData

RECURSIVE SUBROUTINE PDElimitervalue(limiter_value,xx,numberOfObservables, observablesMin, observablesMax)
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar, nDim
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: xx(nDim)        ! 
	INTEGER, INTENT(IN)					:: numberOfObservables
	INTEGER, INTENT(OUT)              :: limiter_value        !
	REAL, INTENT(IN)					:: observablesMin(numberOfObservables), observablesMax(numberOfObservables)
	real	:: rr	

	limiter_value=0
	return
	if(observablesMin(1).ge.1.1 .and. observablesMax(1).le. 2.6) then
		limiter_value=1
	end if
END SUBROUTINE PDElimitervalue

RECURSIVE SUBROUTINE CurvedGaussian(x, t, Q)
    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim
	USE GRGPRmod
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 

	INTEGER	:: iErr
    REAL    :: xp,yp,r,phi,Jac(3,3),kAB(3,3),rho,p,VV(3),ICA, up(nVar), Pi = ACOS(-1.0)
	REAL	:: invJac(3,3),ICsigv(nDim),ICx0(nDim)

	! Parameters for the test case
	ICA=-1.e-3
	ICsigv(1)=0.25
	ICsigv(2)=0.25
	ICx0(1)=3.0
	ICx0(2)=0.0
    ! Initialize parameters
	xp = x(1)*COS(x(2)) 
	yp = x(1)*SIN(x(2)) 
	r   = x(1)                      ! radius 
	phi = x(2)                      ! angle
	Jac(1,:) = (/ COS(phi), -r*SIN(phi), 0.0 /) 
	Jac(2,:) = (/ SIN(phi), +r*COS(phi), 0.0 /)  
	Jac(3,:) = (/ 0.0,       0.0       , 1.0 /)
	CALL MatrixInverse3x3(Jac,invJac,p) 
	kAB = MATMUL( TRANSPOSE(Jac), Jac )
	rho = 1.0
	p   = 1.0
	VV(1:3)  = 0.0

	!IF(ICA.GT.0.0) THEN
	!	VV(3)    = ICA*EXP( - 0.5*( (xp-ICx0(1))**2/ICsigv(1)**2 + (yp-ICx0(2))**2/ICsigv(2)**2  ) ) 
	!ELSE
	   r =  (xp-ICx0(1))**2/ICsigv(1)**2 + (yp-ICx0(2))**2/ICsigv(2)**2 - 1.0  
	   IF( r.LT.0.0) THEN
		   VV(3) = ABS(ICA) 
	   ELSE
		   VV(3) = 0.0 
	   ENDIF           
	!ENDIF  	

       up(1:5)  = (/ rho, VV(1:3), p /) 
       up(6:14) = (/ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 /)    ! A 
       up(15)   = 1.0                                                  ! lapse 
       IF(ICA.GT.0.0) THEN
          up(16:18) = 0.0                                              ! shift 
       ELSE
          up(16:18) = 0.0  
          up(17)    = 2*Pi/10.0                                    ! angular shift vector to rotate the coordinate system 
       ENDIF       
       up(19:24) = (/ kAB(1,1), kAB(1,2), kAB(1,3), kAB(2,2), kAB(2,3), kAB(3,3) /)     ! spatial part of the space-time metric 
       up(25:30) = (/ kAB(1,1), kAB(1,2), kAB(1,3), kAB(2,2), kAB(2,3), kAB(3,3) /)     ! matter metric 
       !
       CALL PDEPrim2Cons(Q,up) 
	
END SUBROUTINE CurvedGaussian

! RECURSIVE SUBROUTINE ShuVortex2D(x, t, Q)
    ! USE, INTRINSIC :: ISO_C_BINDING
    ! USE Parameters, ONLY : nVar, nDim,gamma
    ! IMPLICIT NONE 
    ! ! Argument list 
    ! REAL, INTENT(IN)               :: t
    ! REAL, INTENT(IN)               :: x(nDim)        ! 
    ! REAL, INTENT(OUT)              :: Q(nVar)        ! 

	! INTEGER	:: iErr
    ! REAL    :: up(nVar), Pi = ACOS(-1.0), epsilon, r, du, dv, dT, drho, dp

    ! ! Initialize parameters
       ! epsilon = 5.0
       ! r = SQRT((x(1)-t-5.)**2+(x(2)-t-5.)**2)
       ! du = epsilon/2./Pi*exp(0.5*(1.-r*r))*(5. - x(2) + t)
       ! dv = epsilon/2./Pi*exp(0.5*(1.-r*r))*(x(1)  - 5.- t)
       ! dT = -(gamma-1.)*epsilon**2/8./gamma/Pi**2*exp(1.-r*r)
       ! drho = (1.+dT)**(1./(gamma-1.))-1.
       ! dp   = (1.+dT)**(gamma/(gamma-1.))-1.
       ! !
       ! up = 0.0 
       ! ! Density 
       ! up(1) = 1. + drho 
       ! ! Velocity    
       ! up(2) = 1.  + du
       ! up(3) = 1.  + dv
       ! up(4) = 0.
       ! ! Pressure  
       ! up(5) = 1.  + dp
       ! ! distortion tensor  
       ! up(6:14) = up(1)**(1./3.)*(/ 1., 0., 0., 0., 1., 0., 0., 0., 1. /)   
       ! ! thermal impulse
       ! up(15:17) = 0.0 
       ! !
       ! CALL PDEPrim2Cons(Q,up)
    ! !Q=up
! END SUBROUTINE ShuVortex2D

! RECURSIVE SUBROUTINE SmoothShock(x, t, Q)
    ! USE, INTRINSIC :: ISO_C_BINDING
    ! USE Parameters, ONLY : nVar, nDim,rho0,cs,tau1,gamma
    ! IMPLICIT NONE 
    ! ! Argument list 
    ! REAL, INTENT(IN)               :: t
    ! REAL, INTENT(IN)               :: x(nDim)        ! 
    ! REAL, INTENT(OUT)              :: Q(nVar)        ! 

	! INTEGER	:: iErr
    ! REAL    :: up(nVar), Pi = ACOS(-1.0)
    ! REAL    :: xi, ICQR(2), ICxd, ICsig
    ! REAL    :: Re0,Ms,mu,x0(nDim),nv(nDim),u0,p0, MyAref, vx
    ! ! Initialize parameters
		! mu=1.0/6.0*rho0*cs**2*tau1
		! ms=2.0
		! x0=0.0 
		! nv=0.0
		! nv(1) = 1.0 
        ! up(:) = 0. 
		! u0=0.0
		! p0=1.0/gamma
        ! Re0 = rho0 * Ms / mu 
        ! CALL NSShock(up(1),vx,up(5),DOT_PRODUCT(x(:)-ms*t-x0(:),nv),gamma,mu,Re0,Ms,rho0,u0,p0)
        ! up(4) = 0.0
		! up(2:1+nDim) = vx*nv  
        ! MyAref = up(1)**(1./3.) 
		! up(6)=MyAref
		! up(7)=0.0
		! up(8)=0.0
		! up(9)=0.0
		! up(10)=MyAref
		! up(11)=0.0
		! up(12)=0.0
		! up(13)=0.0
		! up(14)=MyAref 
        ! CALL PDEPrim2Cons(Q,up)
    ! !Q=up
! END SUBROUTINE SmoothShock

! RECURSIVE SUBROUTINE ShearLayer(x, t, Q)
    ! USE, INTRINSIC :: ISO_C_BINDING
    ! USE Parameters, ONLY : nVar, nDim
    ! IMPLICIT NONE 
    ! ! Argument list 
    ! REAL, INTENT(IN)               :: t
    ! REAL, INTENT(IN)               :: x(nDim)        ! 
    ! REAL, INTENT(OUT)              :: Q(nVar)        ! 

	! INTEGER	:: iErr
    ! REAL    :: up(nVar), Pi = ACOS(-1.0)
    ! REAL    :: xi, ICQR(2), ICxd, ICsig
    ! REAL    :: ICuL(nVar), ICuR(nVar),ICA,r, ICx0(3),nv(3)
    ! ! Initialize parameters

        ! up = 0.0
        ! ! Define the material properties
        ! if(x(1) .lt. 0.0) then
            ! up(1) = 1.0  
            ! up(2:4) = (/0.0, 0.1, 0.0/)  
            ! up(5) = 1.0
			! up(6:14) = (/1.0, 0.0, 0.0, 0.0, 1.0, 0.0,0.0, 0.0, 1.0/)
			! up(15:17) = (/0.0, 0.0, 0.0/) 
        ! else
            ! up(1) = 1.0  
            ! up(2:4) = (/0.0, -0.1, 0.0/)  
            ! up(5) = 1.0
			! up(6:14) = (/1.0, 0.0, 0.0, 0.0, 1.0, 0.0,0.0, 0.0, 1.0/)
			! up(15:17) = (/0.0, 0.0, 0.0/) 
        ! end if
        
		
    ! CALL PDEPrim2Cons(Q,up)
    ! !Q=up
! END SUBROUTINE ShearLayer

! RECURSIVE SUBROUTINE NSShock(rho,u,p,x,gamma,mu,Re0,M,rho0,u0,p0) 
      ! IMPLICIT NONE      
      ! ! Argument list declaration
      ! REAL :: x, gamma, mu, Re0, M, rho0, u0, p0
      ! REAL :: rho, u, p 
      ! ! Local variable declarations
      ! INTEGER :: i 
      ! REAL :: lambda,M2,xf 
      ! REAL :: u1,u2,f1,f2,ui,fi,ubar,pbar,mf,dudx 
      ! REAL, PARAMETER :: tol = 1e-12 
	  ! ! Shock parameters
	  ! INTEGER :: ShockType=1
      ! !
      ! M2     = M*M 
      ! lambda = ( 1 + 0.5*(gamma-1)*M2 ) / (0.5*(gamma+1)*M2 ) 
      ! !
      ! ! Use the regula falsi in order to find the root of the transcendental velocity equation 
      ! !
      ! u1 = lambda
      ! u2 = 1 
      ! !
      ! IF(ShockType.EQ.0) THEN
        ! xf = x
      ! ELSE
        ! xf = -x
      ! ENDIF      
      ! !
      ! CALL NSSFunction(f1,u1,xf,gamma,lambda,M2,Re0) 
      ! CALL NSSFunction(f2,u2,xf,gamma,lambda,M2,Re0) 
      ! !
      ! IF(f1*f2.GT.0.) THEN 
        ! PRINT *, ' ERROR: Sign does not change in NSSFunction. '
        ! PRINT *, ' Info: ', u0,u1,x,lambda,M2 
        ! PRINT *, '       ', f1,f2 
        ! STOP  
      ! ENDIF
      ! !mf = MIN( ABS(f1), ABS(f2) )
      ! !
      ! DO i = 1, 100   
          ! IF(ABS(f1).LT.tol) THEN
              ! ubar  = u1 
              ! pbar  = 1.
              ! EXIT 
          ! ENDIF
          ! IF(ABS(f2).LT.tol) THEN
              ! ubar  = u2 
              ! pbar  = 1. 
              ! EXIT
          ! ENDIF          
          ! ui = (f1*u2 - f2*u1)/(f1-f2)
          ! CALL NSSFunction(fi,ui,xf,gamma,lambda,M2,Re0) 
          ! IF(f1*fi.GT.0.) THEN
            ! u1 = ui
            ! f1 = fi 
          ! ELSE
            ! u2 = ui 
            ! f2 = fi
          ! ENDIF          
      ! ENDDO 
      ! ! 
      ! ! The function is very stiff. If 100 iterations do not lead to convergence, then
      ! ! simply take the best values you can get... 
      ! !
      ! IF(ABS(f1).LT.ABS(f2)) THEN
          ! ubar  = u1           
      ! ELSE
          ! ubar  = u2 
      ! ENDIF          
      ! !
      ! dudx = 0.5*(gamma+1)*(ubar-1.)*(ubar-lambda)/( 4./3./Re0*gamma*ubar )
      ! pbar = 1 - ubar +  4./3./Re0*dudx 
      ! !
      ! rho   = rho0 / ubar    
      ! IF(ShockType.EQ.0) THEN
         ! u     = M*SQRT(gamma*p0/rho0) * ubar 
      ! ELSE
         ! u     = M*SQRT(gamma*p0/rho0) * ( 1. - ubar ) 
      ! ENDIF
      ! p     = p0 + pbar*M**2*gamma*p0 
      ! !
  ! END SUBROUTINE NSShock
  ! !
! RECURSIVE SUBROUTINE NSSFunction(f,u,x,gamma,lambda,M2,Re0)
      ! IMPLICIT NONE
      ! ! Argument list declaration
      ! REAL :: f,u,x,Re0,M2,gamma,lambda  
      ! ! Local variable declaration
      ! REAL :: RHS, MYexponent 
      ! !
      ! MYexponent = MIN( 0.75*Re0*(M2-1)/gamma/M2 * x, 10. ) 
      ! RHS      = ABS( 0.5*(1-lambda) )**(1-lambda) * EXP( MYexponent ) 
      ! f        = ABS(u-1) - ABS(u-lambda)**lambda * RHS 
      ! !
  ! END SUBROUTINE  
  


