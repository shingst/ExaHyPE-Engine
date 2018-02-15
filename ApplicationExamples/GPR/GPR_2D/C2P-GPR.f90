! The Con2Prim and Prim2Con routines for MHD.
! Should be merged with SRHD's.

SUBROUTINE PDEPrim2Cons(Q,V)
  USE Parameters, ONLY: nVar, nDim, alpha, cs, gamma, p0, cv 
  IMPLICIT NONE
  ! Argument list declaration
  REAL :: Q(nVar), V(nVar)
  INTENT(IN)  :: V
  INTENT(OUT) :: Q 
  ! Local variable declaration
		REAL :: eh, A(3,3), G(3,3), Id(3,3), devG(3,3), temp(3,3), evv, falpha, p, T
  !
	! Peshkov-Romenski model in 3D with heat conduction 
	Q(1)   = V(1)        ! rho 
	Q(2:4) = V(1)*V(2:4) ! rho*u, rho*v, rho*w  
	Q(6:14) = V(6:14)    ! A11, A12, A13, A21, A22, A23, A31, A32, A33 
	eh     = (V(5)+gamma*p0)/(gamma-1)/V(1) 
	!A(1,:) = (/ V(6),  V(7),  V(8)  /) 
	!A(2,:) = (/ V(9),  V(10), V(11) /)
	!A(3,:) = (/ V(12), V(13), V(14) /) 
	A(1,1)=V(6)
	A(1,2)=V(7)
	A(1,3)=V(8)
	A(2,1)=V(9)
	A(2,2)=V(10)
	A(2,3)=V(11)
	A(3,1)=V(12)
	A(3,2)=V(13)
	A(3,3)=V(14)
	
	G      = MATMUL( TRANSPOSE(A), A ) 
	Id     = 0.0 
	Id(1,1) = 1.0
	Id(2,2) = 1.0
	Id(3,3) = 1.0 
	devG   = G - (G(1,1)+G(2,2)+G(3,3))/3.*Id  
	temp   = MATMUL( TRANSPOSE(devG), devG ) 
	evv     = cs**2/4.*(temp(1,1)+temp(2,2)+temp(3,3)) 
	! Compute the temperature from the ideal gas law 
	p = V(5) 
	T = p/V(1)/cv/(gamma-1)  
	falpha  = alpha**2 

	Q(5)     = V(1)*(eh+evv) + 0.5*V(1)*(V(2)**2+V(3)**2+V(4)**2) + 0.5*V(1)*falpha*(V(15)**2 + V(16)**2 + V(17)**2) ! total energy rhoE     
	Q(15:17) = V(15:17)*V(1) ! rho j
END SUBROUTINE PDEPrim2Cons

SUBROUTINE PDECons2Prim(V,Q)
  USE Parameters, ONLY: nVar, nDim, alpha, cs, gamma, p0
  IMPLICIT NONE
  !--------------------------------------------!
  ! Argument list declaration
  REAL :: Q(nVar), V(nVar)
  INTENT(IN)  :: Q 
  INTENT(OUT) :: V 
  ! Local variables
  REAL :: e, A(3,3), G(3,3), Id(3,3), devG(3,3), tempp(3,3), evv ,ehh
  ! Peshkov-Romenski model with heat conduction 
	V(1)   = Q(1)        ! rho 
	V(2:4) = Q(2:4)/Q(1) ! u, v 
	V(6:14) = Q(6:14)      ! A11, A12, A13, A21, A22, A23, A31, A32, A33 
	V(15:17) = Q(15:17)/Q(1) ! j1, j2, j3  
	e      = Q(5)/Q(1) - 0.5*(V(2)**2+V(3)**2+V(4)**2) -0.5*(V(15)**2+V(16)**2+V(17)**2)*alpha**2     ! e = eh + ev 
	!A(1,:) = (/ V( 6), V( 7),  V( 8) /) 
	!A(2,:) = (/ V( 9), V(10),  V(11) /)
	!A(3,:) = (/ V(12), V(13),  V(14) /) 
	A(1,1)=V(6)
	A(1,2)=V(7)
	A(1,3)=V(8)
	A(2,1)=V(9)
	A(2,2)=V(10)
	A(2,3)=V(11)
	A(3,1)=V(12)
	A(3,2)=V(13)
	A(3,3)=V(14)	
	G      = MATMUL( TRANSPOSE(A), A ) 
	Id     = 0. 
	Id(1,1) = 1.0
	Id(2,2) = 1.0
	Id(3,3) = 1.0 
	devG   = G - (G(1,1)+G(2,2)+G(3,3))/3.*Id  
	tempp  = MATMUL( TRANSPOSE(devG), devG ) 
	evv    = cs**2/4.*(tempp(1,1)+tempp(2,2)+tempp(3,3))         
	ehh    = e-evv  
	V(5)   = ehh*(gamma-1.0)*V(1) - gamma*p0 
END SUBROUTINE PDECons2Prim

