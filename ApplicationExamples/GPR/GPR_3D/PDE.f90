! GRMHD PDE.f90
! Trento (EQNTYPE4)

RECURSIVE SUBROUTINE PDEFlux(f,g,hz,Q)
  USE Parameters, ONLY : nVar, nDim, rho0, cs, cv, gamma, alpha
  USE iso_c_binding
  IMPLICIT NONE
  REAL :: f(nVar), g(nVar), h(nVar), hz(nVar), Q(nVar), V(nVar)
  REAL, PARAMETER :: epsilon = 1e-14 
  INTENT(IN)  :: Q
  INTENT(OUT) :: f, g, hz
  ! Local varialbes
  REAL :: p, A(3,3), AU(3), detA, GT(3,3), devG(3,3), TT(3,3), Id(3,3), T,falpha
	INTEGER :: iErr

  
  !FTensDim = 0
  !RETURN
  
 f = 0. 
  g = 0. 
  h = 0. 
  !
  !  1   2    3    4    5    6   7   8   9   10  11  12  13  14    15    16     17  
  ! rho rhou rhov rhow rhoE A11 A12 A13 A21 A22 A23 A31 A32 A33  rho j1 rho j2 rho j3 
  ! rho u    v    w    p    A11 A12 A13 A21 A22 A23 A31 A32 A33      j1     j2     j3 
  CALL PDECons2Prim(V,Q)  
  ! Compute the pressure (hydrodynamic part) 
  p = V(5) 
  ! Compute the viscous stress tensor 
  !A(1,:) = (/ V( 6), V( 7), V( 8) /) 
  !A(2,:) = (/ V( 9), V(10), V(11) /)
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
  AU = MATMUL( A, V(2:4) ) 
  detA   = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(2,1)*A(1,2)*A(3,3)+A(2,1)*A(1,3)*A(3,2)+A(3,1)*A(1,2)*A(2,3)-A(3,1)*A(1,3)*A(2,2) 
  ! 
  !        
  GT     = MATMUL( TRANSPOSE(A), A ) 
  Id     = 0. 
  Id(1,1) = 1.0
  Id(2,2) = 1.0
  Id(3,3) = 1.0 
  devG   = GT - (GT(1,1)+GT(2,2)+GT(3,3))/3.*Id
  TT      = -rho0*detA*cs**2*MATMUL(GT,devG) 
  ! Compute the temperature from the ideal gas law 
  T = p/V(1)/cv/(gamma-1)  
  falpha  = alpha**2 
  ! 
  f(1) = Q(2) 
  f(2) = Q(2)*V(2) + p - TT(1,1) 
  f(3) = Q(3)*V(2)     - TT(2,1) 
  f(4) = Q(4)*V(2)     - TT(3,1) 
  f(5) = V(2)*(Q(5)+p) - V(2)*TT(1,1) - V(3)*TT(2,1) - V(4)*TT(3,1) + falpha*V(15)*T   
  f(6)  = AU(1) 
  f(9)  = AU(2) 
  f(12) = AU(3) 
  f(15) = Q(15)*V(2) + T
  f(16) = Q(16)*V(2) 
  f(17) = Q(17)*V(2) 
  !
  g(1) = Q(3) 
  g(2) = Q(2)*V(3)      - TT(1,2) 
  g(3) = Q(3)*V(3) + p  - TT(2,2) 
  g(4) = Q(4)*V(3)      - TT(3,2) 
  g(5) = V(3)*(Q(5)+p) - V(2)*TT(1,2) - V(3)*TT(2,2) - V(4)*TT(3,2) + falpha*V(16)*T 
  g(7)  = AU(1) 
  g(10) = AU(2) 
  g(13) = AU(3) 
  g(15) = Q(15)*V(3) 
  g(16) = Q(16)*V(3) + T 
  g(17) = Q(17)*V(3)   
  !
  h(1) = Q(4) 
  h(2) = Q(2)*V(4)      - TT(1,3) 
  h(3) = Q(3)*V(4)      - TT(2,3) 
  h(4) = Q(4)*V(4) + p  - TT(3,3) 
  h(5) = V(4)*(Q(5)+p) - V(2)*TT(1,3) - V(3)*TT(2,3) - V(4)*TT(3,3) + falpha*V(17)*T   
  h(8)  = AU(1) 
  h(11) = AU(2) 
  h(14) = AU(3) 
  h(15) = Q(15)*V(4) 
  h(16) = Q(16)*V(4) 
  h(17) = Q(17)*V(4) + T 
  
  IF(nDim==3) THEN
    hz=h
  END IF  
  END SUBROUTINE PDEFlux


RECURSIVE SUBROUTINE PDENCP(BgradQ,Q,gradQ) 
   USE Parameters, ONLY :  nVar, nDim
   IMPLICIT NONE
   ! 11. Oct 21:40: This was a matrix BGradQ(nVar, nDim) but is a vector in spaceTimePredictorNonlinear
   REAL, INTENT(OUT) :: BgradQ(nVar)
   REAL, INTENT(IN)  :: gradQ(nVar, 3)
   REAL, INTENT(IN)  :: Q(nVar)
  ! Linear elasticity variables
   REAL :: u(3),VP(nVar) 
   REAL :: AQx(nVar), BQy(nVar), CQz(nVar) , Qx(nVar), Qy(nVar), Qz(nVar)
   
   
  ! BgradQ = 0
  ! RETURN
  Qx = gradQ(:,1)
  Qy = gradQ(:,2)
  IF(nDim==3) THEN
	Qz = gradQ(:,3)
  ELSE
	Qz = 0.0 
  ENDIF 
  !
    ! GPR model
   AQx = 0.   
   BQy = 0. 
   CQz = 0. 
      
    u = Q(2:4)/Q(1)  
    ! 
    AQx(6) =  -u(2)*Qx(7) - u(3)*Qx(8) 
    BQy(6) =  +u(2)*Qy(6) 
    CQz(6) =  +u(3)*Qz(6)  
    ! 
    AQx(7) =  +u(1)*Qx(7) 
    BQy(7) =  -u(1)*Qy(6) - u(3)*Qy(8) 
    CQz(7) =  +u(3)*Qz(7) 
    !
    AQx(8) =  +u(1)*Qx(8) 
    BQy(8) =  +u(2)*Qy(8) 
    CQz(8) =  -u(1)*Qz(6) - u(2)*Qz(7) 
    ! 
    AQx(9) =  -u(2)*Qx(10) - u(3)*Qx(11) 
    BQy(9) =  +u(2)*Qy(9) 
    CQz(9) =  +u(3)*Qz(9) 
    !
    AQx(10) = +u(1)*Qx(10)     
    BQy(10) = -u(1)*Qy(9) - u(3)*Qy(11) 
    CQz(10) = +u(3)*Qz(10) 
    !
    AQx(11) = +u(1)*Qx(11) 
    BQy(11) = +u(2)*Qy(11) 
    CQz(11) = -u(1)*Qz(9) - u(2)*Qz(10) 
    !
    AQx(12) = -u(2)*Qx(13) - u(3)*Qx(14)  
    BQy(12) = +u(2)*Qy(12)
    CQz(12) = +u(3)*Qz(12) 
    !
    AQx(13) = +u(1)*Qx(13)  
    BQy(13) = -u(1)*Qy(12) - u(3)*Qy(14) 
    CQz(13) = +u(3)*Qz(13) 
    !
    AQx(14) = +u(1)*Qx(14) 
    BQy(14) = +u(2)*Qy(14) 
    CQz(14) = -u(1)*Qz(12) - u(2)*Qz(13) 

    if( nDim .eq. 2) then
        BgradQ = AQx + BQy         
    else
        BgradQ = AQx + BQy + CQz     
    end if
END SUBROUTINE PDENCP


RECURSIVE SUBROUTINE PDEEigenvalues(L,Q,n)
  USE Parameters, ONLY :  nVar, nDim, p0, cs, gamma, alpha
  USE iso_c_binding
  IMPLICIT NONE
  REAL :: L(nVar), n(3), Q(nVar), Vp(nVar)
  INTENT(IN)  :: Q,n
  INTENT(OUT) :: L 
  ! Local Variables 
  REAL :: p,c, u

  L(:) = 0
  
	CALL PDECons2Prim(Vp,Q)  
    
	p = Vp(5) 
    c = SQRT(gamma*(p+p0)/Vp(1) + 4./3.*cs**2+alpha**2)  
    IF(SUM(n**2).EQ.0.) THEN  
        u = SQRT(Vp(2)**2 + Vp(3)**2 + Vp(4)**2 ) 
    ELSE
        u = Vp(2)*n(1) + Vp(3)*n(2) + Vp(4)*n(3)  
    ENDIF   
    ! 
    L = 0. 
    L(1)  = u-c 
    L(2)  = u 
    L(3)  = u+c  

END SUBROUTINE PDEEigenvalues

RECURSIVE SUBROUTINE PDESource(S,Q) 
  USE Parameters, ONLY:  nVar, nDim, rho0, tau1, tau2, cv, gamma 
  USE iso_c_binding
  IMPLICIT NONE
  ! --------------------------------------------
  ! Argument list declaration
  REAL :: S(nvar), Q(nvar)
  INTENT(IN)  :: Q 
  INTENT(OUT) :: S
  ! --------------------------------------------
  ! Local variable declaration 
  REAL :: AM(3,3), Id(3,3), G(3,3) ,devG(3,3), detA, detA2,psiM(3,3), temp2, T
  REAL :: V(nvar)

	S(1:5)= 0. 
	CALL PDECons2Prim(V,Q) 
	!AM(1,:) = (/ V( 6), V( 7), V( 8) /) 
	!AM(2,:) = (/ V( 9), V(10), V(11) /)
	!AM(3,:) = (/ V(12), V(13), V(14) /) 
	AM(1,1)=V(6)
	AM(1,2)=V(7)
	AM(1,3)=V(8)
	AM(2,1)=V(9)
	AM(2,2)=V(10)
	AM(2,3)=V(11)
	AM(3,1)=V(12)
	AM(3,2)=V(13)
	AM(3,3)=V(14)
	G      = MATMUL( TRANSPOSE(AM), AM ) 
	Id     = 0. 
	Id(1,1) = 1.0; Id(2,2) = 1.0; Id(3,3) = 1.0 
	devG   = G - (G(1,1)+G(2,2)+G(3,3))/3.*Id  
	detA   = AM(1,1)*AM(2,2)*AM(3,3)-AM(1,1)*AM(2,3)*AM(3,2)-AM(2,1)*AM(1,2)*AM(3,3)+AM(2,1)*AM(1,3)*AM(3,2)+AM(3,1)*AM(1,2)*AM(2,3)-AM(3,1)*AM(1,3)*AM(2,2)     ! this is the determinant we have     
	detA2  = Q(1)/rho0                                ! this is the determinant we should have from the compatibility relation 
	psiM   = 3./(detA)*MATMUL(AM,devG)                    ! in the coefficient of the source term, we use the detA2 from the compatibility relation, to reduce nonlinearities in A 
	temp2  = detA2**(1./3.) 
	!S(6:14) = -(/ psiM(1,1), psiM(1,2), psiM(1,3), & 
	!			  psiM(2,1), psiM(2,2), psiM(2,3), & 
	!			  psiM(3,1), psiM(3,2), psiM(3,3) /) /tau1*detA**(8./3.)   &      ! real relaxation term 
	!		  -(detA-detA2)/tau1*Q(6:14)                                          ! artificial relaxation term to correct the errors of the ODE integrator 
    S(6)=-psiM(1,1)/tau1*detA**(8./3.)-(detA-detA2)/tau1*Q(6)
	S(7)=-psiM(1,2)/tau1*detA**(8./3.)-(detA-detA2)/tau1*Q(7)
	S(8)=-psiM(1,3)/tau1*detA**(8./3.)-(detA-detA2)/tau1*Q(8)
	S(9)=-psiM(2,1)/tau1*detA**(8./3.)-(detA-detA2)/tau1*Q(9)
	S(10)=-psiM(2,2)/tau1*detA**(8./3.)-(detA-detA2)/tau1*Q(10)
	S(11)=-psiM(2,3)/tau1*detA**(8./3.)-(detA-detA2)/tau1*Q(11)
	S(12)=-psiM(3,1)/tau1*detA**(8./3.)-(detA-detA2)/tau1*Q(12)
	S(13)=-psiM(3,2)/tau1*detA**(8./3.)-(detA-detA2)/tau1*Q(13)
	S(14)=-psiM(3,3)/tau1*detA**(8./3.)-(detA-detA2)/tau1*Q(14)
	T = V(5)/V(1)/cv/(gamma-1)
	S(15:17) = -Q(15:17)/tau2 * T/V(1) ! the source term for the heat flux is very nice and simple    
      
END SUBROUTINE PDESource

RECURSIVE SUBROUTINE PDEVarName(Name) 
  USE Parameters, ONLY: nVar  
  IMPLICIT NONE     
  CHARACTER(LEN=10):: Name(nVar)

  ! EQNTYPE93
    Name(1)  = 'rho'
    Name(2)  = 'u'
    Name(3)  = 'v'
    Name(4)  = 'w'
    Name(5)  = 'p'
    Name(6)  = 'A11'
    Name(7)  = 'A12'
    Name(8)  = 'A13'
    Name(9)  = 'A21'
    Name(10) = 'A22'
    Name(11) = 'A23'
    Name(12) = 'A31'
    Name(13) = 'A32'
    Name(14) = 'A33'
    Name(15) = 'j1'
    Name(16) = 'j2'
    Name(17) = 'j3'
    END SUBROUTINE PDEVarName

RECURSIVE SUBROUTINE PDEMatrixB(An,Q,nv) 
  USE Parameters, ONLY : nVar, nDim
  USE iso_c_binding
  IMPLICIT NONE
  ! Argument list 
  REAL :: An(nVar,nVar)
  REAL :: Q(nVar), nv(3) 
  INTENT(IN)  :: Q,nv
  INTENT(OUT) :: An  
  ! Local variables
    ! Linear elasticity variables
   REAL :: A(nVar,nVar), B(nVar,nVar), C(nVar,nVar), Vp(nVar)
   
   PRINT *, ' Impossible error! ' 
   
   !An = 0
  !RETURN
    
	!print *, maxval(nv),lam,mu,irho

	!CALL PDECons2Prim(Vp,Q)
	VP = 0.0
	VP(2:4) = Q(2:4)/Q(1) 
	
    A = 0.0
    B = 0.0
    C = 0.0 

    A(6,7) = -VP(3)
    A(6,8) = -VP(4)
    B(6,6) = +VP(3)
    C(6,6) = +VP(4) 
    ! 
    A(7,7) = +VP(2)
    B(7,6) = -VP(2)
    B(7,8) = -VP(4)
    C(7,7) = +VP(4)
    !
    A(8,8) = +VP(2)
    B(8,8) = +VP(3)
    C(8,6) = -VP(2)
    C(8,7) = -VP(3)
    ! 
    A(9,10) = -VP(3)
    A(9,11) = -VP(4)
    B(9,9)  = +VP(3)
    C(9,9)  = +VP(4) 
    !
    A(10,10) = +VP(2)
    B(10, 9) = -VP(2)
    B(10,11) = -VP(4)
    C(10,10) = +VP(4)
    !
    A(11,11) = +VP(2)
    B(11,11) = +VP(3)
    C(11, 9) = -VP(2)
    C(11,10) = -VP(3)
    !
    A(12,13) = -VP(3)
    A(12,14) = -VP(4)
    B(12,12) = +VP(3)
    C(12,12) = +VP(4)
    !
    A(13,13) = +VP(2)
    B(13,12) = -VP(2)
    B(13,14) = -VP(4)
    C(13,13) = +VP(4)
    !
    A(14,14) = +VP(2)
    B(14,14) = +VP(3)
    C(14,12) = -VP(2)
    C(14,13) = -VP(3)
  
    if( nDim .eq. 2) then
        An = A*nv(1) + B*nv(2)         
    else
        An = A*nv(1) + B*nv(2) + C*nv(3)     
    end if
    
    

  
END SUBROUTINE PDEMatrixB


