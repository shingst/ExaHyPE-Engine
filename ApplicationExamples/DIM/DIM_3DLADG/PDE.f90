! GRMHD PDE.f90
! Trento (EQNTYPE4)

RECURSIVE SUBROUTINE PDEFlux(f,g,hz,Q)
  USE Parameters, ONLY : nVar, nDim
  USE iso_c_binding
  IMPLICIT NONE
  REAL :: f(nVar), g(nVar), h(nVar), hz(nVar), Q(nVar), V(nVar)
  REAL, PARAMETER :: epsilon = 1e-14 
  INTENT(IN)  :: Q
  INTENT(OUT) :: f, g, hz

  
  !FTensDim = 0
  !RETURN
  
  f = 0
  g = 0
  h = 0
  END SUBROUTINE PDEFlux


RECURSIVE SUBROUTINE PDENCP(BgradQ,Q,gradQ) 
   USE Parameters, ONLY :  nVar, nDim
   IMPLICIT NONE
   REAL, PARAMETER :: epsilon = 1e-14
   ! 11. Oct 21:40: This was a matrix BGradQ(nVar, nDim) but is a vector in spaceTimePredictorNonlinear
   REAL, INTENT(OUT) :: BgradQ(nVar)
   REAL, INTENT(IN)  :: gradQ(nVar, 3)
   REAL, INTENT(IN)  :: Q(nVar)
  ! Linear elasticity variables
   REAL :: lam,mu,irho, ialpha, u(3)
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
    ! Linear elasticity part
   AQx = 0.   
   BQy = 0. 
   CQz = 0. 

    lam  = Q(10)
    mu   = Q(11) 
    irho = 1./Q(12)
    IF(Q(13)<=1e-3) THEN
        BgradQ = 0.0  
        RETURN 
        !
    ELSE
        ialpha = 1./Q(13)
        u      = Q(7:9)*ialpha 
    ENDIF 
    !
    AQx(1) = - (lam+2*mu)*Qx(7) + (lam+2*mu)*u(1)*Qx(13) 
    AQx(2) = - lam*Qx(7)        + lam*u(1)*Qx(13)
    AQx(3) = - lam*Qx(7)        + lam*u(1)*Qx(13) 
    AQx(4) = - mu *Qx(8)        + mu *u(2)*Qx(13) 
    AQx(5) =   0.0 
    AQx(6) = - mu *Qx(9)        + mu *u(3)*Qx(13) 
    AQx(7) = - irho * Qx(1) - 2*Q(1)*irho*Qx(13)  
    AQx(8) = - irho * Qx(4) - 2*Q(4)*irho*Qx(13)   
    AQx(9) = - irho * Qx(6) - 2*Q(6)*irho*Qx(13)       
    !
    BQy(1) = - lam*Qy(8)        + lam*u(2)*Qy(13) 
    BQy(2) = - (lam+2*mu)*Qy(8) + (lam+2*mu)*u(2)*Qy(13) 
    BQy(3) = - lam*Qy(8)        + lam*u(2)*Qy(13)  
    BQy(4) = - mu *Qy(7)        + mu *u(1)*Qy(13) 
    BQy(5) = - mu *Qy(9)        + mu *u(3)*Qy(13) 
    BQy(6) =   0.0  
    BQy(7) = - irho * Qy(4) - 2*Q(4)*irho*Qy(13)  
    BQy(8) = - irho * Qy(2) - 2*Q(2)*irho*Qy(13)   
    BQy(9) = - irho * Qy(5) - 2*Q(5)*irho*Qy(13)      
    !
    CQz(1) = - lam*Qz(9)        + lam*u(3)*Qz(13) 
    CQz(2) = - lam*Qz(9)        + lam*u(3)*Qz(13) 
    CQz(3) = - (lam+2*mu)*Qz(9) + (lam+2*mu)*u(3)*Qz(13) 
    CQz(4) =  0.0  
    CQz(5) = - mu *Qz(8)        + mu *u(2)*Qz(13)  
    CQz(6) = - mu *Qz(7)        + mu *u(1)*Qz(13) 
    CQz(7) = - irho * Qz(6) - 2*Q(6)*irho*Qz(13)  
    CQz(8) = - irho * Qz(5) - 2*Q(5)*irho*Qz(13)  
    CQz(9) = - irho * Qz(3) - 2*Q(3)*irho*Qz(13)     

    
    BgradQ = AQx + BQy + CQz
END SUBROUTINE PDENCP


RECURSIVE SUBROUTINE PDEEigenvalues(L,Q,n)
  USE Parameters, ONLY :  nVar, nDim
  USE iso_c_binding
  IMPLICIT NONE
  REAL :: L(nVar), n(3), Q(nVar), V(nVar)
  INTENT(IN)  :: Q,n
  INTENT(OUT) :: L 
  REAL  :: lam,mu,irho
  ! Local Variables 
  REAL :: Pi
  Pi = ACOS(-1.0)

  L(:) = 0
  
    lam  = Q(10)   
    mu   = Q(11) 
    irho = 1./Q(12)
    !
    L(1) = -SQRT((lam+2.0*mu)*irho)     
    L(2) = +SQRT((lam+2.0*mu)*irho) 
    L(3) = -SQRT(mu*irho) 
    L(4) = +SQRT(mu*irho) 
    L(5) = 0. 

END SUBROUTINE PDEEigenvalues

RECURSIVE SUBROUTINE PDESource(S,Q) 
  USE Parameters, ONLY:  nVar, nDim
  USE iso_c_binding
  IMPLICIT NONE
  ! --------------------------------------------
  ! Argument list declaration
  REAL :: S(nvar), Q(nvar)
  INTENT(IN)  :: Q 
  INTENT(OUT) :: S
  ! --------------------------------------------
  ! Local variable declaration 

  
  S = 0
      
END SUBROUTINE PDESource

RECURSIVE SUBROUTINE PDEVarName(Name) 
  USE Parameters, ONLY: nVar  
  IMPLICIT NONE     
  CHARACTER(LEN=10):: Name(nVar)

  ! EQNTYPE99
    Name(1)  = 'sxx' 
    Name(2)  = 'syy' 
    Name(3)  = 'szz'
    Name(4)  = 'sxy' 
    Name(5)  = 'syz' 
    Name(6)  = 'sxz' 
    Name(7)  = 'u' 
    Name(8)  = 'v' 
    Name(9)  = 'w' 
    Name(10) = 'lambda'
    Name(11) = 'mu'
    Name(12) = 'rho' 
    Name(13) = 'alpha'
    Name(14) = 'xi' 
    END SUBROUTINE PDEVarName

RECURSIVE SUBROUTINE PDEJacobian(An,Q,gradQ,nv)
  USE Parameters, ONLY : nVar, nDim
  USE iso_c_binding
  IMPLICIT NONE
  ! Argument list 
  REAL :: An(nVar,nVar)
  REAL :: Q(nVar), gradQ(nVar,nDim), nv(3) 
  INTENT(IN)  :: Q, gradQ, nv
  INTENT(OUT) :: An 
  
  CALL PDEMatrixB(An,Q,nv) 
  
END SUBROUTINE PDEJacobian

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
   REAL :: lam,mu,irho, ialpha, uv(3)
   REAL :: A(nVar,nVar), B(nVar,nVar), C(nVar,nVar)
   !An = 0
  !RETURN
	
    lam  = Q(10) 
    mu   = Q(11) 
    irho = 1./Q(12)
    !print *, lam,mu,irho,Q(13)
    IF(Q(13)<=1e-3) THEN        
        ialpha = 0.0
        uv = 0.0 
        An = 0.0
        RETURN 
    ELSE
        ialpha = 1./Q(13)
        uv     = Q(7:9)*ialpha 
    ENDIF 
    
	!print *, maxval(nv),lam,mu,irho
	
    A = 0.0
    B = 0.0
    C = 0.0 

    A(1,7)  = - (lam+2*mu) 
    A(1,13) = + (lam+2*mu)*uv(1)  
    A(2,7)  = - lam 
    A(2,13) = + lam*uv(1) 
    A(3,7)  = - lam 
    A(3,13) = + lam*uv(1)  
    A(4,8)  = - mu
    A(4,13) = + mu*uv(2)
    A(6,9)  = - mu
    A(6,13) = + mu*uv(3)
    A(7,1)  = - irho
    A(7,13) = - 2*Q(1)*irho  
    A(8,4)  = - irho
    A(8,13) = - 2*Q(4)*irho    
    A(9,6)  = - irho 
    A(9,13) = - 2*Q(6)*irho        
    !
    B(1,8)  = - lam 
    B(1,13) = + lam*uv(2) 
    B(2,8)  = - (lam+2*mu) 
    B(2,13) = + (lam+2*mu)*uv(2) 
    B(3,8)  = - lam 
    B(3,13) = + lam*uv(2) 
    B(4,7)  = - mu  
    B(4,13) = + mu*uv(1)  
    B(5,9)  = - mu 
    B(5,13) = + mu*uv(3) 
    B(7,4)  = - irho 
    B(7,13) = - 2*Q(4)*irho 
    B(8,2)  = - irho 
    B(8,13) = - 2*Q(2)*irho 
    B(9,5)  = - irho 
    B(9,13) = - 2*Q(5)*irho       
    !
    C(1,9)  = - lam 
    C(1,13) = + lam*uv(3) 
    C(2,9)  = - lam 
    C(2,13) = + lam*uv(3) 
    C(3,9)  = - (lam+2*mu) 
    C(3,13) = + (lam+2*mu)*uv(3)
    C(5,8)  = - mu  
    C(5,13) = + mu*uv(2)  
    C(6,7)  = - mu 
    C(6,13) = + mu*uv(1) 
    C(7,6)  = - irho 
    C(7,13) = - 2*Q(6)*irho 
    C(8,5)  = - irho 
    C(8,13) = - 2*Q(5)*irho 
    C(9,3)  = - irho 
    C(9,13) = - 2*Q(3)*irho     
    if( nDim .eq. 2) then
        An = A*nv(1) + B*nv(2)         
    else
        An = A*nv(1) + B*nv(2) + C*nv(3)     
    end if
    
    

  
END SUBROUTINE PDEMatrixB
    
RECURSIVE SUBROUTINE PDEIntermediateFields(RL,LL,iRL,Q,nv) 
  USE Parameters, ONLY : nVar, nDim, nLin
  USE iso_c_binding 
  IMPLICIT NONE
  ! Argument list 
  REAL :: RL(nVar,nLin), LL(nLin,nLin), iRL(nLin,nVar)
  REAL :: Q(nVar), nv(3)
  REAL :: x(nDim), time
  INTENT(IN)  :: Q,nv
  INTENT(OUT) :: RL,LL,iRL 
  ! Local variables
  INTEGER :: i,j,k, zero,iErr, minl(1), maxl(1) 
  REAL :: R(nVar,nVar), L(nVar,nVar), Lv(nVar), iR(nVar,nVar)
  !   
  CALL PDEEigenvectors(R,L,iR,Q,nv)
  RL(:,1:4)  = R(:,7:10) 
  iRL(1:4,:) = iR(7:10,:) 
  RL(:,5)  = R(:,14) 
  iRL(5,:) = iR(14,:) 
  LL = 0. 
  LL(1,1) = L(7,7) 
  LL(2,2) = L(8,8) 
  LL(3,3) = L(9,9) 
  LL(4,4) = L(10,10) 
  LL(5,5) = L(14,14)   
  !
END SUBROUTINE PDEIntermediateFields
    
    
RECURSIVE SUBROUTINE PDEEigenvectors(R,L,iR,Q,nv) 
  USE Parameters, ONLY : nVar, nDim
  USE iso_c_binding  
  IMPLICIT NONE
  ! Argument list 
  REAL :: R(nVar,nVar), L(nVar,nVar), iR(nVar,nVar)
  REAL :: Q(nVar), nv(3)
  REAL :: x(nDim), time
  INTENT(IN)  :: Q,nv
  INTENT(OUT) :: R,L,iR 
  ! Local variables
  INTEGER :: i,j,k, zero,iErr, itemp(nVar)    
  REAL    :: rho,u,v,w,p,c,c2,H,v2,M,r2c,Pi,gmo,den,lf,lf2,kk,hrho
  REAL    :: Cplus, Cminus, Aplus, Aminus, Lplus, Lminus, Delta, ftr, vs
  REAL    :: sv(3),tv(3),Lambda(nVar)  
  REAL    :: VPR(3),BVR(3),kappa(2),beta(2)
  REAL    :: TM(nVar,nVar),iTM(nVar,nVar),A(nVar,nVar),RM(nVar,nVar),iRM(nVar,nVar) 
  REAL    :: TestMatrix(nVar,nVar), gradQ(nVar,nDim)  
  REAL    :: R7(7,7), iR7(7,7)
  REAL    :: VP(nVar), QPR(nVar),dudw(nVar,nVar),dwdu(nVar,nVar)
  REAL    :: dfdQ(nVar,nVar), ImLambda(nVar), rtemp(nVar)    
    REAL    :: rho0, k0, uu, vv, ww, eps, dist, alpha, cp, cs, nx, ny, uv(3), lam, mu, irho, ialpha, s11, s12, s13 
  REAL    ::  Qrot(nVar)      
  

  
  REAL, PARAMETER :: epsilon = 1e-11   
  !
  Pi = ACOS(-1.0)
  R = 0. 
  L = 0. 
  iR = 0. 

        !
        IF(ABS(ABS(nv(1))-1.0).LE.1e-14) THEN
            sv = (/ 0., 1., 0. /) 
            tv = (/ 0., 0., 1. /) 
        ENDIF
        IF(ABS(ABS(nv(2))-1.0).LE.1e-14) THEN
            sv = (/ 1., 0., 0. /) 
            tv = (/ 0., 0., 1. /) 
        ENDIF
        IF(ABS(ABS(nv(3))-1.0).LE.1e-14) THEN
            sv = (/ 1., 0., 0. /) 
            tv = (/ 0., 1., 0. /) 
        ENDIF
        !
        ! The exact Riemann solver of Godunov 
        !
        TM = 0.0    ! rotation matrix 
        TM(1,1:9) = (/ nv(1)**2   , sv(1)**2   , tv(1)**2   , 2*nv(1)*sv(1)          , 2*sv(1)*tv(1)          , 2*nv(1)*tv(1)          , 0., 0., 0. /) 
        TM(2,1:9) = (/ nv(2)**2   , sv(2)**2   , tv(2)**2   , 2*nv(2)*sv(2)          , 2*sv(2)*tv(2)          , 2*nv(2)*tv(2)          , 0., 0., 0. /) 
        TM(3,1:9) = (/ nv(3)**2   , sv(3)**2   , tv(3)**2   , 2*nv(3)*sv(3)          , 2*sv(3)*tv(3)          , 2*nv(3)*tv(3)          , 0., 0., 0. /) 
        TM(4,1:9) = (/ nv(2)*nv(1), sv(2)*sv(1), tv(2)*tv(1), nv(2)*sv(1)+nv(1)*sv(2), sv(2)*tv(1)+sv(1)*tv(2), nv(2)*tv(1)+nv(1)*tv(2), 0., 0., 0. /) 
        TM(5,1:9) = (/ nv(3)*nv(2), sv(3)*sv(2), tv(3)*tv(2), nv(3)*sv(2)+nv(2)*sv(3), sv(3)*tv(2)+sv(2)*tv(3), nv(3)*tv(2)+nv(2)*tv(3), 0., 0., 0. /) 
        TM(6,1:9) = (/ nv(3)*nv(1), sv(3)*sv(1), tv(3)*tv(1), nv(3)*sv(1)+nv(1)*sv(3), sv(3)*tv(1)+sv(1)*tv(3), nv(3)*tv(1)+nv(1)*tv(3), 0., 0., 0. /) 
        TM(7,1:9) = (/ 0.,0., 0., 0., 0., 0., nv(1), sv(1), tv(1) /)     
        TM(8,1:9) = (/ 0.,0., 0., 0., 0., 0., nv(2), sv(2), tv(2) /)     
        TM(9,1:9) = (/ 0.,0., 0., 0., 0., 0., nv(3), sv(3), tv(3) /)  
        TM(10,10) = 1.0 
        TM(11,11) = 1.0 
        TM(12,12) = 1.0 
        TM(13,13) = 1.0 
        TM(14,14) = 1.0 
        !
        iTM = 0.0   ! inverse of the rotation matrix 
        iTM(1,1:9) = (/ nv(1)**2   , nv(2)**2   , nv(3)**2   , 2*nv(1)*nv(2)          , 2*nv(3)*nv(2)          , 2*nv(1)*nv(3)          , 0., 0., 0. /) 
        iTM(2,1:9) = (/ sv(1)**2   , sv(2)**2   , sv(3)**2   , 2*sv(1)*sv(2)          , 2*sv(3)*sv(2)          , 2*sv(1)*sv(3)          , 0., 0., 0. /) 
        iTM(3,1:9) = (/ tv(1)**2   , tv(2)**2   , tv(3)**2   , 2*tv(1)*sv(2)          , 2*tv(3)*tv(2)          , 2*tv(1)*tv(3)          , 0., 0., 0. /) 
        iTM(4,1:9) = (/ nv(1)*sv(1), nv(2)*sv(2), nv(3)*sv(3), nv(1)*sv(2)+sv(1)*nv(2), nv(3)*sv(2)+sv(3)*nv(2), nv(1)*sv(3)+sv(1)*nv(3), 0., 0., 0. /) 
        iTM(5,1:9) = (/ sv(1)*tv(1), sv(2)*tv(2), sv(3)*tv(3), sv(1)*tv(2)+tv(1)*sv(2), sv(3)*tv(2)+tv(3)*sv(2), sv(1)*tv(3)+tv(1)*sv(3), 0., 0., 0. /) 
        iTM(6,1:9) = (/ nv(1)*tv(1), nv(2)*tv(2), nv(3)*tv(3), nv(1)*tv(2)+tv(1)*nv(2), nv(3)*tv(2)+tv(3)*nv(2), nv(1)*tv(3)+tv(1)*nv(3), 0., 0., 0. /) 
        iTM(7,1:9) = (/ 0.,0., 0., 0., 0., 0., nv(1), nv(2), nv(3) /)     
        iTM(8,1:9) = (/ 0.,0., 0., 0., 0., 0., sv(1), sv(2), sv(3) /)     
        iTM(9,1:9) = (/ 0.,0., 0., 0., 0., 0., tv(1), tv(2), tv(3) /)   
        iTM(10,10) = 1.0 
        iTM(11,11) = 1.0 
        iTM(12,12) = 1.0         
        iTM(13,13) = 1.0         
        iTM(14,14) = 1.0 
        !
        lam  = Q(10) 
        mu   = Q(11) 
        rho  = Q(12) 
        irho = 1./rho 
        Qrot = MATMUL( iTM, Q) 
        alpha = Q(13) 
        IF(Q(13)<=1e-3) THEN
            ialpha = 0.0
            uv = 0.0 
        ELSE
            ialpha = 1./Q(13)
            uv     = Qrot(7:9)*ialpha 
        ENDIF 
        ! 
        s11 = Qrot(1) 
        s12 = Qrot(4) 
        s13 = Qrot(6) 
        
        !
        cp = SQRT((lam+2*mu)/rho) 
        cs = SQRT(mu/rho) 
        !
        RM = 0.0 
        RM(1,1)   = rho*cp**2
        RM(1,10)  = -2*s11
        RM(1,13)  = rho*cp**2
        RM(2,1)   = -rho*(-cp**2+2*cs**2)
        RM(2,4)   = 1
        RM(2,13)  = -rho*(-cp**2+2*cs**2)
        RM(3,1)   = -rho*(-cp**2+2*cs**2)
        RM(3,5)   = 1
        RM(3,13)  = -rho*(-cp**2+2*cs**2)
        RM(4,2)   = rho*cs**2
        RM(4,10)  = -2*s12
        RM(4,12)  = rho*cs**2
        RM(5,6)   = 1
        RM(6,3)   = rho*cs**2
        RM(6,10)  = -2*s13
        RM(6,11)  = rho*cs**2
        RM(7,1)   = cp
        RM(7,10)  = uv(1) 
        RM(7,13)  = -cp
        RM(8,2)   = cs
        RM(8,10)  = uv(2) 
        RM(8,12)  = -cs
        RM(9,3)   = cs
        RM(9,10)  = uv(3) 
        RM(9,11)  = -cs
        RM(10,9)  = 1
        RM(11,8)  = 1
        RM(12,7)  = 1
        RM(13,10) = 1
        RM(14,14) = 1 
        !
        iRM = 0.0 
        iRM(1,1)   = irho/cp**2/2
        iRM(1,7)   = 1/cp/2
        iRM(1,13)  = -1/cp**2/rho*(cp*rho*uv(1)-2*s11)/2
        iRM(2,4)   = irho/cs**2/2
        iRM(2,8)   = 1/cs/2
        iRM(2,13)  = -1/cs**2/rho*(uv(2)*cs*rho-2*s12)/2
        iRM(3,6)   = irho/cs**2/2
        iRM(3,9)   = 1/cs/2
        iRM(3,13)  = -1/cs**2/rho*(uv(3)*cs*rho-2*s13)/2
        iRM(4,1)   = 1/cp**2*(-cp**2+2*cs**2)
        iRM(4,2)   = 1
        iRM(4,13)  = 1/cp**2*2*s11*(-cp**2+2*cs**2) 
        iRM(5,1)   = 1/cp**2*(-cp**2+2*cs**2)
        iRM(5,3)   = 1
        iRM(5,13)  = 1/cp**2*2*s11*(-cp**2+2*cs**2) 
        iRM(6,5)   = 1
        iRM(7,12)  = 1
        iRM(8,11)  = 1
        iRM(9,10)  = 1
        iRM(10,13) = 1. 
        iRM(11,6)  = irho/cs**2/2
        iRM(11,9)  = -1/cs/2
        iRM(11,13) = 1/cs**2/rho*(uv(3)*cs*rho+2*s13)/2
        iRM(12,4)  = irho/cs**2/2
        iRM(12,8)  = -1/cs/2
        iRM(12,13) = 1/cs**2/rho*(uv(2)*cs*rho+2*s12)/2
        iRM(13,1)  = irho/cp**2/2
        iRM(13,7)  = -1/cp/2
        iRM(13,13) = 1/cp**2/rho*(cp*rho*uv(1)+2*s11)/2
        iRM(14,14) = 1 
        !
        R = MATMUL(   TM,  RM ) 
        iR = MATMUL( iRM, iTM ) 
        !
        L = 0. 
        L(1,1)   = -cp
        L(2,2)   = -cs
        L(3,3)   = -cs
        L(11,11) = +cs 
        L(12,12) = +cs
        L(13,13) = +cp
        L(14,14) = 0.0 
        ! 

    END SUBROUTINE PDEEigenvectors 
    
RECURSIVE SUBROUTINE HLLEMRiemannSolver(basisSize,NormalNonZero,lFbndL,lFbndR,lQbndL,lQbndR,QavL,QavR) 
  USE Parameters, ONLY : nVar, nDim, nLin
  USE iso_c_binding 
  ! Local variables
  INTEGER, INTENT(IN)   :: NormalNonZero, basisSize
  REAL, INTENT(IN)     :: lQbndL(nVar,basisSize,basisSize)
  REAL, INTENT(IN)     :: lQbndR(nVar,basisSize,basisSize)
  REAL, INTENT(INOUT)  :: lFbndL(nVar,basisSize,basisSize)
  REAL, INTENT(INOUT)  :: lFbndR(nVar,basisSize,basisSize)
    ! Local variables 
INTEGER           :: i,j,k,l, ml(1)  
REAL              :: aux(nDim), Id(nVar,nVar), smax, Qav(nVar),QavL(nVar), QavR(nVar) 
REAL              ::  xGP, yGP, xdeb, ydeb  
REAL              :: Bn(nVar,nVar), DM(nVar,nVar), ncp(nVar), nv(3)
REAL              :: gradQ(nVar,nDim), src(nVar),flattener(nLin)
  REAL    :: absA(nVar,nVar), amax  
  REAL    :: QM(nVar),LL(nVar),LR(nVar),LM(nVar)
  REAL    :: deltaL(nLin,nLin),Lam(nLin,nLin),Lap(nLin,nLin) 
  REAL    :: RL(nVar,nLin),iRL(nLin,nVar),LLin(nLin,nLin) 
  
  nv(:)=0.
  nv(NormalNonZero+1)=1.
  !print *, "Normal non zero in fortran=" NormalNonZero
  !print *, "basisSize=", basisSize
  !print *, "NormalNonZero=", NormalNonZero
  !print *, "QavR=",QavR(1)
  !return
  !nv(NormalNonZero)=1.;
  !FL=0.
  !FR=0.
	flattener=1.
	lFbndL=0.
	lFbndR=0.
    CALL PDEEigenvalues(LL,QavL,nv) 
    CALL PDEEigenvalues(LR,QavR,nv) 
    smax = MAX( MAXVAL(ABS(LL)), MAXVAL(ABS(LR)) )
    ! Identity matrix 
    Id = 0.0 
    DO i=1,nVar
        Id(i,i)=1.0
    ENDDO
    gradQ = 0.0      
    ! HLLEM
    Qav = 0.5*(QavL+QavR) 
    CALL PDEIntermediateFields(RL,LLin,iRL,Qav,nv) 
    Lam = 0.5*(LLin-ABS(LLin))
    Lap = 0.5*(LLin+ABS(LLin)) 
    deltaL = 0.0
    DO i = 1, nLin
        deltaL(i,i) = ( 1. - Lam(i,i)/(-smax-1e-14) - Lap(i,i)/(smax+1e-14) )*flattener(i)  
    ENDDO    
    absA = 0. 
    DO i = 1, nVar
        absA(i,i) = -0.5*smax  ! regular Rusanov diffusion, only on the diagonal 
    ENDDO  
    absA = absA + 0.5*smax*MATMUL( RL, MATMUL(deltaL, iRL) )  ! HLLEM anti-diffusion

        DO k = 1, basisSize
          DO j = 1, basisSize
                Qav = 0.5*(lQbndR(:,j,k)+lQbndL(:,j,k)) 
                gradQ(:,NormalNonZero+1) = lQbndR(:,j,k) - lQbndL(:,j,k) 
                CALL PDENCP(ncp,Qav,gradQ)
				lFbndL(:,j,k) = 0.5*( lFbndR(:,j,k) + lFbndL(:,j,k) ) + MATMUL(absA, lQbndR(:,j,k) - lQbndL(:,j,k) )    ! purely conservative flux 
				!lFbndL(:,j,k) = 0.5*( lFbndR(:,j,k) + lFbndL(:,j,k) ) - 0.5*smax*( lQbndR(:,j,k) - lQbndL(:,j,k) )      ! purely conservative flux 
				lFbndR(:,j,k) = lFbndL(:,j,k) - 0.5*ncp(:)                                                              ! subtract the non-conservative product 
                lFbndL(:,j,k) = lFbndL(:,j,k) + 0.5*ncp(:)
            ENDDO
        ENDDO				
END SUBROUTINE HLLEMRiemannSolver


RECURSIVE SUBROUTINE TESTRiemannSolver(basisSize,NormalNonZero,lFbndL,lFbndR,lQbndL,lQbndR,QavL,QavR) 
  USE Parameters, ONLY : nVar, nDim, nLin
  USE iso_c_binding 
  ! Local variables
  INTEGER, INTENT(IN)   :: NormalNonZero, basisSize
  REAL, INTENT(IN)     :: lQbndL(nVar,basisSize,basisSize)
  REAL, INTENT(IN)     :: lQbndR(nVar,basisSize,basisSize)
  REAL, INTENT(INOUT)  :: lFbndL(nVar,basisSize,basisSize)
  REAL, INTENT(INOUT)  :: lFbndR(nVar,basisSize,basisSize)
    ! Local variables 
INTEGER           :: i,j,k,l, ml(1)  
REAL              :: aux(nDim), Id(nVar,nVar), smax, Qav(nVar),QavL(nVar), QavR(nVar) 
REAL              ::  xGP, yGP, xdeb, ydeb  
REAL              :: Bn(nVar,nVar), DM(nVar,nVar), ncp(nVar), nv(3)
REAL              :: gradQ(nVar,3), src(nVar),flattener(nLin)
  REAL    :: absA(nVar,nVar), amax  
  REAL    :: QM(nVar),LL(nVar),LR(nVar),LM(nVar)
  REAL    :: deltaL(nLin,nLin),Lam(nLin,nLin),Lap(nLin,nLin) 
  REAL    :: RL(nVar,nLin),iRL(nLin,nVar),LLin(nLin,nLin) 
  nv(:)=0.
  !if(NormalNonZero .gt. 0) then
	nv(NormalNonZero+1)=1.
  !end if
	Qav = 0.5*(QavL+QavR) 
	!CALL PDEMatrixB(Bn,0.5*(QavL+QavR),nv)   ! evaluate the system matrix just once in the averaged state 
    !CALL DissipationMatrix(DM,QavL,QavR,nv)        ! according to the RIemann solver, compute the dissipation matrix 
    !DO k = 1, basisSize
    !  DO j = 1, basisSize
    !        !CALL PDEMatrixB(Bn,0.5*(lQbndL(:,j,k)+lQbndR(:,j,k)),nv,0.5*(lparbndL(:,j,k)+lparbndR(:,j,k))) ! evaluate the system matrix in each Gaussian point again (slow) 
    !        lFbndL(:,j,k) = 0.5*MATMUL( Bn + DM, lQbndR(:,j,k) - lQbndL(:,j,k) ) ! WAS -
    !        lFbndR(:,j,k) = 0.5*MATMUL( Bn - DM, lQbndR(:,j,k) - lQbndL(:,j,k) )! WAS + 	
    !    ENDDO
    !ENDDO
	!print *, NormalNonZero
	!print *, maxval(abs(lFbndL(:,1,1)-lFbndR(:,1,1)))
	!if(maxval(abs(Bn)) .gt. 0) then
	!	print *, 'OK!', maxval(abs(Bn))
	!end if
	
	
	
	CALL PDEEigenvalues(LL,QavL,nv) 
	CALL PDEEigenvalues(LR,QavR,nv) 
	CALL PDEEigenvalues(LM,Qav,nv) 
	sL = MIN(0.0, MINVAL(LL), MINVAL(LM) ) 
	sR = MAX(0.0, MAXVAL(LR), MAXVAL(LM) )
    ! Identity matrix 
    Id = 0.0 
    DO i=1,nVar
        Id(i,i)=1.0
    ENDDO
	lFbndL=0;
	lFbndR=0;
	gradQ = 0.0  
        DO k = 1, basisSize
          DO j = 1, basisSize
                Qav = 0.5*(lQbndR(:,j,k)+lQbndL(:,j,k)) 
                gradQ(:,NormalNonZero+1) = lQbndR(:,j,k) - lQbndL(:,j,k) 
                CALL PDENCP(ncp,Qav,gradQ)
				!lFbndL(:,j,k) = sL*sR/(sR-sL)*( lQbndR(:,j,k) - lQbndL(:,j,k) )
                lFbndL(:,j,k) = ( sR*lFbndL(:,j,k) - sL*lFbndR(:,j,k) )/(sR-sL) + sR*sL/(sR-sL)*( lQbndR(:,j,k) - lQbndL(:,j,k) )   ! purely conservative flux 
                lFbndR(:,j,k) = lFbndL(:,j,k) - sR/(sR-sL)*ncp(:)                                                                   ! subtract the non-conservative product 
                lFbndL(:,j,k) = lFbndL(:,j,k) - sL/(sR-sL)*ncp(:)                                                                   ! add the non-conservative product 
            ENDDO
        ENDDO	
END SUBROUTINE TESTRiemannSolver

SUBROUTINE DissipationMatrix(DM,QL,QR,nv)
    USE Parameters, ONLY : nVar
    USE iso_c_binding 	
    IMPLICIT NONE 
    ! Argument list 
	integer, parameter :: d=3
    REAL, INTENT(IN)     :: QL(nVar), QR(nVar)                        ! left state, right state 
    REAL, INTENT(IN)     :: nv(d)                                     ! normal vector  
    REAL, INTENT(INOUT)  :: DM(nVar, nVar)                            ! Dissipation Matrix >= 0 
    ! Local variables 
    INTEGER           :: i,j,k,l
	
    REAL              :: aux(d), QavL(nVar), QavR(nVar), Id(nVar,nVar), smax, tv(d), sv(d)
    REAL              :: LL(nVar), LR(nVar), iRM(nVar,nVar), RM(nVar,nVar), TM(nVar,nVar), iTM(nVar,nVar), test(nVar,nVar), LM(nVar,nVar), lambda, mu, cp, cs, rho   
    REAL              :: lam, irho, alpha, ialpha, Qrot(nVar), s11, s12, s13, Q(nVar), uv(3), R(nVar,nVar), iR(nVar,nVar)     

    IF(ABS(nv(1)-1.0) < 1e-14) THEN
        sv = (/ 0., 1., 0. /)
        tv = (/ 0., 0., 1. /) 
    ENDIF    
    IF(ABS(nv(2)-1.0) < 1e-14) THEN
        sv = (/ 0., 0., 1. /)
        tv = (/ 1., 0., 0.  /) 
    ENDIF    
    IF(ABS(nv(3)-1.0) < 1e-14) THEN
        sv = (/ 1., 0., 0. /)
        tv = (/ 0., 1., 0.  /) 
    ENDIF  
    !
    ! The exact Riemann solver of Godunov 
    !
    TM = 0.0    ! rotation matrix 
    TM(1,1:9) = (/ nv(1)**2   , sv(1)**2   , tv(1)**2   , 2*nv(1)*sv(1)          , 2*sv(1)*tv(1)          , 2*nv(1)*tv(1)          , 0., 0., 0. /) 
    TM(2,1:9) = (/ nv(2)**2   , sv(2)**2   , tv(2)**2   , 2*nv(2)*sv(2)          , 2*sv(2)*tv(2)          , 2*nv(2)*tv(2)          , 0., 0., 0. /) 
    TM(3,1:9) = (/ nv(3)**2   , sv(3)**2   , tv(3)**2   , 2*nv(3)*sv(3)          , 2*sv(3)*tv(3)          , 2*nv(3)*tv(3)          , 0., 0., 0. /) 
    TM(4,1:9) = (/ nv(2)*nv(1), sv(2)*sv(1), tv(2)*tv(1), nv(2)*sv(1)+nv(1)*sv(2), sv(2)*tv(1)+sv(1)*tv(2), nv(2)*tv(1)+nv(1)*tv(2), 0., 0., 0. /) 
    TM(5,1:9) = (/ nv(3)*nv(2), sv(3)*sv(2), tv(3)*tv(2), nv(3)*sv(2)+nv(2)*sv(3), sv(3)*tv(2)+sv(2)*tv(3), nv(3)*tv(2)+nv(2)*tv(3), 0., 0., 0. /) 
    TM(6,1:9) = (/ nv(3)*nv(1), sv(3)*sv(1), tv(3)*tv(1), nv(3)*sv(1)+nv(1)*sv(3), sv(3)*tv(1)+sv(1)*tv(3), nv(3)*tv(1)+nv(1)*tv(3), 0., 0., 0. /) 
    TM(7,1:9) = (/ 0.,0., 0., 0., 0., 0., nv(1), sv(1), tv(1) /)     
    TM(8,1:9) = (/ 0.,0., 0., 0., 0., 0., nv(2), sv(2), tv(2) /)     
    TM(9,1:9) = (/ 0.,0., 0., 0., 0., 0., nv(3), sv(3), tv(3) /)  
    TM(10,10) = 1.0 
    TM(11,11) = 1.0 
    TM(12,12) = 1.0 
    TM(13,13) = 1.0 
    TM(14,14) = 1.0 
    !
    iTM = 0.0   ! inverse of the rotation matrix 
    iTM(1,1:9) = (/ nv(1)**2   , nv(2)**2   , nv(3)**2   , 2*nv(1)*nv(2)          , 2*nv(3)*nv(2)          , 2*nv(1)*nv(3)          , 0., 0., 0. /) 
    iTM(2,1:9) = (/ sv(1)**2   , sv(2)**2   , sv(3)**2   , 2*sv(1)*sv(2)          , 2*sv(3)*sv(2)          , 2*sv(1)*sv(3)          , 0., 0., 0. /) 
    iTM(3,1:9) = (/ tv(1)**2   , tv(2)**2   , tv(3)**2   , 2*tv(1)*sv(2)          , 2*tv(3)*tv(2)          , 2*tv(1)*tv(3)          , 0., 0., 0. /) 
    iTM(4,1:9) = (/ nv(1)*sv(1), nv(2)*sv(2), nv(3)*sv(3), nv(1)*sv(2)+sv(1)*nv(2), nv(3)*sv(2)+sv(3)*nv(2), nv(1)*sv(3)+sv(1)*nv(3), 0., 0., 0. /) 
    iTM(5,1:9) = (/ sv(1)*tv(1), sv(2)*tv(2), sv(3)*tv(3), sv(1)*tv(2)+tv(1)*sv(2), sv(3)*tv(2)+tv(3)*sv(2), sv(1)*tv(3)+tv(1)*sv(3), 0., 0., 0. /) 
    iTM(6,1:9) = (/ nv(1)*tv(1), nv(2)*tv(2), nv(3)*tv(3), nv(1)*tv(2)+tv(1)*nv(2), nv(3)*tv(2)+tv(3)*nv(2), nv(1)*tv(3)+tv(1)*nv(3), 0., 0., 0. /) 
    iTM(7,1:9) = (/ 0.,0., 0., 0., 0., 0., nv(1), nv(2), nv(3) /)     
    iTM(8,1:9) = (/ 0.,0., 0., 0., 0., 0., sv(1), sv(2), sv(3) /)     
    iTM(9,1:9) = (/ 0.,0., 0., 0., 0., 0., tv(1), tv(2), tv(3) /)   
    iTM(10,10) = 1.0 
    iTM(11,11) = 1.0 
    iTM(12,12) = 1.0         
    iTM(13,13) = 1.0         
    iTM(14,14) = 1.0     
    !
    !test = MATMUL( TM, TRANSPOSE(TM) )  
    !
    Q = 0.5*(QL+QR) 
    lam  = Q(10) 
    mu   = Q(11) 
    rho  = Q(12) 
    irho = 1./rho 
    Qrot = MATMUL( iTM, Q) 
    alpha = Q(13) 
    IF(Q(13)<=1e-3) THEN
        ialpha = 0.0
        uv = 0.0 
    ELSE
        ialpha = 1./Q(13)
        uv     = Qrot(7:9)*ialpha 
    ENDIF 
    ! 
    s11 = Qrot(1) 
    s12 = Qrot(4) 
    s13 = Qrot(6) 
        
    !
    cp = SQRT((lam+2*mu)/rho) 
    cs = SQRT(mu/rho) 
    !
    RM = 0.0 
    RM(1,1)   = rho*cp**2
    RM(1,10)  = -2*s11
    RM(1,13)  = rho*cp**2
    RM(2,1)   = -rho*(-cp**2+2*cs**2)
    RM(2,4)   = 1
    RM(2,13)  = -rho*(-cp**2+2*cs**2)
    RM(3,1)   = -rho*(-cp**2+2*cs**2)
    RM(3,5)   = 1
    RM(3,13)  = -rho*(-cp**2+2*cs**2)
    RM(4,2)   = rho*cs**2
    RM(4,10)  = -2*s12
    RM(4,12)  = rho*cs**2
    RM(5,6)   = 1
    RM(6,3)   = rho*cs**2
    RM(6,10)  = -2*s13
    RM(6,11)  = rho*cs**2
    RM(7,1)   = cp
    RM(7,10)  = uv(1) 
    RM(7,13)  = -cp
    RM(8,2)   = cs
    RM(8,10)  = uv(2) 
    RM(8,12)  = -cs
    RM(9,3)   = cs
    RM(9,10)  = uv(3) 
    RM(9,11)  = -cs
    RM(10,9)  = 1
    RM(11,8)  = 1
    RM(12,7)  = 1
    RM(13,10) = 1
    RM(14,14) = 1 
    !
    iRM = 0.0 
    iRM(1,1)   = irho/cp**2/2
    iRM(1,7)   = 1/cp/2
    iRM(1,13)  = -1/cp**2/rho*(cp*rho*uv(1)-2*s11)/2
    iRM(2,4)   = irho/cs**2/2
    iRM(2,8)   = 1/cs/2
    iRM(2,13)  = -1/cs**2/rho*(uv(2)*cs*rho-2*s12)/2
    iRM(3,6)   = irho/cs**2/2
    iRM(3,9)   = 1/cs/2
    iRM(3,13)  = -1/cs**2/rho*(uv(3)*cs*rho-2*s13)/2
    iRM(4,1)   = 1/cp**2*(-cp**2+2*cs**2)
    iRM(4,2)   = 1
    iRM(4,13)  = 1/cp**2*2*s11*(-cp**2+2*cs**2) 
    iRM(5,1)   = 1/cp**2*(-cp**2+2*cs**2)
    iRM(5,3)   = 1
    iRM(5,13)  = 1/cp**2*2*s11*(-cp**2+2*cs**2) 
    iRM(6,5)   = 1
    iRM(7,12)  = 1
    iRM(8,11)  = 1
    iRM(9,10)  = 1
    iRM(10,13) = 1. 
    iRM(11,6)  = irho/cs**2/2
    iRM(11,9)  = -1/cs/2
    iRM(11,13) = 1/cs**2/rho*(uv(3)*cs*rho+2*s13)/2
    iRM(12,4)  = irho/cs**2/2
    iRM(12,8)  = -1/cs/2
    iRM(12,13) = 1/cs**2/rho*(uv(2)*cs*rho+2*s12)/2
    iRM(13,1)  = irho/cp**2/2
    iRM(13,7)  = -1/cp/2
    iRM(13,13) = 1/cp**2/rho*(cp*rho*uv(1)+2*s11)/2
    iRM(14,14) = 1 
    !
    R = MATMUL(   TM,  RM ) 
    iR = MATMUL( iRM, iTM ) 
    !
    LM = 0. 
    LM(1,1)   = -cp
    LM(2,2)   = -cs
    LM(3,3)   = -cs
    LM(11,11) = +cs 
    LM(12,12) = +cs
    LM(13,13) = +cp
    LM(14,14) = 0.0 

    !
    DM = MATMUL( R, MATMUL( ABS(LM), iR )  )  
      
    !
END SUBROUTINE DissipationMatrix



