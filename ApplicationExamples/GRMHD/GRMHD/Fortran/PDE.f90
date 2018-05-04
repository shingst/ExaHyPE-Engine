! GRMHD PDE.f90
! Trento (EQNTYPE4)

RECURSIVE SUBROUTINE PDEFlux(f,g,hz,Q)
  USE Parameters, ONLY : nVar, nDim, gamma, DivCleaning_a
  USE iso_c_binding
  IMPLICIT NONE
  REAL :: f(nVar), g(nVar), h(nVar), hz(nVar), Q(nVar), V(nVar)
  REAL, PARAMETER :: epsilon = 1e-14 
  INTENT(IN)  :: Q
  INTENT(OUT) :: f, g, hz
  ! Local Variables 
  REAL :: alpha
  INTEGER :: i, iErr, ii, jj, mm, kk, ll, iDim
  REAL :: rho,vx,vy,vz,p,bx,by,bz,ex,ey,ez
  REAL :: v2,b2,e2,lf,LF2,w,ww,wwx,wwy,wwz,uem,gamma1
  REAL :: irho,ps,pg,rho0,k0,kappa
  REAL :: BQ(3), g_cov(3,3), g_contr(3,3), Fij(3,3), Vtr(3)
  REAL :: S_contr(3),psi,gammaij(6),lapse,shift(3),gv(3),gv_contr(3),gp,gm,ComputeDet
  REAL :: vf_cov(3),vf(3),BV(3),delta(3,3),Vc(nVar),Wim,W_ij
  REAL :: Qv_contr(3),QB_contr(3),vxB(3),vxB_contr(3),BV_contr(3)
  
  INTEGER :: m
  
  !FTensDim = 0
  !RETURN
  
  f = 0
  g = 0
  h = 0
  
  CALL PDECons2Prim(V,Q,iErr)
  !
  gamma1 = gamma/(gamma-1.0)
  rho    = V(1)
  vf_cov = V(2:4)
  p      = V(5)
  !
  BV(1:3) = V(6:8)
  BQ(1:3) = Q(6:8)
  psi = V(9)
  lapse = V(10)
  shift = V(11:13)
  !
  !gammaij = V(14:19) 
  g_cov(1,1) = V(14)
  g_cov(1,2) = V(15)
  g_cov(1,3) = V(16)
  g_cov(2,2) = V(17)
  g_cov(2,3) = V(18)
  g_cov(3,3) = V(19)
  g_cov(2,1) = V(15)
  g_cov(3,1) = V(16)
  g_cov(3,2) = V(18)
  !
  CALL MatrixInverse3x3(g_cov,g_contr,gp)
  gp = SQRT(gp)
  gm = 1./gp
  ! 
  !CALL METRIC(x, lapse, gp, gm, shift, g_cov, g_contr)
  !
  vf     = MATMUL(g_contr,vf_cov)
  Qv_contr = MATMUL(g_contr,Q(2:4))
  QB_contr = MATMUL(g_contr,Q(6:8))
  vxB(1) = vf_cov(2)*BV(3) - vf_cov(3)*BV(2)
  vxB(2) = vf_cov(3)*BV(1) - vf_cov(1)*BV(3)
  vxB(3) = vf_cov(1)*BV(2) - vf_cov(2)*BV(1)
  vxB_contr = MATMUL(g_contr,vxB(1:3))
  BV_contr = MATMUL(g_contr,BV(1:3))
  !
  !v2     = vx**2 + vy**2 + vz**2
  !b2     = bx**2 + by**2 + bz**2
  !e2     = ex**2 + ey**2 + ez**2 
  v2     = vf(1)*vf_cov(1) + vf(2)*vf_cov(2) + vf(3)*vf_cov(3)
  e2     = vxB_contr(1)*vxB(1) + vxB_contr(2)*vxB(2) + vxB_contr(3)*vxB(3)
  b2     = BV_contr(1)*BV(1) + BV_contr(2)*BV(2) + BV_contr(3)*BV(3)
  !
  uem    = 0.5*(b2 + e2) 
  !
  lf     = 1.0/sqrt(1.0 - v2)
  w      = rho + gamma1*p   ! rho*hentalpy
  ww     = w*lf**2
  wwx    = ww*vf(1)
  wwy    = ww*vf(2)
  wwz    = ww*vf(3) 
  !
  ! transport velocity
  Vtr(1:3) = lapse*vf(1:3)-shift(1:3)
  !
  !    Fij(1,1:3) =  (/ f1, g1, h1 /) for Q(6)  
  !    Fij(2,1:3) =  (/ f2, g2, h2 /) for Q(7)  ... without divergence cleaning
  !    Fij(3,1:3) =  (/ f3, g3, h3 /) for Q(8)  
  !
  DO m=1,3
      DO i=1,3
          Fij(i,m) = -Vtr(i)*BQ(m)+Vtr(m)*BQ(i)  ! Fij(i,i) = 0 !!!!
      ENDDO
  ENDDO
  !
  f(1)   = vf(1)*Q(1) !rho*lf   !Q(1) ! rho*lf 
  f(2)   = wwx*vf_cov(1) - vxB_contr(1)*vxB(1) - BV_contr(1)*BV(1) + p + uem
  f(3)   = wwx*vf_cov(2) - vxB_contr(2)*vxB(1) - BV_contr(2)*BV(1) 
  f(4)   = wwx*vf_cov(3) - vxB_contr(3)*vxB(1) - BV_contr(3)*BV(1) 
  f(5)   = Qv_contr(1)-f(1) !wwx - f(1)         ! ADD MAGNETIC COMPONENT
  ! ADD MAGNETIC FIELD and DIV. CLEANING
  f(6)   = Fij(1,1) + V(9) !V(9)
  f(7)   = Fij(2,1) !-ez
  f(8)   = Fij(3,1) !ey  
  f(9)   = DivCleaning_a**2*BQ(1)
  !  lapse&shift&metric fluxes 
  f(10:19) = 0.
  !
  !
  g(1)   = vf(2)*Q(1) !rho*lf   ! rho*lf
  g(2)   = wwy*vf_cov(1) - vxB_contr(1)*vxB(2) - BV_contr(1)*BV(2) 
  g(3)   = wwy*vf_cov(2) - vxB_contr(2)*vxB(2) - BV_contr(2)*BV(2) + p + uem
  g(4)   = wwy*vf_cov(3) - vxB_contr(3)*vxB(2) - BV_contr(3)*BV(2) 
  g(5)   = Qv_contr(2)-g(1) !wwy - g(1)   ! ADD MAGNETIC COMPONENT
  ! ADD MAGNETIC FIELD and DIV. CLEANING
  g(6)   = Fij(1,2) !ez 
  g(7)   = Fij(2,2) + V(9) !V(9) 
  g(8)   = Fij(3,2) !-ex   
  g(9)   = DivCleaning_a**2*BQ(2)
  !  lapse&shift&metric fluxes 
  g(10:19) = 0.
  !
  !
  !
  h(1)   = vf(3)*Q(1) !rho*lf   ! rho*lf   !
  h(2)   = wwz*vf_cov(1) - vxB_contr(1)*vxB(3) - BV_contr(1)*BV(3) 
  h(3)   = wwz*vf_cov(2) - vxB_contr(2)*vxB(3) - BV_contr(2)*BV(3)  
  h(4)   = wwz*vf_cov(3) - vxB_contr(3)*vxB(3) - BV_contr(3)*BV(3) + p + uem
  h(5)   = Qv_contr(3)-h(1) ! wwz - h(1)   !ADD MAGNETIC COMPONENT
  ! ADD MAGNETIC FIELD and DIV. CLEANING
  h(6)   = Fij(1,3)  !-ey  
  h(7)   = Fij(2,3) !ex   
  h(8)   = Fij(3,3) + V(9) !V(9)   
  h(9)   = DivCleaning_a**2*BQ(3)
  !  lapse&shift&metric fluxes 
  h(10:19) = 0.
  ! 
  ! - - - - - - - - - REWRITE THE FOLLOWING 
  f(2:4)   = f(2:4)*gp
  g(2:4)   = g(2:4)*gp
  h(2:4)   = h(2:4)*gp
  ! Remember that Q(:) below contains already the factor gp, which is ok!
  f(1:5)   = lapse*f(1:5) - shift(1)*Q(1:5)
  g(1:5)   = lapse*g(1:5) - shift(2)*Q(1:5)
  h(1:5)   = lapse*h(1:5) - shift(3)*Q(1:5)
  !
  
  IF (nDim == 3) THEN
    hz = h
  ENDIF
  
END SUBROUTINE PDEFlux


RECURSIVE SUBROUTINE PDENCP(BgradQ,Q,gradQ) 
   USE Parameters, ONLY :  nVar, nDim, gamma
   IMPLICIT NONE
   REAL :: AQx(nVar), BQy(nVar), CQz(nVar), Qx(nVar), Qy(nVar), Qz(nVar) 
   REAL :: uu, vv, ww, rho0, kappa, k0, rho, Vp(nVar), k, sigma     
   INTEGER :: i, j, l, m, n, iErr, qq, ii, jj, kk, ll, mm, nn, count
   REAL, PARAMETER :: epsilon = 1e-14
   ! 11. Oct 21:40: This was a matrix BGradQ(nVar, nDim) but is a vector in spaceTimePredictorNonlinear
   REAL, INTENT(OUT) :: BgradQ(nVar)
   REAL, INTENT(IN)  :: gradQ(nVar, nDim)
   REAL, INTENT(IN)  :: Q(nVar)
   REAL :: p, ih,h,pg, A(3,3), u(3)   
   REAL :: vx,vy,cs,cl,k3,ff,fa 
   REAL :: BQ(3), g_cov(3,3), g_contr(3,3), Fij(3,3), Vtr(3)
   REAL :: S_contr(3),psi,gammaij(6),lapse,shift(3),gv(3),gv_contr(3),gp,gm,ComputeDet
   REAL :: vf_cov(3),vf(3),BV(3),delta(3,3),Vc(nVar),Wim,W_ij
   REAL :: Qv_contr(3),QB_contr(3),vxB(3),vxB_contr(3),BV_contr(3)
   REAL :: v2,b2,e2,lf,LF2,w,wwx,wwy,wwz,uem,gamma1

   
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
  !psi = Q(9)
  lapse = Q(10)
  !gS(1:3) = Q(2:4) 
  shift = Q(11:13)
  !
  gammaij = Q(14:19) 
  g_cov(1,1) = Q(14)
  g_cov(1,2) = Q(15)
  g_cov(1,3) = Q(16)
  g_cov(2,2) = Q(17)
  g_cov(2,3) = Q(18)
  g_cov(3,3) = Q(19)
  g_cov(2,1) = Q(15)
  g_cov(3,1) = Q(16)
  g_cov(3,2) = Q(18)
  !
  delta = 0.
  DO i=1,3
      DO j=1,3
          IF(i.eq.j) delta(i,j) = 1.0
      ENDDO
  ENDDO 
  !
  CALL MatrixInverse3x3(g_cov,g_contr,gp)
  gp = SQRT(gp)
  gm = 1./gp
  ! 
  CALL PDECons2Prim(Vc,Q,iErr)
  gamma1 = gamma/(gamma-1.0)
  rho    = Vc(1)
  vf_cov = Vc(2:4)
  p      = Vc(5)
  !
  BV(1:3) = Vc(6:8)
  psi = Vc(9) 
  !
  Qv_contr = MATMUL(g_contr,Q(2:4))
  QB_contr = MATMUL(g_contr,Q(6:8))
  vxB(1) = vf_cov(2)*BV(3) - vf_cov(3)*BV(2)
  vxB(2) = vf_cov(3)*BV(1) - vf_cov(1)*BV(3)
  vxB(3) = vf_cov(1)*BV(2) - vf_cov(2)*BV(1)
  vxB_contr = MATMUL(g_contr,vxB(1:3))
  BV_contr = MATMUL(g_contr,BV(1:3))
  !
  !v2     = vx**2 + vy**2 + vz**2
  !b2     = bx**2 + by**2 + bz**2
  !e2     = ex**2 + ey**2 + ez**2 
  v2     = vf(1)*vf_cov(1) + vf(2)*vf_cov(2) + vf(3)*vf_cov(3)
  e2     = vxB_contr(1)*vxB(1) + vxB_contr(2)*vxB(2) + vxB_contr(3)*vxB(3)
  b2     = BV_contr(1)*BV(1) + BV_contr(2)*BV(2) + BV_contr(3)*BV(3)
  !
  uem    = 0.5*(b2 + e2) 
  !
  !
  vf     = MATMUL(g_contr,vf_cov)
  S_contr = MATMUL(g_contr,Q(2:4))
  !gv_contr = MATMUL(g_contr,gv)
  lf     = 1.0/sqrt(1.0 - v2)
  w      = rho + gamma1*p   ! rho*hentalpy
  ww     = w*lf**2          ! rho*hentalpy*Lorentz^2 
  !
  !DO j=1,3
  AQx = 0.
  BQy = 0.
  CQz = 0.
      count=0
      DO i=1,3 
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=1
                !------ 
                W_ij = ww*vf_cov(i)*vf(j)-vxB(i)*vxB_contr(j)-BV(i)*BV_contr(j)+(p+uem)*delta(i,j)
                !
                AQx(1+j) = AQx(1+j) - Q(1+i)*Qx(10+i)  ! Q(11:13)  shift(i) or shift_contr(i)
                AQx(5) = AQx(5) - gp*W_ij*Qx(10+i)     ! Q(11:13)  shift(i) or shift_contr(i)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=2
                !------ 
                W_ij = ww*vf_cov(i)*vf(j)-vxB(i)*vxB_contr(j)-BV(i)*BV_contr(j)+(p+uem)*delta(i,j)
                !
                BQy(1+j) = BQy(1+j) - Q(1+i)*Qy(10+i)   ! Q(11:13)  shift(i) or shift_contr(i)
                BQy(5) = BQy(5) - gp*W_ij*Qy(10+i)     ! Q(11:13)  shift(i) or shift_contr(i)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=3
                !------
                W_ij = ww*vf_cov(i)*vf(j)-vxB(i)*vxB_contr(j)-BV(i)*BV_contr(j)+(p+uem)*delta(i,j)
                !
                CQz(1+j) = CQz(1+j) - Q(1+i)*Qz(10+i)   ! Q(11:13)  shift(i) or shift_contr(i)
                CQz(5) = CQz(5) - gp*W_ij*Qz(10+i)     ! Q(11:13)  shift(i) or shift_contr(i)
                !
          DO m=1,3
            IF(m.GE.i) THEN  
                count=count+1
                !
                Wim = ww*vf(i)*vf(m)-vxB_contr(i)*vxB_contr(m)-BV_contr(i)*BV_contr(m)+(p+uem)*g_contr(i,m)
                Wim = Wim + (1.0 - delta(i,m))*(ww*vf(m)*vf(i)-vxB_contr(m)*vxB_contr(i)- &
		      BV_contr(m)*BV_contr(i)+(p+uem)*g_contr(m,i))
		      ! account also of the remaining symmetric components of gamma for i.NE.m.
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=1
                !------ 
                AQx(1+j) = AQx(1+j) - 0.5*gp*lapse*Wim*Qx(13+count)  ! Q(14:19) gammaij(count) or  g_cov(i,m)
                AQx(5) = AQx(5) - 0.5*gp*Wim*shift(j)*Qx(13+count)   ! Q(14:19) gammaij(count) or  g_cov(i,m)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=2
                !------ 
                BQy(1+j) = BQy(1+j) - 0.5*gp*lapse*Wim*Qy(13+count)  ! Q(14:19) gammaij(count) or  g_cov(i,m)
                BQy(5) = BQy(5) - 0.5*gp*Wim*shift(j)*Qy(13+count)   ! Q(14:19) gammaij(count) or  g_cov(i,m)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=3
                !------
                CQz(1+j) = CQz(1+j) - 0.5*gp*lapse*Wim*Qz(13+count)  ! Q(14:19) gammaij(count) or  g_cov(i,m)
                CQz(5) = CQz(5) - 0.5*gp*Wim*shift(j)*Qz(13+count)   ! Q(14:19) gammaij(count) or  g_cov(i,m)
                !
            ENDIF
          ENDDO
      ENDDO
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=1
    !------ 
    AQx(1+j) = AQx(1+j) + (Q(5)+Q(1))*Qx(10)    ! Q(10) or lapse
    AQx(5) = AQx(5) + S_contr(j)*Qx(10)         !  Q(10) or lapse
    !
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=2
    !------ 
    BQy(1+j) = BQy(1+j) + (Q(5)+Q(1))*Qy(10)    ! Q(10) or lapse
    BQy(5) = BQy(5) + S_contr(j)*Qy(10)         !  Q(10) or lapse
    !
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=3
    !------
    CQz(1+j) = CQz(1+j) + (Q(5)+Q(1))*Qz(10)    ! Q(10) or lapse
    CQz(5) = CQz(5) + S_contr(j)*Qz(10)         !  Q(10) or lapse
    !

    BgradQ = AQx + BQy + CQz

    ! Check for NaN. Being a quantum computer, the presence of such a
    ! check changes the behaviour of the computer.
!    DO i=1,nVar
!	IF(.not.(BgradQ(i).eq.BgradQ(i))) THEN
!		PRINT *, "NaN: ", BgradQ
!	END IF
 !   END DO
END SUBROUTINE PDENCP


RECURSIVE SUBROUTINE PDEEigenvalues(L,Q,n)
  USE Parameters, ONLY :  nVar, nDim
  USE iso_c_binding
  IMPLICIT NONE
  REAL :: L(nVar), n(nDim), Q(nVar), V(nVar)
  INTENT(IN)  :: Q,n
  INTENT(OUT) :: L 
  ! Local Variables 
  REAL :: Pi
  Pi = ACOS(-1.0)

  L = 1

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

  ! EQNTYPE4
  Name(1) = 'rho' 
  Name(2) = 'u' 
  Name(3) = 'v'
  Name(4) = 'w'
  Name(5) = 'p' 
  Name(6) = 'Bx' 
  Name(7) = 'By' 
  Name(8) = 'Bz' 
  Name(9) = 'psi' 
  Name(10) = 'lapse' 
  Name(11) = 'Shift^1' 
  Name(12) = 'Shift^2' 
  Name(13) = 'Shift^3' 
  Name(14) = 'Gamma_11' 
  Name(15) = 'Gamma_12' 
  Name(16) = 'Gamma_13' 
  Name(17) = 'Gamma_22' 
  Name(18) = 'Gamma_23' 
  Name(19) = 'Gamma_33' 

END SUBROUTINE PDEVarName

RECURSIVE SUBROUTINE Kreuzprodukt(res,vec_1,vec_2)
  !--------------------------------------------------------------------------
  IMPLICIT NONE                                            
  !--------------------------------------------------------------------------
  REAL :: res(3)
  REAL :: vec_1(3)
  REAL :: vec_2(3)
  !--------------------------------------------------------------------------
  INTENT(IN) :: vec_1,vec_2
  !--------------------------------------------------------------------------
  res(1) = vec_1(2)*vec_2(3) - vec_1(3)*vec_2(2)
  res(2) = vec_1(3)*vec_2(1) - vec_1(1)*vec_2(3)
  res(3) = vec_1(1)*vec_2(2) - vec_1(2)*vec_2(1)
END SUBROUTINE Kreuzprodukt


RECURSIVE SUBROUTINE PDEMatrixB(An,Q,nv) 
  USE Parameters, ONLY : nVar, nDim, gamma
  USE iso_c_binding
  IMPLICIT NONE
  ! Argument list 
  REAL :: An(nVar,nVar)
  REAL :: Q(nVar), nv(nDim) 
  INTENT(IN)  :: Q,nv
  INTENT(OUT) :: An  
  ! Local variables
  INTEGER :: i,j,m, count, iErr
  REAL :: A(nVar,nVar), B(nVar,nVar), C(nVar,nVar), Vp(nVar) 
REAL :: BQ(3), g_cov(3,3), g_contr(3,3), Fij(3,3), Vtr(3)
   REAL :: S_contr(3),psi,gammaij(6),lapse,shift(3),gv(3),gv_contr(3),gp,gm,ComputeDet
   REAL :: vf_cov(3),vf(3),BV(3),delta(3,3),Vc(nVar),Wim,W_ij
   REAL :: Qv_contr(3),QB_contr(3),vxB(3),vxB_contr(3),BV_contr(3)
   REAL :: v2,b2,e2,lf,LF2,w,ww,wwx,wwy,wwz,uem,gamma1
   REAL :: rho, p

  !An = 0
  !RETURN

  !psi = Q(9)
  lapse = Q(10)
  !gS(1:3) = Q(2:4) 
  shift = Q(11:13)
  !
  gammaij = Q(14:19) 
  g_cov(1,1) = Q(14)
  g_cov(1,2) = Q(15)
  g_cov(1,3) = Q(16)
  g_cov(2,2) = Q(17)
  g_cov(2,3) = Q(18)
  g_cov(3,3) = Q(19)
  g_cov(2,1) = Q(15)
  g_cov(3,1) = Q(16)
  g_cov(3,2) = Q(18)
  !
  delta = 0.
  DO i=1,3
      DO j=1,3
          IF(i.eq.j) delta(i,j) = 1.0
      ENDDO
  ENDDO 
  !
  CALL MatrixInverse3x3(g_cov,g_contr,gp)
  gp = SQRT(gp)
  gm = 1./gp
  ! 
  CALL PDECons2Prim(Vc,Q,iErr)
  gamma1 = gamma/(gamma-1.0)
  rho    = Vc(1)
  vf_cov = Vc(2:4)
  p      = Vc(5)
  !
  BV(1:3) = Vc(6:8)
  psi = Vc(9) 
  !
  vf     = MATMUL(g_contr,vf_cov)
  S_contr = MATMUL(g_contr,Q(2:4))
  QB_contr = MATMUL(g_contr,Q(6:8))
  vxB(1) = vf_cov(2)*BV(3) - vf_cov(3)*BV(2)
  vxB(2) = vf_cov(3)*BV(1) - vf_cov(1)*BV(3)
  vxB(3) = vf_cov(1)*BV(2) - vf_cov(2)*BV(1)
  vxB_contr = MATMUL(g_contr,vxB(1:3))
  BV_contr = MATMUL(g_contr,BV(1:3))
  !
  !v2     = vx**2 + vy**2 + vz**2
  !b2     = bx**2 + by**2 + bz**2
  !e2     = ex**2 + ey**2 + ez**2 
  v2     = vf(1)*vf_cov(1) + vf(2)*vf_cov(2) + vf(3)*vf_cov(3)
  e2     = vxB_contr(1)*vxB(1) + vxB_contr(2)*vxB(2) + vxB_contr(3)*vxB(3)
  b2     = BV_contr(1)*BV(1) + BV_contr(2)*BV(2) + BV_contr(3)*BV(3)
  !
  uem    = 0.5*(b2 + e2) 
  !
  !gv_contr = MATMUL(g_contr,gv)
  lf     = 1.0/sqrt(1.0 - v2)
  w      = rho + gamma1*p   ! rho*hentalpy
  ww     = w*lf**2          ! rho*hentalpy*Lorentz^2 
  !
  !DO j=1,3
  A = 0.
  B = 0.
  C = 0.
    !lapse
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=1
    !------ 
    A(1+j,10) = + (Q(5)+Q(1))   ! Q(10) or lapse
    A(5,10) =  S_contr(j)     !  Q(10) or lapse
    !
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=2
    !------ 
    B(1+j,10) = + (Q(5)+Q(1))   ! Q(10) or lapse
    B(5,10) =  S_contr(j)     !  Q(10) or lapse
    ! 
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=3
    !------
    C(1+j,10) = + (Q(5)+Q(1))   ! Q(10) or lapse
    C(5,10) =  S_contr(j)     !  Q(10) or lapse
    ! 
    count=0
    DO i=1,3
        ! shift
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        j=1
        !------ 
        W_ij = ww*vf_cov(i)*vf(j)-vxB(i)*vxB_contr(j)-BV(i)*BV_contr(j)+(p+uem)*delta(i,j)
        !
        A(1+j,10+i) = - Q(1+i)  ! Q(11:13)  shift(i) or shift_contr(i)
        A(5,10+i) = - gp*W_ij    ! Q(11:13)  shift(i) or shift_contr(i)
        !
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        j=2
        !------ 
        W_ij = ww*vf_cov(i)*vf(j)-vxB(i)*vxB_contr(j)-BV(i)*BV_contr(j)+(p+uem)*delta(i,j)
        !
        B(1+j,10+i) = - Q(1+i)  ! Q(11:13)  shift(i) or shift_contr(i)
        B(5,10+i) = - gp*W_ij    ! Q(11:13)  shift(i) or shift_contr(i)
        ! 
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        j=3
        !------
        W_ij = ww*vf_cov(i)*vf(j)-vxB(i)*vxB_contr(j)-BV(i)*BV_contr(j)+(p+uem)*delta(i,j)
        !
        C(1+j,10+i) = - Q(1+i)  ! Q(11:13)  shift(i) or shift_contr(i)
        C(5,10+i) = - gp*W_ij    ! Q(11:13)  shift(i) or shift_contr(i)
        ! 
          DO m=1,3
            IF(m.GE.i) THEN  
                !metric
                count=count+1
                !
                Wim = ww*vf(i)*vf(m)-vxB_contr(i)*vxB_contr(m)-BV_contr(i)*BV_contr(m)+(p+uem)*g_contr(i,m)
                Wim = Wim + (1.0 - delta(i,m))*(ww*vf(m)*vf(i)-vxB_contr(m)*vxB_contr(i)-&
                   BV_contr(m)*BV_contr(i)+(p+uem)*g_contr(m,i)) 
                   ! account also of the remaining symmetric components of gamma for i.NE.m.
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=1
                !------ 
                A(1+j,13+count) = - 0.5*gp*lapse*Wim  ! Q(14:19) gammaij(count) or  g_cov(i,m)
                A(5,13+count) = - 0.5*gp*Wim*shift(j)   ! Q(14:19) gammaij(count) or  g_cov(i,m)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=2
                !------ 
                B(1+j,13+count) = - 0.5*gp*lapse*Wim  ! Q(14:19) gammaij(count) or  g_cov(i,m)
                B(5,13+count) = - 0.5*gp*Wim*shift(j)   ! Q(14:19) gammaij(count) or  g_cov(i,m)
                ! 
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=3
                !------
                C(1+j,13+count) = - 0.5*gp*lapse*Wim  ! Q(14:19) gammaij(count) or  g_cov(i,m)
                C(5,13+count) = - 0.5*gp*Wim*shift(j)   ! Q(14:19) gammaij(count) or  g_cov(i,m)
                ! 
            ENDIF
          ENDDO
      ENDDO
  !ENDDO 
  !
  IF(nDim == 3) THEN
     An = A*nv(1) + B*nv(2) + C*nv(3) 
  ELSE
     An = A*nv(1) + B*nv(2)
  ENDIF
  
END SUBROUTINE PDEMatrixB


