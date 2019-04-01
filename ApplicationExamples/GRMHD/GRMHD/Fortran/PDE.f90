! GRMHD PDE.f90
! Trento (EQNTYPE4)

RECURSIVE SUBROUTINE PDEFlux(ff,g,h,Q)
  USE Parameters, ONLY : d, nDim, gamma, DivCleaning_a
  USE iso_c_binding
  IMPLICIT NONE
#ifdef VECTOR    
#ifdef AVX512 
  INTEGER, PARAMETER :: nVar = 24                           ! The number of variables of the PDE system 
#else   
  INTEGER, PARAMETER :: nVar = 20                           ! The number of variables of the PDE system 
#endif 
#else
  INTEGER, PARAMETER :: nVar = 19                           ! The number of variables of the PDE system 
#endif
  REAL :: ff(nVar), g(nVar), h(nVar), Q(nVar) 
  INTENT(IN)  :: Q
  INTENT(OUT) :: ff, g, h
  ! Local Variables 
  REAL :: F(nVar,d) 
    INTEGER :: iErr 
    INTEGER :: i,j,k,m
    REAL :: V(nVar) 
    REAL :: p, irho, lam, mu 
    REAL :: k1, k2,  ggg, e, ds, c, xi, sk, sknl, alpha, fa, k0, beta0(3)  
    REAL :: gamma1, rho, vx, vy, vz, bx, by, bz, ex, ey, ez, v2, b2, e2, lf, w, ww, uem, wwx, wwy, wwz 
    REAL :: gp, gm, g_cov(3,3), g_contr(3,3), vxB_cov(3), vxB_contr(3), B_cov(3), BQ(3), vtr(3), v_contr(3), lapse, shift(3)    
    REAL :: Fij(3,3), v_cov(3), S_contr(3), QB_contr(3), B_contr(3), psi 
    !
  
  !FTensDim = 0
  !RETURN
  
  ff = 0
  g = 0
  h = 0
  !
  F = 0.0 
  !
  CALL PDECons2Prim(V,Q,iErr)
  !
  gamma1 = gamma/(gamma-1.0)
  rho    = V(1)
  !
  DO i=1,3
	  v_cov(i) = V(1+i)
	  B_contr(i) = V(5+i)  		! B is contravariant
	  QB_contr(i) = Q(5+i)  		! B is contravariant
	  shift(i) = V(10+i)
  ENDDO
  p      = V(5)
  psi = V(9)
  lapse = V(10)
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
  v_contr = MATMUL(g_contr,v_cov)
  S_contr = MATMUL(g_contr,Q(2:4))
  BQ = MATMUL(g_cov,QB_contr)
  B_cov = MATMUL(g_cov,B_contr)
  ! COVARIANT =>
  vxB_cov(1) = v_contr(2)*B_contr(3) - v_contr(3)*B_contr(2)
  vxB_cov(2) = v_contr(3)*B_contr(1) - v_contr(1)*B_contr(3)
  vxB_cov(3) = v_contr(1)*B_contr(2) - v_contr(2)*B_contr(1) 
  DO i=1,3
	  vxB_cov(i) = gp*vxB_cov(i)
  ENDDO
  !
  ! CONTRAVARIANT =>
  vxB_contr(1) = v_cov(2)*B_cov(3) - v_cov(3)*B_cov(2)
  vxB_contr(2) = v_cov(3)*B_cov(1) - v_cov(1)*B_cov(3)
  vxB_contr(3) = v_cov(1)*B_cov(2) - v_cov(2)*B_cov(1) 
  DO i=1,3
	  vxB_contr(i) = gm*vxB_contr(i)
  ENDDO
  !
  v2     = v_contr(1)*v_cov(1) + v_contr(2)*v_cov(2) + v_contr(3)*v_cov(3)
  e2     = vxB_contr(1)*vxB_cov(1) + vxB_contr(2)*vxB_cov(2) + vxB_contr(3)*vxB_cov(3)
  b2     = B_contr(1)*B_cov(1) + B_contr(2)*B_cov(2) + B_contr(3)*B_cov(3)
  !
  uem    = 0.5*(b2 + e2) 
  !
  lf     = 1.0/sqrt(1.0 - v2)
  w      = rho + gamma1*p   ! rho*hentalpy
  ww     = w*lf**2
  wwx    = ww*v_contr(1)
  wwy    = ww*v_contr(2)
  wwz    = ww*v_contr(3)
  !
  ! transport velocity
  DO i=1,3
	  Vtr(i) = lapse*v_contr(i)-shift(i)
  ENDDO
  !
  !    Fij(1,1:3) =  (/ f1, g1, h1 /) for Q(6)  
  !    Fij(2,1:3) =  (/ f2, g2, h2 /) for Q(7)  ... without divergence cleaning
  !    Fij(3,1:3) =  (/ f3, g3, h3 /) for Q(8)  
  !
  DO m=1,3
      DO i=1,3
          Fij(i,m) = Vtr(i)*QB_contr(m)-Vtr(m)*QB_contr(i)  ! Fij(i,i) = 0 !!!! ! NB: this is contravariant !!!!
      ENDDO
  ENDDO
  !
  F(1,1)   = v_contr(1)*Q(1) 
  F(2,1)   = wwx*v_cov(1) - vxB_contr(1)*vxB_cov(1) - B_contr(1)*B_cov(1) + p + uem
  F(3,1)   = wwx*v_cov(2) - vxB_contr(1)*vxB_cov(2) - B_contr(1)*B_cov(2) 
  F(4,1)   = wwx*v_cov(3) - vxB_contr(1)*vxB_cov(3) - B_contr(1)*B_cov(3) 
  F(5,1)   = S_contr(1)-F(1,1) 
  ! 
  F(6,1)   = Fij(1,1) !+ V(9)  
  F(7,1)   = Fij(1,2) 
  F(8,1)   = Fij(1,3)  
  !  
#ifdef COVCLEAN     
  F(9,1)   = DivCleaning_a**2*QB_contr(1)   ! flat Newtonian, flat coordinates
#else
  F(9,1)   = 0.  
#endif 
  !  lapse&shift&metric fluxes 
  DO i=10,nVar
	F(i,1) = 0.
  ENDDO
  !
  !
  F(1,2)   = v_contr(2)*Q(1) 
  F(2,2)   = wwy*v_cov(1) - vxB_contr(2)*vxB_cov(1) - B_contr(2)*B_cov(1) 
  F(3,2)   = wwy*v_cov(2) - vxB_contr(2)*vxB_cov(2) - B_contr(2)*B_cov(2) + p + uem
  F(4,2)   = wwy*v_cov(3) - vxB_contr(2)*vxB_cov(3) - B_contr(2)*B_cov(3) 
  F(5,2)   = S_contr(2)-F(1,2)  
  ! 
  F(6,2)   = Fij(2,1)  
  F(7,2)   = Fij(2,2) !+ V(9) 
  F(8,2)   = Fij(2,3)  
  ! 
#ifdef COVCLEAN     
  F(9,2)   = DivCleaning_a**2*QB_contr(2)   
#else
  F(9,2)   = 0. !QB_contr(2) 
#endif
  !  lapse&shift&metric fluxes 
  DO i=10,nVar
	F(i,2) = 0.
  ENDDO
  !
  !
  !
  F(1,3)   = v_contr(3)*Q(1) 
  F(2,3)   = wwz*v_cov(1) - vxB_contr(3)*vxB_cov(1) - B_contr(3)*B_cov(1) 
  F(3,3)   = wwz*v_cov(2) - vxB_contr(3)*vxB_cov(2) - B_contr(3)*B_cov(2)   
  F(4,3)   = wwz*v_cov(3) - vxB_contr(3)*vxB_cov(3) - B_contr(3)*B_cov(3) + p + uem
  F(5,3)   = S_contr(3)-F(1,3) 
  ! 
  F(6,3)   = Fij(3,1)
  F(7,3)   = Fij(3,2)
  F(8,3)   = Fij(3,3) !+ V(9) 
  !
#ifdef COVCLEAN     
  F(9,3)   = DivCleaning_a**2*QB_contr(3)   ! flat Newtonian, flat coordinates
#else
  F(9,3)   = 0. !QB_contr(3) 
#endif 
  !  lapse&shift&metric fluxes 
  DO i=10,nVar
	F(i,3) = 0.
  ENDDO
  ! 
  ! - - - - - - - - - 
  DO i=2,4
	  F(i,1)   = F(i,1)*gp
	  F(i,2)   = F(i,2)*gp
	  F(i,3)   = F(i,3)*gp
  ENDDO
  ! Remember that Q(:) below contains already the factor gp, which is ok!
  
  DO i=1,5
	  F(i,1)   = lapse*F(i,1) - shift(1)*Q(i)
	  F(i,2)   = lapse*F(i,2) - shift(2)*Q(i)
	  F(i,3)   = lapse*F(i,3) - shift(3)*Q(i)
  ENDDO 
  !
#ifdef COVCLEAN   
  F(9,1)   = lapse*F(9,1) - shift(1)*Q(9)
  F(9,2)   = lapse*F(9,1) - shift(2)*Q(9)
  F(9,3)   = lapse*F(9,1) - shift(3)*Q(9)
#endif  
  !   
  ff = F(:,1)
  g = F(:,2)
  !
  IF (nDim == 3) THEN
    h = F(:,3)
  ENDIF
  
END SUBROUTINE PDEFlux


RECURSIVE SUBROUTINE PDENCP(BgradQ,Q,gradQ) 
  USE Parameters, ONLY : d, nDim, gamma, DivCleaning_a
   IMPLICIT NONE
#ifdef VECTOR    
#ifdef AVX512 
  INTEGER, PARAMETER :: nVar = 24                           ! The number of variables of the PDE system 
#else   
  INTEGER, PARAMETER :: nVar = 20                           ! The number of variables of the PDE system 
#endif 
#else
  INTEGER, PARAMETER :: nVar = 19                           ! The number of variables of the PDE system 
#endif
    ! Argument list 
    REAL, INTENT(IN)  :: Q(nVar), gradQ(nVar,d)   
    REAL, INTENT(OUT) :: BgradQ(nVar) 
    ! Local variables 
    REAL :: p, irho
    REAL :: Qx(nVar), Qy(nVar), Qz(nVar) 
    REAL :: e,g_cov(3,3),g_contr(3,3)
    REAL :: AQx(nVar), BQy(nVar), CQz(nVar) 
    REAL :: lapse, shift(3), gammaij(6), delta(3,3), B_cov(3), vxB_cov(3), vxb_contr(3), psi, S_contr(3), qb_contr(3), B_contr(3) 
    REAL :: v2,v_contr(3),uem,b2,e2,gp,gm,vc(nVar),lf,w,ww,gamma1,rho,v_cov(3), w_ij, wim    
    INTEGER :: i,j,k,l,m,iErr, ccount    
    !
    !BgradQ = 0.0 
	DO i=1,nVar
		Qx(i) = gradQ(i,1) 
		Qy(i) = gradQ(i,2)
		Qz(i) = gradQ(i,3)
		AQx(i) = 0.
		BQy(i) = 0.
		CQz(i) = 0.
	ENDDO
    !
  lapse = Q(10)
  !
	shift(1) = Q(11)
	shift(2) = Q(12)
	shift(3) = Q(13)
	!
	DO i=1,6
		gammaij(i) = Q(13+i) 
	ENDDO
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
	DO i=1,3
	  DO j=1,3
		  delta(i,j) = 0.0
	  ENDDO 
	ENDDO 
  DO i=1,3
     delta(i,i) = 1.0 
  ENDDO 
  !
  CALL MatrixInverse3x3(g_cov,g_contr,gp)
  gp = SQRT(gp)
  gm = 1./gp
  ! 
  CALL PDECons2Prim(Vc,Q,iErr)
  gamma1 = gamma/(gamma-1.0)
  rho    = Vc(1)
	v_cov = Vc(2:4)
  p      = Vc(5)
  !
	!B_cov(1:3) = Vc(6:8)
	DO i=1,3
	  B_contr(i) = Vc(5+i) ! contravariant!
	  QB_contr(i) = Q(5+i)
	ENDDO

	B_cov = MATMUL(g_cov,B_contr)
  psi = Vc(9) 
  ! 
	v_contr       = MATMUL(g_contr,v_cov)
	S_contr = MATMUL(g_contr,Q(2:4))
  !
	!! COVARIANT =>
	vxB_cov(1) = v_contr(2)*B_contr(3) - v_contr(3)*B_contr(2)
	vxB_cov(2) = v_contr(3)*B_contr(1) - v_contr(1)*B_contr(3)
	vxB_cov(3) = v_contr(1)*B_contr(2) - v_contr(2)*B_contr(1) 
	DO i=1,3
	  vxB_cov(i) = gp*vxB_cov(i)
	ENDDO
  !
	! CONTRAVARIANT =>
	vxB_contr(1) = v_cov(2)*B_cov(3) - v_cov(3)*B_cov(2)
	vxB_contr(2) = v_cov(3)*B_cov(1) - v_cov(1)*B_cov(3)
	vxB_contr(3) = v_cov(1)*B_cov(2) - v_cov(2)*B_cov(1)
	DO i=1,3
		vxB_contr(i) = gm*vxB_contr(i)
	ENDDO
  !v2     = vx**2 + vy**2 + vz**2
  !b2     = bx**2 + by**2 + bz**2
  !e2     = ex**2 + ey**2 + ez**2 
	v2     = v_contr(1)*v_cov(1) + v_contr(2)*v_cov(2) + v_contr(3)*v_cov(3)
	e2     = vxB_contr(1)*vxB_cov(1) + vxB_contr(2)*vxB_cov(2) + vxB_contr(3)*vxB_cov(3)
	b2     = B_contr(1)*B_cov(1) + B_contr(2)*B_cov(2) + B_contr(3)*B_cov(3)
  !
  uem    = 0.5*(b2 + e2) 
  !
	S_contr = MATMUL(g_contr,Q(2:4))
  !
  lf     = 1.0/sqrt(1.0 - v2)
  w      = rho + gamma1*p   ! rho*hentalpy
  ww     = w*lf**2          ! rho*hentalpy*Lorentz^2 
  !
	ccount=0
      DO i=1,3 
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=1
                !------ 
			W_ij = ww*v_cov(i)*v_contr(j)-vxB_cov(i)*vxB_contr(j)-B_cov(i)*B_contr(j)+(p+uem)*delta(i,j)
                !
                AQx(1+j) = AQx(1+j) - Q(1+i)*Qx(10+i)  ! Q(11:13)  shift(i) or shift_contr(i)
                AQx(5) = AQx(5) - gp*W_ij*Qx(10+i)     ! Q(11:13)  shift(i) or shift_contr(i)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=2
                !------ 
			W_ij = ww*v_cov(i)*v_contr(j)-vxB_cov(i)*vxB_contr(j)-B_cov(i)*B_contr(j)+(p+uem)*delta(i,j)
                !
                BQy(1+j) = BQy(1+j) - Q(1+i)*Qy(10+i)   ! Q(11:13)  shift(i) or shift_contr(i)
                BQy(5) = BQy(5) - gp*W_ij*Qy(10+i)     ! Q(11:13)  shift(i) or shift_contr(i)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=3
                !------
			W_ij = ww*v_cov(i)*v_contr(j)-vxB_cov(i)*vxB_contr(j)-B_cov(i)*B_contr(j)+(p+uem)*delta(i,j)
                !
                CQz(1+j) = CQz(1+j) - Q(1+i)*Qz(10+i)   ! Q(11:13)  shift(i) or shift_contr(i)
                CQz(5) = CQz(5) - gp*W_ij*Qz(10+i)     ! Q(11:13)  shift(i) or shift_contr(i)
                !
          DO m=1,3
            IF(m.GE.i) THEN  
			ccount=ccount+1
                !
			Wim = ww*v_contr(i)*v_contr(m)-vxB_contr(i)*vxB_contr(m)-B_contr(i)*B_contr(m)+(p+uem)*g_contr(i,m)
			!Wim = Wim + (1.0 - delta(i,m))*(ww*v_contr(m)*v_contr(i)-vxB_contr(m)*vxB_contr(i)-B_contr(m)*B_contr(i)+(p+uem)*g_contr(m,i))  ! account also of the remaining symmetric components of gamma for i.NE.m.
			IF(i.NE.m) THEN !TWICE! account also of the remaining symmetric components of gamma §(memo: g_contr and Wim are symmetric)
                    Wim = 2.0*Wim   
                ENDIF
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=1
                !------ 
			AQx(1+j) = AQx(1+j) - 0.5*gp*lapse*Wim*Qx(13+ccount)  ! Q(14:19) gammaij(ccount) or  g_cov(i,m)
			AQx(5) = AQx(5) - 0.5*gp*Wim*shift(j)*Qx(13+ccount)   ! Q(14:19) gammaij(ccount) or  g_cov(i,m)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=2
                !------ 
			BQy(1+j) = BQy(1+j) - 0.5*gp*lapse*Wim*Qy(13+ccount)  ! Q(14:19) gammaij(ccount) or  g_cov(i,m)
			BQy(5) = BQy(5) - 0.5*gp*Wim*shift(j)*Qy(13+ccount)   ! Q(14:19) gammaij(ccount) or  g_cov(i,m)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=3
                !------
			CQz(1+j) = CQz(1+j) - 0.5*gp*lapse*Wim*Qz(13+ccount)  ! Q(14:19) gammaij(ccount) or  g_cov(i,m)
			CQz(5) = CQz(5) - 0.5*gp*Wim*shift(j)*Qz(13+ccount)   ! Q(14:19) gammaij(ccount) or  g_cov(i,m)
                !
#ifdef COVCLEAN     
				CQz(9) = CQz(9) + 0.5*Q(9)*gammaij(ccount)*shift(j)*Qz(13+ccount)   ! D.C.
                CQz(6:8) = CQz(6:8) - 0.5*lapse*Q(9)*g_contr(1:3,j)*gammaij(ccount)*Qz(13+ccount)
                !
                IF(i.NE.m) THEN !TWICE!  account  also of the remaining symmetric components of gamma §(memo: g_contr and Wim are symmetric)
                    CQz(9) = CQz(9) + 0.5*Q(9)*gammaij(ccount)*shift(j)*Qz(13+ccount)   ! D.C.
                    CQz(6:8) = CQz(6:8) - 0.5*lapse*Q(9)*g_contr(1:3,j)*gammaij(ccount)*Qz(13+ccount)
                ENDIF
#endif  		
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
#ifdef COVCLEAN     
	AQx(9) = AQx(9) - Q(5+j)*Qx(10)             ! D.C.
    AQx(9) = AQx(9) + Q(9)*Qx(10+j)             ! D.C.   phi * d(shift^x)/dx
    AQx(6:8) = AQx(6:8) - Q(11:13)*Qx(5+j)             ! D.C.   shift * d(gp*B^x)/dx
    AQx(6:8) = AQx(6:8) + lapse*g_contr(1:3,j)*Qx(9)  ! D.C.  
#else
    AQx(6:8) = AQx(6:8) - Q(11:13)*Qx(5+j)             ! D.C.   shift * d(gp*B^x)/dx 
    AQx(6:8) = AQx(6:8) + lapse*gp*g_contr(1:3,j)*Qx(9)  ! D.C.   
    AQx(9) = AQx(9) + gm*lapse*DivCleaning_a**2*Qx(5+j)             ! D.C.   shift * d(gp*B^x)/dx 
    AQx(9) = AQx(9) - Q(10+j)*Qx(9)             ! D.C.   shift * d(gp*B^x)/dx  
#endif
	!
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=2
    !------ 
    BQy(1+j) = BQy(1+j) + (Q(5)+Q(1))*Qy(10)    ! Q(10) or lapse
    BQy(5) = BQy(5) + S_contr(j)*Qy(10)         !  Q(10) or lapse
    !
#ifdef COVCLEAN     
	BQy(9) = BQy(9) - Q(5+j)*Qy(10)             ! D.C.
    BQy(9) = BQy(9) + Q(9)*Qy(10+j)             ! D.C.   phi * d(shift^y)/dy
    BQy(6:8) = BQy(6:8) - Q(11:13)*Qy(5+j)                  ! D.C.   shift * d(gp*B^y)/dy
    BQy(6:8) = BQy(6:8) + lapse*g_contr(1:3,j)*Qy(9)        ! D.C.  
#else
    BQy(6:8) = BQy(6:8) - Q(11:13)*Qy(5+j)                  ! D.C.   shift * d(gp*B^x)/dx 
    BQy(6:8) = BQy(6:8) + lapse*gp*g_contr(1:3,j)*Qy(9)     ! D.C.  
    BQy(9) = BQy(9) + gm*lapse*DivCleaning_a**2*Qy(5+j)             ! D.C.   shift * d(gp*B^x)/dx 
    BQy(9) = BQy(9) - Q(10+j)*Qy(9)             ! D.C.   shift * d(gp*B^x)/dx  
#endif
	!
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=3
    !------
    CQz(1+j) = CQz(1+j) + (Q(5)+Q(1))*Qz(10)    ! Q(10) or lapse
    CQz(5) = CQz(5) + S_contr(j)*Qz(10)         !  Q(10) or lapse
    !
#ifdef COVCLEAN     
	CQz(9) = CQz(9) - Q(5+j)*Qz(10)             ! D.C.
    CQz(9) = CQz(9) + Q(9)*Qz(10+j)             ! D.C.   phi *d(shift^z)/dz
    CQz(6:8) = CQz(6:8) - Q(11:13)*Qz(5+j)                  ! D.C.   shift * d(gp*B^z)/dz
    CQz(6:8) = CQz(6:8) + lapse*g_contr(1:3,j)*Qz(9)        ! D.C.  
#else
    CQz(6:8) = CQz(6:8) - Q(11:13)*Qz(5+j)                  ! D.C.   shift * d(gp*B^x)/dx 
    CQz(6:8) = CQz(6:8) + lapse*gp*g_contr(1:3,j)*Qz(9)     ! D.C.  
    CQz(9) = CQz(9) + gm*lapse*DivCleaning_a**2*Qz(5+j)             ! D.C.   shift * d(gp*B^x)/dx 
    CQz(9) = CQz(9) - Q(10+j)*Qz(9)             ! D.C.   shift * d(gp*B^x)/dx  
#endif
	!
	DO i=1,nVar
		BgradQ(i) = AQx(i) + BQy(i) + CQz(i)  
	ENDDO
	!
		!IF( ANY( ISNAN(BgradQ(:)) )) THEN
		!	PRINT *,'gradQ, PDENCPGRMHD NAN :' 
		!	WRITE(*,*) i
		!	WRITE(*,*) Q(1:5)
		!	WRITE(*,*) Q(6:10)
		!	WRITE(*,*) Q(11:15)
		!	WRITE(*,*) Q(16:19)
		!	PRINT *,'--------gradQ-------' 
		!	WRITE(*,*) gradQ(1:5,1)
		!	WRITE(*,*) gradQ(6:10,1)
		!	WRITE(*,*) gradQ(11:15,1)
		!	WRITE(*,*) gradQ(16:19,1)
		!	PRINT *,'---------------' 
		!	WRITE(*,*) gradQ(1:5,2)
		!	WRITE(*,*) gradQ(6:10,2)
		!	WRITE(*,*) gradQ(11:15,2)
		!	WRITE(*,*) gradQ(16:19,2)
		!	PRINT *,'---------------' 
		!	WRITE(*,*) gradQ(1:5,3)
		!	WRITE(*,*) gradQ(6:10,3)
		!	WRITE(*,*) gradQ(11:15,3)
		!	WRITE(*,*) gradQ(16:19,3)
		!	PRINT *,'--------BgradQ-------' 
		!	WRITE(*,*) BgradQ(1:5)
		!	WRITE(*,*) BgradQ(6:10)
		!	WRITE(*,*) BgradQ(11:15)
		!	WRITE(*,*) BgradQ(16:19)
		!	PRINT *,'--------AQx-------' 
		!	WRITE(*,*) AQx(1:5)
		!	WRITE(*,*) AQx(6:10)
		!	WRITE(*,*) AQx(11:15)
		!	WRITE(*,*) AQx(16:19)
		!	PRINT *,'--------BQy-------' 
		!	WRITE(*,*) BQy(1:5)
		!	WRITE(*,*) BQy(6:10)
		!	WRITE(*,*) BQy(11:15)
		!	WRITE(*,*) BQy(16:19)
		!	PRINT *,'--------CQz-------' 
		!	WRITE(*,*) CQz(1:5)
		!	WRITE(*,*) CQz(6:10)
		!	WRITE(*,*) CQz(11:15)
		!	WRITE(*,*) CQz(16:19)
		!	PRINT *,'---------------' 
		!	STOP
		!ENDIF 
    CONTINUE
    !            
END SUBROUTINE PDENCP      
!

RECURSIVE SUBROUTINE PDEEigenvalues(L,Q,n)
  USE Parameters, ONLY :  nDim, gamma,d, DivCleaning_a
  USE iso_c_binding
  IMPLICIT NONE
#ifdef VECTOR    
#ifdef AVX512 
  INTEGER, PARAMETER :: nVar = 24                           ! The number of variables of the PDE system 
#else   
  INTEGER, PARAMETER :: nVar = 20                           ! The number of variables of the PDE system 
#endif 
#else
  INTEGER, PARAMETER :: nVar = 19                           ! The number of variables of the PDE system 
#endif
  REAL :: L(nVar), n(nDim), Q(nVar), V(nVar)
  INTENT(IN)  :: Q,n
  INTENT(OUT) :: L 
  ! Local Variables 
    INTEGER, PARAMETER :: nParam=0 
    ! Local variables 
    REAL ::  Vp(nVar), gradQ(nVar,d)
    REAL :: p, irho
    REAL :: Qx(nVar), Qy(nVar), Qz(nVar) , Qp(nVar), Qm(nVar) 
    REAL :: e,g_cov(3,3),g_contr(3,3),nv(3)
    REAL :: AQx(nVar), BQy(nVar), CQz(nVar), FFp(nVar,d), FFm(nVar,d) ,dfdQ(nVar,nVar),AA(nVar,nVar)
    REAL :: lapse, shift(3),shift_cov(3), gammaij(6), delta(3,3), B_cov(3), vxB_cov(3), vxb_contr(3), psi, S_contr(3), qb_contr(3), B_contr(3) 
    REAL :: v2,v_contr(3),uem,b2,e2,gp,gm,vc(nVar),lf,w,ww,gamma1,rho,v_cov(3), w_ij, wim,eps,u,xg(3)
	REAL :: b2_4,cs2,sft,VdotB,vn,den,gg,lf2m,ImLambda(nVar), rtemp(nVar) ,R(nVar,nVar), iR(nVar,nVar),fmm(nVar),gmm(nVar),hmm(nVar),fpp(nVar),gpp(nVar),hpp(nVar),gg2,sft2,vn2
    INTEGER :: i,j,k,m,iErr, ccount, itemp(nVar) 
    LOGICAL :: FLAG
  LOGICAL:: hasNaN
  REAL :: Pi
  Pi = ACOS(-1.0)
    !
    L = 0.0
!    L(1) = -1.0    
!    L(2) = +1.0    
!    RETURN 


    xg(1) = 0.
    xg(2) = 0.
    xg(3) = 0.
   !
   !****************************************************
    ! Eigenvalues as given by Olindo 
    CALL PDECons2Prim(V,Q,iErr)
    rho    = V(1)
	DO i=1,3
		v_cov = V(1+i)
		B_contr(1:3) = V(5+i)
		shift = V(10+i)
	ENDDO
    p      = V(5)
    !
    psi = V(9)
    lapse = V(10)
    !
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
    !  evaluate contr. cov. variables
    v_contr      = MATMUL(g_contr,v_cov)
    B_cov      = MATMUL(g_cov,B_contr)
    !  evaluate useful quantities
    v2     = v_contr(1)*v_cov(1) + v_contr(2)*v_cov(2) + v_contr(3)*v_cov(3)
    lf     = 1.0/sqrt(1.0 - v2)
    lf2m   = 1.0 - v2
    b2     = B_contr(1)*B_cov(1) + B_contr(2)*B_cov(2) + B_contr(3)*B_cov(3)
    VdotB     = v_contr(1)*B_cov(1) + v_contr(2)*B_cov(2) + v_contr(3)*B_cov(3)
    b2_4 = b2*lf2m + VdotB  ! this is b^2
    gamma1 = gamma/(gamma-1.0) 
    w      = rho + gamma1*p + b2_4 ! this is rho*h + b^2
    cs2    = (gamma*p + b2_4)/w
    !
    vn     = v_contr(1)*n(1) + v_contr(2)*n(2) + v_contr(3)*n(3)
    sft    = shift(1)*n(1) + shift(2)*n(2) + shift(3)*n(3) 
    gg     = g_contr(1,1)*ABS(n(1)) + g_contr(2,2)*ABS(n(2)) + g_contr(3,3)*ABS(n(3))
    den    = 1.0/(1.0 - v2*cs2)

    IF(SUM(n**2).EQ.0.) THEN  
        u = SQRT( v2) 
        WRITE(*,*)'Impossible error!'
        STOP
    ELSE
        u = vn
    ENDIF
    L(1)   = ( u*(1.0-cs2) - SQRT( cs2*lf2m*( (1.0-v2*cs2)*gg - u**2*(1.0-cs2) )) )*den
    L(2:4) = u
    L(5)   = ( u*(1.0-cs2) + SQRT( cs2*lf2m*( (1.0-v2*cs2)*gg - u**2*(1.0-cs2) )) )*den 
	DO i=1,5
		L(i)   = lapse*L(i) - sft
	ENDDO
    !
	DO i=6,nVar
		L(i)   = 0.
	ENDDO

  
  !if any abs(lambda) > 1 for any, set to +-1 
  DO i=1,5
    if( abs(L(i)) > 1) then
      L(1)= -1.0
      L(2)= +1.0
    endif
  enddo

  ! if atmo is small, set the eigenvalues to +-1 
  if(Q(1) < 1.e-9) then 
    L(1) = -1.0    
    L(2) = +1.0    
    do i=3,nVar
      L(i)=0.0
    enddo     
  endif

  ! check for NaNs
  hasNaN=.FALSE.
  do i=1,nVar
    if (isnan(L(i))) then
     hasNaN = .TRUE. 
    endif 
  enddo     

  ! if there are any NaNs reset to the default
  if (hasNaN) then
    L(1) = -1.0    
    L(2) = +1.0    
    do i=3,nVar
      L(i)=0.0
    enddo     
  endif
  ! check for NaNs
  hasNaN=.FALSE.
  do i=1,nVar
    if (isnan(L(i))) then
     hasNaN = .TRUE. 
     print *, "lambda(",i,") is a NaN = ", L(i)
    endif 
  enddo     

 

  return 
!     !
!     ! SAFE MODE: define also 'covariant' eigenvalues! (we may use the remaining free slots in L, L(6:8), L(10:19)
!     !
!     shift_cov = MATMUL(g_cov,shift)
!     !
!     vn2  = v_cov(1)*n(1) + v_cov(2)*n(2) + v_cov(3)*n(3)
!     sft2 = shift_cov(1)*n(1) + shift_cov(2)*n(2) + shift_cov(3)*n(3) 
!     gg2  = g_cov(1,1)*ABS(n(1)) + g_cov(2,2)*ABS(n(2)) + g_cov(3,3)*ABS(n(3))
!     den  = 1.0/(1.0 - v2*cs2)
!     !
!     u = vn2 
!     !
!     L(10)  = ( u*(1.0-cs2) - SQRT( cs2*lf2m*( (1.0-v2*cs2)*gg2 - u**2*(1.0-cs2) )) )*den
!     L(6:8) = u
!     L(11)  = ( u*(1.0-cs2) + SQRT( cs2*lf2m*( (1.0-v2*cs2)*gg2 - u**2*(1.0-cs2) )) )*den 
!     L(6:11)   = lapse*L(6:11) - sft2
!     ! 
!     L(9) =  DivCleaning_a   ! 1.  !ch
!     !
!     FLAG = .FALSE.
!     IF(MAXVAL(ABS(L)).GT.1.01) THEN
!         FLAG = .TRUE.
!         continue
!     ENDIF
!     !  
	continue
	!
END SUBROUTINE PDEEigenvalues

RECURSIVE SUBROUTINE PDESource(S,Q) 
  USE Parameters, ONLY:  nDim, d, gamma, DivCleaning_a
  USE iso_c_binding
  IMPLICIT NONE
    ! Argument list 
#ifdef VECTOR    
#ifdef AVX512 
  INTEGER, PARAMETER :: nVar = 24                           ! The number of variables of the PDE system 
#else   
  INTEGER, PARAMETER :: nVar = 20                           ! The number of variables of the PDE system 
#endif 
#else
  INTEGER, PARAMETER :: nVar = 19                           ! The number of variables of the PDE system 
#endif  
  ! --------------------------------------------
  ! Argument list declaration
  REAL :: S(nvar), Q(nvar)
  INTENT(IN)  :: Q 
  INTENT(OUT) :: S
  ! --------------------------------------------
  ! Local variable declaration  
  REAL :: p, irho, lam, mu 
  REAL :: k1, k2, fff, ggg, e, ds, c, xi, sk, sknl, alpha, fa, k0, beta0(3), eta 
  REAL :: g_cov(3,3), det, g_contr(3,3), Christoffel(3,3,3), dgup(3,3,3), faa, b0(3), dChristoffelSrc(3,3,3,3), RiemannSrc(3,3,3,3)
  REAL :: RicciSrc(3,3), gammaup(3), dtmetric(3,3), dtgup(3,3), dtChristoffelSrc(3,3,3), dtGammaUpSrc(3)
  REAL :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s25
    !
    S = 0.0     
    !
END SUBROUTINE PDESource


RECURSIVE SUBROUTINE PDEVarName(Name)  
  IMPLICIT NONE    
#ifdef VECTOR    
#ifdef AVX512 
  INTEGER, PARAMETER :: nVar = 24                           ! The number of variables of the PDE system 
#else   
  INTEGER, PARAMETER :: nVar = 20                           ! The number of variables of the PDE system 
#endif 
#else
  INTEGER, PARAMETER :: nVar = 19                           ! The number of variables of the PDE system 
#endif 
  CHARACTER(LEN=10):: Name(nVar)

  ! EQNTYPE4
#ifdef GRMHD
  Name(1) = 'rho' 
  Name(2) = 'u' 
  Name(3) = 'v'
  Name(4) = 'w'
  Name(5) = 'p' 
  Name(6) = 'Bx' 
  Name(7) = 'By' 
  Name(8) = 'Bz' 
  Name(9) = 'psi' 
  Name(10) = 'alpha' 
  Name(11) = 'beta^1' 
  Name(12) = 'beta^2' 
  Name(13) = 'beta^3' 
  Name(14) = 'gamma_11' 
  Name(15) = 'gamma_12' 
  Name(16) = 'gamma_13' 
  Name(17) = 'gamma_22' 
  Name(18) = 'gamma_23' 
  Name(19) = 'gamma_33' 
#ifdef VECTOR
#ifdef AVX512
  Name(20) = 'pad1' 
  Name(21) = 'pad2' 
  Name(22) = 'pad3' 
  Name(23) = 'pad4' 
  Name(24) = 'pad5' 
#else
  Name(20) = 'pad1' 
#endif 
#endif 
#endif 
    !

END SUBROUTINE PDEVarName
