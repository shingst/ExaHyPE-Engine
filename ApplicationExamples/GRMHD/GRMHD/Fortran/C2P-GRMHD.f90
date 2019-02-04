! C2P for GRMHD
! C2P for GRMHD

RECURSIVE SUBROUTINE PDEPrim2Cons(Q,V)
  USE Parameters, ONLY: gamma,  nDim 
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
  ! Argument list declaration
  REAL :: Q(nVar), V(nVar)
  INTENT(IN)  :: V
  INTENT(OUT) :: Q 
  ! Local variable declaration
  REAL :: rho,p,psi,lapse,gp,gm
  REAL :: shift(3),v_cov(3),v_contr(3),vxB_cov(3),vxB_contr(3),B_cov(3),B_contr(3)
  REAL :: g_cov(3,3),g_contr(3,3)
  REAL :: vb,v2,e2,b2,uem,LF,gamma1,w,ww
  INTEGER :: i


  rho     = V(1)
  v_cov  = V(2:4)
  p       = V(5)
  !
  psi = V(9)
  lapse = V(10)
  DO i=1,3
	  B_contr(i) = V(5+i)  !wrong 
	  shift(i) = V(10+i)          ! NB: we choose V() and Q() being the shift_controvariant!
  ENDDO
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
  v_contr      = MATMUL(g_contr,v_cov)
  B_cov = MATMUL(g_cov,B_contr(1:3))
  !
  ! COVARIANT =>
  vxB_cov(1) = v_contr(2)*B_contr(3) - v_contr(3)*B_contr(2)
  vxB_cov(2) = v_contr(3)*B_contr(1) - v_contr(1)*B_contr(3)
  vxB_cov(3) = v_contr(1)*B_contr(2) - v_contr(2)*B_contr(1)
  vxB_cov(1) = gp*vxB_cov(1)
  vxB_cov(2) = gp*vxB_cov(2)
  vxB_cov(3) = gp*vxB_cov(3)
  !
  ! CONTRAVARIANT =>
  vxB_contr(1) = v_cov(2)*B_cov(3) - v_cov(3)*B_cov(2)
  vxB_contr(2) = v_cov(3)*B_cov(1) - v_cov(1)*B_cov(3)
  vxB_contr(3) = v_cov(1)*B_cov(2) - v_cov(2)*B_cov(1)
  vxB_contr(1) = gm*vxB_contr(1)
  vxB_contr(2) = gm*vxB_contr(2)
  vxB_contr(3) = gm*vxB_contr(3)
  !
  vb      = v_cov(1)*B_contr(1) + v_cov(2)*B_contr(2) + v_cov(3)*B_contr(3)      
  v2     = v_contr(1)*v_cov(1) + v_contr(2)*v_cov(2) + v_contr(3)*v_cov(3)
  e2     = vxB_contr(1)*vxB_cov(1) + vxB_contr(2)*vxB_cov(2) + vxB_contr(3)*vxB_cov(3)
  b2     = B_contr(1)*B_cov(1) + B_contr(2)*B_cov(2) + B_contr(3)*B_cov(3)
  !
  uem    = 0.5*(b2 + e2) 
  !
  !IF (v2 > 1.0) THEN
  !     ! IF(MAXVAL(ABS(xGP)).LE.Exc_radius) THEN  ! RETURN
  !   WRITE(*,*)'Superluminal velocity in PDEPrim2Cons!!'
  !   STOP
  !ENDIF
  lf     = 1.0 / sqrt(1.0 - v2)
  gamma1 = gamma/(gamma-1.0)
  w      = rho + gamma1*p
  ww     = w*lf**2
  !
  Q(1)    = rho*lf
  Q(2)  = ww*v_cov(1) + b2*v_cov(1) - vb*B_cov(1)
  Q(3)  = ww*v_cov(2) + b2*v_cov(2) - vb*B_cov(2)
  Q(4)  = ww*v_cov(3) + b2*v_cov(3) - vb*B_cov(3) 
  Q(5)    = ww - p + uem - Q(1)     !!!!! we subtract PDE(Q(1))!!!!
  Q(6) = V(6)
  Q(7) = V(7)
  Q(8) = V(8)
  !
  DO i=1,8
		Q(i) = gp*Q(i)  
  ENDDO
  DO i=9,nVar
		Q(i) = V(i)
  ENDDO
  ! 
  continue
    !
END SUBROUTINE PDEPrim2Cons     


RECURSIVE SUBROUTINE PDECons2Prim(V,Q,iErr)
  USE Parameters, ONLY: gamma, nDim,p_floor,rho_floor,NSTOV_rho_atmo,NSTOV_p_atmo
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
    REAL, INTENT(IN)     :: Q(nVar)     ! vector of conserved quantities 
    REAL, INTENT(OUT)    :: V(nVar)     ! primitive variables 
    INTEGER, INTENT(OUT) :: iErr        ! error flag 
    ! Local variables 
	REAL :: rho,p,psi,lapse,gp,gm
	REAL :: shift(3),v_cov(3),v_contr(3),vxB_cov(3),vxB_contr(3),B_cov(3),B_contr(3)
	REAL :: sm(3),sm_cov(3),Qtmp(nVar),Qloc(nVar)
	REAL :: g_cov(3,3),g_contr(3,3),gammaij(6)
	REAL :: vb,v2,e2,b2,s2,sb,sb2,e,x1,x2,eps,uem,LF,gamma1,gam,w,ww
	REAL :: vx,vy,vz,den,dd
	REAL, PARAMETER :: tol = 1.0e-18
	INTEGER :: i
    LOGICAL              :: FAILED
    !    
    iErr = 0     
    !
	!V = Q
	!RETURN
	
  psi = Q(9)
  lapse = Q(10)
	DO i=1,3
		shift(i) = Q(10+i)
	ENDDO
	DO i=1,6
		gammaij(i) = Q(13+i)
	ENDDO
  !
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
	IF(SUM(g_cov**2).LT.1e-9) THEN
        print *, 'FATAL ERROR: g_cov = 0'
		WRITE(*,*) g_cov(1,1:3)
		WRITE(*,*) g_cov(2,1:3)
		WRITE(*,*) g_cov(3,1:3)
	ENDIF
		
		
  CALL MatrixInverse3x3(g_cov,g_contr,gp)
  !
  gp = SQRT(gp)
  gm = 1./gp
  ! 
	DO i=1,8
		Qloc(i)  = gm*Q(i)
	ENDDO
  !
	dd = Qloc(1)
	DO i=1,3
		B_contr(i) = Qloc(5+i) ! !wrong: magnetic field is contravariant
		sm_cov(i)  = Qloc(1+i)
	ENDDO
  !
  gamma1 = gamma/(gamma - 1.0)
  gam    = 1.0/gamma1
  ! Solve for p
  FAILED  = .FALSE.
  !
  sm   = MATMUL (g_contr, sm_cov)
	!
	B_cov   = MATMUL (g_cov, B_contr)
  s2   = sm_cov(1)*sm(1) + sm_cov(2)*sm(2) + sm_cov(3)*sm(3)
	b2   = B_contr(1)*B_cov(1) + B_contr(2)*B_cov(2) + B_contr(3)*B_cov(3)
	sb   = sm_cov(1)*B_contr(1) + sm_cov(2)*B_contr(2) + sm_cov(3)*B_contr(3)
  sb2  = sb**2
	!
  ! First option [Del Zanna et al. (2007) A&A, 473, 11-30 (method 3)]
	e    = Qloc(5) + dd  ! Q(5) = gamma^1/2 ( U - dd )
  x1   = 0.      ! 
	eps = 1.0e-10
	x2   = 1.0-eps !tol ! 
  w=0 
  !
	CALL RTSAFE_C2P_RMHD1(v2,x1,x2,tol,gam,dd,e,s2,b2,sb2,w,FAILED) 
  !
  IF (FAILED) THEN
		 ! 
     iErr = -1
     p    = p_floor
     rho  = rho_floor
     vx   = 0.0
     vy   = 0.0
     vz   = 0.0
		 v_cov(1) = 0. 
		 v_cov(2) = 0. 
		 v_cov(3) = 0. 
		 !
  ELSE
     den  = 1.0/(w+b2)
     vb   = sb/w
     !
     rho  = dd*sqrt(1.-v2)
	 v_cov(1) = (sm_cov(1) + vb*B_cov(1))*den
	 v_cov(2) = (sm_cov(2) + vb*B_cov(2))*den
	 v_cov(3) = (sm_cov(3) + vb*B_cov(3))*den
     !p = max(1.e-15, gam*(w*(1.-v2)-rho))
     p = gam*(w*(1.-v2)-rho)
  ENDIF
  !
    IF(rho<NSTOV_rho_atmo*1.01) THEN   
        v_cov(:)= 0.0
        !rho = 1e-12 
        !p   = 1e-12 
    ENDIF
    IF(rho<NSTOV_rho_atmo) THEN 
        rho = NSTOV_rho_atmo
    ENDIF 
    !
    IF(p<NSTOV_p_atmo*1.01) THEN 
        v_cov(:)= 0.0 
    ENDIF   
    IF(p<NSTOV_p_atmo) THEN 
        p=NSTOV_p_atmo
    ENDIF   
  !
	V(1) = rho
	V(2) = v_cov(1)
	V(3) = v_cov(2)
	V(4) = v_cov(3)
	V(5) = p
	V(6) = B_contr(1)
	V(7) = B_contr(2)
	V(8) = B_contr(3)
	V(9) = psi
	V(10) = lapse
	V(11) = shift(1)
	V(12) = shift(2)
	V(13) = shift(3)
	DO i=1,6
		V(13+i) = gammaij(i)
    ENDDO
    V(19:nVar)=Q(19:nVar)
	!
	!CALL PDEPrim2ConsGRMHD(Qtmp,V)
	!
	!IF(MAXVAL(ABS(Qtmp-Q)).GT.1e-8) THEN
	!  iErr = +10   
	!		print *, 'FATAL ERROR: C2P o P2C'
	!		WRITE(*,*) MAXVAL(ABS(Qtmp-Q))  
	!		STOP
	!  continue
	!  !
	!ENDIF
	!  
  !  
END SUBROUTINE PDECons2Prim  

    
    
    
    
