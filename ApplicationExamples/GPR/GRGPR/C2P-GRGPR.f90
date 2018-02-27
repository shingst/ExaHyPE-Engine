! The Con2Prim and Prim2Con routines for MHD.
! Should be merged with SRHD's.

SUBROUTINE PDEPrim2Cons(Q,V)
  USE Parameters, ONLY: nVar, nDim, cs, p_floor, rho_floor, gamma,tol
  USE GRGPRmod
  IMPLICIT NONE
  ! Argument list declaration
  REAL :: Q(nVar), V(nVar)
  INTENT(IN)  :: V
  INTENT(OUT) :: Q 
  ! Local variable declaration
		INTEGER :: iter, Ierr,ii,jj,ll,mm,nn
		REAL	:: kAB(3,3),shift(3),gammaij(6),g_cov(3,3),g_contr(3,3),v2
		REAL    :: lapse, gp, gm,Qloc(nVar), gamma1, gam, A(3,3),detA, vf_cov(3),vf(3),Vtr(3),LF,shift_cov(3)
		REAL 	:: rho,pi_work(3),ep,pi_scalar,A_st(1:3,0:3),g_cov_st(0:3,0:3),Gmunu(0:3,0:3),hmunu(0:3,0:3),devG_st(0:3,0:3),g_contr_st(0:3,0:3)
		REAL	:: pi_mix(0:3,0:3), pi_st(0:3,0:3)
		LOGICAL	:: FAILED
		REAL    :: sm_cov(3), sm(3),Bv(3),s2,b2,sb,sb2,eps,e,x1,x2,w,p,vx,vy,vz,bx,by,bz,den,vb,vtr_old(3),res, ww
  !
  ! For the moment, we assume a flat metric in matter space 
  kAB(1,1) = V(25) 
  kAB(1,2) = V(26) 
  kAB(1,3) = V(27) 
  kAB(2,1) = V(26) 
  kAB(2,2) = V(28) 
  kAB(2,3) = V(29) 
  kAB(3,1) = V(27) 
  kAB(3,2) = V(29) 
  kAB(3,3) = V(30) 
  !
  rho     = V(1)
  vf_cov  = V(2:4)
  p       = V(5)
  !
  lapse = V(15)
  shift = V(16:18)          ! NB: we choose V() and Q() being the shift_controvariant!
  !
  g_cov(1,1) = V(19)
  g_cov(1,2) = V(20)
  g_cov(1,3) = V(21)
  g_cov(2,2) = V(22)
  g_cov(2,3) = V(23)
  g_cov(3,3) = V(24)
  g_cov(2,1) = g_cov(1,2) 
  g_cov(3,1) = g_cov(1,3) 
  g_cov(3,2) = g_cov(2,3) 
  !
  ! The deformation gradient tensor 
  A(1,:) = (/ V( 6), V( 7),  V( 8) /) 
  A(2,:) = (/ V( 9), V(10),  V(11) /)
  A(3,:) = (/ V(12), V(13),  V(14) /)         
  detA   = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(2,1)*A(1,2)*A(3,3)+A(2,1)*A(1,3)*A(3,2)+A(3,1)*A(1,2)*A(2,3)-A(3,1)*A(1,3)*A(2,2) 
  !
  CALL MatrixInverse3x3(g_cov,g_contr,gp)
  gp = SQRT(gp)
  gm = 1./gp
  !
  vf       = MATMUL(g_contr,vf_cov) 
  v2     = vf(1)*vf_cov(1) + vf(2)*vf_cov(2) + vf(3)*vf_cov(3)
  !
  IF (v2 > 1.0) THEN
     WRITE(*,*)'Superluminal velocity in PDEPrim2Cons!!'
     STOP
  ENDIF
  lf     = 1.0 / sqrt(1.0 - v2)
  gamma1 = gamma/(gamma-1.0)
  !
  ! transport velocity
  Vtr(1:3) = lapse*vf(1:3)-shift(1:3)   ! this is contravariant !!!!
  !
  ! compute the missing time components of the 3x4 deformation gradient 
  A_st(1:3,1:3) = A(1:3,1:3) 
  DO ii = 1, 3
        A_st(ii,0) = -DOT_PRODUCT(Vtr(1:3),A(ii,:)) 
  ENDDO
  ! compute the covariant 4D metric using the 3+1 split 
  g_cov_st(1:3,1:3) = g_cov
  shift_cov = MATMUL(g_cov,shift) 
  g_cov_st(0,0)   = -lapse**2 + DOT_PRODUCT(shift,shift_cov) 
  g_cov_st(0,1:3) = shift_cov 
  g_cov_st(1:3,0) = shift_cov 
  ! compute the inverse of the 4D covariant metric 
  DO ii = 1, 3
   DO jj = 1, 3
     g_contr_st(ii,jj) = g_contr(ii,jj) - 1./lapse**2*shift(ii)*shift(jj) 
   ENDDO
  ENDDO
  g_contr_st(0,0)   = -1./lapse**2 
  g_contr_st(0,1:3) =  1./lapse**2*shift(1:3) 
  g_contr_st(1:3,0) =  1./lapse**2*shift(1:3) 
  hmunu = 0.0 
  DO mm = 0, 3 
   hmunu(mm,mm) = 1.0
  ENDDO  
  Gmunu = 0.0 
  DO mm = 0, 3 
   DO nn = 0, 3 
    DO ll = 0, 3 
     DO jj = 1, 3 
      DO ii = 1, 3 
         Gmunu(mm,nn) = Gmunu(mm,nn) + g_contr_st(mm,ll)*kAB(ii,jj)*A_st(ii,ll)*A_st(jj,nn)  
      ENDDO        
     ENDDO 
    ENDDO
   ENDDO 
  ENDDO    
  !
  devG_st   = Gmunu - (Gmunu(0,0)+Gmunu(1,1)+Gmunu(2,2)+Gmunu(3,3))/3.0 * hmunu  
  pi_mix    = rho*cs**2*MATMUL(Gmunu,devG_st)    
  pi_st     = MATMUL( g_cov_st, pi_mix )  
  !
  CALL GRTepsilonA(g_contr_st,kAB,A_st,ep)   
  !
  pi_scalar = SUM( g_contr(1:3,1:3)*pi_st(1:3,1:3) ) 
  ! Compute the work of the shear stress tensor 
  pi_work   = MATMUL( pi_st(1:3,1:3), vf )  
  !
  w      = rho + rho*ep + gamma1*p 
  ww     = w*lf**2    ! = rho*h*W**2
  !
  Q(1)    = rho*lf 
  Q(2:4)  = ww*vf_cov(1:3) + pi_work    ! covariant!!!!
  Q(5)    = ww - p + pi_scalar - Q(1)     !!!!! we subtract PDE(Q(1))!!!!
  ! 
  Q(6:nvar) = V(6:nvar)
  !
  Q(1:5)  = gp*Q(1:5)
  !
END SUBROUTINE PDEPrim2Cons

SUBROUTINE PDECons2Prim(V,Q)
  USE Parameters, ONLY: nVar, nDim, cs, p_floor, rho_floor, gamma,tol
  USE GRGPRmod
  IMPLICIT NONE
  ! Argument list declaration
  REAL :: Q(nVar), V(nVar)
  INTENT(OUT)  :: V
  INTENT(IN) :: Q 
  ! Local variable declaration
		INTEGER :: iter, Ierr,ii,jj,ll,mm,nn
		REAL	:: kAB(3,3),shift(3),gammaij(6),g_cov(3,3),g_contr(3,3),v2
		REAL    :: lapse, gp, gm,Qloc(nVar), gamma1, gam, A(3,3),detA, vf_cov(3),vf(3),Vtr(3),LF,shift_cov(3)
		REAL 	:: rho,pi_work(3),ep,pi_scalar,A_st(1:3,0:3),g_cov_st(0:3,0:3),Gmunu(0:3,0:3),hmunu(0:3,0:3),devG_st(0:3,0:3),g_contr_st(0:3,0:3)
		REAL	:: pi_mix(0:3,0:3), pi_st(0:3,0:3)
		LOGICAL	:: FAILED
		REAL    :: sm_cov(3), sm(3),Bv(3),s2,b2,sb,sb2,eps,e,x1,x2,w,p,vx,vy,vz,bx,by,bz,den,vb,vtr_old(3),res
  !
  !
  iter = 0   
  !
  ! For the moment, we assume a flat metric in matter space 
  kAB(1,1) = Q(25) 
  kAB(1,2) = Q(26) 
  kAB(1,3) = Q(27) 
  kAB(2,1) = Q(26) 
  kAB(2,2) = Q(28) 
  kAB(2,3) = Q(29) 
  kAB(3,1) = Q(27) 
  kAB(3,2) = Q(29) 
  kAB(3,3) = Q(30) 
  !
  lapse = Q(15)
  shift = Q(16:18)
  !
  gammaij = Q(19:24) 
  g_cov(1,1) = Q(19)
  g_cov(1,2) = Q(20)
  g_cov(1,3) = Q(21)
  g_cov(2,2) = Q(22)
  g_cov(2,3) = Q(23)
  g_cov(3,3) = Q(24)
  g_cov(2,1) = g_cov(1,2) 
  g_cov(3,1) = g_cov(1,3) 
  g_cov(3,2) = g_cov(2,3) 
  !
  CALL MatrixInverse3x3(g_cov,g_contr,gp)
  !     
  gp = SQRT(gp)
  gm = 1./gp
  ! 
  Qloc(1:5)  = gm*Q(1:5)
  !
  gamma1 = gamma/(gamma - 1.0)
  gam    = 1.0/gamma1
  !
  ! The deformation gradient tensor 
  A(1,:) = (/ Q( 6), Q( 7),  Q( 8) /) 
  A(2,:) = (/ Q( 9), Q(10),  Q(11) /)
  A(3,:) = (/ Q(12), Q(13),  Q(14) /)         
  detA   = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(2,1)*A(1,2)*A(3,3)+A(2,1)*A(1,3)*A(3,2)+A(3,1)*A(1,2)*A(2,3)-A(3,1)*A(1,3)*A(2,2) 
  !
  ! initial guess for the transport velocity
  !Vtr(1:3) = lapse*vf(1:3)-shift(1:3)   ! this is contravariant !!!!
  !
  vf_cov   = Qloc(2:4)/(Qloc(5)+Qloc(1)) 
  vf       = MATMUL(g_contr, vf_cov) 
  v2       = vf(1)*vf_cov(1) + vf(2)*vf_cov(2) + vf(3)*vf_cov(3) 
  IF( v2 > 1.0 ) THEN
      vf_cov   = 0.5* vf_cov / SQRT( v2 )  
      vf       = MATMUL(g_contr, vf_cov) 
      v2       = vf(1)*vf_cov(1) + vf(2)*vf_cov(2) + vf(3)*vf_cov(3) 
      IF( v2 > 1.0 ) THEN
          PRINT *, ' impossible error! v2 still > 1 even after speed reduction. ' 
          STOP 
      ENDIF      
  ENDIF  
  Vtr(1:3)      = lapse*vf(1:3)-shift(1:3)   ! this is contravariant !!!! 
  LF            = 1.0 / sqrt(1.0 - v2) 
  rho           = Qloc(1)/LF 
  pi_work       = 0.0 
  ep            = 0.0 
  pi_scalar     = 0.0 
  !
  !GOTO 22345 
  ! 
  !do while(iter<10)
12345 CONTINUE
  ! compute the missing time components of the 3x4 deformation gradient 
  A_st(1:3,1:3) = A(1:3,1:3) 
  DO ii = 1, 3
        A_st(ii,0) = -DOT_PRODUCT(Vtr(1:3),A(ii,:)) 
  ENDDO
  ! compute the covariant 4D metric using the 3+1 split 
  g_cov_st(1:3,1:3) = g_cov
  shift_cov = MATMUL(g_cov,shift) 
  g_cov_st(0,0)   = -lapse**2 + DOT_PRODUCT(shift,shift_cov) 
  g_cov_st(0,1:3) = shift_cov 
  g_cov_st(1:3,0) = shift_cov 
  ! compute the inverse of the 4D covariant metric 
  DO ii = 1, 3
   DO jj = 1, 3
     g_contr_st(ii,jj) = g_contr(ii,jj) - 1./lapse**2*shift(ii)*shift(jj) 
   ENDDO
  ENDDO
  g_contr_st(0,0)   = -1./lapse**2 
  g_contr_st(0,1:3) =  1./lapse**2*shift(1:3) 
  g_contr_st(1:3,0) =  1./lapse**2*shift(1:3) 

  hmunu = 0.0 
  DO mm = 0, 3 
   hmunu(mm,mm) = 1.0
  ENDDO     
  Gmunu = 0.0 
  DO mm = 0, 3 
   DO nn = 0, 3 
    DO ll = 0, 3 
     DO jj = 1, 3 
      DO ii = 1, 3 
         Gmunu(mm,nn) = Gmunu(mm,nn) + g_contr_st(mm,ll)*kAB(ii,jj)*A_st(ii,ll)*A_st(jj,nn)  
      ENDDO        
     ENDDO 
    ENDDO
   ENDDO 
  ENDDO    
  devG_st   = Gmunu - (Gmunu(0,0)+Gmunu(1,1)+Gmunu(2,2)+Gmunu(3,3))/3.0 * hmunu 
  pi_mix    = rho*cs**2*MATMUL(Gmunu,devG_st)        
  pi_st  = MATMUL( g_cov_st, pi_mix )  
  !
  CALL GRTepsilonA(g_contr_st,kAB,A_st,ep) 
  ! 
  pi_scalar = SUM( g_contr(1:3,1:3)*pi_st(1:3,1:3) ) 
  ! Compute the work of the shear stress tensor 
  pi_work   = MATMUL( pi_st(1:3,1:3), vf )  
  FAILED  = .FALSE.
  sm_cov  = Qloc(2:4) - pi_work - rho*ep*LF**2*vf_cov 
  sm   = MATMUL (g_contr, sm_cov)
  BV   = 0.0 
  s2   = sm_cov(1)*sm(1) + sm_cov(2)*sm(2) + sm_cov(3)*sm(3)
  b2   = 0.0 
  sb   = 0.0 
  sb2  = 0.0 
  eps  = 1.e-10  ! 1e-8 
  e    = Qloc(5) + Qloc(1)  - pi_scalar  - rho*LF**2*ep    ! Q(5) = gamma^1/2 ( U - D )
  x1   = 0.      ! 
  x2   = 1.0-eps ! 
  w=0
  IF(ANY(ISNAN((/x1,x2,tol,gam,Qloc(1),e,s2,b2,sb2,w/)))) THEN
      continue
  ENDIF
  !
  CALL RTSAFE_C2P_RMHD1(v2,x1,x2,tol,gam,Qloc(1),e,s2,b2,sb2,w,FAILED) 
  !
  IF (FAILED) THEN
     IF(iter>1) THEN 
         CALL GRGPRNewton(V,Q,iErr,iter) 
         RETURN 
     ENDIF 
     iErr = -1
     p    = p_floor
     rho  = rho_floor
     vx   = 0.0
     vy   = 0.0
     vz   = 0.0
     bx   = 0.0 
     by   = 0.0 
     bz   = 0.0 
  ELSE
     den  = 1.0/(w+b2)
     vb   = sb/w
     !
     rho  = Qloc(1)*SQRT(1.0-v2)
     vf_cov(1) = (sm_cov(1) + vb*BV(1))*den
     vf_cov(2) = (sm_cov(2) + vb*BV(2))*den
     vf_cov(3) = (sm_cov(3) + vb*BV(3))*den
     p = max(1.e-15, gam*(w*(1.-v2)-rho))
     !
     vtr_old = vtr 
     vf = MATMUL( g_contr, vf_cov )   
     v2     = vf(1)*vf_cov(1) + vf(2)*vf_cov(2) + vf(3)*vf_cov(3)
     LF     = 1.0 / sqrt(1.0 - v2)
     Vtr(1:3) = lapse*vf(1:3)-shift(1:3)
     res = SQRT( SUM( (vtr-vtr_old)**2 ) ) 
     IF(res<1e-9) THEN 
        CONTINUE
     ELSE
        iter = iter + 1
        IF(iter<10) THEN 
            GOTO 12345    
        ELSE
             CALL GRGPRNewton(V,Q,iErr,iter) 
             RETURN 
        ENDIF
     ENDIF 
     !
  ENDIF
  !end do
  !
  V(1:5) = (/ rho, vf_cov(1:3), p /) 
  V(6:nVar) = Q(6:nVar) 
  !
END SUBROUTINE PDECons2Prim

SUBROUTINE GRGPRNewton(V,Q,iErr,iNewton)
    USE Parameters, ONLY : nVar, gamma 
	USE GRGPRmod
    IMPLICIT NONE 
    INTEGER :: ii,jj,kk,ll,mm,nn,iNewton,MaxNewton = 100 
    INTEGER :: k, iErr, inner  
    REAL    :: Q(nVar), V(nVar)  
    REAL    :: g_cov(3,3), g_contr(3,3), Aex(3,3), traceA, phi, DD(3,3,3), traceDk, tol = 1e-11         
    REAL    :: rho, kAB(3,3), detkAB, vf_cov(3), vf(3), lapse, shift(3), shift_cov(3), A(3,3)
    REAL    :: A_st(3,0:3), Gmunu(0:3,0:3), gp, gm, lf, vtr(3), g_cov_st(0:3,0:3), g_contr_st(0:3,0:3), ep
    REAL    :: I1, I2, devG_st(0:3,0:3), rhoA, detA, detA2, v2, Qloc(nVar), gam, gamma1, p     
    REAL    :: f(5), df(5,5), idf(5,5), dv(5), res, newres, delta   
    ! 
  ! For the moment, we assume a flat metric in matter space 
  kAB(1,1) = Q(25) 
  kAB(1,2) = Q(26) 
  kAB(1,3) = Q(27) 
  kAB(2,1) = Q(26) 
  kAB(2,2) = Q(28) 
  kAB(2,3) = Q(29) 
  kAB(3,1) = Q(27) 
  kAB(3,2) = Q(29) 
  kAB(3,3) = Q(30) 
  !
  lapse = Q(15)
  shift = Q(16:18)
  !
  g_cov(1,1) = Q(19)
  g_cov(1,2) = Q(20)
  g_cov(1,3) = Q(21)
  g_cov(2,2) = Q(22)
  g_cov(2,3) = Q(23)
  g_cov(3,3) = Q(24)
  g_cov(2,1) = g_cov(1,2) 
  g_cov(3,1) = g_cov(1,3) 
  g_cov(3,2) = g_cov(2,3) 
  !
  CALL MatrixInverse3x3(g_cov,g_contr,gp)
  !     
  gp = SQRT(gp)
  gm = 1./gp
  ! 
  Qloc(1:5)  = gm*Q(1:5)
  !
  gamma1 = gamma/(gamma - 1.0)
  gam    = 1.0/gamma1
  !
  ! The deformation gradient tensor 
  A(1,:) = (/ Q( 6), Q( 7),  Q( 8) /) 
  A(2,:) = (/ Q( 9), Q(10),  Q(11) /)
  A(3,:) = (/ Q(12), Q(13),  Q(14) /)         
  detA   = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(2,1)*A(1,2)*A(3,3)+A(2,1)*A(1,3)*A(3,2)+A(3,1)*A(1,2)*A(2,3)-A(3,1)*A(1,3)*A(2,2) 
  !
  ! initial guess for the transport velocity
  !Vtr(1:3) = lapse*vf(1:3)-shift(1:3)   ! this is contravariant !!!!
  !
  vf_cov   = Qloc(2:4)/(Qloc(5)+Qloc(1)) 
  vf       = MATMUL(g_contr, vf_cov) 
  v2       = vf(1)*vf_cov(1) + vf(2)*vf_cov(2) + vf(3)*vf_cov(3) 
  IF( v2 > 1.0 ) THEN
      vf_cov   = 0.5* vf_cov / SQRT( v2 )  
      vf       = MATMUL(g_contr, vf_cov) 
      v2       = vf(1)*vf_cov(1) + vf(2)*vf_cov(2) + vf(3)*vf_cov(3) 
      IF( v2 > 1.0 ) THEN
          PRINT *, ' impossible error! v2 still > 1 even after speed reduction. ' 
          STOP 
      ENDIF      
  ENDIF  
  LF   = 1.0 / sqrt(1.0 - v2) 
  rho  = Qloc(1)/lf 
  p    = ( Qloc(5)+Qloc(1) - LF**2*rho )/( LF**2 * gamma1 - 1.0 ) 
  !
  V(1:5) = (/ rho, vf_cov(1:3), p /) 
  V(6:30) = Q(6:30) 
  !
  iErr = 0 
  DO iNewton = 1, MaxNewton
      CALL func1(f,V(1:5),Q(:)) 
      res = SQRT(SUM(f**2)) 
      IF(res<tol) THEN 
          RETURN 
      ENDIF
      CALL dfunc1(df,V(1:5),Q(:)) 
      CALL MatrixInverse(5,df,idf) 
      dv = MATMUL( idf, f ) 
      delta = 1.0 
      DO inner = 1, 50
          IF( (SUM( (V(2:4)-delta*dv(2:4))**2 )) > 0.999 .OR. ( V(1)-delta*dv(1) < 1e-6 ) .OR. ( V(5)-delta*dv(5) < 1e-6 )   ) THEN
              delta = delta / 2.0 
          ELSE
              CALL func1(f,V(1:5)-delta*dv,Q(:))  
              newres = SQRT(SUM(f**2)) 
              IF(newres >= res ) THEN
                  delta = delta / 2.0
              ELSE
                  EXIT 
              ENDIF 
          ENDIF          
      ENDDO      
      V(1:5) = V(1:5) - delta*dv 
  ENDDO  
  ! Newton did not converge 
  PRINT *, ' Newton did not converge ', res  
  iErr = -10 
  !
  CONTINUE
  !
END SUBROUTINE GRGPRNewton  

SUBROUTINE GRTepsilonA(g_contr_st,kAB,A_st,myepsilon)
    USE Parameters, ONLY : cs
    IMPLICIT NONE
    ! Argument list 
    REAL :: g_contr_st(0:3,0:3), kAB(3,3), A_st(3,0:3), myepsilon 
    ! Local variables 
  INTEGER :: i, iErr, ii, jj, mm, kk, ll, nn, m
  REAL    :: Gmunu(0:3,0:3), temp_st(0:3,0:3)   
  REAL    :: I1, I2 

  Gmunu = 0.0 
  DO mm = 0, 3 
   DO nn = 0, 3 
    DO ll = 0, 3 
     DO jj = 1, 3 
      DO ii = 1, 3 
         Gmunu(mm,nn) = Gmunu(mm,nn) + g_contr_st(mm,ll)*kAB(ii,jj)*A_st(ii,ll)*A_st(jj,nn)  
      ENDDO        
     ENDDO 
    ENDDO
   ENDDO 
  ENDDO  
  !
  I1 = Gmunu(0,0) + Gmunu(1,1) + Gmunu(2,2) + Gmunu(3,3) 
  temp_st = MATMUL(Gmunu,Gmunu) 
  I2 = temp_st(0,0) + temp_st(1,1) + temp_st(2,2) + temp_st(3,3)   
  !
  myepsilon = cs**2/4.0*( I2 - 1./3.*I1**2 )    
  !
END SUBROUTINE GRTepsilonA 


SUBROUTINE func1(f,Vin,Qbar) 
    USE Parameters, ONLY : nVar 
    IMPLICIT NONE 
    REAL    :: f(5), Vin(5), Qbar(nVar), Q(nVar), V(nVar) 
    ! 
    V(1:5)    = Vin(1:5) 
    V(6:nVar) = Qbar(6:nVar) 
    !CALL PDEPrim2Cons(Q,V,(/0.,0.,0./),0.0)
	CALL PDEPrim2Cons(Q,V)	
    f = Qbar(1:5) - Q(1:5) 
    ! 
    CONTINUE
    !
END SUBROUTINE func1 


SUBROUTINE dfunc1(df,V,Qbar) 
    USE Parameters, ONLY : nVar
    IMPLICIT NONE 
    INTEGER :: i 
    REAL    :: df(5,5), Qbar(nVar), V(5) 
    REAL    :: Vp(5), Vm(5), fp(5), fm(5), eps = 1e-7  
    ! 
    DO i = 1, 5
        Vp = V
        Vp(i) = Vp(i) + eps 
        Vm = V
        Vm(i) = Vm(i) - eps 
        CALL func1(fp,Vp,Qbar) 
        CALL func1(fm,Vm,Qbar) 
        df(:,i) = (fp-fm)/(2*eps) 
    ENDDO    
    !
END SUBROUTINE dfunc1 
    

