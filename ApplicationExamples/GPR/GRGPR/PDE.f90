! GRMHD PDE.f90
! Trento (EQNTYPE4)

RECURSIVE SUBROUTINE PDEFlux(f,g,hz,Q)
  USE Parameters, ONLY: nVar, nDim, cs, gamma
  USE GRGPRmod
  USE iso_c_binding
  IMPLICIT NONE
  REAL :: f(nVar), g(nVar), h(nVar), hz(nVar), Q(nVar), V(nVar)
  REAL, PARAMETER :: epsilon = 1e-14 
  INTENT(IN)  :: Q
  INTENT(OUT) :: f, g, hz
  ! Local varialbes
INTEGER :: iter, Ierr,ii,jj,ll,mm,nn
REAL	:: kAB(3,3),shift(3),gammaij(6),g_cov(3,3),g_contr(3,3),v2
REAL    :: lapse, gp, gm,Qloc(nVar), gamma1, gam, A(3,3),detA, vf_cov(3),vf(3),Vtr(3),LF,shift_cov(3)
REAL 	:: rho,pi_work(3),ep,pi_scalar,A_st(1:3,0:3),g_cov_st(0:3,0:3),Gmunu(0:3,0:3),hmunu(0:3,0:3),devG_st(0:3,0:3),g_contr_st(0:3,0:3)
REAL	:: pi_mix(0:3,0:3), pi_st(0:3,0:3)
LOGICAL	:: FAILED
REAL    :: sm_cov(3), sm(3),Bv(3),s2,b2,sb,sb2,eps,e,x1,x2,w,p,vx,vy,vz,bx,by,bz,den,vb,vtr_old(3),res, ww, Spi_mix_st(0:3,0:3)
REAL 	::gamma_cov_st(0:4,0:4),gamma_mix_st(0:3,0:3),pi_contr_st(0:3,0:3),temp_st(0:3,0:3),Spi_contr_st(0:3,0:3), AU(3),detkAB,Qv_contr(3),rhoA,wwx,wwy,wwz

  ! Metric in matter space 
  kAB(1,1) = Q(25) 
  kAB(1,2) = Q(26) 
  kAB(1,3) = Q(27) 
  kAB(2,1) = Q(26) 
  kAB(2,2) = Q(28) 
  kAB(2,3) = Q(29) 
  kAB(3,1) = Q(27) 
  kAB(3,2) = Q(29) 
  kAB(3,3) = Q(30) 
  detkAB   = kAB(1,1)*kAB(2,2)*kAB(3,3)-kAB(1,1)*kAB(2,3)*kAB(3,2)-kAB(2,1)*kAB(1,2)*kAB(3,3)+kAB(2,1)*kAB(1,3)*kAB(3,2)+kAB(3,1)*kAB(1,2)*kAB(2,3)-kAB(3,1)*kAB(1,3)*kAB(2,2)  

  CALL PDECons2Prim(V,Q)
  !
  gamma1 = gamma/(gamma-1.0)
  rho    = V(1)
  vf_cov = V(2:4)
  p      = V(5)
  ! The 4 metric (decomposed in 3+1) 
  lapse = V(15)
  shift = V(16:18)
  g_cov(1,1) = V(19)
  g_cov(1,2) = V(20)
  g_cov(1,3) = V(21)
  g_cov(2,2) = V(22)
  g_cov(2,3) = V(23)
  g_cov(3,3) = V(24)
  g_cov(2,1) = g_cov(1,2)
  g_cov(3,1) = g_cov(1,3)
  g_cov(3,2) = g_cov(2,3) 
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
  vf = MATMUL(g_contr,vf_cov)
  Qv_contr = MATMUL(g_contr,Q(2:4))
  !
  v2     = vf(1)*vf_cov(1) + vf(2)*vf_cov(2) + vf(3)*vf_cov(3)
  !
  ! transport velocity
  Vtr(1:3) = lapse*vf(1:3)-shift(1:3)   ! this is contravariant !!!!
  ! conservative contribution to the transport of A 
  AU = MATMUL( A, Vtr(1:3) ) 
  ! compute the missing time components of the 3x4 deformation gradient (eq. 10 Gundlach2012)
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
  !
  ! compute the covariant 4D SPATIAL metric using the 3+1 split 
  gamma_cov_st(1:3,1:3) = g_cov
  shift_cov = MATMUL(g_cov,shift) 
  gamma_cov_st(0,0)   = DOT_PRODUCT(shift,shift_cov) 
  gamma_cov_st(0,1:3) = shift_cov 
  gamma_cov_st(1:3,0) = shift_cov 
  gamma_mix_st = MATMUL(g_contr_st(0:3,0:3),gamma_cov_st(0:3,0:3)) 
  !

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
  devG_st   = Gmunu - (Gmunu(0,0)+Gmunu(1,1)+Gmunu(2,2)+Gmunu(3,3))/3.0*hmunu 
  pi_mix    = rho*cs**2*MATMUL(Gmunu,devG_st)  
  !
  pi_contr_st = MATMUL(pi_mix, g_contr_st)     ! contravariant pi 
  pi_st       = MATMUL( g_cov_st, pi_mix )     ! covariant pi 
  !

  !
  CALL GRTepsilonA(g_contr_st,kAB,A_st,ep)   
  !
  lf     = 1.0/sqrt(1.0 - v2)
  w      = rho + rho*ep + gamma1*p   ! rho*enthalpy
  ww     = w*lf**2
  wwx    = ww*vf(1)
  wwy    = ww*vf(2)
  wwz    = ww*vf(3) 
  !
  rhoA = SQRT(detkAB)*detA/LF/gp                ! just for curiosity 
  !
  !!test = pi_contr_st(0,0)  
  !!testv = pi_mix(0,1:3) 
  !!DO ii = 1, 3 
  !! testv(ii) = testv(ii) - 1./lapse*SUM(pi_st(ii,1:3)*vf(1:3))        
  !! DO jj = 1, 3 
  !!    test = test - 1./lapse**2*pi_st(ii,jj)*vf(ii)*vf(jj) 
  !! ENDDO
  !!ENDDO  
  !
  !
  ! compute the stress tensor due to the anisotropic stresses:  S^ab = gamma^a_c gamma^b_d pi^cd
  temp_st = MATMUL(gamma_mix_st(0:3,0:3),pi_contr_st(0:3,0:3))   ! project columns
  Spi_contr_st = MATMUL(temp_st(0:3,0:3),gamma_mix_st(0:3,0:3))    ! project rows
  !
  ! compute the corresponding mixed stress tensor 
  Spi_mix_st = MATMUL(Spi_contr_st(0:3,0:3),g_cov_st(0:3,0:3))    ! contract second index (rows)
  ! 
  f(1)    = vf(1)*Q(1) !rho*lf   !Q(1) ! rho*lf 
  f(2)    = wwx*vf_cov(1) + p + Spi_mix_st(1,1) 
  f(3)    = wwx*vf_cov(2)     + Spi_mix_st(1,2) 
  f(4)    = wwx*vf_cov(3)     + Spi_mix_st(1,3) 
  f(5)    = Qv_contr(1) - f(1)  
  ! fluxes for the deformation gradient tensor A 
  f(6:14) = 0.0 
  f(6)    = AU(1) 
  f(9)    = AU(2) 
  f(12)   = AU(3) 
  !  lapse&shift&metric fluxes 
  f(15:nVar) = 0.
  !
  g(1)    = vf(2)*Q(1) !rho*lf   ! rho*lf
  g(2)    = wwy*vf_cov(1)      + Spi_mix_st(2,1)  
  g(3)    = wwy*vf_cov(2) + p  + Spi_mix_st(2,2) 
  g(4)    = wwy*vf_cov(3)      + Spi_mix_st(2,3)  
  g(5)    = Qv_contr(2) - g(1)  
  ! fluxes for the deformation gradient tensor A 
  g(6:14) = 0.0 
  g(7)    = AU(1) 
  g(10)   = AU(2) 
  g(13)   = AU(3) 
  !  lapse&shift&metric fluxes 
  g(15:nVar) = 0.
  !
  h(1)    = vf(3)*Q(1) !rho*lf   ! rho*lf   !
  h(2)    = wwz*vf_cov(1)      + Spi_mix_st(3,1)  
  h(3)    = wwz*vf_cov(2)      + Spi_mix_st(3,2)   
  h(4)    = wwz*vf_cov(3) + p  + Spi_mix_st(3,3)  
  h(5)    = Qv_contr(3) - h(1)  
  ! fluxes for the deformation gradient tensor A 
  h(6:14) = 0.0 
  h(8)    = AU(1) 
  h(11)   = AU(2) 
  h(14)   = AU(3) 
  !  lapse&shift&metric fluxes 
  h(15:nVar) = 0.
  ! 
  f(2:4)   = f(2:4)*gp
  g(2:4)   = g(2:4)*gp
  h(2:4)   = h(2:4)*gp
  ! Remember that Q(:) below contains already the factor gp, which is ok!
  f(1:5)   = lapse*f(1:5) - shift(1)*Q(1:5)
  g(1:5)   = lapse*g(1:5) - shift(2)*Q(1:5)
  h(1:5)   = lapse*h(1:5) - shift(3)*Q(1:5)	
	
	
	
	
  IF(nDim==3) THEN
    hz=h
  END IF  
  END SUBROUTINE PDEFlux


RECURSIVE SUBROUTINE PDENCP(BgradQ,Q,gradQ) 
   USE Parameters, ONLY :  nVar, nDim, gamma, cs
   USE GRGPRmod
   IMPLICIT NONE
   ! 11. Oct 21:40: This was a matrix BGradQ(nVar, nDim) but is a vector in spaceTimePredictorNonlinear
   REAL, INTENT(OUT) :: BgradQ(nVar)
   REAL, INTENT(IN)  :: gradQ(nVar, 3)
   REAL, INTENT(IN)  :: Q(nVar)
  ! Linear elasticity variables
   REAL :: u(3),VP(nVar) 
   REAL :: AQx(nVar), BQy(nVar), CQz(nVar) , Qx(nVar), Qy(nVar), Qz(nVar)

	REAL	:: kAB(3,3),shift(3),gammaij(6),g_cov(3,3),v2
	REAL    :: lapse, gp, gm,Qloc(nVar), gamma1, gam, A(3,3),detA, vf_cov(3),vf(3),Vtr(3),LF,shift_cov(3)
	REAL 	:: rho,pi_work(3),ep,pi_scalar,A_st(1:3,0:3),g_cov_st(0:3,0:3),Gmunu(0:3,0:3),hmunu(0:3,0:3),devG_st(0:3,0:3),g_contr_st(0:3,0:3)
	REAL	:: pi_mix(0:3,0:3), pi_st(0:3,0:3)
	LOGICAL	:: FAILED
	REAL    :: sm_cov(3), sm(3),Bv(3),s2,b2,sb,sb2,eps,e,x1,x2,w,p,vx,vy,vz,bx,by,bz,den,vb,vtr_old(3),res, ww, Spi_mix_st(0:3,0:3)
	REAL 	::gamma_cov_st(0:4,0:4),gamma_mix_st(0:3,0:3),pi_contr_st(0:3,0:3),temp_st(0:3,0:3),Spi_contr_st(0:3,0:3), AU(3),detkAB,Qv_contr(3),rhoA,wwx,wwy,wwz
   
	real	:: delta(3,3), W_ij,Wim,g_contr(3,3),S_contr(3),BV_contr(3),vxB_contr(3),vxB(3)
	real	:: Vc(nVar),uem	
	integer	:: i,j,m,nn,mm,ll, count,ii, jj
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
 
! =====================================================================================================================

  kAB(1,1) = Q(25) 
  kAB(1,2) = Q(26) 
  kAB(1,3) = Q(27) 
  kAB(2,1) = Q(26) 
  kAB(2,2) = Q(28) 
  kAB(2,3) = Q(29) 
  kAB(3,1) = Q(27) 
  kAB(3,2) = Q(29) 
  kAB(3,3) = Q(30) 
  detkAB   = kAB(1,1)*kAB(2,2)*kAB(3,3)-kAB(1,1)*kAB(2,3)*kAB(3,2)-kAB(2,1)*kAB(1,2)*kAB(3,3)+kAB(2,1)*kAB(1,3)*kAB(3,2)+kAB(3,1)*kAB(1,2)*kAB(2,3)-kAB(3,1)*kAB(1,3)*kAB(2,2) 
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
  delta = 0.
  DO i=1,3
      delta(i,i) = 1.0
  ENDDO 
  !
  CALL MatrixInverse3x3(g_cov,g_contr,gp)
  gp = SQRT(gp)
  gm = 1./gp
  ! 
  CALL PDECons2Prim(Vc,Q)

  gamma1 = gamma/(gamma-1.0)
  rho    = Vc(1)
  vf_cov = Vc(2:4)
  vf     = MATMUL(g_contr,vf_cov)
  ! transport velocity
  Vtr(1:3) = Q(15)*vf(1:3)-Q(16:18)   ! this is contravariant
  !
  ! The deformation gradient tensor 
  A(1,:) = (/ Vc( 6), Vc( 7),  Vc( 8) /) 
  A(2,:) = (/ Vc( 9), Vc(10),  Vc(11) /)
  A(3,:) = (/ Vc(12), Vc(13),  Vc(14) /)         
  detA   = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(2,1)*A(1,2)*A(3,3)+A(2,1)*A(1,3)*A(3,2)+A(3,1)*A(1,2)*A(2,3)-A(3,1)*A(1,3)*A(2,2) 
  ! 
  ! conservative contribution to the transport of A 
  AU = MATMUL( A, Vtr(1:3) ) 
  ! compute the missing time components of the 3x4 deformation gradient 
  A_st(1:3,1:3) = A(1:3,1:3) 
  DO ii = 1, 3
        A_st(ii,0) = -DOT_PRODUCT(Vtr(1:3),A(ii,:)) 
  ENDDO
  !
  p      = Vc(5)
  !
  Qv_contr = MATMUL(g_contr,Q(2:4))  
  !
  v2     = vf(1)*vf_cov(1) + vf(2)*vf_cov(2) + vf(3)*vf_cov(3)
  !
  S_contr = MATMUL(g_contr,Q(2:4))
  lf     = 1.0/sqrt(1.0 - v2)
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
  !
  ! compute the covariant 4D SPATIAL metric using the 3+1 split 
  gamma_cov_st(1:3,1:3) = g_cov
  shift_cov = MATMUL(g_cov,shift) 
  gamma_cov_st(0,0)   = DOT_PRODUCT(shift,shift_cov) 
  gamma_cov_st(0,1:3) = shift_cov 
  gamma_cov_st(1:3,0) = shift_cov 

  gamma_mix_st = MATMUL(g_contr_st(0:3,0:3),gamma_cov_st(0:3,0:3)) 
  !
  rhoA = SQRT(detkAB)*detA/LF/gp                ! just for curiosity 

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
  devG_st   = Gmunu - (Gmunu(0,0)+Gmunu(1,1)+Gmunu(2,2)+Gmunu(3,3))/3.0*hmunu 
  pi_mix    = rho*cs**2*MATMUL(Gmunu,devG_st)      
  pi_contr_st = MATMUL(pi_mix, g_contr_st)     
  !
  CALL GRTepsilonA(g_contr_st,kAB,A_st,ep) 
  !
  w      = rho + rho*ep + gamma1*p   ! rho*enthalpy
  ww     = w*lf**2          ! rho*enthalpy*Lorentz^2 
  !
  ! compute the stress tensor due to the anisotropic stresses:  S^ab = gamma^a_c gamma^b_d pi^cd
  temp_st = MATMUL(gamma_mix_st(0:3,0:3),pi_contr_st(0:3,0:3))     ! project columns
  Spi_contr_st = MATMUL(temp_st(0:3,0:3),gamma_mix_st(0:3,0:3))    ! project rows
  !
  ! compute the corresponding mixed stress tensor 
  Spi_mix_st = MATMUL(Spi_contr_st(0:3,0:3),g_cov_st(0:3,0:3))    ! contract second index (rows)
  !
  uem = 0.0
  BV  = 0.0
  BV_contr  = 0.0
  vxB       = 0.0
  vxB_contr = 0.0
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
                ! this is the stress tensor S^i_j
                W_ij = ww*vf_cov(i)*vf(j)-vxB(i)*vxB_contr(j)-BV(i)*BV_contr(j)+(p+uem)*delta(i,j) + Spi_mix_st(i,j)  
                !
                AQx(1+j) = AQx(1+j) - Q(1+i)*Qx(15+i)   ! Q(16:18)  shift(i) or shift_contr(i)
                AQx(5) = AQx(5) - gp*W_ij*Qx(15+i)      ! Q(16:18)  shift(i) or shift_contr(i)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=2
                !------ 
                W_ij = ww*vf_cov(i)*vf(j)-vxB(i)*vxB_contr(j)-BV(i)*BV_contr(j)+(p+uem)*delta(i,j) + Spi_mix_st(i,j)
                !
                BQy(1+j) = BQy(1+j) - Q(1+i)*Qy(15+i)   ! Q(16:18)  shift(i) or shift_contr(i)
                BQy(5) = BQy(5) - gp*W_ij*Qy(15+i)      ! Q(16:18)  shift(i) or shift_contr(i)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=3
                !------
                W_ij = ww*vf_cov(i)*vf(j)-vxB(i)*vxB_contr(j)-BV(i)*BV_contr(j)+(p+uem)*delta(i,j) + Spi_mix_st(i,j)
                !
                CQz(1+j) = CQz(1+j) - Q(1+i)*Qz(15+i)   ! Q(16:18)  shift(i) or shift_contr(i)
                CQz(5) = CQz(5) - gp*W_ij*Qz(15+i)      ! Q(16:18)  shift(i) or shift_contr(i)
                !
          DO m=1,3
            IF(m.GE.i) THEN  
                count=count+1
                !
                ! this is the contravariant stress tensor S^ij
                Wim = ww*vf(i)*vf(m)-vxB_contr(i)*vxB_contr(m)-BV_contr(i)*BV_contr(m)+(p+uem)*g_contr(i,m) + Spi_contr_st(i,m)
                Wim = Wim + (1.0 - delta(i,m))*(ww*vf(m)*vf(i)-vxB_contr(m)*vxB_contr(i)-BV_contr(m)*BV_contr(i)+(p+uem)*g_contr(m,i) + Spi_contr_st(m,i))  ! account also of the remaining symmetric components of gamma for i.NE.m.
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=1
                !------ 
                AQx(1+j) = AQx(1+j) - 0.5*gp*lapse*Wim*Qx(18+count)  ! Q(19:24) gammaij(count) or  g_cov(i,m)
                AQx(5) = AQx(5) - 0.5*gp*Wim*shift(j)*Qx(18+count)   ! Q(19:24) gammaij(count) or  g_cov(i,m)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=2
                !------ 
                BQy(1+j) = BQy(1+j) - 0.5*gp*lapse*Wim*Qy(18+count)  ! Q(19:24) gammaij(count) or  g_cov(i,m)
                BQy(5) = BQy(5) - 0.5*gp*Wim*shift(j)*Qy(18+count)   ! Q(19:24) gammaij(count) or  g_cov(i,m)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=3
                !------
                CQz(1+j) = CQz(1+j) - 0.5*gp*lapse*Wim*Qz(18+count)  ! Q(19:24) gammaij(count) or  g_cov(i,m)
                CQz(5) = CQz(5) - 0.5*gp*Wim*shift(j)*Qz(18+count)   ! Q(19:24) gammaij(count) or  g_cov(i,m)
                !
            ENDIF
          ENDDO
      ENDDO
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=1
    !------ 
    AQx(1+j) = AQx(1+j) + (Q(5)+Q(1))*Qx(15)    ! Q(15) or lapse
    AQx(5) = AQx(5) + S_contr(j)*Qx(15)         ! Q(15) or lapse
    !
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=2
    !------ 
    BQy(1+j) = BQy(1+j) + (Q(5)+Q(1))*Qy(15)    ! Q(15) or lapse
    BQy(5) = BQy(5) + S_contr(j)*Qy(15)         ! Q(15) or lapse
    !
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=3
    !------
    CQz(1+j) = CQz(1+j) + (Q(5)+Q(1))*Qz(15)    ! Q(15) or lapse
    CQz(5) = CQz(5) + S_contr(j)*Qz(15)         ! Q(15) or lapse
    !
    ! Non conservative product for the evolution equations of A
    ! transport velocity in general relativity is fluid velocity times lapse minus shift 
    ! 
    AQx(6) =  -vtr(2)*Qx(7) - vtr(3)*Qx(8) 
    BQy(6) =  +vtr(2)*Qy(6) 
    CQz(6) =  +vtr(3)*Qz(6)  
    ! 
    AQx(7) =  +vtr(1)*Qx(7) 
    BQy(7) =  -vtr(1)*Qy(6) - vtr(3)*Qy(8) 
    CQz(7) =  +vtr(3)*Qz(7) 
    !
    AQx(8) =  +vtr(1)*Qx(8) 
    BQy(8) =  +vtr(2)*Qy(8) 
    CQz(8) =  -vtr(1)*Qz(6) - vtr(2)*Qz(7) 
    ! 
    AQx(9) =  -vtr(2)*Qx(10) - vtr(3)*Qx(11) 
    BQy(9) =  +vtr(2)*Qy(9) 
    CQz(9) =  +vtr(3)*Qz(9) 
    !
    AQx(10) = +vtr(1)*Qx(10)     
    BQy(10) = -vtr(1)*Qy(9) - vtr(3)*Qy(11) 
    CQz(10) = +vtr(3)*Qz(10) 
    !
    AQx(11) = +vtr(1)*Qx(11) 
    BQy(11) = +vtr(2)*Qy(11) 
    CQz(11) = -vtr(1)*Qz(9) - vtr(2)*Qz(10) 
    !
    AQx(12) = -vtr(2)*Qx(13) - vtr(3)*Qx(14)  
    BQy(12) = +vtr(2)*Qy(12)
    CQz(12) = +vtr(3)*Qz(12) 
    !
    AQx(13) = +vtr(1)*Qx(13)  
    BQy(13) = -vtr(1)*Qy(12) - vtr(3)*Qy(14) 
    CQz(13) = +vtr(3)*Qz(13) 
    !
    AQx(14) = +vtr(1)*Qx(14) 
    BQy(14) = +vtr(2)*Qy(14) 
    CQz(14) = -vtr(1)*Qz(12) - vtr(2)*Qz(13) 
    !
    ! Advection of the matter metric 
    AQx(25:30) = vtr(1)*Qx(25:30) 
    BQy(25:30) = vtr(2)*Qy(25:30) 
    CQz(25:30) = vtr(3)*Qz(25:30) 


! =====================================================================================================================


 
 

    if( nDim .eq. 2) then
        BgradQ = AQx + BQy         
    else
        BgradQ = AQx + BQy + CQz     
    end if
END SUBROUTINE PDENCP


RECURSIVE SUBROUTINE PDEEigenvalues(L,Q,n)
  USE Parameters, ONLY :  nVar, nDim, gamma
  USE recipies_mod
  USE GRGPRmod
  USE iso_c_binding
  IMPLICIT NONE
  REAL :: L(nVar), n(3), Q(nVar), V(nVar),Qp(nVar),Qm(nVar)
  INTENT(IN)  :: Q,n
  INTENT(OUT) :: L 
  ! Local Variables 
  REAL :: fpp(nVar),gpp(nVar),hpp(nVar),fmm(nVar),gmm(nVar),hmm(nVar),dfdQ(nVar,nVar),AA(nVar,nVar)
  REAL :: shift(3),rho, vf_cov(3),p, lapse, g_cov(3,3),g_contr(3,3),gp,gm,vf(3),gamma1,w,cs2,vn,sft,gg,den,v2,eps
  REAL	:: u,nv(3),t1,rtemp(nVar), R(nVar,nVar),ml(2),ImLambda(nVar), GradQ(nVar,3), iR(nVar,nVar)
  INTEGER  :: nParam=0
  integer  :: ierr,itemp(nVar),j

  lapse = Q(15)
  shift = Q(16:18)
  ! 
  L = 0.0
  L = lapse*2.0 - DOT_PRODUCT(shift,n) 
	return
CALL PDECons2Prim(V,Q)
  rho    = V(1)
  vf_cov = V(2:4)
  p      = V(5)
  !
  lapse = V(15)
  shift = V(16:18)
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
  CALL MatrixInverse3x3(g_cov,g_contr,gp)
  gp = SQRT(gp)
  gm = 1./gp
  !  
  !CALL METRIC(xg, lapse, gp, gm, shift, g_cov, g_contr)
  !
  vf      = MATMUL(g_contr,vf_cov)
  v2     = vf(1)*vf_cov(1) + vf(2)*vf_cov(2) + vf(3)*vf_cov(3)
  !
  gamma1 = gamma/(gamma-1.0)
  !
  w      = rho + gamma1*p
  cs2    = gamma*p/w
  vn     = vf(1)*n(1) + vf(2)*n(2) + vf(3)*n(3)
  sft    = shift(1)*(n(1)) + shift(2)*(n(2)) + shift(3)*(n(3))
  gg     = g_contr(1,1)*ABS(n(1)) + g_contr(2,2)*ABS(n(2)) + g_contr(3,3)*ABS(n(3))
  den    = 1.0/(1.0 - v2*cs2)
  IF(SUM(n**2).EQ.0.) THEN  
     u = SQRT( v2) 
     WRITE(*,*)'Impossible error!'
     STOP
  ELSE
     u = vn 
  ENDIF
  L(1)   = ( u*(1.0-cs2) - SQRT( cs2*(1.0-v2)*( (1.0-v2*cs2)*gg - u**2*(1.0-cs2) )) )*den
  L(2:4) = u
  L(5)   = ( u*(1.0-cs2) + SQRT( cs2*(1.0-v2)*( (1.0-v2*cs2)*gg - u**2*(1.0-cs2) )) )*den 
  L(1:5)   = lapse*L(1:5) - sft
  !
  L(6:) = 0.
  
 ! RETURN 
  !
   IF(SUM(n**2).EQ.0.) THEN
       nv = (/ 1.0, 0.0, 0.0 /) 
   ELSE
       nv = n 
   ENDIF   
   gradQ = 0.0  
   dfdQ  = 0.0 
   DO j = 1, nVar-nParam  
       Qp = Q 
       Qm = Q 
       eps = MAX( 1e-7, 1e-7*ABS(Q(j)) ) 
       Qp(j) = Qp(j) + eps
       Qm(j) = Qm(j) - eps 
       CALL PDEFlux(fpp,gpp,hpp,Qp) 
       CALL PDEFlux(fmm,gmm,hmm,Qm)
       dfdQ(1:nVar-nParam,j) = (fpp(1:nVar-nParam)-fmm(1:nVar-nParam))/(2*eps)*nv(1) + (gpp(1:nVar-nParam)-gmm(1:nVar-nParam))/(2*eps)*nv(2) + (hpp(1:nVar-nParam)-hmm(1:nVar-nParam))/(2*eps)*nv(3)  
   ENDDO 
   CALL PDEJacobian(AA,Q,gradQ,nv) 
   AA = dfdQ + AA 
   !
   
   !  Compute the eigenstructure numerically 
   !where(abs(AA)<1e-14) AA = 0.0  
   CALL RG(nVar,nVar,AA,L,ImLambda,0,R,itemp,rtemp,ierr)
!	print*, L 
	!pause
   !where(abs(R)<1e-10) R = 0.0 
   !CALL MatrixInverse(nVar,R,iR) 
   ! 
   !t1 = MAXVAL(ABS(iR)) 
   !ml = MAXLOC(ABS(iR)) 
   !
  
END SUBROUTINE PDEEigenvalues

RECURSIVE SUBROUTINE PDESource(S,Q) 
  USE Parameters, ONLY:  nVar, nDim, gamma, cs,tau
  USE GRGPRmod
  USE iso_c_binding
  IMPLICIT NONE
  ! --------------------------------------------
  ! Argument list declaration
  REAL :: S(nvar), Q(nvar)
  INTENT(IN)  :: Q 
  INTENT(OUT) :: S
  ! --------------------------------------------
  ! Local variable declaration 
	REAL	:: kAB(3,3),shift(3),gammaij(6),g_cov(3,3),v2
	REAL    :: lapse, gp, gm,Qloc(nVar), gamma1, gam, A(3,3),detA, vf_cov(3),vf(3),Vtr(3),LF,shift_cov(3)
	REAL 	:: rho,pi_work(3),ep,pi_scalar,A_st(3,0:3),g_cov_st(0:3,0:3),Gmunu(0:3,0:3),hmunu(0:3,0:3),devG_st(0:3,0:3),g_contr_st(0:3,0:3)
	REAL	:: pi_mix(0:3,0:3), pi_st(0:3,0:3)
	LOGICAL	:: FAILED
	REAL    :: sm_cov(3), sm(3),Bv(3),s2,b2,sb,sb2,eps,e,x1,x2,w,p,vx,vy,vz,bx,by,bz,den,vb,vtr_old(3),res, ww, Spi_mix_st(0:3,0:3)
  REAL :: V(nvar),eAst(1:3,0:3),detA2, detkAB, I1, g_contr(3,3)
	INTEGER	:: ii,jj,nn,mm,ll

	S= 0. 
  S = 0.
  !
  CALL PDECons2Prim(V,Q)  
  rho = V(1) 
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
  detkAB   = kAB(1,1)*kAB(2,2)*kAB(3,3)-kAB(1,1)*kAB(2,3)*kAB(3,2)-kAB(2,1)*kAB(1,2)*kAB(3,3)+kAB(2,1)*kAB(1,3)*kAB(3,2)+kAB(3,1)*kAB(1,2)*kAB(2,3)-kAB(3,1)*kAB(1,3)*kAB(2,2) 
  
  ! 
  vf_cov = V(2:4)
  lapse  = V(15)
  shift  = V(16:18)
  g_cov(1,1) = V(19)
  g_cov(1,2) = V(20)
  g_cov(1,3) = V(21)
  g_cov(2,2) = V(22)
  g_cov(2,3) = V(23)
  g_cov(3,3) = V(24)
  g_cov(2,1) = g_cov(1,2)
  g_cov(3,1) = g_cov(1,3)
  g_cov(3,2) = g_cov(2,3) 
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
  vf = MATMUL(g_contr,vf_cov)
  v2     = vf(1)*vf_cov(1) + vf(2)*vf_cov(2) + vf(3)*vf_cov(3)
  lf     = 1.0/sqrt(1.0 - v2) 
  !
  ! transport velocity
  Vtr(1:3) = lapse*vf(1:3)-shift(1:3)   ! this is contravariant !!!!
  ! compute the missing time components of the 3x4 deformation gradient (eq. 10 Gundlach2012)
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

  
  !eA = 0.0 
  !eps = 1e-7 
  !DO ii = 1, 3   
  ! DO jj = 1, 3
  !      Ap_st = A_st   
  !      Ap_st(ii,jj) = Ap_st(ii,jj) + eps 
  !      Am_st = A_st   
  !      Am_st(ii,jj) = Am_st(ii,jj) - eps 
  !      CALL GRTepsilonA(g_contr_st,kAB,Ap_st,ep) 
  !      CALL GRTepsilonA(g_contr_st,kAB,Am_st,em)
  !      eA(ii,jj) = (ep-em)/(2.0*eps)  
  ! ENDDO
  !ENDDO 
  !
  ! in the derivative e_\psi^A_i the matter indices go down and the spatial indices come up, we therefore must multiply with kAB_contr and g_ij 
  !
  !eA = MATMUL(  MATMUL( kAB_contr, eA ), g_cov_st(1:3,1:3) ) 
  !
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
  I1 = ( Gmunu(0,0)+Gmunu(1,1)+Gmunu(2,2)+Gmunu(3,3) ) 
  devG_st   = Gmunu - I1/3.*MATMUL(g_contr_st,g_cov_st) 
  
  !
  !eAst = EQN%cs**2*MATMUL(kAB,MATMUL(A_st,MATMUL(devG_st,g_contr_st)))  
  !eAst = MATMUL(  MATMUL( kAB_contr, eAst ), g_cov_st ) 
  eAst = MATMUL(A_st,devG_st)  
  !
  S(6)  = eAst(1,1) 
  S(7)  = eAst(1,2) 
  S(8)  = eAst(1,3) 
  S(9)  = eAst(2,1) 
  S(10) = eAst(2,2) 
  S(11) = eAst(2,3) 
  S(12) = eAst(3,1) 
  S(13) = eAst(3,2) 
  S(14) = eAst(3,3) 
  ! 
  IF(cs>0.0) THEN
      detA2 = SQRT(detkAB)*detA/LF/gp     ! the density we want 
      S(6:14) = - rho/SQRT(detkAB)*I1/tau * S(6:14) ! - (detA2-rho)/detA2/EQN%tau*Q(6:14) 
  ENDIF 
	  
	  
END SUBROUTINE PDESource

RECURSIVE SUBROUTINE PDEVarName(MyNameOUT,ind) 
  USE Parameters, ONLY: nVar  
  IMPLICIT NONE     
  CHARACTER(LEN=10):: MyName(nVar),MyNameOUT
  INTEGER			:: ind

  ! EQNTYPE93
  MyName(1) = 'rho' 
  MyName(2) = 'u' 
  MyName(3) = 'v'
  MyName(4) = 'w'
  MyName(5) = 'p' 
  MyName(6)  = 'A11' 
  MyName(7)  = 'A12' 
  MyName(8)  = 'A13' 
  MyName(9)  = 'A21' 
  MyName(10) = 'A22' 
  MyName(11) = 'A23' 
  MyName(12) = 'A31' 
  MyName(13) = 'A32' 
  MyName(14) = 'A33' 
  ! 
  MyName(15) = 'lapse' 
  MyName(16) = 'Shift^1' 
  MyName(17) = 'Shift^2' 
  MyName(18) = 'Shift^3' 
  MyName(19) = 'Gamma_11' 
  MyName(20) = 'Gamma_12' 
  MyName(21) = 'Gamma_13' 
  MyName(22) = 'Gamma_22' 
  MyName(23) = 'Gamma_23' 
  MyName(24) = 'Gamma_33' 
  MyName(25) = 'k_11' 
  MyName(26) = 'k_12' 
  MyName(27) = 'k_13' 
  MyName(28) = 'k_22' 
  MyName(29) = 'k_23' 
  MyName(30) = 'k_33' 
	
	MyNameOUT=MyName(ind+1)
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

    if( nDim .eq. 2) then
        An = A*nv(1) + B*nv(2)         
    else
        An = A*nv(1) + B*nv(2) + C*nv(3)     
    end if
    
    

  
END SUBROUTINE PDEMatrixB


RECURSIVE SUBROUTINE PDEJacobian(An,Q,gradQ,nv) 
  USE Parameters, ONLY : nVar, nDim, gamma, cs
  USE GRGPRmod
  IMPLICIT NONE
  ! Argument list 
  REAL :: An(nVar,nVar)
  REAL :: Q(nVar), gradQ(nVar,3), nv(3) 
  INTENT(IN)  :: Q, gradQ, nv
  INTENT(OUT) :: An  
  ! Local variables
  INTEGER :: i,j,ii,jj,ll,kk,mm,nn,iErr,count,m
  REAL :: A(nVar,nVar), B(nVar,nVar), C(nVar,nVar), Vp(nVar) 
  REAL :: AQx(nVar), BQy(nVar), CQz(nVar), AQx0(nVar), eu(nVar), eu0(nVar), ev(nVar,3), uv(3)  
  REAL :: AQxp(nVar), AQxm(nVar), BQyp(nVar), BQym(nVar), temp(nVar)  
  REAL :: vx, vy, vz, cl 
  REAL :: u,v,h,pg,k,sigma, vtr(3) 
  REAL :: nx, ny, uu, vv, ww, MD12, rho   
  REAL :: t1,t2,t3,t5,t8,t11,t13,t17,t25,t29
  REAL :: g_cov(3,3), g_contr(3,3), g_covx(3,3), g_covy(3,3), g_covz(3,3)
  REAL :: k1,k2,k3,ff,e,det, alpha, fa, s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,k0,beta0(3),b0(3) 
  REAL :: Gtilde(3), Christoffel(3,3,3), Z(3), Zup(3), traceA, Atildeup(3,3), ggg, fff, QG(3), dgup(3,3,3), xi, ds, sk, sknl
  REAL :: lam, mu, irho, ialpha 
  REAL, PARAMETER :: epsilon = 1e-7  
  REAL :: S_contr(3),psi,gammaij(6),lapse,shift(3),gp,gm,ComputeDet,delta(3,3),Wim,W_ij,gamma1,Vc(nVar),vf(3),vf_cov(3),BV(3),p,v2,lf,w !,gv(3),gv_contr(3)
  REAL :: Qv_contr(3),vxB(3),vxB_contr(3),BV_contr(3),e2,b2,uem ! QB_contr(3),
  REAL :: alp1, alp2, rho1a1sq, rho2a2sq, kappaorho, ibnorm
  REAL :: A_st(3,0:3), rhoA   , g_cov_st(0:3,0:3),g_contr_st(0:3,0:3), gamma_cov_st(0:3,0:3), gamma_contr_st(0:3,0:3), gamma_mix_st(0:3,0:3), Spi_contr_st(0:3,0:3), Spi_mix_st(0:3,0:3),temp_st(0:3,0:3)
  REAL :: I1,I2,detgAB, piAB(3,3), piAB1(3,3), piAB2(3,3), piAB_st(0:3,0:3),piAB_contr_st(0:3,0:3), pi_mix_st(0:3,0:3), pi_scalar, f1, f2   
  REAL :: gAB_contr(3,3), gAB(3,3), kAB(3,3), kAB_contr(3,3), eta_mix(3,3), eta_cov(3,3)    
  REAL :: detA,  AU(3), AM(3,3), shift_cov(3), eps, dedA(1:3,0:3), App(1:3,0:3), Amm(1:3,0:3), ep, em, pi_mix(0:3,0:3), pi_contr_st(0:3,0:3)     
  REAL :: u_contr(0:3), u_cov(0:3), gmunu(0:3,0:3), hmunu(0:3,0:3), devG_st(0:3,0:3), pi_st(0:3,0:3), chi    
  ! 
  ! 

  kAB(1,1) = Q(25) 
  kAB(1,2) = Q(26) 
  kAB(1,3) = Q(27) 
  kAB(2,1) = Q(26) 
  kAB(2,2) = Q(28) 
  kAB(2,3) = Q(29) 
  kAB(3,1) = Q(27) 
  kAB(3,2) = Q(29) 
  kAB(3,3) = Q(30) 


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
  delta = 0.
  DO i=1,3
      delta(i,i) = 1.0
  ENDDO 
  !
  CALL MatrixInverse3x3(g_cov,g_contr,gp)
  gp = SQRT(gp)
  gm = 1./gp
  ! 
  CALL PDECons2Prim(Vc,Q)
  gamma1 = gamma/(gamma-1.0)
  rho    = Vc(1)
  vf_cov = Vc(2:4)
  p      = Vc(5)
  !
  vf     = MATMUL(g_contr,vf_cov)
  S_contr = MATMUL(g_contr,Q(2:4))
  !
  v2     = vf(1)*vf_cov(1) + vf(2)*vf_cov(2) + vf(3)*vf_cov(3)
  !
  ! transport velocity
  Vtr(1:3) = Q(15)*vf(1:3)-Q(16:18)   ! this is contravariant !!!!
  !
  ! The deformation gradient tensor 
  AM(1,:) = (/ Vc( 6), Vc( 7),  Vc( 8) /) 
  AM(2,:) = (/ Vc( 9), Vc(10),  Vc(11) /)
  AM(3,:) = (/ Vc(12), Vc(13),  Vc(14) /)         
  detA   = AM(1,1)*AM(2,2)*AM(3,3)-AM(1,1)*AM(2,3)*AM(3,2)-AM(2,1)*AM(1,2)*AM(3,3)+AM(2,1)*AM(1,3)*AM(3,2)+AM(3,1)*AM(1,2)*AM(2,3)-AM(3,1)*AM(1,3)*AM(2,2) 
  ! 
  ! conservative contribution to the transport of A 
  AU = MATMUL( AM, Vtr(1:3) ) 
  ! compute the missing time components of the 3x4 deformation gradient 
  A_st(1:3,1:3) = AM(1:3,1:3) 
  DO ii = 1, 3
        A_st(ii,0) = -DOT_PRODUCT(Vtr(1:3),AM(ii,:)) 
  ENDDO
  !
  !
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
  !
  ! compute the covariant 4D SPATIAL metric using the 3+1 split 
  gamma_cov_st(1:3,1:3) = g_cov
  shift_cov = MATMUL(g_cov,shift) 
  gamma_cov_st(0,0)   = DOT_PRODUCT(shift,shift_cov) 
  gamma_cov_st(0,1:3) = shift_cov 
  gamma_cov_st(1:3,0) = shift_cov 
  !! NOT USED: compute the inverse of the 4D covariant SPATIAL metric 
  !gamma_cov_st(1:3,1:3) = g_cov
  !gamma_contr_st(0,0)   = 0.
  !gamma_contr_st(0,1:3) = 0. 
  !gamma_contr_st(1:3,0) = 0. 
  ! compute the PROJECTION SPATIAL OPERATOR i.e. mixed SPATIAL metric
  gamma_mix_st = MATMUL(g_contr_st(0:3,0:3),gamma_cov_st(0:3,0:3)) 

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
  pi_st  = MATMUL( g_cov_st, pi_mix )     
  pi_contr_st = MATMUL(pi_mix, g_contr_st)   
  !
  CALL GRTepsilonA(g_contr_st,kAB,A_st,ep) 
  !
  lf     = 1.0/sqrt(1.0 - v2)
  w      = rho + rho*ep + gamma1*p   ! rho*enthalpy
  ww     = w*lf**2          ! rho*enthalpy*Lorentz^2 
  !
  ! compute the stress tensor due to the anisotropic stresses:  S^ab = gamma^a_c gamma^b_d pi^cd
  temp_st = MATMUL(gamma_mix_st(0:3,0:3),pi_contr_st(0:3,0:3))     ! project columns
  Spi_contr_st = MATMUL(temp_st(0:3,0:3),gamma_mix_st(0:3,0:3))    ! project rows
  !
  ! compute the corresponding mixed stress tensor 
  Spi_mix_st = MATMUL(Spi_contr_st(0:3,0:3),g_cov_st(0:3,0:3))    ! contract second index (rows)
  !
  uem = 0.0
  BV  = 0.0
  BV_contr  = 0.0
  vxB       = 0.0 
  vxB_contr = 0.0 
  !
  !DO j=1,3
  A = 0.
  B = 0.
  C = 0.
    !lapse
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=1
    !------ 
    A(1+j,15) = + (Q(5)+Q(1))   ! Q(15) or lapse
    A(5,15) =  S_contr(j)     !  Q(15) or lapse
    !
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=2
    !------ 
    B(1+j,15) = + (Q(5)+Q(1))   ! Q(15) or lapse
    B(5,15) =  S_contr(j)     !  Q(15) or lapse
    ! 
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=3
    !------
    C(1+j,15) = + (Q(5)+Q(1))   ! Q(15) or lapse
    C(5,15) =  S_contr(j)     !  Q(15) or lapse
    ! 
    count=0
    DO i=1,3
        ! shift
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        j=1
        !------ 
        W_ij = ww*vf_cov(i)*vf(j)-vxB(i)*vxB_contr(j)-BV(i)*BV_contr(j)+(p+uem)*delta(i,j) + Spi_mix_st(i,j)
        !
        A(1+j,15+i) = - Q(1+i)  ! Q(16:18)  shift(i) or shift_contr(i)
        A(5,15+i) = - gp*W_ij    ! Q(16:18)  shift(i) or shift_contr(i)
        !
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        j=2
        !------ 
        W_ij = ww*vf_cov(i)*vf(j)-vxB(i)*vxB_contr(j)-BV(i)*BV_contr(j)+(p+uem)*delta(i,j) + Spi_mix_st(i,j)
        !
        B(1+j,15+i) = - Q(1+i)  ! Q(16:18)  shift(i) or shift_contr(i)
        B(5,15+i) = - gp*W_ij    ! Q(16:18)  shift(i) or shift_contr(i)
        ! 
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        j=3
        !------
        W_ij = ww*vf_cov(i)*vf(j)-vxB(i)*vxB_contr(j)-BV(i)*BV_contr(j)+(p+uem)*delta(i,j) + Spi_mix_st(i,j)
        !
        C(1+j,15+i) = - Q(1+i)  ! Q(16:18)  shift(i) or shift_contr(i)
        C(5,15+i) = - gp*W_ij    ! Q(16:18)  shift(i) or shift_contr(i)
        ! 
          DO m=1,3
            IF(m.GE.i) THEN  
                !metric
                count=count+1
                !
                Wim = ww*vf(i)*vf(m)-vxB_contr(i)*vxB_contr(m)-BV_contr(i)*BV_contr(m)+(p+uem)*g_contr(i,m)+Spi_contr_st(i,m)
                Wim = Wim + (1.0 - delta(i,m))*(ww*vf(m)*vf(i)-vxB_contr(m)*vxB_contr(i)-BV_contr(m)*BV_contr(i)+(p+uem)*g_contr(m,i)+Spi_contr_st(m,i))  ! account also of the remaining symmetric components of gamma for i.NE.m.
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=1
                !------ 
                A(1+j,18+count) = - 0.5*gp*lapse*Wim  ! Q(19:24) gammaij(count) or  g_cov(i,m)
                A(5,18+count) = - 0.5*gp*Wim*shift(j)   ! Q(19:24) gammaij(count) or  g_cov(i,m)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=2
                !------ 
                B(1+j,18+count) = - 0.5*gp*lapse*Wim  ! Q(19:24) gammaij(count) or  g_cov(i,m)
                B(5,18+count) = - 0.5*gp*Wim*shift(j)   ! Q(19:24) gammaij(count) or  g_cov(i,m)
                ! 
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=3
                !------
                C(1+j,18+count) = - 0.5*gp*lapse*Wim  ! Q(19:24) gammaij(count) or  g_cov(i,m)
                C(5,18+count) = - 0.5*gp*Wim*shift(j)   ! Q(19:24) gammaij(count) or  g_cov(i,m)
                ! 
            ENDIF
          ENDDO
      ENDDO
  !ENDDO 
  !
    A(6,7) = -vtr(2)
    A(6,8) = -vtr(3)
    B(6,6) = +vtr(2)
    C(6,6) = +vtr(3) 
    ! 
    A(7,7) = +vtr(1)
    B(7,6) = -vtr(1)
    B(7,8) = -vtr(3)
    C(7,7) = +vtr(3)
    !
    A(8,8) = +vtr(1)
    B(8,8) = +vtr(2)
    C(8,6) = -vtr(1)
    C(8,7) = -vtr(2)
    ! 
    A(9,10) = -vtr(2)
    A(9,11) = -vtr(3)
    B(9,9)  = +vtr(2)
    C(9,9)  = +vtr(3) 
    !
    A(10,10) = +vtr(1)
    B(10, 9) = -vtr(1)
    B(10,11) = -vtr(3)
    C(10,10) = +vtr(3)
    !
    A(11,11) = +vtr(1)
    B(11,11) = +vtr(2)
    C(11, 9) = -vtr(1)
    C(11,10) = -vtr(2)
    !
    A(12,13) = -vtr(2)
    A(12,14) = -vtr(3)
    B(12,12) = +vtr(2)
    C(12,12) = +vtr(3)
    !
    A(13,13) = +vtr(1)
    B(13,12) = -vtr(1)
    B(13,14) = -vtr(3)
    C(13,13) = +vtr(3)
    !
    A(14,14) = +vtr(1)
    B(14,14) = +vtr(2)
    C(14,12) = -vtr(1)
    C(14,13) = -vtr(2)
    !
    DO i = 25, 30 
      A(i,i) = vtr(1) 
      B(i,i) = vtr(2) 
      C(i,i) = vtr(3)
    ENDDO     
    !
    An = A*nv(1) + B*nv(2) + C*nv(3) 
    !
 

END SUBROUTINE PDEJacobian

RECURSIVE SUBROUTINE getNumericalSolution(V,Q) 
  USE Parameters, ONLY: nVar  
  IMPLICIT NONE     
  REAL				:: V(nVar), Q(nVar)
  V=0.
END SUBROUTINE getNumericalSolution

RECURSIVE SUBROUTINE getExactSolution(V,pos, timeStamp) 
  USE Parameters, ONLY: nVar , nDim  
  IMPLICIT NONE     
  REAL				:: V(nVar), Q(nVar), pos(nDim), timeStamp
  call InitialData(pos, timeStamp, Q)
  V=1
END SUBROUTINE getExactSolution


