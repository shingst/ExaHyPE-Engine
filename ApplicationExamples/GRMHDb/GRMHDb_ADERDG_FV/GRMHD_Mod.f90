
#define C2PFF    

MODULE GRMHD_Mod
	USE Parameters, ONLY : nDim, EQN, nVar
	IMPLICIT NONE
	PRIVATE
	
	! - PRIVATE
    !INTERFACE FUNC_C2P_RMHD1_VECTORF
    !    MODULE PROCEDURE FUNC_C2P_RMHD1_VECTORF
    !END INTERFACE 
    INTERFACE RTSAFE_C2P_RMHD1
        MODULE PROCEDURE RTSAFE_C2P_RMHD1
    END INTERFACE 
	! - PUBLIC
    INTERFACE PDEFluxGRMHD
        MODULE PROCEDURE PDEFluxGRMHD
    END INTERFACE
    INTERFACE PDESourceGRMHD
        MODULE PROCEDURE PDESourceGRMHD
    END INTERFACE
    INTERFACE PDENCPGRMHD
        MODULE PROCEDURE PDENCPGRMHD
    END INTERFACE
    INTERFACE PDEMatrixBGRMHD
        MODULE PROCEDURE PDEMatrixBGRMHD
    END INTERFACE
    INTERFACE PDECons2PrimGRMHD
        MODULE PROCEDURE PDECons2PrimGRMHD
    END INTERFACE
    INTERFACE PDEPrim2ConsGRMHD
        MODULE PROCEDURE PDEPrim2ConsGRMHD
    END INTERFACE
    INTERFACE PDEEigenvaluesGRMHD
        MODULE PROCEDURE PDEEigenvaluesGRMHD
    END INTERFACE 
    INTERFACE PDEFluxPrimGRMHD
    	MODULE PROCEDURE PDEFluxPrimGRMHD
    END INTERFACE	
    INTERFACE 	PDEFluxPrimVectorGRMHD
    	MODULE PROCEDURE PDEFluxPrimVectorGRMHD
    END INTERFACE
    INTERFACE  PDENCPPrimGRMHD
        MODULE PROCEDURE PDENCPPrimGRMHD
    END INTERFACE 
	INTERFACE  PDENCPPrimVectorGRMHD
		MODULE PROCEDURE PDENCPPrimVectorGRMHD
	END INTERFACE 
	INTERFACE  PDECons2PrimVectorGRMHD
		MODULE PROCEDURE PDECons2PrimVectorGRMHD
	END INTERFACE 
	INTERFACE  PDEPrim2ConsGRMHDvector
		MODULE PROCEDURE PDEPrim2ConsGRMHDvector
	END INTERFACE 
	INTERFACE  PDESourcePrimGRMHD
		MODULE PROCEDURE PDESourcePrimGRMHD
	END INTERFACE 
	INTERFACE  PDEVarNameGRMHD
		MODULE PROCEDURE PDEVarNameGRMHD
	END INTERFACE 
	INTERFACE  PDEAuxNameGRMHD
		MODULE PROCEDURE PDEAuxNameGRMHD
	END INTERFACE 
	INTERFACE  PDEAuxVarGRMHD
		MODULE PROCEDURE PDEAuxVarGRMHD
	END INTERFACE 
	
	
	PUBLIC :: PDEFluxGRMHD, PDESourceGRMHD,PDENCPGRMHD,PDEMatrixBGRMHD,PDECons2PrimGRMHD,PDEPrim2ConsGRMHD,  PDEEigenvaluesGRMHD, &
 PDEFluxPrimGRMHD, PDEFluxPrimVectorGRMHD, PDENCPPrimGRMHD, PDENCPPrimVectorGRMHD,   &
PDECons2PrimVectorGRMHD, PDEPrim2ConsGRMHDvector,PDESourcePrimGRMHD, &
RTSAFE_C2P_RMHD1,PDEVarNameGRMHD,PDEAuxVarGRMHD,PDEAuxNameGRMHD
				

CONTAINS


!
! Conservative part of the PDE ( flux tensor F(Q) ) 
!
RECURSIVE SUBROUTINE PDEFluxGRMHD(F,Q)
    USE Parameters, ONLY : d, EQN
    USE AstroMod
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
    REAL, INTENT(IN)  :: Q(nVar) 
    REAL, INTENT(OUT) :: F(nVar,d) 
    ! Local variables 
    INTEGER :: iErr 
    INTEGER :: i,j,k,m
    REAL :: V(nVar) 
    REAL :: p, irho, lam, mu 
    REAL :: k1, k2, fff, ggg, e, ds, c, xi, sk, sknl, alpha, fa, k0, beta0(3)  
    REAL :: gamma1, rho, vx, vy, vz, bx, by, bz, ex, ey, ez, v2, b2, e2, lf, w, ww, uem, wwx, wwy, wwz 
    REAL :: gp, gm, g_cov(3,3), g_contr(3,3), vxB_cov(3), vxB_contr(3), B_cov(3), BQ(3), vtr(3), v_contr(3), lapse, shift(3)    
    REAL :: Fij(3,3), v_cov(3), S_contr(3), QB_contr(3), B_contr(3), psi
    REAL :: rhoeps,dd,tau,s2,sm_cov(3),sm(3)
    !
    F = 0.0 
    !
	!IF( ANY( ISNAN(Q) )) THEN
	!	PRINT *,' Q NAN ',nVar
 !       PRINT *, Q
	!	STOP
 !   ENDIF 
 !   !
	!IF(ANY(Q(6:8).NE.0.0)) THEN
 !       PRINT *, Q(6:8)
 !       PRINT *, "I feel magnetized :-( PDEFluxGRMHD in"
 !       ERROR STOP
 !   ENDIF
    !
  CALL PDECons2PrimGRMHD(V,Q,iErr)
  !
	!IF(ANY(Q(6:8).NE.0.0)) THEN
 !       PRINT *, Q(6:8)
 !       PRINT *, "I feel magnetized :-( PDEFluxGRMHD in 2"
 !       ERROR STOP
 !   ENDIF
    !
  gamma1 = EQN%gamma/(EQN%gamma-1.0)
  rho    = V(1)
  !
	!IF( ANY( ISNAN(Q) )) THEN
	!	PRINT *,' Q NAN ',iErr,nVar
 !       PRINT *, Q
	!	ERROR STOP
 !   ENDIF 
 !   !
	!IF( ANY( ISNAN(V) )) THEN
	!	PRINT *,' V NAN ',iErr ,nVar
 !       PRINT *, V
	!	PRINT *,' Q:' 
 !       PRINT *, Q
 !       !i=20
 !       !V(i) = Q(1)/(Q(2)-Q(2))
	!	ERROR STOP
 !   ENDIF 
    
  DO i=1,3
	  v_cov(i) = V(1+i)
	  B_contr(i) = V(5+i)  		! B is contravariant
	  QB_contr(i) = Q(5+i)  		! B is contravariant
	  shift(i) = V(10+i)
  ENDDO
  p = V(5)
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
  !
  dd = gm*Q(1)
  tau = gm*Q(5)
  DO i=1,3 
      sm_cov(i)  = Q(1+i)
  ENDDO
  sm_cov = gm*sm_cov
  sm   = MATMUL (g_contr, sm_cov)  
  s2   = sm_cov(1)*sm(1) + sm_cov(2)*sm(2) + sm_cov(3)*sm(3)
  !
#ifdef GEOS  
  CALL EOS_eps(rhoeps,rho,p,dd,tau,s2) 
  w      = rho + rhoeps + p   ! rho*enthalpy
#else
  w    = rho + gamma1*p  ! this is rho*h
#endif
  !IF(ABS(w-(rho + gamma1*p)).GT.1e-13) THEN
  !    continue
  !ENDIF
  !
  !w      = rho + gamma1*p
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
	!IF(ANY(F(6:8,1:3).NE.0.0)) THEN
 !       PRINT *,F(6:8,1)
 !       PRINT *,F(6:8,2)
 !       PRINT *,F(6:8,3)
 !       PRINT *, "I feel magnetized :-( PDEFluxGRMHD "
 !       ERROR STOP
 !   ENDIF
    !
  CONTINUE      
  !
END SUBROUTINE PDEFluxGRMHD
!
    
!
! Conservative part of the PDE ( flux tensor F(Q) ) 
!
RECURSIVE  SUBROUTINE PDEFluxPrimGRMHD(F,V,Q)
    USE Parameters, ONLY : d, EQN
    USE AstroMod
    USE EOS_mod
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
    REAL, INTENT(IN)  :: Q(nVar), V(nVar)  
    REAL, INTENT(OUT) :: F(nVar,d) 
    ! Local variables 
    INTEGER :: iErr 
    INTEGER :: i,j,k,m 
    REAL :: p, irho, lam, mu 
    REAL :: k1, k2, fff, ggg, e, ds, c, xi, sk, sknl, alpha, fa, k0, beta0(3)  
    REAL :: gamma1, rho, vx, vy, vz, bx, by, bz, ex, ey, ez, v2, b2, e2, lf, w, ww, uem, wwx, wwy, wwz 
    REAL :: gp, gm, g_cov(3,3), g_contr(3,3), vxB_cov(3), vxB_contr(3), B_cov(3), BQ(3), vtr(3), v_contr(3), lapse, shift(3)    
    REAL :: Fij(3,3), v_cov(3), S_contr(3), QB_contr(3), B_contr(3), psi
    REAL :: rhoeps,dd,tau,s2,sm_cov(3),sm(3)
    !
    F = 0.0 
    !
  gamma1 = EQN%gamma/(EQN%gamma-1.0)
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
  !
  dd = gm*Q(1)
  tau = gm*Q(5)
  DO i=1,3 
      sm_cov(i)  = Q(1+i)
  ENDDO
  sm_cov = gm*sm_cov
  sm   = MATMUL (g_contr, sm_cov)  
  s2   = sm_cov(1)*sm(1) + sm_cov(2)*sm(2) + sm_cov(3)*sm(3)
  !
#ifdef GEOS  
  CALL EOS_eps(rhoeps,rho,p,dd,tau,s2) 
  w      = rho + rhoeps + p   ! rho*enthalpy
#else
  w    = rho + gamma1*p  ! this is rho*h
#endif
  !IF(ABS(w-(rho + gamma1*p)).GT.1e-13) THEN
  !    continue
  !ENDIF
  !
  !w      = rho + gamma1*p
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
  CONTINUE      
  !   
END SUBROUTINE PDEFluxPrimGRMHD
!
RECURSIVE SUBROUTINE PDEFluxPrimVectorGRMHD(F,V,Q) 
    USE Parameters, ONLY : d, EQN,VECTORLENGTH
    USE AstroMod
    USE EOS_mod
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
    REAL(8), INTENT(IN)  :: V(VECTORLENGTH,nVar),  Q(VECTORLENGTH,nVar)
    REAL(8), INTENT(OUT) :: F(VECTORLENGTH,nVar,d) 
    ! Local variables 
    INTEGER :: iErr 
    INTEGER :: i,j,k,m 
    REAL(8) :: p(VECTORLENGTH),rhoeps(VECTORLENGTH),DD(VECTORLENGTH),tau(VECTORLENGTH),s2(VECTORLENGTH)
    REAL(8) :: e(VECTORLENGTH), ds(VECTORLENGTH), c(VECTORLENGTH), xi(VECTORLENGTH), sk(VECTORLENGTH), sknl(VECTORLENGTH)  
    REAL(8) :: alpha(VECTORLENGTH), fa(VECTORLENGTH), k0(VECTORLENGTH), beta0(VECTORLENGTH,3), det(VECTORLENGTH), idet(VECTORLENGTH)  
    REAL(8) :: gamma1(VECTORLENGTH), rho(VECTORLENGTH), vx(VECTORLENGTH), vy(VECTORLENGTH), vz(VECTORLENGTH), bx(VECTORLENGTH)
    REAL(8) :: by(VECTORLENGTH), bz(VECTORLENGTH), ex(VECTORLENGTH), ey(VECTORLENGTH), ez(VECTORLENGTH), v2(VECTORLENGTH), b2(VECTORLENGTH)
    REAL(8) :: e2(VECTORLENGTH), lf(VECTORLENGTH), w(VECTORLENGTH), ww(VECTORLENGTH), uem(VECTORLENGTH), wwx(VECTORLENGTH), wwy(VECTORLENGTH), wwz(VECTORLENGTH)  
    REAL(8) :: gp(VECTORLENGTH), gm(VECTORLENGTH), g_cov(VECTORLENGTH,3,3), g_contr(VECTORLENGTH,3,3), vxB_cov(VECTORLENGTH,3), vxB_contr(VECTORLENGTH,3), B_cov(VECTORLENGTH,3),vtr(VECTORLENGTH,3), v_contr(VECTORLENGTH,3), lapse(VECTORLENGTH), shift(VECTORLENGTH,3)    
    REAL(8) :: Fij(VECTORLENGTH,3,3),v_cov(VECTORLENGTH,3), S_contr(VECTORLENGTH,3), gS_contr(VECTORLENGTH,3), S_cov(VECTORLENGTH,3), QB_contr(VECTORLENGTH,3), B_contr(VECTORLENGTH,3), psi(VECTORLENGTH)  
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN:  64 :: p,e,ds,c,xi,sk,sknl,alpha,fa,k0,beta0,det,idet,gamma1,rho,vx,vy,vz,bx,by,bz,ex,ey,ez,v2,b2,e2,lf,w,ww,uem,wwx,wwy,wwz
    !DIR$ ATTRIBUTES ALIGN:  64 :: gp,gm,g_cov,g_contr,vxB_cov,vxB_contr, B_cov, vtr, v_contr, lapse, shift, Fij, v_cov, S_contr,  gS_contr,  S_cov, QB_contr, B_contr, psi,rhoeps,DD,tau,s2  
    !DIR$ ASSUME_ALIGNED F : 64  
    !DIR$ ASSUME_ALIGNED V : 64  
    !DIR$ ASSUME_ALIGNED Q : 64    
#else
    !DIR$ ATTRIBUTES ALIGN:  32 :: p,e,ds,c,xi,sk,sknl,alpha,fa,k0,beta0,det,idet,gamma1,rho,vx,vy,vz,bx,by,bz,ex,ey,ez,v2,b2,e2,lf,w,ww,uem,wwx,wwy,wwz
    !DIR$ ATTRIBUTES ALIGN:  32 :: gp,gm,g_cov,g_contr,vxB_cov,vxB_contr,B_cov,vtr,v_contr,lapse,shift,Fij,v_cov,S_contr, gS_contr, S_cov,QB_contr,B_contr,psi,rhoeps,DD,tau,s2 
    !DIR$ ASSUME_ALIGNED F : 32 
    !DIR$ ASSUME_ALIGNED V : 32 
    !DIR$ ASSUME_ALIGNED Q : 32 
#endif 
  !
  gamma1 = EQN%gamma/(EQN%gamma-1.0)
  rho    = V(:,1)
  !
  DO i=1,3
	  v_cov(:,i)    = V(:,1+i)
	  B_contr(:,i)  = V(:,5+i)  		! B is contravariant
	  QB_contr(:,i) = Q(:,5+i)  		! B is contravariant
	  shift(:,i)    = V(:,10+i)
  ENDDO
  !
  p(:)     = V(:,5)
  psi(:)   = V(:,9)
  lapse(:) = V(:,10)
  !
  !gammaij = V(14:19) 
  g_cov(:,1,1) = V(:,14)
  g_cov(:,1,2) = V(:,15)
  g_cov(:,1,3) = V(:,16)
  g_cov(:,2,2) = V(:,17)
  g_cov(:,2,3) = V(:,18)
  g_cov(:,3,3) = V(:,19)
  g_cov(:,2,1) = V(:,15)
  g_cov(:,3,1) = V(:,16)
  g_cov(:,3,2) = V(:,18)
  !
  det(:) = (g_cov(:,1,1)*g_cov(:,2,2)*g_cov(:,3,3)-g_cov(:,1,1)*g_cov(:,2,3)**2-g_cov(:,1,2)**2*g_cov(:,3,3)+2*g_cov(:,1,2)*g_cov(:,1,3)*g_cov(:,2,3)-g_cov(:,1,3)**2*g_cov(:,2,2))  
  idet(:) = 1./det(:)  
  !
  g_contr(:,1,1) =  ( g_cov(:,2,2)*g_cov(:,3,3)-g_cov(:,2,3)*g_cov(:,3,2))*idet(:) 
  g_contr(:,1,2) = -( g_cov(:,1,2)*g_cov(:,3,3)-g_cov(:,1,3)*g_cov(:,3,2))*idet(:)
  g_contr(:,1,3) = -(-g_cov(:,1,2)*g_cov(:,2,3)+g_cov(:,1,3)*g_cov(:,2,2))*idet(:) 
  g_contr(:,2,1) = -( g_cov(:,2,1)*g_cov(:,3,3)-g_cov(:,2,3)*g_cov(:,3,1))*idet(:) 
  g_contr(:,2,2) =  ( g_cov(:,1,1)*g_cov(:,3,3)-g_cov(:,1,3)*g_cov(:,3,1))*idet(:) 
  g_contr(:,2,3) = -( g_cov(:,1,1)*g_cov(:,2,3)-g_cov(:,1,3)*g_cov(:,2,1))*idet(:) 
  g_contr(:,3,1) = -(-g_cov(:,2,1)*g_cov(:,3,2)+g_cov(:,2,2)*g_cov(:,3,1))*idet(:) 
  g_contr(:,3,2) = -( g_cov(:,1,1)*g_cov(:,3,2)-g_cov(:,1,2)*g_cov(:,3,1))*idet(:) 
  g_contr(:,3,3) =  ( g_cov(:,1,1)*g_cov(:,2,2)-g_cov(:,1,2)*g_cov(:,2,1))*idet(:)       
  !
  gp(:) = SQRT(det)
  gm(:) = 1./gp(:) 
  !
  DD(:)         = Q(:,1)*gm(:) 
  S_cov(:,1)    = Q(:,2)*gm(:)
  S_cov(:,2)    = Q(:,3)*gm(:)
  S_cov(:,3)    = Q(:,4)*gm(:)
  tau(:)        = Q(:,5)*gm(:)
  !
  !v_contr = MATMUL(g_contr,v_cov)
  !S_contr = MATMUL(g_contr,S_cov)
  !B_cov = MATMUL(g_cov,B_contr)
  ! 
  v_contr(:,1) = g_contr(:,1,1)*v_cov(:,1) + g_contr(:,1,2)*v_cov(:,2) + g_contr(:,1,3)*v_cov(:,3) 
  v_contr(:,2) = g_contr(:,2,1)*v_cov(:,1) + g_contr(:,2,2)*v_cov(:,2) + g_contr(:,2,3)*v_cov(:,3) 
  v_contr(:,3) = g_contr(:,3,1)*v_cov(:,1) + g_contr(:,3,2)*v_cov(:,2) + g_contr(:,3,3)*v_cov(:,3) 
  gS_contr(:,1) = g_contr(:,1,1)*Q(:,2) + g_contr(:,1,2)*Q(:,3) + g_contr(:,1,3)*Q(:,4) 
  gS_contr(:,2) = g_contr(:,2,1)*Q(:,2) + g_contr(:,2,2)*Q(:,3) + g_contr(:,2,3)*Q(:,4) 
  gS_contr(:,3) = g_contr(:,3,1)*Q(:,2) + g_contr(:,3,2)*Q(:,3) + g_contr(:,3,3)*Q(:,4) 
  S_contr(:,1) = g_contr(:,1,1)*S_cov(:,1) + g_contr(:,1,2)*S_cov(:,2) + g_contr(:,1,3)*S_cov(:,3) 
  S_contr(:,2) = g_contr(:,2,1)*S_cov(:,1) + g_contr(:,2,2)*S_cov(:,2) + g_contr(:,2,3)*S_cov(:,3) 
  S_contr(:,3) = g_contr(:,3,1)*S_cov(:,1) + g_contr(:,3,2)*S_cov(:,2) + g_contr(:,3,3)*S_cov(:,3) 
  !   
  ! COVARIANT =>
  vxB_cov(:,1) = gp(:)*( v_contr(:,2)*B_contr(:,3) - v_contr(:,3)*B_contr(:,2) )
  vxB_cov(:,2) = gp(:)*( v_contr(:,3)*B_contr(:,1) - v_contr(:,1)*B_contr(:,3) )
  vxB_cov(:,3) = gp(:)*( v_contr(:,1)*B_contr(:,2) - v_contr(:,2)*B_contr(:,1) )
  !
  ! CONTRAVARIANT =>
  vxB_contr(:,1) = gm(:)*( v_cov(:,2)*B_cov(:,3) - v_cov(:,3)*B_cov(:,2) ) 
  vxB_contr(:,2) = gm(:)*( v_cov(:,3)*B_cov(:,1) - v_cov(:,1)*B_cov(:,3) ) 
  vxB_contr(:,3) = gm(:)*( v_cov(:,1)*B_cov(:,2) - v_cov(:,2)*B_cov(:,1) ) 
  !
  v2(:) = v_contr(:,1)*v_cov(:,1) + v_contr(:,2)*v_cov(:,2) + v_contr(:,3)*v_cov(:,3)
  e2(:) = vxB_contr(:,1)*vxB_cov(:,1) + vxB_contr(:,2)*vxB_cov(:,2) + vxB_contr(:,3)*vxB_cov(:,3)
  b2(:) = B_contr(:,1)*B_cov(:,1) + B_contr(:,2)*B_cov(:,2) + B_contr(:,3)*B_cov(:,3)
  s2(:) = S_contr(:,1)*S_cov(:,1) + S_contr(:,2)*S_cov(:,2) + S_contr(:,3)*S_cov(:,3)
  !
  uem(:) = 0.5*(b2(:)+e2(:))
  lf(:)  = 1.0/SQRT(1.0-v2(:))
  ! 
#ifdef GEOS  
  CALL EOS_eps_vector(rhoeps(:),rho(:),p(:),DD(:),tau(:),s2(:))
  w(:) = rho(:) + rhoeps(:) + p(:)          ! rho*enthalpy
#else
  w(:)   = rho(:) + gamma1(:)*p(:)
#endif
!
  !IF(SQRT(SUM((w-(rho + gamma1*p))**2)).GT.1e-10) THEN
  !    continue
  !ENDIF
  !
  ! 
  ww(:)  = w(:)*lf(:)**2
  wwx(:) = ww(:)*v_contr(:,1)
  wwy(:) = ww(:)*v_contr(:,2)
  wwz(:) = ww(:)*v_contr(:,3)
  !
  ! transport velocity
  DO i=1,3
	  Vtr(:,i) = lapse(:)*v_contr(:,i)-shift(:,i)
  ENDDO
  !
  !    Fij(1,1:3) =  (/ f1, g1, h1 /) for Q(6)  
  !    Fij(2,1:3) =  (/ f2, g2, h2 /) for Q(7)  ... without divergence cleaning
  !    Fij(3,1:3) =  (/ f3, g3, h3 /) for Q(8)  
  !
  DO m=1,3
      DO i=1,3
          Fij(:,i,m) = Vtr(:,i)*QB_contr(:,m)-Vtr(:,m)*QB_contr(:,i)  ! Fij(i,i) = 0 !!!! ! NB: this is contravariant !!!!
      ENDDO
  ENDDO
  !
  F(:,1,1)   = v_contr(:,1)*Q(:,1) 
  F(:,2,1)   = wwx(:)*v_cov(:,1) - vxB_contr(:,1)*vxB_cov(:,1) - B_contr(:,1)*B_cov(:,1) + p(:) + uem(:) 
  F(:,3,1)   = wwx(:)*v_cov(:,2) - vxB_contr(:,1)*vxB_cov(:,2) - B_contr(:,1)*B_cov(:,2) 
  F(:,4,1)   = wwx(:)*v_cov(:,3) - vxB_contr(:,1)*vxB_cov(:,3) - B_contr(:,1)*B_cov(:,3) 
  F(:,5,1)   = gS_contr(:,1)-F(:,1,1) 
  ! 
  F(:,6,1)   = Fij(:,1,1) !+ V(9)  
  F(:,7,1)   = Fij(:,1,2) 
  F(:,8,1)   = Fij(:,1,3)  
  !  
#ifdef COVCLEAN     
  F(:,9,1)   = DivCleaning_a**2*QB_contr(:,1)   ! flat Newtonian, flat coordinates
#else
  F(:,9,1)   = 0.  
#endif 
  !  lapse&shift&metric fluxes 
  DO i=10,nVar
	F(:,i,1) = 0.
  ENDDO
  !
  !
  F(:,1,2)   = v_contr(:,2)*Q(:,1) 
  F(:,2,2)   = wwy(:)*v_cov(:,1) - vxB_contr(:,2)*vxB_cov(:,1) - B_contr(:,2)*B_cov(:,1) 
  F(:,3,2)   = wwy(:)*v_cov(:,2) - vxB_contr(:,2)*vxB_cov(:,2) - B_contr(:,2)*B_cov(:,2) + p(:) + uem(:) 
  F(:,4,2)   = wwy(:)*v_cov(:,3) - vxB_contr(:,2)*vxB_cov(:,3) - B_contr(:,2)*B_cov(:,3) 
  F(:,5,2)   = gS_contr(:,2)-F(:,1,2)  
  ! 
  F(:,6,2)   = Fij(:,2,1)  
  F(:,7,2)   = Fij(:,2,2) !+ V(9) 
  F(:,8,2)   = Fij(:,2,3)  
  ! 
#ifdef COVCLEAN     
  F(:,9,2)   = DivCleaning_a**2*QB_contr(:,2)   
#else
  F(:,9,2)   = 0. !QB_contr(2) 
#endif
  !  lapse&shift&metric fluxes 
  DO i=10,nVar
	F(:,i,2) = 0.
  ENDDO
  !
  !
  !
  F(:,1,3)   = v_contr(:,3)*Q(:,1) 
  F(:,2,3)   = wwz(:)*v_cov(:,1) - vxB_contr(:,3)*vxB_cov(:,1) - B_contr(:,3)*B_cov(:,1) 
  F(:,3,3)   = wwz(:)*v_cov(:,2) - vxB_contr(:,3)*vxB_cov(:,2) - B_contr(:,3)*B_cov(:,2)   
  F(:,4,3)   = wwz(:)*v_cov(:,3) - vxB_contr(:,3)*vxB_cov(:,3) - B_contr(:,3)*B_cov(:,3) + p(:) + uem(:) 
  F(:,5,3)   = gS_contr(:,3)-F(:,1,3) 
  ! 
  F(:,6,3)   = Fij(:,3,1)
  F(:,7,3)   = Fij(:,3,2)
  F(:,8,3)   = Fij(:,3,3) !+ V(9) 
  !
#ifdef COVCLEAN     
  F(:,9,3)   = DivCleaning_a**2*QB_contr(:,3)   ! flat Newtonian, flat coordinates
#else
  F(:,9,3)   = 0. !QB_contr(3) 
#endif 
  !  lapse&shift&metric fluxes 
  DO i=10,nVar
	F(:,i,3) = 0.
  ENDDO
  ! 
  ! - - - - - - - - - 
  DO i=2,4
	  F(:,i,1)   = F(:,i,1)*gp(:) 
	  F(:,i,2)   = F(:,i,2)*gp(:) 
	  F(:,i,3)   = F(:,i,3)*gp(:) 
  ENDDO
  ! Remember that Q(:) below contains already the factor gp, which is ok!
  
  DO i=1,5
	  F(:,i,1)   = lapse(:)*F(:,i,1) - shift(:,1)*Q(:,i)
	  F(:,i,2)   = lapse(:)*F(:,i,2) - shift(:,2)*Q(:,i)
	  F(:,i,3)   = lapse(:)*F(:,i,3) - shift(:,3)*Q(:,i)
  ENDDO 
  !
#ifdef COVCLEAN   
  F(:,9,1)   = lapse(:)*F(:,9,1) - shift(:,1)*Q(:,9)
  F(:,9,2)   = lapse(:)*F(:,9,1) - shift(:,2)*Q(:,9)
  F(:,9,3)   = lapse(:)*F(:,9,1) - shift(:,3)*Q(:,9)
#endif 
  !
  CONTINUE    
    !  
  !
END SUBROUTINE PDEFluxPrimVectorGRMHD           
    
!
! Nonconservative part of the PDE ( B(Q) * gradQ ) 
!    
RECURSIVE SUBROUTINE PDENCPGRMHD(BgradQ,Q,gradQ)
    USE Parameters, ONLY : d, EQN
    USE EOS_mod
    USE AstroMod
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
    REAL, INTENT(OUT) :: BgradQ(nVar) 
    REAL, INTENT(IN)  :: Q(nVar), gradQ(nVar,d)
    ! Local variables 
    REAL :: p, irho
  REAL :: rhoeps,dd,tau,s2,sm_cov(3),sm(3)
    REAL :: Qx(nVar), Qy(nVar), Qz(nVar) 
    REAL :: e,g_cov(3,3),g_contr(3,3)
    REAL :: AQx(nVar), BQy(nVar), CQz(nVar) 
  REAL :: lapse, shift(3), gammaij(6), delta(3,3), B_cov(3), vxB_cov(3), vxb_contr(3), psi, S_contr(3), qb_contr(3), B_contr(3) 
  REAL :: v2,v_contr(3),uem,b2,e2,gp,gm,vc(nVar),lf,w,ww,gamma1,rho,v_cov(3), w_ij, wim    
    INTEGER :: i,j,k,l,m,iErr, ccount    
    !
    !BgradQ = 0.0 
	DO i=1,nVar
        BgradQ(i) = 0.
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
	CALL PDECons2PrimGRMHD(Vc,Q,iErr)
	gamma1 = EQN%gamma/(EQN%gamma-1.0)
	rho    = Vc(1)
	v_cov = Vc(2:4)
	p      = Vc(5)
	!
	!B_cov(1:3) = Vc(6:8)
	DO i=1,3
		B_contr(i) = Vc(5+i) ! contravariant!
	  QB_contr(i) = Q(5+i)
	ENDDO
	!
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
	!
	dd = gm*Q(1)
	tau = gm*Q(5)
	DO i=1,3 
		sm_cov(i)  = Q(1+i)
	ENDDO
	sm_cov = gm*sm_cov
	sm   = MATMUL (g_contr, sm_cov)  
	s2   = sm_cov(1)*sm(1) + sm_cov(2)*sm(2) + sm_cov(3)*sm(3)
	!
#ifdef GEOS  
	CALL EOS_eps(rhoeps,rho,p,dd,tau,s2) 
	w      = rho + rhoeps + p   ! rho*enthalpy
#else
	w   = rho + gamma1*p
#endif
	!
	!IF(ABS(w-(rho + gamma1*p)).GT.1e-13) THEN
	!    continue
	!ENDIF
	!
	!w      = rho + gamma1*p   ! rho*hentalpy
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
    AQx(9) = AQx(9) + gm*lapse*EQN%DivCleaning_a**2*Qx(5+j)             ! D.C.   shift * d(gp*B^x)/dx 
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
    BQy(9) = BQy(9) + gm*lapse*EQN%DivCleaning_a**2*Qy(5+j)             ! D.C.   shift * d(gp*B^x)/dx 
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
    CQz(9) = CQz(9) + gm*lapse*EQN%DivCleaning_a**2*Qz(5+j)             ! D.C.   shift * d(gp*B^x)/dx 
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
END SUBROUTINE PDENCPGRMHD  
!
!
!    
RECURSIVE SUBROUTINE PDENCPPrimGRMHD(BgradQ,Vc,Q,gradQ) 
    USE Parameters, ONLY : d, EQN
    USE AstroMod
    USE EOS_mod
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
	REAL, INTENT(IN)  :: Q(nVar), Vc(nVar), gradQ(nVar,d) 
	REAL, INTENT(OUT) :: BgradQ(nVar) 
	! Local variables 
	REAL :: p, irho
	REAL :: rhoeps,dd,tau,s2,sm_cov(3),sm(3)
	REAL :: Qx(nVar), Qy(nVar), Qz(nVar) 
	REAL :: e,g_cov(3,3),g_contr(3,3)
	REAL :: AQx(nVar), BQy(nVar), CQz(nVar) 
	REAL :: lapse, shift(3), gammaij(6), delta(3,3), B_cov(3), vxB_cov(3), vxb_contr(3), psi, S_contr(3), qb_contr(3), B_contr(3) 
	REAL :: v2,v_contr(3),uem,b2,e2,gp,gm,lf,w,ww,gamma1,rho,v_cov(3), w_ij, wim    
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
	gamma1 = EQN%gamma/(EQN%gamma-1.0)
	rho    = Vc(1)
	v_cov = Vc(2:4)
	p      = Vc(5)
	!
	!B_cov(1:3) = Vc(6:8)
	DO i=1,3
		B_contr(i) = Vc(5+i) ! contravariant!
		QB_contr(i) = Q(5+i)
	ENDDO
	!
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
	!
	dd = gm*Q(1)
	tau = gm*Q(5)
	DO i=1,3 
		sm_cov(i)  = Q(1+i)
	ENDDO
	sm_cov = gm*sm_cov
	sm   = MATMUL (g_contr, sm_cov)  
	s2   = sm_cov(1)*sm(1) + sm_cov(2)*sm(2) + sm_cov(3)*sm(3)
	!
#ifdef GEOS  
	CALL EOS_eps(rhoeps,rho,p,dd,tau,s2) 
	w      = rho + rhoeps + p   ! rho*enthalpy
#else
	w   = rho + gamma1*p
#endif
	!
	!IF(ABS(w-(rho + gamma1*p)).GT.1e-13) THEN
	!    continue
	!ENDIF
	!
	!w      = rho + gamma1*p   ! rho*hentalpy
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
    AQx(9) = AQx(9) + gm*lapse*EQN%DivCleaning_a**2*Qx(5+j)             ! D.C.   shift * d(gp*B^x)/dx 
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
    BQy(9) = BQy(9) + gm*lapse*EQN%DivCleaning_a**2*Qy(5+j)             ! D.C.   shift * d(gp*B^x)/dx 
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
    CQz(9) = CQz(9) + gm*lapse*EQN%DivCleaning_a**2*Qz(5+j)             ! D.C.   shift * d(gp*B^x)/dx 
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
END SUBROUTINE PDENCPPrimGRMHD      
!
RECURSIVE SUBROUTINE PDENCPPrimVectorGRMHD(BgradQ,Vc,Q,Qx,Qy,Qz) 
    USE Parameters, ONLY : d, VECTORLENGTH,EQN
    USE AstroMod
    USE EOS_mod 
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
    REAL(8), INTENT(IN)  :: Vc(VECTORLENGTH,nVar), Q(VECTORLENGTH,nVar), Qx(VECTORLENGTH,nVar), Qy(VECTORLENGTH,nVar), Qz(VECTORLENGTH,nVar) 
    REAL(8), INTENT(OUT) :: BgradQ(VECTORLENGTH,nVar) 
    ! Local variables 
    INTEGER :: i,j,k,l,m,n,ip,iq,ii,jj,kk,ll,iErr,count,ccount  
    REAL(8) :: rhoeps(VECTORLENGTH)
    REAL(8) :: p(VECTORLENGTH), irho(VECTORLENGTH), det(VECTORLENGTH), idet(VECTORLENGTH) , DD(VECTORLENGTH) , tau(VECTORLENGTH) , s2(VECTORLENGTH) 
    REAL(8) :: e(VECTORLENGTH),g_cov(VECTORLENGTH,3,3),g_contr(VECTORLENGTH,3,3)
    REAL(8) :: AQx(VECTORLENGTH,nVar), BQy(VECTORLENGTH,nVar), CQz(VECTORLENGTH,nVar) 
    REAL(8) :: lapse(VECTORLENGTH), shift(VECTORLENGTH,3), gammaij(VECTORLENGTH,6), delta(VECTORLENGTH,3,3), B_cov(VECTORLENGTH,3), vxB_cov(VECTORLENGTH,3), vxb_contr(VECTORLENGTH,3), psi(VECTORLENGTH), gS_contr(VECTORLENGTH,3),S_contr(VECTORLENGTH,3), qb_contr(VECTORLENGTH,3), B_contr(VECTORLENGTH,3) 
    REAL(8) :: v2(VECTORLENGTH),v_contr(VECTORLENGTH,3),S_cov(VECTORLENGTH,3),uem(VECTORLENGTH),b2(VECTORLENGTH),e2(VECTORLENGTH),gp(VECTORLENGTH),gm(VECTORLENGTH),lf(VECTORLENGTH),w(VECTORLENGTH),ww(VECTORLENGTH),gamma1(VECTORLENGTH),rho(VECTORLENGTH),v_cov(VECTORLENGTH,3), w_ij(VECTORLENGTH), wim(VECTORLENGTH)    
    !
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN:  64 :: rhoeps,p,irho,det,idet,DD,tau,s2,e,g_cov,g_contr,AQx,BQy,CQz,lapse,shift,gammaij,delta,B_cov,vxB_cov,vxb_contr,psi,S_contr,gS_contr,qb_contr,B_contr
    !DIR$ ATTRIBUTES ALIGN:  64 :: v2,v_contr,S_cov,uem,b2,e2,gp,gm,lf,w,ww,gamma1,rho,v_cov,w_ij,wim
    !DIR$ ASSUME_ALIGNED BgradQ : 64  
    !DIR$ ASSUME_ALIGNED Vc : 64  
    !DIR$ ASSUME_ALIGNED Q  : 64  
    !DIR$ ASSUME_ALIGNED Qx : 64  
    !DIR$ ASSUME_ALIGNED Qy : 64  
    !DIR$ ASSUME_ALIGNED Qz : 64   
#else
    !DIR$ ATTRIBUTES ALIGN:  32 :: rhoeps,p,irho,det,idet,DD,tau,s2,e,g_cov,g_contr,AQx,BQy,CQz,lapse,shift,gammaij,delta,B_cov,vxB_cov,vxb_contr,psi,S_contr,gS_contr,qb_contr,B_contr
    !DIR$ ATTRIBUTES ALIGN:  32 :: v2,v_contr,S_cov,uem,b2,e2,gp,gm,lf,w,ww,gamma1,rho,v_cov,w_ij,wim
    !DIR$ ASSUME_ALIGNED BgradQ : 32 
    !DIR$ ASSUME_ALIGNED Vc : 32 
    !DIR$ ASSUME_ALIGNED Q  : 32  
    !DIR$ ASSUME_ALIGNED Qx : 32  
    !DIR$ ASSUME_ALIGNED Qy : 32  
    !DIR$ ASSUME_ALIGNED Qz : 32   
#endif 
    !
    !BgradQ = 0.0 
    AQx = 0.0
    BQy = 0.0
    CQz = 0.0 
    !
    lapse(:) = Q(:,10)
    !
	shift(:,1) = Q(:,11)
	shift(:,2) = Q(:,12)
	shift(:,3) = Q(:,13)
	!
	DO i=1,6
		gammaij(:,i) = Q(:,13+i) 
	ENDDO
    g_cov(:,1,1) = Q(:,14)
    g_cov(:,1,2) = Q(:,15)
    g_cov(:,1,3) = Q(:,16)
    g_cov(:,2,2) = Q(:,17)
    g_cov(:,2,3) = Q(:,18)
    g_cov(:,3,3) = Q(:,19)
    g_cov(:,2,1) = Q(:,15)
    g_cov(:,3,1) = Q(:,16)
    g_cov(:,3,2) = Q(:,18)
    !
    delta = 0.0
    DO i=1,3
        delta(:,i,i) = 1.0 
    ENDDO 
    !
    det(:) = (g_cov(:,1,1)*g_cov(:,2,2)*g_cov(:,3,3)-g_cov(:,1,1)*g_cov(:,2,3)**2-g_cov(:,1,2)**2*g_cov(:,3,3)+2*g_cov(:,1,2)*g_cov(:,1,3)*g_cov(:,2,3)-g_cov(:,1,3)**2*g_cov(:,2,2))  
    idet(:) = 1./det(:)  
    !
    g_contr(:,1,1) =  ( g_cov(:,2,2)*g_cov(:,3,3)-g_cov(:,2,3)*g_cov(:,3,2)) * idet(:) 
    g_contr(:,1,2) = -( g_cov(:,1,2)*g_cov(:,3,3)-g_cov(:,1,3)*g_cov(:,3,2)) * idet(:)
    g_contr(:,1,3) = -(-g_cov(:,1,2)*g_cov(:,2,3)+g_cov(:,1,3)*g_cov(:,2,2)) * idet(:) 
    g_contr(:,2,1) = -( g_cov(:,2,1)*g_cov(:,3,3)-g_cov(:,2,3)*g_cov(:,3,1)) * idet(:) 
    g_contr(:,2,2) =  ( g_cov(:,1,1)*g_cov(:,3,3)-g_cov(:,1,3)*g_cov(:,3,1)) * idet(:) 
    g_contr(:,2,3) = -( g_cov(:,1,1)*g_cov(:,2,3)-g_cov(:,1,3)*g_cov(:,2,1)) * idet(:) 
    g_contr(:,3,1) = -(-g_cov(:,2,1)*g_cov(:,3,2)+g_cov(:,2,2)*g_cov(:,3,1)) * idet(:) 
    g_contr(:,3,2) = -( g_cov(:,1,1)*g_cov(:,3,2)-g_cov(:,1,2)*g_cov(:,3,1)) * idet(:) 
    g_contr(:,3,3) =  ( g_cov(:,1,1)*g_cov(:,2,2)-g_cov(:,1,2)*g_cov(:,2,1)) * idet(:)       

    gp(:) = SQRT(det(:))
    gm(:) = 1./gp(:) 
    !  
    DD(:)         = Q(:,1)*gm(:) 
    S_cov(:,1)    = Q(:,2)*gm(:)
    S_cov(:,2)    = Q(:,3)*gm(:)
    S_cov(:,3)    = Q(:,4)*gm(:)
    tau(:)        = Q(:,5)*gm(:)
    !
    gamma1 = EQN%gamma/(EQN%gamma-1.0)
    rho(:) = Vc(:,1)
    v_cov(:,:) = Vc(:,2:4)
    p(:)       = Vc(:,5)
    !
    DO i=1,3
	    B_contr(:,i) = Vc(:,5+i) ! contravariant!
	    QB_contr(:,i) = Q(:,5+i)
    ENDDO
    !
    !B_cov = MATMUL(g_cov,B_contr)
    B_cov(:,1) = g_cov(:,1,1)*B_contr(:,1) + g_cov(:,1,2)*B_contr(:,2) + g_cov(:,1,3)*B_contr(:,3) 
    B_cov(:,2) = g_cov(:,2,1)*B_contr(:,1) + g_cov(:,2,2)*B_contr(:,2) + g_cov(:,2,3)*B_contr(:,3) 
    B_cov(:,3) = g_cov(:,3,1)*B_contr(:,1) + g_cov(:,3,2)*B_contr(:,2) + g_cov(:,3,3)*B_contr(:,3) 
    !
    psi(:) = Vc(:,9) 
    ! 
    !v_contr       = MATMUL(g_contr,v_cov)
    !gS_contr = MATMUL(g_contr,Q(2:4))
    !
    v_contr(:,1) = g_contr(:,1,1)*v_cov(:,1) + g_contr(:,1,2)*v_cov(:,2) + g_contr(:,1,3)*v_cov(:,3) 
    v_contr(:,2) = g_contr(:,2,1)*v_cov(:,1) + g_contr(:,2,2)*v_cov(:,2) + g_contr(:,2,3)*v_cov(:,3) 
    v_contr(:,3) = g_contr(:,3,1)*v_cov(:,1) + g_contr(:,3,2)*v_cov(:,2) + g_contr(:,3,3)*v_cov(:,3) 
    S_contr(:,1) = g_contr(:,1,1)*S_cov(:,1) + g_contr(:,1,2)*S_cov(:,2) + g_contr(:,1,3)*S_cov(:,3) 
    S_contr(:,2) = g_contr(:,2,1)*S_cov(:,1) + g_contr(:,2,2)*S_cov(:,2) + g_contr(:,2,3)*S_cov(:,3) 
    S_contr(:,3) = g_contr(:,3,1)*S_cov(:,1) + g_contr(:,3,2)*S_cov(:,2) + g_contr(:,3,3)*S_cov(:,3)   
    gS_contr(:,1) = g_contr(:,1,1)*Q(:,2) + g_contr(:,1,2)*Q(:,3) + g_contr(:,1,3)*Q(:,4) 
    gS_contr(:,2) = g_contr(:,2,1)*Q(:,2) + g_contr(:,2,2)*Q(:,3) + g_contr(:,2,3)*Q(:,4) 
    gS_contr(:,3) = g_contr(:,3,1)*Q(:,2) + g_contr(:,3,2)*Q(:,3) + g_contr(:,3,3)*Q(:,4)   
    !
    !! COVARIANT =>
    vxB_cov(:,1) = gp(:)*( v_contr(:,2)*B_contr(:,3) - v_contr(:,3)*B_contr(:,2) ) 
    vxB_cov(:,2) = gp(:)*( v_contr(:,3)*B_contr(:,1) - v_contr(:,1)*B_contr(:,3) ) 
    vxB_cov(:,3) = gp(:)*( v_contr(:,1)*B_contr(:,2) - v_contr(:,2)*B_contr(:,1) ) 
    !
    ! CONTRAVARIANT   
    vxB_contr(:,1) = gm(:)* ( v_cov(:,2)*B_cov(:,3) - v_cov(:,3)*B_cov(:,2) ) 
    vxB_contr(:,2) = gm(:)* ( v_cov(:,3)*B_cov(:,1) - v_cov(:,1)*B_cov(:,3) ) 
    vxB_contr(:,3) = gm(:)* ( v_cov(:,1)*B_cov(:,2) - v_cov(:,2)*B_cov(:,1) ) 
    !v2     = vx**2 + vy**2 + vz**2
    !b2     = bx**2 + by**2 + bz**2
    !e2     = ex**2 + ey**2 + ez**2 
    v2     = v_contr(:,1)*v_cov(:,1) + v_contr(:,2)*v_cov(:,2) + v_contr(:,3)*v_cov(:,3)
    e2     = vxB_contr(:,1)*vxB_cov(:,1) + vxB_contr(:,2)*vxB_cov(:,2) + vxB_contr(:,3)*vxB_cov(:,3)
    b2     = B_contr(:,1)*B_cov(:,1) + B_contr(:,2)*B_cov(:,2) + B_contr(:,3)*B_cov(:,3)
    s2     = S_contr(:,1)*S_cov(:,1) + S_contr(:,2)*S_cov(:,2) + S_contr(:,3)*S_cov(:,3)
    !
    uem(:)    = 0.5*(b2(:) + e2(:) ) 
    !
    lf(:)     = 1.0/sqrt(1.0 - v2(:) ) 
    !
#ifdef GEOS  
    CALL EOS_eps_vector(rhoeps(:),rho(:),p(:),DD(:),tau(:),s2(:))
    w(:)      = rho(:) + rhoeps(:) + p(:)   ! rho*enthalpy
#else
    w(:)     = rho(:) + gamma1(:)*p(:)   ! rho*hentalpy
#endif
    !IF(SQRT(SUM((w-(rho + gamma1*p))**2)).GT.1e-10) THEN
    !    continue
    !ENDIF
    !
    !w(:)     = rho(:) + gamma1(:)*p(:)   ! rho*hentalpy
    ww(:)     = w(:)*lf(:)**2          ! rho*hentalpy*Lorentz^2 
    !
	ccount=0
    DO i=1,3 
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=1
                !------ 
			W_ij(:) = ww(:)*v_cov(:,i)*v_contr(:,j)-vxB_cov(:,i)*vxB_contr(:,j)-B_cov(:,i)*B_contr(:,j)+(p(:)+uem(:))*delta(:,i,j)
                !
                AQx(:,1+j) = AQx(:,1+j) - Q(:,1+i)*Qx(:,10+i)          ! Q(11:13)  shift(i) or shift_contr(i)
                AQx(:,5)   = AQx(:,5)   - gp(:)*W_ij(:)*Qx(:,10+i)     ! Q(11:13)  shift(i) or shift_contr(i)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=2
                !------ 
			W_ij(:) = ww(:)*v_cov(:,i)*v_contr(:,j)-vxB_cov(:,i)*vxB_contr(:,j)-B_cov(:,i)*B_contr(:,j)+(p(:)+uem(:))*delta(:,i,j)
                !
                BQy(:,1+j) = BQy(:,1+j) - Q(:,1+i)*Qy(:,10+i)   ! Q(11:13)  shift(i) or shift_contr(i)
                BQy(:,5) = BQy(:,5) - gp(:)*W_ij(:)*Qy(:,10+i)     ! Q(11:13)  shift(i) or shift_contr(i)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=3
                !------
			W_ij(:) = ww(:)*v_cov(:,i)*v_contr(:,j)-vxB_cov(:,i)*vxB_contr(:,j)-B_cov(:,i)*B_contr(:,j)+(p(:)+uem(:))*delta(:,i,j)
                !
                CQz(:,1+j) = CQz(:,1+j) - Q(:,1+i)*Qz(:,10+i)   ! Q(11:13)  shift(i) or shift_contr(i)
                CQz(:,5) = CQz(:,5) - gp(:)*W_ij(:)*Qz(:,10+i)     ! Q(11:13)  shift(i) or shift_contr(i)
                !
        DO m=1,3
            IF(m.GE.i) THEN  
			ccount=ccount+1
                !
			Wim(:) = ww(:)*v_contr(:,i)*v_contr(:,m)-vxB_contr(:,i)*vxB_contr(:,m)-B_contr(:,i)*B_contr(:,m)+(p(:)+uem(:))*g_contr(:,i,m)
			!Wim = Wim + (1.0 - delta(i,m))*(ww*v_contr(m)*v_contr(i)-vxB_contr(m)*vxB_contr(i)-B_contr(m)*B_contr(i)+(p+uem)*g_contr(m,i))  ! account also of the remaining symmetric components of gamma for i.NE.m.
			IF(i.NE.m) THEN !TWICE! account also of the remaining symmetric components of gamma §(memo: g_contr and Wim are symmetric)
                    Wim(:) = 2.0*Wim(:)    
                ENDIF
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=1
                !------ 
			AQx(:,1+j) = AQx(:,1+j) - 0.5*gp(:)*lapse(:)*Wim(:)*Qx(:,13+ccount)  ! Q(14:19) gammaij(ccount) or  g_cov(i,m)
			AQx(:,5) = AQx(:,5) - 0.5*gp(:)*Wim(:)*shift(:,j)*Qx(:,13+ccount)   ! Q(14:19) gammaij(ccount) or  g_cov(i,m)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=2
                !------ 
			BQy(:,1+j) = BQy(:,1+j) - 0.5*gp(:)*lapse(:)*Wim(:)*Qy(:,13+ccount)  ! Q(14:19) gammaij(ccount) or  g_cov(i,m)
			BQy(:,5)  = BQy(:,5) - 0.5*gp(:)*Wim(:)*shift(:,j)*Qy(:,13+ccount)   ! Q(14:19) gammaij(ccount) or  g_cov(i,m)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=3
                !------
			CQz(:,1+j) = CQz(:,1+j) - 0.5*gp(:)*lapse(:)*Wim(:)*Qz(:,13+ccount)  ! Q(14:19) gammaij(ccount) or  g_cov(i,m)
			CQz(:,5)  = CQz(:,5) - 0.5*gp(:)*Wim(:)*shift(:,j)*Qz(:,13+ccount)   ! Q(14:19) gammaij(ccount) or  g_cov(i,m)
                !
#ifdef COVCLEAN     
				CQz(:,9) = CQz(:,9) + 0.5*Q(:,9)*gammaij(:,ccount)*shift(:,j)*Qz(:,13+ccount)   ! D.C.
                CQz(:,6:8) = CQz(:,6:8) - 0.5*lapse(:)*Q(:,9)*g_contr(:,1:3,j)*gammaij(:,ccount)*Qz(:,13+ccount)
                !
                IF(i.NE.m) THEN !TWICE!  account  also of the remaining symmetric components of gamma §(memo: g_contr and Wim are symmetric)
                    CQz(:,9) = CQz(:,9) + 0.5*Q(:,9)*gammaij(:,ccount)*shift(:,j)*Qz(:,13+ccount)   ! D.C.
                    CQz(:,6:8) = CQz(:,6:8) - 0.5*lapse(:)*Q(:,9)*g_contr(:,1:3,j)*gammaij(:,ccount)*Qz(:,13+ccount)
                ENDIF
#endif  		
                !
            ENDIF
        ENDDO
    ENDDO
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=1
    !------ 
    AQx(:,1+j) = AQx(:,1+j) + (Q(:,5)+Q(:,1))*Qx(:,10)    ! Q(10) or lapse
    AQx(:,5) = AQx(:,5) + gS_contr(:,j)*Qx(:,10)         !  Q(10) or lapse
    !
#ifdef COVCLEAN     
	AQx(:,9) = AQx(:,9) - Q(:,5+j)*Qx(:,10)             ! D.C.
    AQx(:,9) = AQx(:,9) + Q(:,9)*Qx(:,10+j)             ! D.C.   phi * d(shift^x)/dx
    AQx(:,6:8) = AQx(:,6:8) - Q(:,11:13)*Qx(:,5+j)             ! D.C.   shift * d(gp*B^x)/dx
    AQx(:,6:8) = AQx(:,6:8) + lapse(:)*g_contr(:,1:3,j)*Qx(:,9)  ! D.C.  
#else
    AQx(:,6) = AQx(:,6)  - Q(:,11)*Qx(:,5+j) + lapse(:)*gp(:)*g_contr(:,1,j)*Qx(:,9)              ! D.C.   shift * d(gp*B^x)/dx 
    AQx(:,7) = AQx(:,7)  - Q(:,12)*Qx(:,5+j) + lapse(:)*gp(:)*g_contr(:,2,j)*Qx(:,9)              ! D.C.   shift * d(gp*B^x)/dx 
    AQx(:,8) = AQx(:,8)  - Q(:,13)*Qx(:,5+j) + lapse(:)*gp(:)*g_contr(:,3,j)*Qx(:,9)              ! D.C.   shift * d(gp*B^x)/dx 
    AQx(:,9) = AQx(:,9)  + gm(:)*lapse(:)*EQN%DivCleaning_a**2*Qx(:,5+j) - Q(:,10+j)*Qx(:,9)             ! D.C.   shift * d(gp*B^x)/dx 
#endif
	!
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=2
    !------ 
    BQy(:,1+j) = BQy(:,1+j) + (Q(:,5)+Q(:,1))*Qy(:,10)    ! Q(10) or lapse
    BQy(:,5)   = BQy(:,5) + gS_contr(:,j)*Qy(:,10)         !  Q(10) or lapse
    !
#ifdef COVCLEAN     
	BQy(:,9) = BQy(:,9) - Q(:,5+j)*Qy(:,10)             ! D.C.
    BQy(:,9) = BQy(:,9) + Q(:,9)*Qy(:,10+j)             ! D.C.   phi * d(shift^y)/dy
    BQy(:,6:8) = BQy(:,6:8) - Q(:,11:13)*Qy(:,5+j)                  ! D.C.   shift * d(gp*B^y)/dy
    BQy(:,6:8) = BQy(:,6:8) + lapse(:)*g_contr(:,1:3,j)*Qy(:,9)        ! D.C.  
#else
    BQy(:,6) = BQy(:,6) - Q(:,11)*Qy(:,5+j)  + lapse(:)*gp(:)*g_contr(:,1,j)*Qy(:,9)                 ! D.C.   shift * d(gp*B^x)/dx 
    BQy(:,7) = BQy(:,7) - Q(:,12)*Qy(:,5+j)  + lapse(:)*gp(:)*g_contr(:,2,j)*Qy(:,9)                 ! D.C.   shift * d(gp*B^x)/dx 
    BQy(:,8) = BQy(:,8) - Q(:,13)*Qy(:,5+j)  + lapse(:)*gp(:)*g_contr(:,3,j)*Qy(:,9)                 ! D.C.   shift * d(gp*B^x)/dx 
    BQy(:,9) = BQy(:,9) + gm(:)*lapse(:)*EQN%DivCleaning_a**2*Qy(:,5+j) - Q(:,10+j)*Qy(:,9)          ! D.C.   shift * d(gp*B^x)/dx 
#endif
	!
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=3
    !------
    CQz(:,1+j) = CQz(:,1+j) + (Q(:,5)+Q(:,1))*Qz(:,10)    ! Q(10) or lapse
    CQz(:,5) = CQz(:,5) + gS_contr(:,j)*Qz(:,10)         !  Q(10) or lapse
    !
#ifdef COVCLEAN     
	CQz(:,9) = CQz(:,9) - Q(:,5+j)*Qz(:,10)             ! D.C.
    CQz(:,9) = CQz(:,9) + Q(:,9)*Qz(:,10+j)             ! D.C.   phi *d(shift^z)/dz
    CQz(:,6:8) = CQz(:,6:8) - Q(:,11:13)*Qz(:,5+j)                  ! D.C.   shift * d(gp*B^z)/dz
    CQz(:,6:8) = CQz(:,6:8) + lapse(:)*g_contr(:,1:3,j)*Qz(:,9)        ! D.C.  
#else
    CQz(:,6) = CQz(:,6) - Q(:,11)*Qz(:,5+j) + lapse(:)*gp(:)*g_contr(:,1,j)*Qz(:,9)                  ! D.C.   shift * d(gp*B^x)/dx  
    CQz(:,7) = CQz(:,7) - Q(:,12)*Qz(:,5+j) + lapse(:)*gp(:)*g_contr(:,2,j)*Qz(:,9)                  ! D.C.   shift * d(gp*B^x)/dx  
    CQz(:,8) = CQz(:,8) - Q(:,13)*Qz(:,5+j) + lapse(:)*gp(:)*g_contr(:,3,j)*Qz(:,9)                  ! D.C.   shift * d(gp*B^x)/dx  
    CQz(:,9) = CQz(:,9) + gm(:)*lapse(:)*EQN%DivCleaning_a**2*Qz(:,5+j)  - Q(:,10+j)*Qz(:,9)         ! D.C.   shift * d(gp*B^x)/dx 
#endif
	!
	DO i=1,nVar
		BgradQ(:,i) = AQx(:,i) + BQy(:,i) + CQz(:,i)  
    ENDDO
    !    
    !WHERE(Vc(:,1) < 1e-7 ) 
    !    BgradQ(:,1) = 0.0
    !    BgradQ(:,2) = 0.0
    !    BgradQ(:,3) = 0.0
    !    BgradQ(:,4) = 0.0
    !    BgradQ(:,5) = 0.0
    !    BgradQ(:,6) = 0.0
    !    BgradQ(:,7) = 0.0
    !    BgradQ(:,8) = 0.0
    !    BgradQ(:,9) = 0.0
    !ENDWHERE    
    !
    continue
    !
END SUBROUTINE PDENCPPrimVectorGRMHD           
!
RECURSIVE SUBROUTINE PDECons2PrimVectorGRMHD(V,Q,iErr)  
    USE Parameters, ONLY : d, EQN, VECTORLENGTH, p_floor, rho_floor ,NSTOV_rho_atmo,NSTOV_p_atmo
    USE AstroMod
    USE EOS_mod
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
    REAL(8), INTENT(INOUT) :: V(VECTORLENGTH,nVar) 
    REAL(8), INTENT(IN)    :: Q(VECTORLENGTH,nVar) 
    ! Local variables 
    INTEGER :: iErr 
    INTEGER :: i,j,k,m
    INTEGER(8) :: iErrV(VECTORLENGTH)
    LOGICAL(8)    :: FAILED(VECTORLENGTH) 
    REAL(8)    :: p0(VECTORLENGTH), tau(VECTORLENGTH),dprho(VECTORLENGTH),dpeps(VECTORLENGTH), Lor(VECTORLENGTH), h(VECTORLENGTH)
    REAL(8)    :: rhoeps(VECTORLENGTH)
    REAL(8)    :: gp(VECTORLENGTH), gm(VECTORLENGTH), g_cov(VECTORLENGTH,3,3), g_contr(VECTORLENGTH,3,3), det(VECTORLENGTH), idet(VECTORLENGTH) 
    REAL(8)    :: psi(VECTORLENGTH), lapse(VECTORLENGTH), shift(VECTORLENGTH,3), Qloc(VECTORLENGTH,nVar), dd(VECTORLENGTH)
    REAL(8)    :: B_contr(VECTORLENGTH,3), B_cov(VECTORLENGTH,3), sm_cov(VECTORLENGTH,3), sm(VECTORLENGTH,3),gamma1(VECTORLENGTH),gam(VECTORLENGTH) 
    REAL(8)    :: s2(VECTORLENGTH), b2(VECTORLENGTH), sb(VECTORLENGTH), sb2(VECTORLENGTH), e(VECTORLENGTH), x1(VECTORLENGTH), x2(VECTORLENGTH), w(VECTORLENGTH)   
    REAL(8)    :: eps(VECTORLENGTH), v2(VECTORLENGTH), p(VECTORLENGTH), rho(VECTORLENGTH), vx(VECTORLENGTH), vy(VECTORLENGTH), vz(VECTORLENGTH), v_cov(VECTORLENGTH,3)  
    REAL(8)    :: den(VECTORLENGTH), vb(VECTORLENGTH), v_contr(VECTORLENGTH,3), ww(VECTORLENGTH), ptest(VECTORLENGTH), tol(VECTORLENGTH), sqrtX(VECTORLENGTH)
    !
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64:: gp,gm,g_cov,g_contr,det,idet,psi,lapse,shift,Qloc,B_contr,B_cov,sm_cov,sm,dd,gamma1,gam,s2,b2,sb,sb2,e
    !DIR$ ATTRIBUTES ALIGN: 64:: x1,x2,w,eps,v2,p,rho,vx,vy,vz,v_cov,v_contr,den,vb,p0,tau,dprho,dpeps,Lor,h,rhoeps,iErrV,FAILED,ptest,tol,sqrtX
    !DIR$ ASSUME_ALIGNED V : 64 
    !DIR$ ASSUME_ALIGNED Q : 64 
#else
    !DIR$ ATTRIBUTES ALIGN: 32:: gp,gm,g_cov,g_contr,det,idet,psi,lapse,shift,Qloc,B_contr,B_cov,sm_cov,sm,dd,gamma1,gam,s2,b2,sb,sb2,e 
    !DIR$ ATTRIBUTES ALIGN: 32:: x1,x2,w,eps,v2,p,rho,vx,vy,vz,v_cov,v_contr,den,vb,p0,tau,dprho,dpeps,Lor,h,rhoeps,iErrV,FAILED,ptest,tol,sqrtX
    !DIR$ ASSUME_ALIGNED V : 32 
    !DIR$ ASSUME_ALIGNED Q : 32 
#endif 
    !
    !RETURN 
    tol = 0.
    !
    psi(:)   = Q(:,9)
    lapse(:) = Q(:,10)
	DO i=1,3
		shift(:,i) = Q(:,10+i)
	ENDDO
	DO i=1,6
		V(:,13+i) = Q(:,13+i)
	ENDDO
    !
    g_cov(:,1,1) = Q(:,14)
    g_cov(:,1,2) = Q(:,15)
    g_cov(:,1,3) = Q(:,16)
    g_cov(:,2,2) = Q(:,17)
    g_cov(:,2,3) = Q(:,18)
    g_cov(:,3,3) = Q(:,19)
    g_cov(:,2,1) = Q(:,15)
    g_cov(:,3,1) = Q(:,16)
    g_cov(:,3,2) = Q(:,18)
    !
    det(:) = (g_cov(:,1,1)*g_cov(:,2,2)*g_cov(:,3,3)-g_cov(:,1,1)*g_cov(:,2,3)**2-g_cov(:,1,2)**2*g_cov(:,3,3)+2*g_cov(:,1,2)*g_cov(:,1,3)*g_cov(:,2,3)-g_cov(:,1,3)**2*g_cov(:,2,2))  
    idet(:) = 1./det(:)  
    !
    g_contr(:,1,1) =  ( g_cov(:,2,2)*g_cov(:,3,3)-g_cov(:,2,3)*g_cov(:,3,2)) * idet(:) 
    g_contr(:,1,2) = -( g_cov(:,1,2)*g_cov(:,3,3)-g_cov(:,1,3)*g_cov(:,3,2)) * idet(:)
    g_contr(:,1,3) = -(-g_cov(:,1,2)*g_cov(:,2,3)+g_cov(:,1,3)*g_cov(:,2,2)) * idet(:) 
    g_contr(:,2,1) = -( g_cov(:,2,1)*g_cov(:,3,3)-g_cov(:,2,3)*g_cov(:,3,1)) * idet(:) 
    g_contr(:,2,2) =  ( g_cov(:,1,1)*g_cov(:,3,3)-g_cov(:,1,3)*g_cov(:,3,1)) * idet(:) 
    g_contr(:,2,3) = -( g_cov(:,1,1)*g_cov(:,2,3)-g_cov(:,1,3)*g_cov(:,2,1)) * idet(:) 
    g_contr(:,3,1) = -(-g_cov(:,2,1)*g_cov(:,3,2)+g_cov(:,2,2)*g_cov(:,3,1)) * idet(:) 
    g_contr(:,3,2) = -( g_cov(:,1,1)*g_cov(:,3,2)-g_cov(:,1,2)*g_cov(:,3,1)) * idet(:) 
    g_contr(:,3,3) =  ( g_cov(:,1,1)*g_cov(:,2,2)-g_cov(:,1,2)*g_cov(:,2,1)) * idet(:)       
    ! 
    gp(:) = SQRT(det)
    gm(:) = 1./gp(:) 
    ! 
	DO i=1,8
		Qloc(:,i)  = gm(:)*Q(:,i)
	ENDDO
    ! 
	dd(:) = Qloc(:,1)
	tau(:) = Qloc(:,5)
	DO i=1,3
		B_contr(:,i) = Qloc(:,5+i)  
		sm_cov(:,i)  = Qloc(:,1+i)
	ENDDO
    !
    gamma1 = EQN%gamma/(EQN%gamma - 1.0)
    gam    = 1.0/gamma1
    ! Solve for p
    FAILED  = .FALSE.
    !
    !sm   = MATMUL (g_contr, sm_cov)
	!B_cov   = MATMUL (g_cov, B_contr)
    ! 
    sm(:,1) = g_contr(:,1,1)*sm_cov(:,1) + g_contr(:,1,2)*sm_cov(:,2) + g_contr(:,1,3)*sm_cov(:,3) 
    sm(:,2) = g_contr(:,2,1)*sm_cov(:,1) + g_contr(:,2,2)*sm_cov(:,2) + g_contr(:,2,3)*sm_cov(:,3) 
    sm(:,3) = g_contr(:,3,1)*sm_cov(:,1) + g_contr(:,3,2)*sm_cov(:,2) + g_contr(:,3,3)*sm_cov(:,3) 
	!
    B_cov(:,1) = g_cov(:,1,1)*B_contr(:,1) + g_cov(:,1,2)*B_contr(:,2) + g_cov(:,1,3)*B_contr(:,3) 
    B_cov(:,2) = g_cov(:,2,1)*B_contr(:,1) + g_cov(:,2,2)*B_contr(:,2) + g_cov(:,2,3)*B_contr(:,3) 
    B_cov(:,3) = g_cov(:,3,1)*B_contr(:,1) + g_cov(:,3,2)*B_contr(:,2) + g_cov(:,3,3)*B_contr(:,3) 
    
    s2(:)  = sm_cov(:,1)*sm(:,1) + sm_cov(:,2)*sm(:,2) + sm_cov(:,3)*sm(:,3)
	b2(:)  = B_contr(:,1)*B_cov(:,1) + B_contr(:,2)*B_cov(:,2) + B_contr(:,3)*B_cov(:,3)
	sb(:)  = sm_cov(:,1)*B_contr(:,1) + sm_cov(:,2)*B_contr(:,2) + sm_cov(:,3)*B_contr(:,3)
    sb2(:) = sb(:)*sb(:)  
	!
    ! First option [Del Zanna et al. (2007) A&A, 473, 11-30 (method 3)]
	e(:)    = Qloc(:,5) + dd(:)      ! Q(5) = gamma^1/2 ( U - dd )
    x1(:)   = 0.         ! 
	eps(:)  = 1.0e-18
	x2(:)   = 1.0-1e-10    !tol ! 
    tol(:)  = 1e-18
    w(:)    = 0 
    !
    ! THIS IS ONLY FOR THE GRHD.
#ifdef GEOS
  IF(ANY(v(:,5).GT.0)) THEN
    p = v(:,5)
  ELSE
    p = NSTOV_p_atmo  ! initial guess 
  ENDIF
  !p0(:) = NSTOV_p_atmo  ! initial guess 
  !CALL C2P_GEOS_vector(p,rho,eps,p0,tau,DD,s2,FAILED)  !  (tau,D,s2) => (p,rho,eps)
  !
  x1(:)   = 1.0e-18         ! 
  x2(:)   = 1.0e7    !tol ! 
  tol(:)  = 1e-18 
  ! 
  CALL RTSAFE_C2P_RHD1_vector(p,X1,X2,tol,dd,tau,s2,FAILED)
  !
  !
  den = tau+p+DD 
  sqrtX = SQRT(den**2-s2)
  !
  rho =DD*sqrtX/den
  eps =sqrtX-p*den/sqrtX-DD
  !
  continue
  !
#ifdef GEOS  
        CALL EOS_eps_VECTOR(rhoeps,rho,p,dd,tau,s2) 
        w      = rho + rhoeps + p   ! rho*enthalpy
#else
        w    = rho + gamma1*p  ! this is rho*h
#endif
  !
  IF(ANY(FAILED)) THEN
    WHERE(FAILED(:)) 
	    ! 
        iErrV(:)    = -1
#ifdef C2PFF
        p(:)    = NSTOV_p_atmo
        rho(:)  = NSTOV_rho_atmo
#else
        p(:)    = p_floor
        rho (:) = rho_floor
#endif 
        vx(:)   = 0.0
        vy(:)   = 0.0
        vz(:)   = 0.0
	    v_cov(:,1) = 0. 
	    v_cov(:,2) = 0. 
	    v_cov(:,3) = 0. 
	    !
    ELSEWHERE 
        Lor(:) = DD(:)/rho(:)
        h(:) = (tau(:) + DD(:) + p(:))  ! rho*enthampy*Lor**2 
        v_cov(:,1) = sm_cov(:,1)/h(:) 
        v_cov(:,2) = sm_cov(:,2)/h(:) 
        v_cov(:,3) = sm_cov(:,3)/h(:) 
        !
#ifdef C2PFF    
        WHERE(Lor.LE.1.)  
            v2 = 0.
        ELSEWHERE 
            v2 =  SQRT( 1.0 - 1./LOR(:)**2 )  
        ENDWHERE  
        !   
        WHERE(rho(:)<NSTOV_rho_atmo*(1.0+0.02))     
            v_cov(:,1)=0.0
            v_cov(:,2)=0.0
            v_cov(:,3)=0.0
            v2 = 0.
        ENDWHERE
        !  
        WHERE( rho(:) < NSTOV_rho_atmo ) 
            rho(:) = NSTOV_rho_atmo
        END WHERE                    
        !
        p(:) = gam(:)*(w(:)*(1.-v2(:))-rho(:))  
        !
        !WHERE(p(:)<1e-12)
        WHERE(p(:)<NSTOV_p_atmo*(1.0+0.02))
            v_cov(:,1)=0.0
            v_cov(:,2)=0.0
            v_cov(:,3)=0.0
        ENDWHERE        
        WHERE(p(:)<NSTOV_p_atmo)
            p(:)   = NSTOV_p_atmo
        ENDWHERE      
#else
        WHERE(rho(:)<1e-5)    ! for the first successfull TOV star simulations, this stuff was computed before p. 
        !WHERE(rho(:)<NSTOV_rho_atmo*(1.0+0.02))     
            v_cov(:,1)=0.0
            v_cov(:,2)=0.0
            v_cov(:,3)=0.0
            v2 = 0.
        ENDWHERE
        ! 
        
        WHERE( rho(:) < 1e-10 ) 
            rho(:) = 1e-10 
        !WHERE( rho(:) < NSTOV_rho_atmo ) 
        !    rho(:) = NSTOV_rho_atmo
        END WHERE                    
        !
        p(:) = gam(:)*(w(:)*(1.-v2(:))-rho(:)) 
        !p(:) = gam(:)*(w(:)*(1.-v2(:))-rho(:))
        !
        WHERE(p(:)<1e-12)
        !WHERE(p(:)<NSTOV_p_atmo*(1.0+0.02))
            v_cov(:,1)=0.0
            v_cov(:,2)=0.0
            v_cov(:,3)=0.0
        ENDWHERE        
        !WHERE(p(:)<NSTOV_p_atmo)
        !    p(:)   = NSTOV_p_atmo
        !ENDWHERE      
#endif           
        !
    ENDWHERE
  ELSE
        Lor(:) = DD(:)/rho(:)
        h(:) = (tau(:) + DD(:) + p(:))  ! rho*enthampy*Lor**2 
        v_cov(:,1) = sm_cov(:,1)/h(:) 
        v_cov(:,2) = sm_cov(:,2)/h(:) 
        v_cov(:,3) = sm_cov(:,3)/h(:) 
        !
#ifdef C2PFF   
        WHERE(Lor.LE.1.)  
            v2 = 0.
        ELSEWHERE 
            v2 =  SQRT( 1.0 - 1./LOR(:)**2 )  
        ENDWHERE  
        !     
        WHERE(rho(:)<NSTOV_rho_atmo*(1.0+0.02))     
            v_cov(:,1)=0.0
            v_cov(:,2)=0.0
            v_cov(:,3)=0.0
            v2 = 0.
        ENDWHERE
        !  
        WHERE( rho(:) < NSTOV_rho_atmo ) 
            rho(:) = NSTOV_rho_atmo
        END WHERE                    
        !
        p(:) = gam(:)*(w(:)*(1.-v2(:))-rho(:))  
        !
        !WHERE(p(:)<1e-12)
        WHERE(p(:)<NSTOV_p_atmo*(1.0+0.02))
            v_cov(:,1)=0.0
            v_cov(:,2)=0.0
            v_cov(:,3)=0.0
        ENDWHERE        
        WHERE(p(:)<NSTOV_p_atmo)
            p(:)   = NSTOV_p_atmo
        ENDWHERE      
#else
        WHERE(rho(:)<1e-5)    ! for the first successfull TOV star simulations, this stuff was computed before p. 
        !WHERE(rho(:)<NSTOV_rho_atmo*(1.0+0.02))     
            v_cov(:,1)=0.0
            v_cov(:,2)=0.0
            v_cov(:,3)=0.0
            v2 = 0.
        ENDWHERE
        ! 
        
        WHERE( rho(:) < 1e-10 ) 
            rho(:) = 1e-10 
        !WHERE( rho(:) < NSTOV_rho_atmo ) 
        !    rho(:) = NSTOV_rho_atmo
        END WHERE                    
        !
        p(:) = gam(:)*(w(:)*(1.-v2(:))-rho(:)) 
        !p(:) = gam(:)*(w(:)*(1.-v2(:))-rho(:))
        !
        WHERE(p(:)<1e-12)
        !WHERE(p(:)<NSTOV_p_atmo*(1.0+0.02))
            v_cov(:,1)=0.0
            v_cov(:,2)=0.0
            v_cov(:,3)=0.0
        ENDWHERE        
        !WHERE(p(:)<NSTOV_p_atmo)
        !    p(:)   = NSTOV_p_atmo
        !ENDWHERE      
#endif            
        !
  ENDIF
  !
  !
    IF(ANY(p.LT.NSTOV_p_atmo)) THEN
        continue
    ENDIF
    !
#else
    v2 = 0.0 
    !
    v_cov(:,1) = V(:,2) 
    v_cov(:,2) = V(:,3) 
    v_cov(:,3) = V(:,4) 
    !
    v_contr(:,1) = g_contr(:,1,1)*V(:,2) + g_contr(:,1,2)*V(:,3) + g_contr(:,1,3)*V(:,4) 
    v_contr(:,2) = g_contr(:,2,1)*V(:,2) + g_contr(:,2,2)*V(:,3) + g_contr(:,2,3)*V(:,4) 
    v_contr(:,3) = g_contr(:,3,1)*V(:,2) + g_contr(:,3,2)*V(:,3) + g_contr(:,3,3)*V(:,4) 
    !
    v2 = v_cov(:,1)*v_contr(:,1) + v_cov(:,2)*v_contr(:,2) + v_cov(:,3)*v_contr(:,3) 
    !
	CALL RTSAFE_C2P_RMHD1_VECTOR(v2,x1,x2,tol,gam,dd,e,s2,b2,sb2,w,FAILED,VECTORLENGTH) 
    !
    IF(ANY(FAILED)) THEN
        WHERE(FAILED(:)) 
	        ! 
            iErrV(:)    = -1 
#ifdef C2PFF
            p(:)    = NSTOV_p_atmo
            rho(:)  = NSTOV_rho_atmo
#else
            p(:)    = p_floor
            rho (:) = rho_floor
#endif 
            vx(:)   = 0.0
            vy(:)   = 0.0
            vz(:)   = 0.0
            v_cov(:,1) = 0. 
            v_cov(:,2) = 0. 
            v_cov(:,3) = 0. 
            !
        ELSEWHERE 
            
            iErrV(:)    = 0 

            den(:)  = 1.0/(w(:)+b2(:))
            vb(:)   = sb(:)/w(:)
            !
            rho(:)  = dd(:)*sqrt(1.-v2(:))
            v_cov(:,1) = (sm_cov(:,1) + vb(:)*B_cov(:,1))*den(:)
            v_cov(:,2) = (sm_cov(:,2) + vb(:)*B_cov(:,2))*den(:)
            v_cov(:,3) = (sm_cov(:,3) + vb(:)*B_cov(:,3))*den(:)
            ! 
#ifdef C2PFF   
            WHERE(rho(:)<NSTOV_rho_atmo*(1.0+0.02))     
                v_cov(:,1)=0.0
                v_cov(:,2)=0.0
                v_cov(:,3)=0.0
                v2(:) = 0.
            ENDWHERE
            !  
            WHERE( rho(:) < NSTOV_rho_atmo ) 
                rho(:) = NSTOV_rho_atmo
            END WHERE                    
            !
            p(:) = gam(:)*(w(:)*(1.-v2(:))-rho(:))  
            !
            !WHERE(p(:)<1e-12)
            WHERE(p(:)<NSTOV_p_atmo*(1.0+0.02))
                v_cov(:,1)=0.0
                v_cov(:,2)=0.0
                v_cov(:,3)=0.0
            ENDWHERE        
            WHERE(p(:)<NSTOV_p_atmo)
                p(:)   = NSTOV_p_atmo
            ENDWHERE      
#else
            WHERE(rho(:)<1e-5)    ! for the first successfull TOV star simulations, this stuff was computed before p. 
            !WHERE(rho(:)<NSTOV_rho_atmo*(1.0+0.02))     
                v_cov(:,1)=0.0
                v_cov(:,2)=0.0
                v_cov(:,3)=0.0
                v2(:) = 0.
            ENDWHERE
            ! 
            
            WHERE( rho(:) < 1e-10 ) 
                rho(:) = 1e-10 
            !WHERE( rho(:) < NSTOV_rho_atmo ) 
            !    rho(:) = NSTOV_rho_atmo
            END WHERE                    
            !
            p(:) = gam(:)*(w(:)*(1.-v2(:))-rho(:)) 
            !p(:) = gam(:)*(w(:)*(1.-v2(:))-rho(:))
            !
            WHERE(p(:)<1e-12)
            !WHERE(p(:)<NSTOV_p_atmo*(1.0+0.02))
                v_cov(:,1)=0.0
                v_cov(:,2)=0.0
                v_cov(:,3)=0.0
            ENDWHERE        
            !WHERE(p(:)<NSTOV_p_atmo)
            !    p(:)   = NSTOV_p_atmo
            !ENDWHERE      
#endif  
            !
            
            
        ENDWHERE 
    ELSE
            iErrV(:)    = 0 
            den(:)  = 1.0/(w(:)+b2(:))
            vb(:)   = sb(:)/w(:)
            !
            rho(:)  = dd(:)*sqrt(1.-v2(:))
            v_cov(:,1) = (sm_cov(:,1) + vb(:)*B_cov(:,1))*den(:)
            v_cov(:,2) = (sm_cov(:,2) + vb(:)*B_cov(:,2))*den(:)
            v_cov(:,3) = (sm_cov(:,3) + vb(:)*B_cov(:,3))*den(:)
            ! 
#ifdef C2PFF   
            WHERE(rho(:)<NSTOV_rho_atmo*(1.0+0.02))     
                v_cov(:,1)=0.0
                v_cov(:,2)=0.0
                v_cov(:,3)=0.0
                v2(:) = 0.
            ENDWHERE
            !  
            WHERE( rho(:) < NSTOV_rho_atmo ) 
                rho(:) = NSTOV_rho_atmo
            END WHERE                    
            !
            p(:) = gam(:)*(w(:)*(1.-v2(:))-rho(:))  
            !
            !WHERE(p(:)<1e-12)
            WHERE(p(:)<NSTOV_p_atmo*(1.0+0.02))
                v_cov(:,1)=0.0
                v_cov(:,2)=0.0
                v_cov(:,3)=0.0
            ENDWHERE        
            WHERE(p(:)<NSTOV_p_atmo)
                p(:)   = NSTOV_p_atmo
            ENDWHERE      
#else
            WHERE(rho(:)<1e-5)    ! for the first successfull TOV star simulations, this stuff was computed before p. 
            !WHERE(rho(:)<NSTOV_rho_atmo*(1.0+0.02))     
                v_cov(:,1)=0.0
                v_cov(:,2)=0.0
                v_cov(:,3)=0.0
                v2(:) = 0.
            ENDWHERE
            ! 
            
            WHERE( rho(:) < 1e-10 ) 
                rho(:) = 1e-10 
            !WHERE( rho(:) < NSTOV_rho_atmo ) 
            !    rho(:) = NSTOV_rho_atmo
            END WHERE                    
            !
            p(:) = gam(:)*(w(:)*(1.-v2(:))-rho(:)) 
            !p(:) = gam(:)*(w(:)*(1.-v2(:))-rho(:))
            !
            WHERE(p(:)<1e-12)
            !WHERE(p(:)<NSTOV_p_atmo*(1.0+0.02))
                v_cov(:,1)=0.0
                v_cov(:,2)=0.0
                v_cov(:,3)=0.0
            ENDWHERE        
            !WHERE(p(:)<NSTOV_p_atmo)
            !    p(:)   = NSTOV_p_atmo
            !ENDWHERE      
#endif  
            !
    ENDIF
    !
#endif
    !
    !WHERE(Q(:,1)<1e-6)  
    !    rho(:) = Q(:,1)   
    !    p(:)   = (EQN%gamma-1.0)*Q(:,5)  
    !    ww(:)  = rho(:) + EQN%gamma/(EQN%gamma-1.0)*p(:)   
    !    v_cov(:,1) = Q(:,2)*ww(:)/(ww(:)*ww(:)+1e-3) 
    !    v_cov(:,2) = Q(:,3)*ww(:)/(ww(:)*ww(:)+1e-3) 
    !    v_cov(:,3) = Q(:,4)*ww(:)/(ww(:)*ww(:)+1e-3) 
    !END WHERE 
    !     
	V(:,1)  = rho(:)
	V(:,2)  = v_cov(:,1)
	V(:,3)  = v_cov(:,2)
	V(:,4)  = v_cov(:,3)
	V(:,5)  = p(:)
	V(:,6)  = B_contr(:,1)
	V(:,7)  = B_contr(:,2)
	V(:,8)  = B_contr(:,3)
	V(:,9)  = psi(:)
	V(:,10) = lapse(:)
	V(:,11) = shift(:,1)
	V(:,12) = shift(:,2)
	V(:,13) = shift(:,3)
    !
    iErr = MINVAL(iErrV)
    !
END SUBROUTINE PDECons2PrimVectorGRMHD        
!    
! Conservative part of the PDE ( flux tensor F(Q) ) 
!
RECURSIVE SUBROUTINE PDESourcePrimGRMHD(S,V,Q,time) 
    USE Parameters, ONLY : d, EQN 
    USE AstroMod
    USE EOS_mod
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
    REAL, INTENT(IN)  :: Q(nVar), V(nVar), time   
    REAL, INTENT(OUT) :: S(nVar) 
    ! Local variables 
    REAL :: p, irho, lam, mu 
    REAL :: k1, k2, fff, ggg, e, ds, c, xi, sk, sknl, alpha, fa, k0, beta0(3), eta 
    REAL :: g_cov(3,3), det, g_contr(3,3), Christoffel(3,3,3), dgup(3,3,3), faa, b0(3), dChristoffelSrc(3,3,3,3), RiemannSrc(3,3,3,3)
    REAL :: RicciSrc(3,3), gammaup(3), dtmetric(3,3), dtgup(3,3), dtChristoffelSrc(3,3,3), dtGammaUpSrc(3)
    REAL :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s25
    !
    S = 0.0     
    !
END SUBROUTINE PDESourcePrimGRMHD
!

RECURSIVE SUBROUTINE PDESourceGRMHD(S,Q) 
    USE Parameters, ONLY : d, EQN
    USE AstroMod
    USE EOS_mod
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
    REAL, INTENT(IN)  :: Q(nVar)   
    REAL, INTENT(OUT) :: S(nVar) 
    ! Local variables 
    REAL :: p, irho, lam, mu 
    REAL :: k1, k2, fff, ggg, e, ds, c, xi, sk, sknl, alpha, fa, k0, beta0(3), eta 
    REAL :: g_cov(3,3), det, g_contr(3,3), Christoffel(3,3,3), dgup(3,3,3), faa, b0(3), dChristoffelSrc(3,3,3,3), RiemannSrc(3,3,3,3)
    REAL :: RicciSrc(3,3), gammaup(3), dtmetric(3,3), dtgup(3,3), dtChristoffelSrc(3,3,3), dtGammaUpSrc(3)
    REAL :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s25
    !
    S = 0.0     
	!
END SUBROUTINE PDESourceGRMHD
!

RECURSIVE SUBROUTINE PDEEigenvaluesGRMHD(L,Q,n)
    USE Parameters, ONLY : d, EQN
    USE EOS_mod
	USE AstroMod
	!USE recipies_mod
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
    INTEGER, PARAMETER :: nParam=0 
    REAL, INTENT(IN)  :: Q(nVar),n(3)
    REAL, INTENT(OUT) :: L(nVar) 
    ! Local variables 
    REAL ::  Vp(nVar), gradQ(nVar,d)
    REAL :: p, irho
    REAL :: V(nVar),Qx(nVar), Qy(nVar), Qz(nVar) , Qp(nVar), Qm(nVar) 
    REAL :: e,g_cov(3,3),g_contr(3,3),nv(3)
    REAL :: AQx(nVar), BQy(nVar), CQz(nVar), FFp(nVar,d), FFm(nVar,d) ,dfdQ(nVar,nVar),AA(nVar,nVar)
    REAL :: lapse, shift(3),shift_cov(3), gammaij(6), delta(3,3), B_cov(3), vxB_cov(3), vxb_contr(3), psi, S_contr(3), qb_contr(3), B_contr(3) 
    REAL :: v2,v_contr(3),uem,b2,e2,gp,gm,vc(nVar),lf,w,ww,gamma1,rho,v_cov(3), w_ij, wim,eps,u,xg(3)
	REAL :: b2_4,cs2,sft,VdotB,vn,den,gg,lf2m,ImLambda(nVar), rtemp(nVar) ,R(nVar,nVar), iR(nVar,nVar),fmm(nVar),gmm(nVar),hmm(nVar),fpp(nVar),gpp(nVar),hpp(nVar),gg2,sft2,vn2
    INTEGER :: i,j,k,m,iErr, ccount, itemp(nVar) 
    LOGICAL :: FLAG
    REAL :: rhoeps,dd,tau,s2,sm_cov(3),sm(3)
    !
   !
    xg(1) = 0.
    xg(2) = 0.
    xg(3) = 0.
   !
    !
    L = 0.0
    !L(1) = -1.0    
    !L(2) = +1.0    
    !   RETURN  
    
   !****************************************************
    CALL PDECons2PrimGRMHD(V,Q,iErr)
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
    gamma1 = EQN%gamma/(EQN%gamma-1.0) 
    !
    !
    dd = gm*Q(1)
    tau = gm*Q(5)
    DO i=1,3 
        sm_cov(i)  = Q(1+i)
    ENDDO
    sm_cov = gm*sm_cov
    sm   = MATMUL (g_contr, sm_cov)  
    s2   = sm_cov(1)*sm(1) + sm_cov(2)*sm(2) + sm_cov(3)*sm(3)
    !
#ifdef GEOS  
    CALL EOS_eps(rhoeps,rho,p,dd,tau,s2) 
    w      = rho + rhoeps + p  + b2_4  ! rho*enthalpy
#else
    w      = rho + gamma1*p + b2_4 ! this is rho*h + b^2
#endif
    !IF(ABS(w-(rho + gamma1*p + b2_4)).GT.1e-13) THEN
    !    continue
    !ENDIF
    !
    !
    !w      = rho + gamma1*p + b2_4 ! this is rho*h + b^2
    cs2    = (EQN%gamma*p + b2_4)/w
    !
    vn     = v_contr(1)*n(1) + v_contr(2)*n(2) + v_contr(3)*n(3)
    sft    = shift(1)*n(1) + shift(2)*n(2) + shift(3)*n(3) 
    gg     = g_contr(1,1)*ABS(n(1)) + g_contr(2,2)*ABS(n(2)) + g_contr(3,3)*ABS(n(3))
    den    = 1.0/(1.0 - v2*cs2)
    !IF(SUM(n**2).EQ.0.) THEN  
    !    u = SQRT( v2) 
    !    WRITE(*,*)'Impossible error!'
    !        STOP
    !ELSE
        u = vn
    !ENDIF
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
    !
    ! SAFE MODE: define also 'covariant' eigenvalues! (we may use the remaining free slots in L, L(6:8), L(10:19)
	!
    shift_cov = MATMUL(g_cov,shift)
    !
    vn2  = v_cov(1)*n(1) + v_cov(2)*n(2) + v_cov(3)*n(3)
    sft2 = shift_cov(1)*n(1) + shift_cov(2)*n(2) + shift_cov(3)*n(3) 
    gg2  = g_cov(1,1)*ABS(n(1)) + g_cov(2,2)*ABS(n(2)) + g_cov(3,3)*ABS(n(3))
    den  = 1.0/(1.0 - v2*cs2)
    !
    u = vn2 
    !
    L(10)  = ( u*(1.0-cs2) - SQRT( cs2*lf2m*( (1.0-v2*cs2)*gg2 - u**2*(1.0-cs2) )) )*den
    L(6:8) = u
    L(11)  = ( u*(1.0-cs2) + SQRT( cs2*lf2m*( (1.0-v2*cs2)*gg2 - u**2*(1.0-cs2) )) )*den 
    L(6:11)   = lapse*L(6:11) - sft2
    !
    RETURN
    ! 
    !L(9) =  EQN%DivCleaning_a   ! 1.  !EQN%ch
    !
    !FLAG = .FALSE.
    !IF(MAXVAL(ABS(L)).GT.1.01) THEN
    !    FLAG = .TRUE.
    !    continue
    !ENDIF
    ! 
    !
    
    !
	!
	continue
	!
END SUBROUTINE PDEEigenvaluesGRMHD

RECURSIVE SUBROUTINE PDEMatrixBGRMHD(Bn,Q,nv)
    USE Parameters, ONLY :  d, EQN
	USE AstroMod
    USE EOS_mod
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
    REAL, INTENT(IN)  :: Q(nVar),nv(d)   
    REAL, INTENT(OUT) :: Bn(nVar,nVar) 
    ! Local variables 
    INTEGER :: i 
    REAL    :: p, irho, lam, mu, ialpha  
    REAL    :: B1(nVar,nVar), B2(nVar,nVar), B3(nVar,nVar)  
    REAL    :: A(nVar,nVar), B(nVar,nVar), C(nVar,nVar)  
    REAL    :: k1, k2, fff, ggg, e, ds, cs, xi, sk, sknl, alpha, fa, k0, beta0(3), b0(3)   
    REAL    :: g_cov(3,3), g_contr(3,3), det, Christoffel(3,3,3), dgup(3,3,3), uv(3), bb2   
  REAL    :: lapse, shift(3), gammaij(6), delta(3,3), B_cov(3), vxB_cov(3), vxb_contr(3), psi, S_contr(3),  B_contr(3) 
  REAL    :: v2,v_contr(3),uem,e2,gp,gm,vc(nVar),lf,w,ww,gamma1,rho,v_cov(3), w_ij, wim    
    INTEGER :: j,k,l,m,iErr, ccount    
  REAL :: rhoeps,dd,tau,s2,sm_cov(3),sm(3)
    !
    Bn = 0.0
    ! 
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
  DO j=1,3
	  DO i=1,3
		delta(i,j) = 0.
	  ENDDO 
  ENDDO
  !  
  DO i=1,3
     delta(i,i) = 1.0
  ENDDO 
  !
  CALL MatrixInverse3x3(g_cov,g_contr,gp)
  gp = SQRT(gp)
  gm = 1./gp
  ! 
  CALL PDECons2PrimGRMHD(Vc,Q,iErr)
  gamma1 = EQN%gamma/(EQN%gamma-1.0)
  rho    = Vc(1)
  DO i=1,3
	v_cov(i) = Vc(1+i)
	B_contr(i) = Vc(5+i)
  ENDDO
  p      = Vc(5)
  !
  B_cov = MATMUL(g_cov,B_contr(1:3))
  psi = Vc(9) 
  !
  v_contr     = MATMUL(g_contr,v_cov)
  S_contr = MATMUL(g_contr,Q(2:4))
  ! 
  ! COVARIANT =>
  vxB_cov(1) = v_contr(2)*B_contr(3) - v_contr(3)*B_contr(2)
  vxB_cov(2) = v_contr(3)*B_contr(1) - v_contr(1)*B_contr(3)
  vxB_cov(3) = v_contr(1)*B_contr(2) - v_contr(2)*B_contr(1)
  vxB_cov(1:3) = gp*vxB_cov(1:3)
  ! 
  ! CONTRAVARIANT =>
  vxB_contr(1) = v_cov(2)*B_cov(3) - v_cov(3)*B_cov(2)
  vxB_contr(2) = v_cov(3)*B_cov(1) - v_cov(1)*B_cov(3)
  vxB_contr(3) = v_cov(1)*B_cov(2) - v_cov(2)*B_cov(1)
  vxB_contr(1:3) = gm*vxB_contr(1:3)
  !
  v2     = v_contr(1)*v_cov(1) + v_contr(2)*v_cov(2) + v_contr(3)*v_cov(3)
  e2     = vxB_contr(1)*vxB_cov(1) + vxB_contr(2)*vxB_cov(2) + vxB_contr(3)*vxB_cov(3)
  bb2    = B_contr(1)*B_cov(1) + B_contr(2)*B_cov(2) + B_contr(3)*B_cov(3)
  !
  uem    = 0.5*(bb2 + e2) 
  !
  lf     = 1.0/sqrt(1.0 - v2)
  !
  !
  dd = gm*Q(1)
  tau = gm*Q(5)
  DO i=1,3 
      sm_cov(i)  = Q(1+i)
  ENDDO
  sm_cov = gm*sm_cov
  sm   = MATMUL (g_contr, sm_cov)  
  s2   = sm_cov(1)*sm(1) + sm_cov(2)*sm(2) + sm_cov(3)*sm(3)
  !
#ifdef GEOS  
  CALL EOS_eps(rhoeps,rho,p,dd,tau,s2) 
  w      = rho + rhoeps + p   ! rho*enthalpy
#else
  w      = rho + gamma1*p  ! this is rho*h
#endif
    !IF(ABS(w-(rho + gamma1*p)).GT.1e-13) THEN
    !    continue
    !ENDIF
    !
  !
  !w      = rho + gamma1*p   ! rho*hentalpy
  ww     = w*lf**2          ! rho*hentalpy*Lorentz^2 
  !
  !DO j=1,3
  A = 0.
  B = 0.
  C = 0.
    !
    ccount=0
    ! 
    DO i=1,3
        ! shift
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        j=1
        !------ 
        W_ij = ww*v_cov(i)*v_contr(j)-vxB_cov(i)*vxB_contr(j)-B_cov(i)*B_contr(j)+(p+uem)*delta(i,j)
        !
        A(1+j,10+i) = - Q(1+i)  ! Q(11:13)  shift(i) or shift_contr(i)
        A(5,10+i) = - gp*W_ij    ! Q(11:13)  shift(i) or shift_contr(i)
        !
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        j=2
        !------ 
        W_ij = ww*v_cov(i)*v_contr(j)-vxB_cov(i)*vxB_contr(j)-B_cov(i)*B_contr(j)+(p+uem)*delta(i,j)
        !
        B(1+j,10+i) = - Q(1+i)  ! Q(11:13)  shift(i) or shift_contr(i)
        B(5,10+i) = - gp*W_ij    ! Q(11:13)  shift(i) or shift_contr(i)
        ! 
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        j=3
        !------
        W_ij = ww*v_cov(i)*v_contr(j)-vxB_cov(i)*vxB_contr(j)-B_cov(i)*B_contr(j)+(p+uem)*delta(i,j)
        !
        C(1+j,10+i) = - Q(1+i)  ! Q(11:13)  shift(i) or shift_contr(i)
        C(5,10+i) = - gp*W_ij    ! Q(11:13)  shift(i) or shift_contr(i)
        ! 
          DO m=1,3
            IF(m.GE.i) THEN  
                !metric
                ccount=ccount+1
                !
                Wim = ww*v_contr(i)*v_contr(m)-vxB_contr(i)*vxB_contr(m)-B_contr(i)*B_contr(m)+(p+uem)*g_contr(i,m)
                !Wim = Wim + (1.0 - delta(i,m))*(ww*v_contr(m)*v_contr(i)-vxB_contr(m)*vxB_contr(i)-B_contr(m)*B_contr(i)+(p+uem)*g_contr(m,i))  ! account also of the remaining symmetric components of gamma for i.NE.m.
                IF(i.NE.m) THEN !TWICE! account also of the remaining symmetric components of gamma §(memo: g_contr and Wim are symmetric)
                    Wim = 2.0*Wim   
                ENDIF
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=1
                !------ 
                A(1+j,13+ccount) = - 0.5*gp*lapse*Wim  ! Q(14:19) gammaij(ccount) or  g_cov(i,m)
                A(5,13+ccount) = - 0.5*gp*Wim*shift(j)   ! Q(14:19) gammaij(ccount) or  g_cov(i,m)

#ifdef COVCLEAN     
				A(9,13+ccount) =  + 0.5*Q(9)*gammaij(ccount)*shift(j)    ! D.C.
                A(6:8,13+ccount) =  - 0.5*lapse*Q(9)*g_contr(1:3,j)*gammaij(ccount) 
                !
                IF(i.NE.m) THEN !TWICE!  account  also of the remaining symmetric components of gamma §(memo: g_contr and Wim are symmetric)
                    A(9,13+ccount) = A(9,13+ccount) + 0.5*Q(9)*gammaij(ccount)*shift(j)  ! D.C.
                    A(6:8,13+ccount) = A(6:8,13+ccount) - 0.5*lapse*Q(9)*g_contr(1:3,j)*gammaij(ccount) 
                ENDIF
#endif  	
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=2
                !------ 
                B(1+j,13+ccount) = - 0.5*gp*lapse*Wim  ! Q(14:19) gammaij(ccount) or  g_cov(i,m)
                B(5,13+ccount) = - 0.5*gp*Wim*shift(j)   ! Q(14:19) gammaij(ccount) or  g_cov(i,m)

#ifdef COVCLEAN     
				B(9,13+ccount) =  + 0.5*Q(9)*gammaij(ccount)*shift(j)   ! D.C.
                B(6:8,13+ccount) =  - 0.5*lapse*Q(9)*g_contr(1:3,j)*gammaij(ccount)
                !
                IF(i.NE.m) THEN !TWICE!  account  also of the remaining symmetric components of gamma §(memo: g_contr and Wim are symmetric)
                    B(9,13+ccount) = B(9,13+ccount) + 0.5*Q(9)*gammaij(ccount)*shift(j)   ! D.C.
                    B(6:8,13+ccount) = B(6:8,13+ccount) - 0.5*lapse*Q(9)*g_contr(1:3,j)*gammaij(ccount)
                ENDIF
#endif  		
                ! 
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=3
                !------
                C(1+j,13+ccount) = - 0.5*gp*lapse*Wim  ! Q(14:19) gammaij(ccount) or  g_cov(i,m)
                C(5,13+ccount) = - 0.5*gp*Wim*shift(j)   ! Q(14:19) gammaij(ccount) or  g_cov(i,m)
#ifdef COVCLEAN     
				C(9,13+ccount) =  + 0.5*Q(9)*gammaij(ccount)*shift(j)   ! D.C.
                C(6:8,13+ccount) =  - 0.5*lapse*Q(9)*g_contr(1:3,j)*gammaij(ccount)
                !
                IF(i.NE.m) THEN !TWICE!  account  also of the remaining symmetric components of gamma §(memo: g_contr and Wim are symmetric)
                    C(9,13+ccount) = C(9,13+ccount) + 0.5*Q(9)*gammaij(ccount)*shift(j)   ! D.C.
                    C(6:8,13+ccount) = C(6:8,13+ccount) - 0.5*lapse*Q(9)*g_contr(1:3,j)*gammaij(ccount)
                ENDIF
#endif  		
                ! 
            ENDIF
          ENDDO
      ENDDO
    !lapse
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=1
    !------ 
    A(1+j,10) = + (Q(5)+Q(1))   ! Q(10) or lapse
    A(5,10) =  S_contr(j)     !  Q(10) or lapse

#ifdef COVCLEAN     
	A(9,10) = - Q(5+j)             ! D.C.
    A(9,10+j) = + Q(9)             ! D.C.   phi * d(shift^x)/dx
    A(6:8,5+j) = - Q(11:13)             ! D.C.   shift * d(gp*B^x)/dx
    A(6:8,9) = lapse*g_contr(1:3,j)  ! D.C.  
#else
    A(6:8,5+j) = - Q(11:13)             ! D.C.   shift * d(gp*B^x)/dx 
    A(6:8,9) = lapse*gp*g_contr(1:3,j)  ! D.C.   
    A(9,5+j) = gm*lapse*EQN%DivCleaning_a**2             ! D.C.   shift * d(gp*B^x)/dx 
    A(9,9) = - Q(10+j)            ! D.C.   shift * d(gp*B^x)/dx  
#endif
    !
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=2
    !------ 
    B(1+j,10) = + Q(5)+Q(1)   ! Q(10) or lapse
    B(5,10) =  S_contr(j)     !  Q(10) or lapse
!
#ifdef COVCLEAN     
	B(9,10) = - Q(5+j)            ! D.C.
    B(9,10+j) = + Q(9)            ! D.C.   phi * d(shift^y)/dy
    B(6:8,5+j) = - Q(11:13)                  ! D.C.   shift * d(gp*B^y)/dy
    B(6:8,9) = + lapse*g_contr(1:3,j)       ! D.C.  
#else
    B(6:8,5+j) =  - Q(11:13)                  ! D.C.   shift * d(gp*B^x)/dx 
    B(6:8,9) =  + lapse*gp*g_contr(1:3,j)    ! D.C.  
    B(9,5+j) = + gm*lapse*EQN%DivCleaning_a**2            ! D.C.   shift * d(gp*B^x)/dx 
    B(9,9) = - Q(10+j)             ! D.C.   shift * d(gp*B^x)/dx  
#endif
    ! 
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=3
    !------
    C(1+j,10) = + Q(5)+Q(1)   ! Q(10) or lapse
    C(5,10) =  S_contr(j)     !  Q(10) or lapse
	!
#ifdef COVCLEAN     
	C(9,10) = - Q(5+j)             ! D.C.
    C(9,10+j) = + Q(9)             ! D.C.   phi *d(shift^z)/dz
    C(6:8,5+j) = - Q(11:13)                  ! D.C.   shift * d(gp*B^z)/dz
    C(6:8,9) = + lapse*g_contr(1:3,j)        ! D.C.  
#else
    C(6:8,5+j) = - Q(11:13)                  ! D.C.   shift * d(gp*B^x)/dx 
    C(6:8,9) = + lapse*gp*g_contr(1:3,j)     ! D.C.  
    C(9,5+j) = + gm*lapse*EQN%DivCleaning_a**2             ! D.C.   shift * d(gp*B^x)/dx 
    C(9,9) = - Q(10+j)             ! D.C.   shift * d(gp*B^x)/dx  
#endif
  !
  Bn = A*nv(1) + B*nv(2) + C*nv(3) 
  !  
END SUBROUTINE PDEMatrixBGRMHD  
    

RECURSIVE SUBROUTINE PDECons2PrimGRMHD(V,Q,iErr)
    USE Parameters, ONLY : EQN,p_floor,rho_floor,NSTOV_rho_atmo,NSTOV_p_atmo
    USE EOS_mod
	USE AstroMod
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
    REAL :: g_cov(3,3),g_contr(3,3),gammaij(6),rhoeps
    REAL :: vb,v2,e2,b2,s2,sb,sb2,e,x1,x2,eps,uem,LF,gamma1,gam,w,ww,depsdrho,depsdp
	REAL :: vx,vy,vz,den,dd
    REAL :: p0, tau,dprho,dpeps, Lor, h,ptest,sqrtX
    REAL :: tol
	INTEGER :: i
	LOGICAL :: FAILED
	!  
	iErr = 0     
        tol = 1.0e-18
        !
	!V = Q
	!RETURN
	V=0.
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
    tau = Qloc(5)
	DO i=1,3
        B_contr(i) = Qloc(5+i) ! !wrong: magnetic field is contravariant
		sm_cov(i)  = Qloc(1+i)
	ENDDO
	!
	gamma1 = EQN%gamma/(EQN%gamma - 1.0)
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
    e = tau + dd  ! Q(5) = gamma^1/2 ( U - dd )
    ! THIS IS ONLY FOR THE GRHD.
#ifdef GEOS
  IF(v(5).GT.0) THEN
    p = v(5)
  ELSE
    !p0 = 0.5*(Q(5)+Q(1)) !1000000000*NSTOV_p_atmo  ! initial guess 
    p = 0.5*(1e-18+1e7)
  ENDIF
  !CALL C2P_GEOS(p,rho,eps,p0,tau,DD,s2,FAILED)  !  (tau,D,s2) => (p,rho,eps)  
  !
  x1   = 1.0e-18      
  x2   = 1.0e7    !tol ! 
  tol = 1e-18
  CALL RTSAFE_C2P_RHD1(p,X1,X2,tol,dd,tau,s2,FAILED)
  !
  IF(FAILED) THEN
      continue
  ENDIF
  !
  den = tau+p+DD 
  sqrtX = SQRT(den**2-s2)
  !
  rho =DD*sqrtX/den
  eps =sqrtX-p*den/sqrtX-DD
  continue
  !
  !rho = V(1)      
  !v_cov(1) = V(2) 
  !v_cov(2) = V(3) 
  !v_cov(3) = V(4) 
  !p = V(5)        
  !!
  !IF(V(1).EQ.0) THEN
  !  p = NSTOV_p_atmo  ! initial guess 
  !  w = 0.            ! initial guess 
  !  v2=0.
  !  rho = dd*SQRT(1.0-v2)            ! initial guess 
  !  w = (rho*(1.0+0.5)+0.01)
  !  den  = 1.0/(w+b2)
  !  vb   = sb/w
  !  v_cov(1) = (sm_cov(1) + vb*B_cov(1))*den
  !  v_cov(2) = (sm_cov(2) + vb*B_cov(2))*den
  !  v_cov(3) = (sm_cov(3) + vb*B_cov(3))*den
  !  p = max(NSTOV_p_atmo, gam*(w*(1.-v2)-rho)) 
  !ENDIF
  !!
  !v2 =      g_contr(1,1)*v_cov(1)*v_cov(1)+g_contr(1,2)*v_cov(2)*v_cov(1)+g_contr(1,3)*v_cov(3)*v_cov(1)
  !v2 = v2 + g_contr(2,1)*v_cov(1)*v_cov(2)+g_contr(2,2)*v_cov(2)*v_cov(2)+g_contr(2,3)*v_cov(3)*v_cov(2)
  !v2 = v2 + g_contr(3,1)*v_cov(1)*v_cov(3)+g_contr(3,2)*v_cov(2)*v_cov(3)+g_contr(3,3)*v_cov(3)*v_cov(3)
  !! 
  !Lor = 1.0/SQRT(1.0-v2)
  !CALL EOS_eps_depsdrho_depsdp(eps,depsdrho,depsdp,rho,p) 
  !w = (rho*(1.0+eps)+p)*Lor**2
  !!  
  !CALL C2P_GEOS_3D(w,p,v2,gam,dd,e,s2,b2,sb2,FAILED)  !  (tau,D,s2) => (p,rho,eps)
  !!
  IF(FAILED) THEN
     IF(Q(1).GT.1e-10.AND.Q(5).GT.1e-10) iErr = +1
     rho = NSTOV_rho_atmo
     CALL EOS_p_dprho_dpeps(p,dprho,dpeps,rho,eps)
     rho  = rho_floor
     vx   = 0.0
     vy   = 0.0
     vz   = 0.0
     v_cov = 0.0  
     V(1:19) = (/ rho, v_cov(1:3), p, B_contr(1:3), psi , lapse, shift(1:3), gammaij(1:6)/)
     !CALL PDEPrim2Cons(Qtmp,V,x0_000,0.) 
     !CALL PDECons2Prim(V,Qtmp,x0_000,0.,iErrTMP)   
     continue
     RETURN
  ELSE
       ! C2P_GEOS pure hydro
      Lor = DD/rho
      h = (tau + DD + p)  ! rho*enthampy*Lor**2
      v_cov(1:3) = sm_cov(1:3)/h
      !
      ! C2P_GEOS_3D
      !den  = 1.0/(w+b2)
      !vb   = sb/w
      !!
      !rho  = dd*sqrt(1.-v2)
      !v_cov(1) = (sm_cov(1) + vb*B_cov(1))*den
      !v_cov(2) = (sm_cov(2) + vb*B_cov(2))*den
      !v_cov(3) = (sm_cov(3) + vb*B_cov(3))*den
      !!p = max(1.e-15, gam*(w*(1.-v2)-rho))
      !!p = gam*(w*(1.-v2)-rho)
  ENDIF
  !

#else
	x1   = 0.      ! 
    eps = 1.0e-18
    x2   = 1.0-1e-10 !tol ! 
    tol = 1e-18
	w=0
	!
	CALL RTSAFE_C2P_RMHD1(v2,x1,x2,tol,gam,dd,e,s2,b2,sb2,w,FAILED) 
	!
	IF(FAILED) THEN
		 ! 
		 iErr = -1
#ifdef C2PFF
        p    = NSTOV_p_atmo
        rho  = NSTOV_rho_atmo
#else
		 p    = p_floor
		 rho  = rho_floor
#endif
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
#endif  
  !
#ifdef C2PFF   
#ifdef GEOS
        IF(Lor.LE.1.) THEN
            v2 = 0.
        ELSE 
            v2 = SQRT(1.0 - 1.0/Lor**2)
        ENDIF
        !   
        CALL EOS_eps(rhoeps,rho,p,dd,tau,s2) 
        w      = rho + rhoeps + p   ! rho*enthalpy  
#endif
        IF(rho<NSTOV_rho_atmo*(1.0+0.02))      THEN
            v_cov(1)=0.0
            v_cov(2)=0.0
            v_cov(3)=0.0
            v2 = 0.
            p = gam*(w*(1.-v2)-rho)  
        ENDIF
        !  
        IF( rho < NSTOV_rho_atmo )  THEN
            rho = NSTOV_rho_atmo
            p = gam*(w*(1.-v2)-rho)  
        END IF                    
        !
        !
        !IF(p(:)<1e-12)
        IF(p<NSTOV_p_atmo*(1.0+0.02)) THEN
            v_cov(1)=0.0
            v_cov(2)=0.0
            v_cov(3)=0.0
            v2 = 0.
            p = gam*(w*(1.-v2)-rho)  
        ENDIF        
        IF(p<NSTOV_p_atmo) THEN
            p   = NSTOV_p_atmo
        ENDIF      
#else
    IF(rho<NSTOV_rho_atmo*1.01) THEN   ! for the first successfull TOV star simulations, this stuff was computed before p. 
        v_cov(:)= 0.0
        !rho = 1e-12 
        !p   = 1e-12 
    ENDIF
    IF(rho<1e-11) THEN 
        rho = 1e-11
    !IF(rho<NSTOV_rho_atmo) THEN 
    !    rho = NSTOV_rho_atmo
	ENDIF
	!
    IF(p<1e-12) THEN 
        v_cov(:)= 0.0
        !PRINT *,'I"m here!'
        !pause
        !rho = 1e-12 
        p   = 1e-12  
    !IF(p<NSTOV_p_atmo*1.01) THEN 
    !    v_cov(:)= 0.0 
    ENDIF   
    !IF(p<NSTOV_p_atmo) THEN 
    !    p=NSTOV_p_atmo
    !ENDIF  
#endif    
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
    !Qtmp =Qtmp
    continue
	!
END SUBROUTINE PDECons2PrimGRMHD  

    
    
    
RECURSIVE SUBROUTINE PDEPrim2ConsGRMHD(Q,V)
    USE Parameters, ONLY : d, EQN 
    USE EOS_mod
	USE AstroMod
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
    REAL, INTENT(IN)      :: V(nVar)     ! primitive variables 
    REAL, INTENT(OUT)     :: Q(nVar)     ! vector of conserved quantities 
    ! Local variables 
  ! Local variable declaration
  REAL :: rho,p,psi,lapse,gp,gm
  REAL :: shift(3),v_cov(3),v_contr(3),vxB_cov(3),vxB_contr(3),B_cov(3),B_contr(3)
  REAL :: g_cov(3,3),g_contr(3,3)
  REAL :: vb,v2,e2,b2,uem,LF,gamma1,w,ww
  REAL :: rhoeps,b2_4
  INTEGER :: i
  !
  REAL :: Dd,gam,sb2,s_cov(3),s_contr(3),e,s2,f(3),df(3,3),x(3)
	!Q=V
	!RETURN
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
  gamma1 = EQN%gamma/(EQN%gamma-1.0)
  !
#ifdef GEOS  
  CALL EOS_eps_P2C(rhoeps,rho,p) 
  w      = rho + rhoeps + p  ! rho*enthalpy
#else
  w      = rho + gamma1*p  ! this is rho*h
#endif
    !IF(ABS(w-(rho + gamma1*p)).GT.1e-13) THEN
    !    continue
    !ENDIF
  !
  !w      = rho + gamma1*p
  ww     = w*lf**2
  !
  Q(1)  = rho*lf
  Q(2)  = ww*v_cov(1) + b2*v_cov(1) - vb*B_cov(1)
  Q(3)  = ww*v_cov(2) + b2*v_cov(2) - vb*B_cov(2)
  Q(4)  = ww*v_cov(3) + b2*v_cov(3) - vb*B_cov(3) 
  Q(5)    = ww - p + uem - Q(1)     !!!!! we subtract PDE(Q(1))!!!!
  Q(6) = V(6)
  Q(7) = V(7)
  Q(8) = V(8)
  !
  !! test functional
  !gam = 1.0/gamma1
  !dd=Q(1)
  !s_cov(1) = Q(2)
  !s_cov(2) = Q(3)
  !s_cov(3) = Q(4)
  !e=(Q(5)+Q(1))
  !s_contr = MATMUL(g_contr,s_cov)
  !s2 = s_cov(1)*s_contr(1) + s_cov(2)*s_contr(2) + s_cov(3)*s_contr(3)
  !sb2= (s_cov(1)*B_contr(1) + s_cov(2)*B_contr(2) + s_cov(3)*B_contr(3))**2
  !x(1) = v2
  !x(2) = ww
  !x(3) = p
  !CALL FUNC_C2P_GEOS_3D(x,f,df,gam,dd,e,s2,b2,sb2)
  !!
  !IF(SUM(ABS(f)).GT.1e-13)THEN
  !      continue
  !ENDIF
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
END SUBROUTINE PDEPrim2ConsGRMHD       
            
 
RECURSIVE SUBROUTINE PDEPrim2ConsGRMHDvector(Q,V)
    USE Parameters, ONLY : d,VECTORLENGTH , EQN
    USE EOS_mod
	USE AstroMod
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
  REAL(8), INTENT(IN)  :: V(VECTORLENGTH,nVar)     ! primitive variables
  REAL(8), INTENT(OUT) :: Q(VECTORLENGTH,nVar)     ! vector of conserved quantities
  ! Local variables
  REAL(8), DIMENSION(VECTORLENGTH)     :: rho,p,psi,lapse,gp,gm,vb,v2,e2,b2,uem,LF,gamma1,w,ww,rhoeps,b2_4,Dd,gam,sb2,e,s2,det,idet
  REAL(8), DIMENSION(VECTORLENGTH,3)   :: shift,v_cov,v_contr,vxB_cov,vxB_contr,B_cov,B_contr,s_cov,s_contr,x,f
  REAL(8), DIMENSION(VECTORLENGTH,3,3) :: g_cov,g_contr,df
  INTEGER :: i
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64:: rho,p,psi,lapse,gp,gm,vb,v2,e2,b2,uem,LF,gamma1,w,ww,rhoeps,b2_4,Dd,gam,sb2,e,s2
    !DIR$ ATTRIBUTES ALIGN: 64:: shift,v_cov,v_contr,vxB_cov,vxB_contr,B_cov,B_contr,s_cov,s_contr,x,f
    !DIR$ ATTRIBUTES ALIGN: 64:: g_cov,g_contr,df
    !DIR$ ASSUME_ALIGNED Q     : 64
    !DIR$ ASSUME_ALIGNED V     : 64
#else
    !DIR$ ATTRIBUTES ALIGN: 32:: rho,p,psi,lapse,gp,gm,vb,v2,e2,b2,uem,LF,gamma1,w,ww,rhoeps,b2_4,Dd,gam,sb2,e,s2
    !DIR$ ATTRIBUTES ALIGN: 32:: shift,v_cov,v_contr,vxB_cov,vxB_contr,B_cov,B_contr,s_cov,s_contr,x,f
    !DIR$ ATTRIBUTES ALIGN: 32:: g_cov,g_contr,df
    !DIR$ ASSUME_ALIGNED Q     : 32
    !DIR$ ASSUME_ALIGNED V     : 32
#endif
  !Q=V
  !RETURN
  DO i=1,3
    v_cov(:,i) = V(:,i+1)
	B_contr(:,i) = V(:,5+i)  !wrong 
	shift(:,i) = V(:,10+i)          ! NB: we choose V() and Q() being the shift_controvariant!
  ENDDO
  rho(:)    = V(:,1)
  p(:)      = V(:,5)
  psi(:)    = V(:,9)
  lapse(:)  = V(:,10)
  !
  !gammaij = V(14:19) 
  g_cov(:,1,1) = V(:,14)
  g_cov(:,1,2) = V(:,15)
  g_cov(:,1,3) = V(:,16)
  g_cov(:,2,2) = V(:,17)
  g_cov(:,2,3) = V(:,18)
  g_cov(:,3,3) = V(:,19)
  g_cov(:,2,1) = V(:,15)
  g_cov(:,3,1) = V(:,16)
  g_cov(:,3,2) = V(:,18)
  !
  det(:) = (g_cov(:,1,1)*g_cov(:,2,2)*g_cov(:,3,3)-g_cov(:,1,1)*g_cov(:,2,3)**2-g_cov(:,1,2)**2*g_cov(:,3,3)+2*g_cov(:,1,2)*g_cov(:,1,3)*g_cov(:,2,3)-g_cov(:,1,3)**2*g_cov(:,2,2))  
  idet(:) = 1./det(:)  
  !
  g_contr(:,1,1) =  ( g_cov(:,2,2)*g_cov(:,3,3)-g_cov(:,2,3)*g_cov(:,3,2)) * idet(:) 
  g_contr(:,1,2) = -( g_cov(:,1,2)*g_cov(:,3,3)-g_cov(:,1,3)*g_cov(:,3,2)) * idet(:)
  g_contr(:,1,3) = -(-g_cov(:,1,2)*g_cov(:,2,3)+g_cov(:,1,3)*g_cov(:,2,2)) * idet(:) 
  g_contr(:,2,1) = -( g_cov(:,2,1)*g_cov(:,3,3)-g_cov(:,2,3)*g_cov(:,3,1)) * idet(:) 
  g_contr(:,2,2) =  ( g_cov(:,1,1)*g_cov(:,3,3)-g_cov(:,1,3)*g_cov(:,3,1)) * idet(:) 
  g_contr(:,2,3) = -( g_cov(:,1,1)*g_cov(:,2,3)-g_cov(:,1,3)*g_cov(:,2,1)) * idet(:) 
  g_contr(:,3,1) = -(-g_cov(:,2,1)*g_cov(:,3,2)+g_cov(:,2,2)*g_cov(:,3,1)) * idet(:) 
  g_contr(:,3,2) = -( g_cov(:,1,1)*g_cov(:,3,2)-g_cov(:,1,2)*g_cov(:,3,1)) * idet(:) 
  g_contr(:,3,3) =  ( g_cov(:,1,1)*g_cov(:,2,2)-g_cov(:,1,2)*g_cov(:,2,1)) * idet(:)       
  ! 
  gp(:) = SQRT(det(:))
  gm(:) = 1./gp(:)
  !  
  v_contr(:,1) = g_contr(:,1,1)*v_cov(:,1) + g_contr(:,1,2)*v_cov(:,2) + g_contr(:,1,3)*v_cov(:,3) 
  v_contr(:,2) = g_contr(:,2,1)*v_cov(:,1) + g_contr(:,2,2)*v_cov(:,2) + g_contr(:,2,3)*v_cov(:,3) 
  v_contr(:,3) = g_contr(:,3,1)*v_cov(:,1) + g_contr(:,3,2)*v_cov(:,2) + g_contr(:,3,3)*v_cov(:,3) 
  B_cov(:,1) = g_cov(:,1,1)*B_contr(:,1) + g_cov(:,1,2)*B_contr(:,2) + g_cov(:,1,3)*B_contr(:,3) 
  B_cov(:,2) = g_cov(:,2,1)*B_contr(:,1) + g_cov(:,2,2)*B_contr(:,2) + g_cov(:,2,3)*B_contr(:,3) 
  B_cov(:,3) = g_cov(:,3,1)*B_contr(:,1) + g_cov(:,3,2)*B_contr(:,2) + g_cov(:,3,3)*B_contr(:,3)   
  !
  ! COVARIANT =>
  vxB_cov(:,1) = v_contr(:,2)*B_contr(:,3) - v_contr(:,3)*B_contr(:,2)
  vxB_cov(:,2) = v_contr(:,3)*B_contr(:,1) - v_contr(:,1)*B_contr(:,3)
  vxB_cov(:,3) = v_contr(:,1)*B_contr(:,2) - v_contr(:,2)*B_contr(:,1)
  vxB_cov(:,1) = gp(:)*vxB_cov(:,1)
  vxB_cov(:,2) = gp(:)*vxB_cov(:,2)
  vxB_cov(:,3) = gp(:)*vxB_cov(:,3)
  !
  ! CONTRAVARIANT =>
  vxB_contr(:,1) = v_cov(:,2)*B_cov(:,3) - v_cov(:,3)*B_cov(:,2)
  vxB_contr(:,2) = v_cov(:,3)*B_cov(:,1) - v_cov(:,1)*B_cov(:,3)
  vxB_contr(:,3) = v_cov(:,1)*B_cov(:,2) - v_cov(:,2)*B_cov(:,1)
  vxB_contr(:,1) = gm(:)*vxB_contr(:,1)
  vxB_contr(:,2) = gm(:)*vxB_contr(:,2)
  vxB_contr(:,3) = gm(:)*vxB_contr(:,3)
  !
  vb(:) = v_cov(:,1)*B_contr(:,1) + v_cov(:,2)*B_contr(:,2) + v_cov(:,3)*B_contr(:,3)      
  v2(:) = v_contr(:,1)*v_cov(:,1) + v_contr(:,2)*v_cov(:,2) + v_contr(:,3)*v_cov(:,3)
  e2(:) = vxB_contr(:,1)*vxB_cov(:,1) + vxB_contr(:,2)*vxB_cov(:,2) + vxB_contr(:,3)*vxB_cov(:,3)
  b2(:) = B_contr(:,1)*B_cov(:,1) + B_contr(:,2)*B_cov(:,2) + B_contr(:,3)*B_cov(:,3)
  !
  uem(:) = 0.5*(b2(:) + e2(:)) 
  !
  lf(:)     = 1.0 / sqrt(1.0 - v2(:))
  gamma1(:) = EQN%gamma/(EQN%gamma-1.0)  
  !
#ifdef GEOS  
  CALL EOS_eps_P2C_vector(rhoeps(:),rho(:),p(:)) 
  w(:) = rho(:) + rhoeps(:) + p(:)  ! rho*enthalpy
#else
  w(:)    = rho(:) + gamma1(:)*p(:)  ! this is rho*h
#endif
  !IF(SUM(ABS(w-(rho + gamma1*p))).GT.1e-13) THEN
  !    continue
  !ENDIF
  !
  !w      = rho + gamma1*p
  ww(:) = w(:)*lf(:)**2
  !
  Q(:,1) = rho(:)*lf(:)
  Q(:,2) = ww(:)*v_cov(:,1) + b2(:)*v_cov(:,1) - vb(:)*B_cov(:,1)
  Q(:,3) = ww(:)*v_cov(:,2) + b2(:)*v_cov(:,2) - vb(:)*B_cov(:,2)
  Q(:,4) = ww(:)*v_cov(:,3) + b2(:)*v_cov(:,3) - vb(:)*B_cov(:,3) 
  Q(:,5) = ww(:) - p(:) + uem(:) - Q(:,1)     !!!!! we subtract PDE(:,Q(:,1))!!!!
  Q(:,6) = V(:,6)
  Q(:,7) = V(:,7)
  Q(:,8) = V(:,8)
  !
  !! test functional
  !gam = 1.0/gamma1
  !dd=Q(:,1)
  !s_cov(:,1) = Q(:,2)
  !s_cov(:,2) = Q(:,3)
  !s_cov(:,3) = Q(:,4)
  !e=(:,Q(:,5)+Q(:,1))
  !s_contr = MATMUL(:,g_contr,s_cov)
  !s2 = s_cov(:,1)*s_contr(:,1) + s_cov(:,2)*s_contr(:,2) + s_cov(:,3)*s_contr(:,3)
  !sb2= (:,s_cov(:,1)*B_contr(:,1) + s_cov(:,2)*B_contr(:,2) + s_cov(:,3)*B_contr(:,3))**2
  !x(:,1) = v2
  !x(:,2) = ww
  !x(:,3) = p
  !CALL FUNC_C2P_GEOS_3D(:,x,f,df,gam,dd,e,s2,b2,sb2)
  !!
  !IF(:,SUM(:,ABS(:,f)).GT.1e-13)THEN
  !      continue
  !ENDIF
  !
  DO i=1,8
		Q(:,i) = gp(:)*Q(:,i)  
  ENDDO
  DO i=9,nVar
		Q(:,i) = V(:,i)
  ENDDO
  ! 
  continue
    !
END SUBROUTINE PDEPrim2ConsGRMHDvector      
    
			
RECURSIVE SUBROUTINE RTSAFE_C2P_RMHD1(v2,X1,X2,XACC,gam,d,e,s2,b2,sb2,w,FAILED)
  USE Parameters, ONLY : nVar
    USE AstroMod
  IMPLICIT NONE
  ! argument variables
  REAL, INTENT(INOUT)       :: v2,X1,X2,XACC,w 
  REAL, INTENT(IN)          :: gam,d,e,s2,b2,sb2
  LOGICAL, INTENT(INOUT)    :: FAILED
  ! auxiliary variables
  INTEGER, PARAMETER   :: MAXIT=400
  INTEGER              :: J,i
  REAL                 :: FL,FH,DF,XH,XL,SWAP,DXOLD,XACC2,DX,F,TEMP,Xtmp,Vtmp(7),Ftmp,dFtmp,dXtmp
  ! 
  !IF(DebugPrint) THEN 
  !    WRITE(DebugFile,'(a)') 'DebugFile.dat' 
  !    OPEN(UNIT=4321,FILE=TRIM(DebugFile), RECL=800) 
  !    WRITE(4321,*) '  VARIABLES = "i"   "Xtmp"   "Ftmp"   "dFtmp"   "gam" "d" "e" "s2" "b2" "sb2" "w" '
  !    WRITE(4321,*)
  !    Vtmp(1:7) = (/ gam,d,e,s2,b2,sb2,w /)
  !    dXtmp = (1.0-1.0e-12)/300
  !    DO i =1, 350
  !        Xtmp = 0.5*dXtmp + real(i)*dXtmp
  !      CALL FUNC_C2P_RMHD1(Xtmp,Ftmp,dFtmp,Vtmp(1),Vtmp(2),Vtmp(3),Vtmp(4),Vtmp(5),Vtmp(6),Vtmp(7))
  !      WRITE(4321,'(i9.4,10(E16.6))') i,Xtmp,Ftmp,dFtmp,Vtmp(1:7)
  !    ENDDO
  !    CLOSE(4321)   
  !ENDIF
  !
  !PRINT *,e
	!
  XACC2 = 1e-60 !XACC(:)**2
  FAILED = .FALSE.
  CALL FUNC_C2P_RMHD1(X1,FL,DF,gam,d,e,s2,b2,sb2,w)
#ifdef C2PFF  
  IF(ABS(FL).LT.XACC2) THEN
      continue
     v2=X1
     RETURN
  ENDIF
#else
  IF(FL.EQ.0.) THEN
     v2=X1
     RETURN
  ENDIF
#endif 
  CALL FUNC_C2P_RMHD1(X2,FH,DF,gam,d,e,s2,b2,sb2,w)
#ifdef C2PFF  
  IF(ABS(FH).LT.XACC2) THEN
      continue
     v2=X2
     RETURN
  ENDIF
#else
  IF(FH.EQ.0.) THEN
     v2=X2
     RETURN
  ENDIF
#endif 
  IF(FL*FH.GT.0.) THEN
     FAILED = .TRUE.
     v2 = 0. ! assign value even if it fails
     RETURN
  ENDIF
  IF(FL.LT.0.)THEN
     XL=X1
     XH=X2
  ELSE
     XH=X1
     XL=X2
     SWAP=FL
     FL=FH
     FH=SWAP
  ENDIF
  v2=.5*(X1+X2)
  DXOLD=ABS(X2-X1)
  DX=DXOLD
  CALL FUNC_C2P_RMHD1(v2,F,DF,gam,d,e,s2,b2,sb2,w)
  DO 11 J=1,MAXIT
     IF(((v2-XH)*DF-F)*((v2-XL)*DF-F).GE.0. &
          .OR. ABS(2.*F).GT.ABS(DXOLD*DF) ) THEN
        DXOLD=DX
        DX=0.5*(XH-XL)
        v2=XL+DX
        IF(XL.EQ.v2) THEN
            continue
            RETURN
        ENDIF
     ELSE
        DXOLD=DX
        DX=F/DF
        TEMP=v2
        v2=v2-DX
        IF(TEMP.EQ.v2) THEN
            continue
            RETURN
        ENDIF
     ENDIF
     IF(ABS(DX).LT.XACC) THEN
            !PRINT *,j,DX
            continue
            RETURN
     ENDIF
     CALL FUNC_C2P_RMHD1(v2,F,DF,gam,d,e,s2,b2,sb2,w)
     IF(F.LT.0.) THEN
        XL=v2
        FL=F
     ELSE
        XH=v2
        FH=F
     ENDIF
#ifdef C2PFF  
     IF(ABS(F).LT.XACC2) THEN
         continue
        RETURN
     ENDIF
#else
     IF(F.EQ.0.) THEN
        RETURN
     ENDIF
#endif  
11   CONTINUE
     FAILED = .TRUE.
     v2 = 0. ! assign value even if it fails
     RETURN
  END SUBROUTINE RTSAFE_C2P_RMHD1
   !
    
RECURSIVE SUBROUTINE FUNC_C2P_RMHD1_VECTORF(x,f,df,gam,d,e,s2,b2,sb2,w,VECTORLENGTH)
  !
  ! This is the CONS2PRIM strategy adopted by Del Zanna et al. (2007) A&A, 473, 11-30
  ! and it corresponds to their choice 3 in Section 3.2
  ! 
  IMPLICIT NONE
  INTEGER       :: iter
  INTEGER, INTENT(IN) :: VECTORLENGTH
  LOGICAL(8)       :: MASK(VECTORLENGTH) 
  REAL(8)       :: x(VECTORLENGTH),f(VECTORLENGTH),df(VECTORLENGTH),v2(VECTORLENGTH),rho(VECTORLENGTH),c0(VECTORLENGTH),tmp(VECTORLENGTH),c2(VECTORLENGTH),c3(VECTORLENGTH),dw(VECTORLENGTH),dc2(VECTORLENGTH),dc3(VECTORLENGTH),dlogw(VECTORLENGTH),wb(VECTORLENGTH),vb2(VECTORLENGTH) 
  REAL(8)       :: gam(VECTORLENGTH),d(VECTORLENGTH),e(VECTORLENGTH),s2(VECTORLENGTH),b2(VECTORLENGTH),sb2(VECTORLENGTH),w(VECTORLENGTH),wtmp(VECTORLENGTH)
  REAL(8), PARAMETER :: tolerance = 1e-10, third=1./3.,tolerance1 = 1e-12
  INTENT(IN) :: x,gam,d,e,s2,b2,sb2  
  INTENT(INOUT):: f,df,w
  REAL(8)   :: sqrtdet(VECTORLENGTH),q(VECTORLENGTH),w0(VECTORLENGTH),p(VECTORLENGTH),c2ic3(VECTORLENGTH),c2ic33(VECTORLENGTH)
  REAL, PARAMETER ::  i27=1./27.
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64:: v2, rho, c0, tmp, c2, c3, dw, dc2, dc3, dlogw, wb, vb2  ,MASK,sqrtdet,q,w0,p,c2ic3,c2ic33,wtmp
    !DIR$ ASSUME_ALIGNED x   : 64 
    !DIR$ ASSUME_ALIGNED f   : 64 
    !DIR$ ASSUME_ALIGNED df  : 64 
    !DIR$ ASSUME_ALIGNED gam : 64 
    !DIR$ ASSUME_ALIGNED d   : 64 
    !DIR$ ASSUME_ALIGNED e   : 64 
    !DIR$ ASSUME_ALIGNED s2  : 64 
    !DIR$ ASSUME_ALIGNED b2  : 64 
    !DIR$ ASSUME_ALIGNED sb2 : 64 
    !DIR$ ASSUME_ALIGNED w   : 64  
#else
    !DIR$ ATTRIBUTES ALIGN: 32:: v2, rho, c0, tmp, c2, c3, dw, dc2, dc3, dlogw, wb, vb2  ,MASK,sqrtdet,q,w0,p,c2ic3,c2ic33,wtmp
    !DIR$ ASSUME_ALIGNED x   : 32 
    !DIR$ ASSUME_ALIGNED f   : 32 
    !DIR$ ASSUME_ALIGNED df  : 32 
    !DIR$ ASSUME_ALIGNED gam : 32 
    !DIR$ ASSUME_ALIGNED d   : 32 
    !DIR$ ASSUME_ALIGNED e   : 32 
    !DIR$ ASSUME_ALIGNED s2  : 32 
    !DIR$ ASSUME_ALIGNED b2  : 32 
    !DIR$ ASSUME_ALIGNED sb2 : 32 
    !DIR$ ASSUME_ALIGNED w   : 32  
#endif

  v2(:)  = x(:)
  rho(:) = d(:)*sqrt(1.-v2(:))
  
  c3(:) = 1.0-gam*(1.-v2(:))
  c2(:) = gam(:)*rho(:) + 0.5*b2(:)*(1.0+v2(:))-e(:) 
  c0(:) = -0.5*sb2(:)
  
  ! For every x=v2, we solve a cubic in W of the form: 
  ! c3*W^3+c2*W^2+c0=0 (c3>0, c0<=0)
  ! W=y of the paper. If sb=0 ( meaning c0 = 0), 
  ! w = -c2/c3 > 0 and dw = 0 in the do loop below. 
  ! If -c2/c3 < 0 when sb=0, which makes w=0, 
  ! this is a signature that something was wrong before.
  !
  
#ifdef C2PFF  
  c2ic3(:) = c2(:)/c3(:) 
  IF(ANY( abs ( c0).GT.1.0d-20 ) ) THEN
       WHERE( abs ( c0) < 1.0d-20) 
          w(:) = -c2ic3(:)
       ELSEWHERE 
          c2ic33(:) = c2ic3(:)**3
          q(:) = c0(:)/c3(:) + 2.0*i27*c2ic33(:)
          p(:) = -third*c2ic3(:)**2
          sqrtdet(:) = SQRT(0.25*q(:)**2+i27*p(:)**3)
          w(:) = -third*c2ic3(:) + (-0.5*q(:)+sqrtdet(:))**third + (-0.5*q(:)-sqrtdet(:))**third 
       ENDWHERE
  ELSE
          w(:) = -c2ic3(:)
  ENDIF
  !
  !
  ! then the NR step, solving the 1D problem F=F(w(v2),v2)=0
  dc3(:)   = gam(:)
  dc2(:)   = 0.5 * ( b2(:) - gam(:) * rho(:) / (1.0 - v2(:)))
  dlogw(:) = ( dc3(:) * w(:) + dc2(:) ) / ( 3.0 * c3(:) * w(:) + 2.0 * c2(:))
  wb(:)    = w(:) + b2(:)
  vb2(:)   = sb2(:) / w(:)**2
  f(:)     = wb(:)**2 * v2(:) - ( 2.0 * w(:) + b2(:)) * vb2(:) - s2(:)
  df(:)    = wb(:) * ( wb(:) + 2.0 * dlogw(:) * ( w(:) * v2(:) + vb2(:)))
#else  
        w(:) = max ( - c2(:) / c3(:), ( -c0(:) / c3(:) )**third )
        ! unroll the first six Newton iterations. Flops are for free.  
        tmp(:) = ((3*c3(:)*w(:) + 2*c2(:))*w(:))/(((3*c3(:)*w(:) + 2*c2(:))*w(:))**2+tolerance1)
        dw(:) = -((c3(:)*w(:) + c2(:))*w(:)**2 + c0(:))*tmp(:)
        w(:) = w(:) + dw(:)
        tmp(:) = ((3*c3(:)*w(:) + 2*c2(:))*w(:))/(((3*c3(:)*w(:) + 2*c2(:))*w(:))**2+tolerance1)
        dw(:) = -((c3(:)*w(:) + c2(:))*w(:)**2 + c0(:))*tmp(:) 
        w(:) = w(:) + dw(:)
        tmp(:) = ((3*c3(:)*w(:) + 2*c2(:))*w(:))/(((3*c3(:)*w(:) + 2*c2(:))*w(:))**2+tolerance1)
        dw(:) = -((c3(:)*w(:) + c2(:))*w(:)**2 + c0(:))*tmp(:)
        w(:) = w(:) + dw(:)
        tmp(:) = ((3*c3(:)*w(:) + 2*c2(:))*w(:))/(((3*c3(:)*w(:) + 2*c2(:))*w(:))**2+tolerance1)
        dw(:) = -((c3(:)*w(:) + c2(:))*w(:)**2 + c0(:))*tmp(:)
        w(:) = w(:) + dw(:)
        tmp(:) = ((3*c3(:)*w(:) + 2*c2(:))*w(:))/(((3*c3(:)*w(:) + 2*c2(:))*w(:))**2+tolerance1)
        dw(:) = -((c3(:)*w(:) + c2(:))*w(:)**2 + c0(:))*tmp(:)
        w(:) = w(:) + dw(:)
        tmp(:) = ((3*c3(:)*w(:) + 2*c2(:))*w(:))/(((3*c3(:)*w(:) + 2*c2(:))*w(:))**2+tolerance1)
        dw(:) = -((c3(:)*w(:) + c2(:))*w(:)**2 + c0(:))*tmp(:) 
        w(:) = w(:) + dw(:)
  !
  WHERE ( abs ( c0(:) ) < 1.0d-20) 
        w(:) = -c2(:) / c3(:)
        wtmp(:) = w(:)    ! save results
  ENDWHERE
  !
  IF(MAXVAL(ABS(dw(:)/w(:)))>tolerance) THEN 
      DO iter = 1, 100
            tmp(:) = ((3*c3(:)*w(:) + 2*c2(:))*w(:))/(((3*c3(:)*w(:) + 2*c2(:))*w(:))**2+tolerance1)
            dw(:) = -((c3(:)*w(:) + c2(:))*w(:)**2 + c0(:))*tmp(:)
            IF(MAXVAL(abs(dw(:)/w(:))).LT.tolerance) THEN
                  EXIT 
            ENDIF
            w(:) = w(:) + dw(:)
      ENDDO   
      WHERE ( abs ( c0(:)) < 1.0d-20) 
            w(:) = wtmp(:)  ! get correct values
      ENDWHERE 
  ENDIF          
  !
  dc3(:)   = gam(:)
  dc2(:)   = 0.5 * ( b2(:) - gam(:) * rho(:) / (1.0 - v2(:)))
  dlogw(:) = -( dc3(:) * w(:) + dc2(:) ) / ( 3.0 * c3(:) * w(:) + 2.0 * c2(:) )
  wb(:)    = w(:) + b2(:)
  vb2(:)   = sb2(:) / w(:)**2
  f(:)     = wb(:)**2 * v2(:) - ( 2.0 * w(:) + b2(:) ) * vb2(:) - s2(:) 
  df(:)    = wb(:) * ( wb(:) + 2.0 * dlogw(:) * ( w(:) * v2(:) + vb2(:) ))
#endif  
  !
END SUBROUTINE FUNC_C2P_RMHD1_VECTORF 


RECURSIVE SUBROUTINE RTSAFE_C2P_RMHD1_VECTOR(v2,X1,X2,XACC,gam,d,e,s2,b2,sb2,w,FAILED,VECTORLENGTH)
  !USE VSMM_mod 
  IMPLICIT NONE
  INTEGER               :: VECTORLENGTH 
  INTEGER, PARAMETER    :: MAXIT=400
  INTEGER               :: i,j 
  REAL(8)               :: v2(VECTORLENGTH), XG(VECTORLENGTH), FG(VECTORLENGTH), Vtmp1(VECTORLENGTH), Vtmp2(VECTORLENGTH) 
  REAL(8)               :: X1(VECTORLENGTH), X2(VECTORLENGTH), XACC(VECTORLENGTH), XACC2(VECTORLENGTH), gam(VECTORLENGTH)
  REAL(8)               :: d(VECTORLENGTH),e(VECTORLENGTH),s2(VECTORLENGTH),b2(VECTORLENGTH),sb2(VECTORLENGTH),w(VECTORLENGTH) 
  REAL(8)               :: FL(VECTORLENGTH),FH(VECTORLENGTH),DF(VECTORLENGTH),XH(VECTORLENGTH),XL(VECTORLENGTH),SWAP(VECTORLENGTH)
  REAL(8)               :: DXOLD(VECTORLENGTH),DX(VECTORLENGTH),F(VECTORLENGTH),TEMP(VECTORLENGTH),Xtmp(VECTORLENGTH),Ftmp(VECTORLENGTH),dFtmp(VECTORLENGTH),dXtmp(VECTORLENGTH) 
  LOGICAL(8)               :: FAILED(VECTORLENGTH) 
  INTENT(INOUT) :: v2,X1,X2,gam,w,FAILED
  INTENT(IN)    :: VECTORLENGTH,d,e,s2,b2,sb2
  !
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64:: FL,FH,DF,XH,XL,SWAP,DXOLD,DX,F,TEMP,Xtmp
    !DIR$ ATTRIBUTES ALIGN: 64:: Ftmp,dFtmp,dXtmp,XACC2,XG,FG
    !DIR$ ASSUME_ALIGNED v2 : 64 
    !DIR$ ASSUME_ALIGNED X1 : 64 
    !DIR$ ASSUME_ALIGNED X2 : 64 
    !DIR$ ASSUME_ALIGNED XACC : 64 
    !DIR$ ASSUME_ALIGNED gam  : 64 
    !DIR$ ASSUME_ALIGNED d    : 64 
    !DIR$ ASSUME_ALIGNED e    : 64 
    !DIR$ ASSUME_ALIGNED s2   : 64 
    !DIR$ ASSUME_ALIGNED b2   : 64 
    !DIR$ ASSUME_ALIGNED sb2  : 64 
    !DIR$ ASSUME_ALIGNED w    : 64 
    !DIR$ ASSUME_ALIGNED FAILED    : 64 
#else
    !DIR$ ATTRIBUTES ALIGN: 32:: FL,FH,DF,XH,XL,SWAP,DXOLD,DX,F,TEMP,Xtmp
    !DIR$ ATTRIBUTES ALIGN: 32:: Ftmp,dFtmp,dXtmp,XACC2,XG,FG
    !DIR$ ASSUME_ALIGNED v2   : 32 
    !DIR$ ASSUME_ALIGNED X1   : 32 
    !DIR$ ASSUME_ALIGNED X2   : 32 
    !DIR$ ASSUME_ALIGNED XACC : 32 
    !DIR$ ASSUME_ALIGNED gam  : 32 
    !DIR$ ASSUME_ALIGNED d    : 32 
    !DIR$ ASSUME_ALIGNED e    : 32 
    !DIR$ ASSUME_ALIGNED s2   : 32 
    !DIR$ ASSUME_ALIGNED b2   : 32 
    !DIR$ ASSUME_ALIGNED sb2  : 32 
    !DIR$ ASSUME_ALIGNED w    : 32 
    !DIR$ ASSUME_ALIGNED FAILED    : 32 
#endif 
  !
  XACC2(:) = 1e-60 !XACC(:)**2
  FAILED(:) = .FALSE.
  !!
  !XACCold=XACC
  WHERE(d.gt.0) XACC=MAX(MIN(XACC,XACC*D),1e-22)
  XG(:) = v2(:) ! initial guess
  ! do bracketing with the initial guess: 
  CALL FUNC_C2P_RMHD1_VECTORF(XG,FG,DF,gam,d,e,s2,b2,sb2,w,VECTORLENGTH)
  IF(ANY(ABS(FG(:)).GT.XACC2(:))) THEN
      continue
  ELSE 
     v2(:)=XG(:) 
     RETURN
  ENDIF
  !
  
  !! First try four Newton iterations. 
  !!   
  !DO j = 1 , 5 
  !    CALL FUNC_C2P_RMHD1_VECTORF(v2,F,DF,gam,d,e,s2,b2,sb2,w,VECTORLENGTH)
  !    WHERE( ISNAN(F) )
  !        F = 1e20
  !    ENDWHERE 
  !    IF( ALL(ABS(F).LT.XACC) ) THEN
  !        RETURN
  !    ENDIF 
  !    v2(:) = v2(:) - F(:)/DF(:)  
  !ENDDO 
  !!
  CALL FUNC_C2P_RMHD1_VECTORF(X1,FL,DF,gam,d,e,s2,b2,sb2,w,VECTORLENGTH)
        !IF(ANY(ABS(FL(:)).GT.XACC2(:))) THEN
        !  WHERE(ABS(FL(:)).LT.XACC2(:))  
        !     v2(:) = X1(:) 
        !     X2(:) = X1(:) 
        !  ENDWHERE
        !ELSE 
        !   v2(:)=X1(:)
        !   RETURN
        !ENDIF
  !
  CALL FUNC_C2P_RMHD1_VECTORF(X2,FH,DF,gam,d,e,s2,b2,sb2,w,VECTORLENGTH)
        !IF(ANY(ABS(FH(:)).GT.XACC2(:))) THEN
        !      WHERE(ABS(FH(:)).LT.XACC2(:))  
        !         v2(:) = X2(:) 
        !         X1(:) = X2(:) 
        !      ENDWHERE 
        !ELSE 
        !   v2(:)=X2(:)
        !   RETURN
        !ENDIF
  WHERE(FL(:)*FH(:).GT.0.)  
     FAILED(:) = .TRUE.
     v2(:)     = 0. ! assign value even if it fails     
  ENDWHERE 
  WHERE(FL(:).LT.0.0)  
        WHERE(FG(:).LT.0.0)  
            XL(:) = XG(:) 
            XH(:) = X2(:)
        ELSEWHERE   
            XL(:) = X1(:) 
            XH(:) = XG(:)
        ENDWHERE
  ELSEWHERE 
     XH(:) = X1(:)
     XL(:) = X2(:) 
     SWAP(:) = FL(:) 
     FL(:) = FH(:) 
     FH(:) = SWAP(:) 
        WHERE(FG(:).LT.0.0)  
            XL(:) = XG(:)  
        ELSEWHERE    
            XH(:) = XG(:)
        ENDWHERE
  ENDWHERE 
  !
  ! First try four Newton iterations. 
  !   
  Vtmp1 = ABS(FL)-ABS(FH)
  WHERE(Vtmp1.GT.0.) 
      v2 = XH
  ELSEWHERE
      v2 = XL
  ENDWHERE
  ! Do 5 Newton iterations
  DO j = 1,5 
      CALL FUNC_C2P_RMHD1_VECTORF(v2,F,DF,gam,d,e,s2,b2,sb2,w,VECTORLENGTH)
      WHERE( ISNAN(F) )
          F = 1e20
      ENDWHERE 
      !IF(ALL(ABS(F).LT.XACC2)) THEN
      !    RETURN
      !ENDIF 
      DX(:) = F(:)/DF(:)
      IF( ALL(ABS(DX).LT.XACC) ) THEN 
        CONTINUE 
        RETURN
      ENDIF
      v2(:) = MIN(v2(:) - DX(:),XH)
  ENDDO
  !
  !v2(:)=0.5*(X1(:)+X2(:))
  !
  DXOLD(:)=DX(:) !ABS(X2(:)-X1(:))
  DX(:)=DXOLD(:)
  CALL FUNC_C2P_RMHD1_VECTORF(v2,F,DF,gam,d,e,s2,b2,sb2,w,VECTORLENGTH)  
  DO j = 1, MAXIT 
      Vtmp1 = ((v2-XH)*DF-F)*((v2-XL)*DF-F)
      Vtmp2 = ABS(2.*F)-ABS(DXOLD*DF) 
     WHERE(Vtmp1.GE.0. .OR. Vtmp2.GT.0. )  
        DXOLD(:) = DX(:)
        DX(:) = 0.5*(XH(:)-XL(:)) 
        v2(:)=XL(:)+DX(:) 
     ELSEWHERE 
        DXOLD(:) = DX(:)
        DX(:)=F(:)/DF(:) 
        TEMP(:)=v2(:)
        v2(:)=v2(:)-DX(:) 
     ENDWHERE 
     IF( ALL(ABS(DX).LT.XACC) ) THEN
        CONTINUE 
        RETURN
     ENDIF
     CALL FUNC_C2P_RMHD1_VECTORF(v2,F,DF,gam,d,e,s2,b2,sb2,w,VECTORLENGTH)
     WHERE(F.LT.0.)  
        XL(:) = v2(:)
        FL(:) = F(:) 
     ELSEWHERE 
        XH(:)=v2(:)
        FH(:)=F(:)
     ENDWHERE 
#ifdef C2PFF  
     WHERE(ABS(F).LT.XACC2)
#else
     WHERE(F.EQ.0.)
#endif 
         DF = 0. 
     ENDWHERE
    ! 
  ENDDO 
  !
  WHERE(ABS(DX).GT.1e-10.OR.d.LT.0.OR.e.LT.0.OR.v2.GT.1.0)
      FAILED = .TRUE.
      v2 = 0. ! assign value even if it fails
  ELSEWHERE
      FAILED = .FALSE.
  ENDWHERE
  ! 
END SUBROUTINE RTSAFE_C2P_RMHD1_VECTOR
 
 RECURSIVE SUBROUTINE FUNC_C2P_RMHD1(x,f,df,gam,d,e,s2,b2,sb2,w)
  !
  ! This is the CONS2PRIM strategy adopted by Del Zanna et al. (2007) A&A, 473, 11-30
  ! and it corresponds to their choice 3 in Section 3.2
  !
  IMPLICIT NONE
  REAL, PARAMETER :: third=1./3.
  INTEGER    :: iter
  REAL       :: x,f,df,v2,rho,c0,c2,c3,dw,dc2,dc3,dlogw,wb,vb2
  REAL       :: gam,d,e,s2,b2,sb2,w
  REAL, PARAMETER :: tolerance = 1e-10
  INTENT(IN) :: x,gam,d,e,s2,b2,sb2
  INTENT(OUT):: f,df,w
  REAL :: sqrtdet,q,w0,p,c2ic3,c2ic33
  REAL, PARAMETER ::  i27=1./27.
  v2=x
  rho=d*sqrt(1.-v2)
  c3=1.-gam*(1.-v2)
  c2=gam*rho+.5*b2*(1.+v2)-e
  c0=-0.5*sb2
  ! For every x=v2, we solve a cubic in W of the form: 
  ! c3*W^3+c2*W^2+c0=0 (c3>0, c0<=0)
  ! W=y of the paper. If sb=0 ( meaning c0 = 0), 
  ! w = -c2/c3 > 0 and dw = 0 in the do loop below. 
  ! If -c2/c3 < 0 when sb=0, which makes w=0, 
  ! this is a signature that something was wrong before.
#ifdef C2PFF  
  c2ic3 = c2/c3 
  if ( abs ( c0) < 1.0d-20) then
     w = -c2ic3
  else 
     c2ic33 = c2ic3**3
     q = c0/c3 + 2.0*i27*c2ic33
     p = -third*c2ic3**2
     sqrtdet = SQRT(0.25*q**2+i27*p**3)
     w = -third*c2ic3 + (-0.5*q+sqrtdet)**third + (-0.5*q-sqrtdet)**third 
  endif
  !
  ! then the NR step, solving the 1D problem F=F(w(v2),v2)=0
  dc3   = gam
  dc2   = 0.5 * ( b2 - gam * rho / (1.0 - v2))
  dlogw = ( dc3 * w + dc2 ) / ( 3.0 * c3 * w + 2.0 * c2)
  wb    = w + b2
  vb2   = sb2 / w**2
  f     = wb**2 * v2 - ( 2.0 * w + b2) * vb2 - s2
  df    = wb * ( wb + 2.0 * dlogw * ( w * v2 + vb2))
#else  
  if ( abs ( c0) < 1.0d-20) then
     w = -c2 / c3
  else
     w = max ( - c2 / c3, ( -c0 / c3)**third)
     do iter = 1,100
        dw = -((c3*w + c2)*w**2 + c0)/((3*c3*w + 2*c2)*w)
        if (abs(dw/w).LT.tolerance) THEN
            exit
        ENDIF
        w = w + dw
     end do
  endif
  !
  dc3   = gam
  dc2   = 0.5 * ( b2 - gam * rho / (1.0 - v2))
  !dlogw = ( dc3 * w + dc2 ) / ( 3.0 * c3 * w + 2.0 * c2) 
  dlogw = -( dc3 * w + dc2 ) / ( 3.0 * c3 * w + 2.0 * c2)   ! where this 'minus' comes from????
  wb    = w + b2
  vb2   = sb2 / w**2
  f     = wb**2 * v2 - ( 2.0 * w + b2) * vb2 - s2
  df    = wb * ( wb + 2.0 * dlogw * ( w * v2 + vb2))
#endif  
  !
END SUBROUTINE FUNC_C2P_RMHD1

RECURSIVE SUBROUTINE RTSolve_V2_P_C2P_RMHD(w,p,v2,gam,d,e,s2,b2,sb2,FAILED)
  !***********************************************************
  ! given w=rho*h*W^2, solve the 2d nonlinear system for x=(w,p)
  !*****
  IMPLICIT NONE
  ! arguments
  REAL :: w,p,v2,gam,d,e,s2,b2,sb2
  LOGICAL :: FAILED
  INTENT (INOUT) :: w,p
  INTENT (IN)    :: v2,gam,d,e,s2,b2,sb2
  ! local variables
  INTEGER:: i
  REAL :: dx(2),x(2),f(2),df(2,2),idf(2,2),det,drhodx,dydx,depsdrho,depsdp,eps
  !
  REAL, PARAMETER :: third=1./3.
  INTEGER, PARAMETER :: MaxNewton=100
  REAL, PARAMETER :: tolerance = 1e-10
  FAILED = .FALSE.
  !v2=x0
  ! check the value fo
  x(1) = w
  x(2) = p
  CALL FUNC_C2P_RMHD_F2F3(x,f,df,v2,gam,d,e,s2,b2,sb2)
  det =df(1,1)*df(2,2)-df(1,2)*df(2,1)
  idf(1,1) = df(2,2)/det
  idf(1,2) =-df(1,2)/det
  idf(2,1) =-df(2,1)/det
  idf(2,2) =-df(1,1)/det
  !
  dx(1) = -idf(1,1)*f(1)-idf(1,2)*f(2)     ! = dv2
  dx(2) = -idf(2,1)*f(1)-idf(2,2)*f(2)     ! = dp
  !
  IF(SUM(ABS(dx)).LT.tolerance) THEN
      continue
      RETURN
  ENDIF
  !
  DO i=1,MaxNewton
        x=x+dx
        CALL FUNC_C2P_RMHD_F2F3(x,f,df,v2,gam,d,e,s2,b2,sb2)
        det =df(1,1)*df(2,2)-df(1,2)*df(2,1)
        idf(1,1) = df(2,2)/det
        idf(1,2) =-df(1,2)/det
        idf(2,1) =-df(2,1)/det
        idf(2,2) =-df(1,1)/det
        !
        dx(1) = -idf(1,1)*f(1)-idf(1,2)*f(2)     ! = dv2
        dx(2) = -idf(2,1)*f(1)-idf(2,2)*f(2)     ! = dp
        !
        IF(SUM(ABS(dx)).LT.tolerance) THEN
            w=x(1)
            p=x(2)
			continue
            RETURN
        ENDIF
  ENDDO
  !
  FAILED = .TRUE.
  w=x(1)
  p=x(2)
  !
  continue
  !
END SUBROUTINE RTSolve_V2_P_C2P_RMHD


RECURSIVE SUBROUTINE FUNC_C2P_RMHD_1D_2D(x,f,df,gam,d,e,s2,b2,sb2,w)
  !
  ! This is the CONS2PRIM strategy adopted by Del Zanna et al. (2007) A&A, 473, 11-30
  ! and it corresponds to their choice 3 in Section 3.2
  !
   USE EOS_mod
  IMPLICIT NONE
  REAL, PARAMETER :: third=1./3.
  INTEGER    :: iter
  REAL       :: x,f,df,v2,rho,c0,c2,c3,dw,dc2,dc3,dlogw,wb,vb2
  REAL       :: gam,d,e,s2,b2,sb2,w,dfdx,ydfdy,dydx,ydfdp,dpdx,dydp
  REAL       :: eps,depsdrho,depsdp,den
  REAL, PARAMETER :: tolerance = 1e-14
  INTENT(IN) :: x,gam,d,e,s2,b2,sb2
  INTENT(OUT):: f,df,w 
  REAL :: sqrtdet,q,w0,p,c2ic3,c2ic33
  REAL, PARAMETER ::  i27=1./27.
  LOGICAL :: FAILED
  v2=x
  rho=d*sqrt(1.-v2) 
  c3=1.-gam*(1.-v2)
  c2=gam*rho+.5*b2*(1.+v2)-e
  c0=-0.5*sb2
  
  ! For every x=v2, we solve a cubic in W of the form:
  ! c3*W^3+c2*W^2+c0=0 (c3>0, c0<=0)
  ! W=y of the paper. If sb=0 ( meaning c0 = 0),
  ! w = -c2/c3 > 0 and dw = 0 in the do loop below.
  ! If -c2/c3 < 0 when sb=0, which makes w=0,
  ! this is a signature that something was wrong before.
  !
  ! given w=rho*h*W^2, solve the 2d nonlinear system for x=(w,p) 
  CALL RTSolve_V2_P_C2P_RMHD(w,p,v2,gam,d,e,s2,b2,sb2,FAILED)
  ! then compute f and df for the 1D NR iteration F=F(w(v2),p(v2),v2)=0
  !
  dc3   = gam
  dc2   = 0.5 * ( b2 - gam * rho / (1.0 - v2))
  dlogw = -( dc3 * w + dc2 ) / ( 3.0 * c3 * w + 2.0 * c2)
  wb    = w + b2
  vb2   = sb2 / w**2
  !
  f     = wb**2 * v2 - ( 2.0 * w + b2) * vb2 - s2
  !
  dfdx = wb*wb
  ydfdy = 2.0*wb*(w*v2+vb2)
  dlogw = ( dc3 * w + dc2 ) / ( 3.0 * c3 * w + 2.0 * c2)
  !
  CALL  EOS_eps_depsdrho_depsdp(eps,depsdrho,depsdp,rho,p)
  !
  ydfdp = w*D*D/rho*depsdp + w/(1.0-x)
  den = (1.0 + rho*depsdp)
  !dp/dx = 
  
  df    = wb * ( wb + 2.0 * dlogw * ( w * v2 + vb2)) 
  !
  ! compute df = df/dx + df/dy*(dy/dx + dy/dp*dp/dx) + df/dp*dp/dx .... ???
  !  df/dp = 0 ...
END SUBROUTINE FUNC_C2P_RMHD_1D_2D



RECURSIVE SUBROUTINE FUNC_C2P_RMHD_F2F3(x,f,df,v2,gam,d,e,s2,b2,sb2)
  USE EOS_mod
  IMPLICIT NONE
  ! arguments
  REAL :: x(2),f(2),df(2,2)
  REAL :: v2,gam,d,e,s2,b2,sb2 
  INTENT (INOUT) :: x,f,df
  INTENT (IN)    :: v2,gam,d,e,s2,b2,sb2
  ! local variables
  INTEGER    :: i
  REAL :: det,drhodx,dydx,depsdrho,depsdp,eps ,w
  REAL, PARAMETER :: third=1./3.
  !
  REAL :: rho,c0,c2,c3,dw,dc2,dc3,dlogw,wb,vb2 
  ! ... we actually know the root of a cubic!!!!
  REAL :: sqrtdet,q,w0,p,c2ic3,c2ic33
  REAL, PARAMETER ::  i27=1./27. 
  w=x(1)
  p=x(2)
  rho=d*sqrt(1.-v2) 
  c3=1. !-gam*(1.-v2)
  c2=.5*b2*(1.+v2)-e-p
  c0=-0.5*sb2
  !
  f(1) = c3*w**3 + c2*w**2 + c0
  !
  CALL EOS_eps_depsdrho_depsdp(eps,depsdrho,depsdp,rho,p) 
  !
  f(2) = p + eps*rho + rho - (1.-v2)*w 
  !
  ! compute the jacobian: (d/dw, d/dp)[ (/ f(1), f(2) /) ]
  dc3   = 1.0
  dc2   = 0.5 * ( b2 - gam * rho / (1.0 - v2))
  dlogw = -( dc3 * w + dc2 ) / ( 3.0 * c3 * w + 2.0 * c2)
  wb    = w + b2
  vb2   = sb2 / w**2
  f     = wb**2 * v2 - ( 2.0 * w + b2) * vb2 - s2
  df    = wb * ( wb + 2.0 * dlogw * ( w * v2 + vb2))
  
  df(1,1) = (3*c3*w + 2*c2)*w
  df(1,2) = -w**2
  !
  df(2,1) = -(1-v2)
  df(2,2) = 1.0+rho*depsdp
  !
  continue
  !
END SUBROUTINE FUNC_C2P_RMHD_F2F3


 
RECURSIVE SUBROUTINE PDEVarNameGRMHD(MyName) 
  USE ISO_C_BINDING 
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
  CHARACTER(LEN=20):: MyName(nVar)
  INTEGER			:: ind
	!
    MyName(1)  = 'rho'      
    MyName(2)  = 'u'        
    MyName(3)  = 'v'        
    MyName(4)  = 'w'        
    MyName(5)  = 'p'        
    MyName(6)  = 'Bx'       
    MyName(7)  = 'By'       
    MyName(8)  = 'Bz'       
    MyName(9)  = 'psi'      
    MyName(10) = 'lapse'    
    MyName(11) = 'shift_x'  
    MyName(12) = 'shift_y'  
    MyName(13) = 'shift_z'  
    MyName(14) = 'gammaij_1'
    MyName(15) = 'gammaij_2'
    MyName(16) = 'gammaij_3'
    MyName(17) = 'gammaij_4'
    MyName(18) = 'gammaij_5'
    MyName(19) = 'gammaij_6'
#ifdef VECTOR
    PRINT *,'nVar = ', nVar
#endif
	!
END SUBROUTINE PDEVarNameGRMHD   


	
	
    

RECURSIVE SUBROUTINE PDEAuxNameGRMHD(AuxName) 
  IMPLICIT NONE  
#if defined(SPHERICAL) || defined(CYLINDRICAL)
  INTEGER, PARAMETER :: nAux = 9                           ! The number of variables of the PDE system 
#else
  INTEGER, PARAMETER :: nAux = 0                           ! The number of variables of the PDE system 
#endif   
  CHARACTER(LEN=10):: AuxName(nAux)
  !
#if defined(SPHERICAL) || defined(CYLINDRICAL)
  AuxName(1) = 'x1' //C_NULL_CHAR
  AuxName(2) = 'x2' //C_NULL_CHAR
  AuxName(3) = 'x3' //C_NULL_CHAR
  AuxName(4) = 'V_x'//C_NULL_CHAR
  AuxName(5) = 'V_y'//C_NULL_CHAR
  AuxName(6) = 'V_z'//C_NULL_CHAR
  AuxName(7) = 'B_x'//C_NULL_CHAR
  AuxName(8) = 'B_y'//C_NULL_CHAR
  AuxName(9) = 'B_z'//C_NULL_CHAR
#endif 
  !
 END SUBROUTINE PDEAuxNameGRMHD     
  
      

RECURSIVE SUBROUTINE PDEAuxVarGRMHD(aux,V,x) 
  USE Parameters, ONLY : EQN, nVar, nAux, aom, d,nDim
  IMPLICIT NONE
  REAL :: aux(nAux), V(nVar) 
  REAL :: time, x(d)
  REAL, PARAMETER :: epsilon = 1e-14  
  ! Local Variables 
  INTEGER :: i,ii,jj,kk,ll,mm,nn,qq,iErr
  INTEGER :: itemp(nDim)
  REAL :: rho,vx,vy,vz,p,bx,by,bz,ex,ey,ez,pl
  REAL :: v2,b2,e2,lf,LF2,w,ww,wwx,wwy,wwz,uem,gamma1
  REAL :: irho,ps,pg,rho0,k0,kappa,sigma,psi,phi,falpha,cs,cl,c0
  REAL :: uu,vv,Bx2,By2,Bz2,B28P,vB,ru2,rv2,rw2,a2,pMag,hh
  REAL :: lapse, gp, gm, rhoE, dcs,dc0, rhoA, Y, eps, xi, c, tau 
  REAL :: g_contr(3,3), g_cov(3,3)
  REAL :: shift(3), v_contr(3), v_cov(3), Bv(3), Ev(3), vxB(3), vxB_cov(3), ExB(3), J(3), B_cov(3)
  REAL :: A(3,3),A_contr(3,3),iA(3,3),iA_contr(3,3),C2S(3,3),iC2S(3,3),C2S_contr(3,3),iC2S_contr(3,3),detC2S, devG(3,3), GT(3,3), Id(3,3), detA, check, TT(3,3), Atilde(3,3)
  REAL :: Fv(nVar), Gv(nVar), Hv(nVar), Qx(nVar), Qy(nVar), Qz(nVar), L(nVar)  
  REAL :: uux,vvx,Tx
  REAL :: uuy,vvy,Ty
  REAL :: uuz,vvz,Tz
  REAL, EXTERNAL :: PDERefVar
  REAL :: gradQT(d,nVar) 
  REAL :: divV,divV23,mu,icv,T 
  REAL :: iRho2,iRho3,w1,w2,w3 
  REAL :: dTdW1,dTdW2,dTdW3,dTdW4
  REAL :: vorx,vory,vorz,Vel,Vor, DetGamma, TrA
  REAL :: lambda(nDim), GradV(nDim,nDim), Sym(nDim,nDim), ASym(nDim,nDim), Matrix(nDim,nDim) , RR(nDim,nDim), temp1(nDim), temp2(nDim)
  REAL :: K(3,3), Ricci(3,3), R, traceK, KK2, H   
  REAL :: DD(3,3,3), PP(3), GG(3), dP(3,3)
  REAL :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10, s11, s12, s13 
  REAL :: dDD(3,3,3,3), Christoffel(3,3,3), ChristoffelNC(3,3,3), dgup(3,3,3), Riemann(3,3,3,3), dChristoffel(3,3,3,3)
  REAL :: Christoffel_diff(3,3,3,3), Ham, Mom(3), dK(3,3,3), dAtilde(3,3,3), dtraceK(3)
  REAL :: dg_cov(3,3,3), g_covx(3,3), g_covy(3,3), g_covz(3,3), det 
  REAL :: a_st(3,0:3), vtr(3), g_cov_st(0:3,0:3), g_contr_st(0:3,0:3), gAB(3,3), gAB_contr(3,3), kAB(3,3), kAB_contr(3,3) 
  REAL :: piAB_st(0:3,0:3), shift_cov(3), piAB1(3,3), piAB2(3,3), piAB(3,3), pi_mix_st(0:3,0:3), I1, I2, f1, f2 
  REAL :: eta_cov(3,3), eta_mix(3,3), detgAB ,eta
  REAL :: Aex(3,3), dAex(3,3,3), AMix(3,3), Aup(3,3), Theta, dTheta(3), Ghat(3), dGhat(3,3), traceA, AA(3), dAA(3,3), b(3), BB(3,3), dBB(3,3,3), dphi(3), dPP(3,3), detkAB 
  REAL :: Kex(3,3), Gtilde(3), Christoffel_tilde(3,3,3), Z(3), Zup(3), u_contr(0:3), u_cov(0:3), gmunu(0:3,0:3), hmunu(0:3,0:3), devG_st(0:3,0:3), pi_mix(0:3,0:3), pi_st(0:3,0:3)      
  ! 
  REAL, PARAMETER :: iPi  = .31830988618379
  REAL, PARAMETER :: iPi4 = .79577471545948e-1
  REAL, PARAMETER :: iPi8 = .39788735772974e-1
  ! 
  !    f = Q  
  !    g = 0. 
  !    h = 0. 
  !    RETURN 
  !
#if defined(SPHERICAL) || defined(CYLINDRICAL)
  aux = 0
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
  CALL Curved2Cartesian(aux(1:3),x)
#ifdef SPHERICAL
  CALL Cart2SphMatrix_cov(A,iA,x)
#else
  !phi = x(3) 
  CALL Cart2CylMatrix_cov(A,iA,x)
  !aux(1:3) = (/ x(1)*DCOS(phi) , x(2), x(1)*DSIN(phi)  /)  
#endif
  ! 
  IF(X(1).EQ.0) THEN    ! r=0
      !
      aux(4:6) = 0.0
      aux(7:9) = 0.0
      !
  ELSE
      !
      !IF(MAXVAL(ABS(A)).LT.1e-14) THEN
      !      print *, 'FATAL ERROR: MAXVAL(ABS(A))'
      !      WRITE(*,*) A(1,1:3)
      !      WRITE(*,*) A(2,1:3)
      !      WRITE(*,*) A(3,1:3)
      !      WRITE(*,*) x
      !      STOP
      !ENDIF 
      !CALL MatrixInverse3x3(A,iA,detA)      ! iA = Sph2CartMatrix_cov    = Cart2SphMatrix_contr
      A_contr(:,:) = TRANSPOSE(iA(:,:))     ! Cart2SphMatrix_contr = TRANSPOSE( Sph2CartMatrix_cov )
      iA_contr(:,:) = TRANSPOSE(A(:,:))    ! Sph2CartMatrix_contr = TRANSPOSE( Cart2SphMatrix_cov )
      !  
      aux(4:6) = MATMUL(iA,V(2:4))  ! Spherical to Cartesian: Velocity is covariant
      ! 
      ! IF we plot COVARIANT B
      B_cov = MATMUL(g_cov,V(6:8)) !  Magnetic field is contravariant: 
      aux(7:9) = MATMUL(iA,B_cov)  !  Spherical to Cartesian:  Magnetic field is controvariant (B_cov is covariant) 
  ENDIF 
  !  
#endif 


  !
END SUBROUTINE PDEAuxVarGRMHD     

	
	
	
END MODULE GRMHD_Mod