! The Con2Prim and Prim2Con routines for MHD.
! Should be merged with SRHD's.

SUBROUTINE PDEPrim2Cons(Q,V)
  USE Parameters, ONLY: nVar, nDim 
  IMPLICIT NONE
  ! Argument list declaration
  REAL :: Q(nVar), V(nVar)
  INTENT(IN)  :: V
  INTENT(OUT) :: Q 
  ! Local variable declaration
  REAL :: Prim(nVar), Buf(nVar)
  REAL :: rho,vx,vy,vz,p,bx,by,bz,ex,ey,ez,cs,c0
  REAL :: v2,b2,e2,lf,w,ww,uem,gamma1
  REAL :: lapse, gp, gm, dcs, dc0, eel 
  REAL :: g_contr(3,3), g_cov(3,3)
  REAl :: shift(3), vf(3), vf_cov(3), Ev(3), Bv(3), ExB(3)
  REAL :: A(3,3), devG(3,3), G(3,3), temp(3,3), Id(3,3), detA, eh, S, evv, T, falpha  
  REAL :: alphas, alphal, rhos, rhol, us, vs, ul, vl, es, el
  REAL :: EE(3),BB(3),vv(3),detvEB, vxE(3), vxB(3)  
  INTEGER :: i

  !
  Q(:)=V(:)  
  Q(7:9)=V(7:9)*Q(13)
END SUBROUTINE PDEPrim2Cons


SUBROUTINE PDECons2Prim(V,Q)
  USE Parameters, ONLY: nVar, nDim, epsilon1
  IMPLICIT NONE
  !--------------------------------------------!
  ! Argument list declaration
  REAL :: Q(nVar), V(nVar)
  INTENT(IN)  :: Q 
  INTENT(OUT) :: V 
  ! Local variable declaration
  REAL        :: iRho
  REAL, PARAMETER :: epsilon = 1e-14
  INTEGER     :: i, iter, indx_nr(2)
  REAL        :: dr, sx, sy, sz, bx, by, bz
  REAL        :: v2, sb, den, vb, zeta, cs    
  REAL        :: gamma1, G1, G12, x1, x2, eps
  REAL        :: rho, vx, vy, vz, p
  REAL        :: lf, lf2, lf3, lf4, gam,e,s2,b2,e2,sb2,w,ww
  REAL        :: RTSAFE_C2P_RMHD1, RTSAFE_C2P_RMHD2, RTSAFE_C2P_RHD1, RTSAFE_C2P_RHD2
  REAL        :: k(3), B(3), vel(3)
  REAL        :: p1,q1,dd,phi,temp1, H, k2, kB, T2
  REAL        :: epsilon0, mu0, f, df, drho
  REAL        :: iPi,B28P, dcs, dc0
  REAL        :: dv2, c0, dw
  REAL        :: ri,qi,kappa,z,kappa_max
  REAL        :: ZBRENT_C2P_RHD2, ZBRENT_LF_POLY
  REAL        :: lapse, gp, gm
  REAL        :: g_contr(3,3), g_cov(3,3)
  REAL        :: shift(3), sm(3), sm_cov(3), vf_cov(3), vf(3), Ev(3), Bv(3)
  REAL        :: Qtest(14), dQ(14)
  REAL        :: ExB(3),ExB_up(3), U2e(3), S_up(3)  
  REAL        :: LFsolutions(4), cLF(5), C2, C3 
  COMPLEX     :: LFsolutionsC(4)
  REAL        :: xi(2), alpha_nr(2,2), beta_nr(2), d_nr, S  
  REAL        :: A(3,3), devG(3,3), G(3,3), tempp(3,3), Id(3,3), detA, ehh, evv
  REAL        :: EE(3), vv(3), BB(3), eel, detvEB, vxE(3), vxB(3)  
  REAL        :: vs,vl,us,ul,rhos,rhol,rhoeh,rhoe,alphas,alphal
  REAL        :: xx(9), temp(3) 
  LOGICAL     :: FAILED
  REAL, PARAMETER    :: tol = 1e-8, third=1.0/3.0, p_floor = 1.0e-5, rho_floor = 1.0e-4
  !
    V=Q
    V(7:9)=Q(7:9)*Q(13)/(Q(13)**2+epsilon1)
END SUBROUTINE PDECons2Prim

