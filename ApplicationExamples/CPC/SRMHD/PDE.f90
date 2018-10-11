!!!
!!! A Special Relativistic ideal Magnetohydrodynamics ExaHyPE Kernel.
!!! Based on files provided by Olindo, adopted by Sven.
!!! Based on srhd3dfortran/PDE.f90.
!!!

! To avoid needless copying, these preprocessor rules allow
! neat naming conventions
#define rho       V(1)
#define vx        V(2)
#define vy        V(3)
#define vz        V(4)
#define p         V(5)
#define bx        V(6)
#define by        V(7)
#define bz        V(8)
#define cleaning  V(9)

SUBROUTINE PDEEigenvalues(Lambda,Q,nv) 
  USE Parameters, ONLY : nVar, nDim , gamma
  USE, INTRINSIC :: ISO_C_BINDING 
  IMPLICIT NONE 
  ! Argument list  
  REAL, INTENT(IN)  :: Q(nVar), nv(nDim)  
  REAL, INTENT(OUT) :: Lambda(nVar)  
  ! Local variables  
  INTEGER :: iErr
  REAL :: cs2, cs, c0, v2, w, gamma1, vn, den, u, c
  REAL :: ex, ey, ez, b2, e2, a2, ca2, vf1, vf2, vel(nDim)
  REAL :: V(nVar)
  
  ! MUST be 1 as Divergence cleaning factor is 1.
  ! they travel with speed 1
  Lambda = 1
  RETURN


  ! These are not the exact eigenvalues, instead of Lambda(1..9)
  ! we compute only two eigenvalues: Approximate magnetosonics

  CALL PDECons2Prim(V,Q,iErr)

  ex     = - (vy*bz - vz*by)
  ey     = - (vz*bx - vx*bz)
  ez     = - (vx*by - vy*bx)
  gamma1 = gamma/(gamma-1.0)
  b2     = bx*bx + by*by + bz*bz
  v2     = vx*vx + vy*vy + vz*vz
  e2     = ex*ex + ey*ey + ez*ez
  w      = rho + gamma1*p
  vn     = vx*nv(1) + vy*nv(2) + vz*nv(3)

  cs2    = gamma * p / w
  ca2    = (b2 - e2) / ( w + b2 - e2)
  a2     = cs2 + ca2 - cs2 * ca2
  den    = 1.0/( 1.0 - v2 * a2)
  vf1    = den * ( 1.0 - a2) * vn
  vf2    = den * sqrt(a2 * ( 1.0 - v2) * ((1.0 - v2 * a2) - (1.0 - a2) * vn**2))

! Svens version:
!  ex     = - (vy*bz - vz*by)
!  ey     = - (vz*bx - vx*bz)
!  ez     = - (vx*by - vy*bx)
!  gamma1 = gamma/(gamma-1.0)
!  b2     = bx*bx + by*by + bz*bz
!  !v2     = vx*vx + vy*vy + vz*vz
!  e2     = ex*ex + ey*ey + ez*ez
!  w      = rho + gamma1*p
!  !vn     = vx*nv(1) + vy*nv(2) + vz*nv(3)
!
!  ! the velocity vector is dimension agnostic
!  vel    = V(2:2+nDim-1)
!  v2     = SUM(vel**2)
!  vn     = SUM(vel * nv)
!
!  cs2    = gamma * p / w
!  ca2    = (b2 - e2) / ( w + b2 - e2)
!  a2     = cs2 + ca2 - cs2 * ca2
!  den    = 1.0/( 1.0 - v2 * a2)
!  vf1    = den * ( 1.0 - a2) * vn
!  vf2    = den * sqrt(a2 * ( 1.0 - v2) * ((1.0 - v2 * a2) - (1.0 - a2) * vn**2))

  Lambda(1) = vf1 + vf2
  Lambda(2) = vf1 - vf2
  Lambda(3) = 0.0
  Lambda(4) = 0.0
  Lambda(5) = 0.0
  Lambda(6) = 0.0
  Lambda(7) = 0.0
  Lambda(8) = 0.0
  Lambda(9) = 0.0

END SUBROUTINE PDEEigenvalues
 


SUBROUTINE PDEFlux(Fx,Fy,Fz,Q) 
  USE Parameters, ONLY : nVar, nDim, gamma, DivCleaning_a
  USE, INTRINSIC :: ISO_C_BINDING 
  IMPLICIT NONE 
  ! Argument list  
  REAL, INTENT(IN)  :: Q(nVar) 
  REAL, INTENT(OUT) :: Fx(nVar), Fy(nvar), Fz(nVar)
  ! Local variables  
  REAL :: v2,lf,w,ww,wwx,wwy,wwz,gamma1
  REAL :: ex,ey,ez,b2,e2,uem
  REAL :: V(nVar)
  INTEGER :: iErr
  
  CALL PDECons2Prim(V,Q,iErr)
  gamma1 = gamma/(gamma-1.0)

  ex     = - (vy*bz - vz*by)
  ey     = - (vz*bx - vx*bz)
  ez     = - (vx*by - vy*bx)

  v2     = vx**2 + vy**2 + vz**2
  b2     = bx**2 + by**2 + bz**2
  e2     = ex**2 + ey**2 + ez**2
  lf     = 1.0/sqrt(1.0 - v2)
  w      = rho + gamma1*p
  ww     = w*lf**2
  uem    = 0.5*(b2 + e2)
  wwx    = ww*vx
  wwy    = ww*vy
  wwz    = ww*vz

  !
  Fx(1) = vx*rho*lf
  Fx(2) = wwx*vx - bx*bx - ex*ex + p + uem
  Fx(3) = wwx*vy - bx*by - ex*ey
  Fx(4) = wwx*vz - bx*bz - ex*ez 
  Fx(5) = wwx + (ey*bz - ez*by) 
  Fx(6) = V(9)
  Fx(7) = -ez
  Fx(8) = ey  
  Fx(9) = DivCleaning_a**2*bx
  !
  Fy(1) = vy*rho*lf
  Fy(2) = wwy*vx - by*bx - ey*ex 
  Fy(3) = wwy*vy - by*by - ey*ey + p + uem
  Fy(4) = wwy*vz - by*bz - ey*ez 
  Fy(5) = wwy + (ez*bx - ex*bz) 
  Fy(6) = ez 
  Fy(7) = V(9) 
  Fy(8) = -ex   
  Fy(9) = DivCleaning_a**2*by
  !
  IF ( nDim > 2 ) then
    Fz(1) = vz*rho*lf
    Fz(2) = wwz*vx - bz*bx - ez*ex 
    Fz(3) = wwz*vy - bz*by - ez*ey 
    Fz(4) = wwz*vz - bz*bz - ez*ez + p + uem
    Fz(5) = wwz + (ex*by - ey*bx) 
    Fz(6) = -ey  
    Fz(7) = ex   
    Fz(8) = V(9)   
    Fz(9) = DivCleaning_a**2*bz
  ENDIF
  !
END SUBROUTINE PDEFlux 
 
SUBROUTINE PDESource(S,Q) 
  USE Parameters, ONLY : nVar, nDim, gamma, DivCleaning_a
  USE, INTRINSIC :: ISO_C_BINDING 
  IMPLICIT NONE 
  ! Argument list  
  REAL, INTENT(IN)  :: Q(nVar) 
  REAL, INTENT(OUT) :: S(nVar)
  ! Local variables  
  
  S = 0.0

  S(9) = - DivCleaning_a * Q(9)

END SUBROUTINE PDESource
