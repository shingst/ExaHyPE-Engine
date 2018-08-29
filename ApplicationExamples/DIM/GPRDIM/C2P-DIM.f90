! The Con2Prim and Prim2Con routines for MHD.
! Should be merged with SRHD's.

SUBROUTINE PDEPrim2Cons(Q,V)
  USE Parameters, ONLY: nVar, nDim , epsilon1
  IMPLICIT NONE
  ! Argument list declaration
  REAL :: Q(nVar), V(nVar)
  INTENT(IN)  :: V
  INTENT(OUT) :: Q 
  ! Local variable declaration
	Q=V
    Q(2:4)=V(2:4)*V(14)*V(1)
END SUBROUTINE PDEPrim2Cons


SUBROUTINE PDECons2Prim(V,Q)
  USE Parameters, ONLY: nVar, nDim, epsilon1
  IMPLICIT NONE
  !--------------------------------------------!
  ! Argument list declaration
  REAL :: Q(nVar), V(nVar),alpha,ialpha
  INTENT(IN)  :: Q 
  INTENT(OUT) :: V 
  ! Local variable declaration
    V=Q
    alpha=Q(14) ! GPR+LE+DI
    if(epsilon1<0) then
        ialpha=1.0/alpha
    else
        ialpha=alpha/(alpha**2+epsilon1*(1-alpha))
    end if
    V(2:4)=Q(2:4)/Q(1)*ialpha

END SUBROUTINE PDECons2Prim

