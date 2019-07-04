! The Con2Prim and Prim2Con routines for MHD.
! Should be merged with SRHD's.

SUBROUTINE PDEPrim2Cons(Q,V)
  USE MainVariables, ONLY: nVar, nDim, EQN
  IMPLICIT NONE
  ! Argument list declaration
  REAL :: Q(nVar), V(nVar)
  INTENT(IN)  :: V
  INTENT(OUT) :: Q 
  ! Local variable declaration
   Q(1) = V(1)*V(6)      
   Q(2) = V(2)*V(1)*V(6)
   Q(3) = V(3)*V(1)*V(6)
   Q(4) = V(4)*V(1)*V(6)
   Q(5) = V(6)*(V(5)/(EQN%gamma-1.0)+V(1)*(V(2)**2+V(3)**2+V(4)**2)/2.0)
   Q(6:9) = V(6:9)
END SUBROUTINE PDEPrim2Cons

SUBROUTINE PDECons2Prim(V,Q)
  USE MainVariables, ONLY: nVar, nDim, EQN
  IMPLICIT NONE
  !--------------------------------------------!
  ! Argument list declaration
  REAL :: Q(nVar), V(nVar)
  INTENT(IN)  :: Q 
  INTENT(OUT) :: V 
  ! Local variables
  REAL :: ialpha
  
  ialpha = Q(6)/(Q(6)**2 + 1e-10) 
  V(1)   = Q(1)*ialpha 
  V(2:4) = Q(2:4)*Q(1)/(Q(1)**2 + 1e-10)  
  V(5) = (EQN%gamma-1.)*( Q(5)*ialpha - 0.5*V(1)*( V(2)**2 + V(3)**2 +V(4)**2 ) ) 
  V(6:9) = Q(6:9) 
END SUBROUTINE PDECons2Prim

