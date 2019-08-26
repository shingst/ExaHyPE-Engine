! Tools.f90

RECURSIVE SUBROUTINE PDECons2PrimFix(V,Q,iErr) 
  USE Parameters, ONLY: nVar  
  IMPLICIT NONE     
  REAL :: V(nVar)
  REAL, INTENT(IN) :: Q(nVar)
  INTEGER          :: iErr
  REAL :: Qaux(nVar)
  Qaux(1:nVar) = Q(1:nVar)
  CALL PDECons2Prim(V,Qaux,iErr)
  !
END SUBROUTINE PDECons2PrimFix
    

RECURSIVE SUBROUTINE getNumericalSolution(V,Q) 
  USE Parameters, ONLY: nVar  
  IMPLICIT NONE     
  REAL				:: V(nVar), Q(nVar)
  INTEGER 			:: iErr
  !
  CALL PDECons2Prim(V,Q,iErr)
  !
END SUBROUTINE getNumericalSolution

RECURSIVE SUBROUTINE getExactSolution(x,timeStamp,V) 
  USE Parameters, ONLY: nAux,nVar,nDim
  USE GRMHD_Mod  
  IMPLICIT NONE     
  REAL				:: V(nVar), Q(nVar), x(nDim), timeStamp
  INTEGER 			:: iErr
  !
  call InitialData(x, timeStamp, Q)
  CALL PDECons2Prim(V,Q,iErr)
  ! 
END SUBROUTINE getExactSolution


