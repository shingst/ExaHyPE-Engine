! DIM Initial Data
#ifndef HEADER_INITALDATA
#define HEADER_INITALDATA
!#include "MainVariables.f90"
!#include "expintegrator_type.f90"

RECURSIVE SUBROUTINE InitParameters(STRLEN,PARSETUP) 
	USE MainVariables, ONLY: nVar , nDim, EQN, ICType
	IMPLICIT NONE  
	integer          :: STRLEN
	character(len=STRLEN) :: PARSETUP

	ICType=trim(parsetup)
	print *, "****************************************************************"
	print *, 'Chosen setup: 	',ICType
	print *, "****************************************************************"
	!stop
	
	select case(ICType)
		case('NLOPRUPTURE')

	end select
END SUBROUTINE InitParameters

RECURSIVE SUBROUTINE InitialData(xGP, tGp, Q)
	USE, INTRINSIC :: ISO_C_BINDING
	USE MainVariables, ONLY : nVar, nDim, EQN, ICType
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: xGP(3), tGp        ! 
	REAL, INTENT(OUT)              :: Q(nVar)        ! 
	REAL :: up(nVar)

	
	select case(ICType)
		case('NLOPRUPTURE')

	CASE DEFAULT
		PRINT *, 'NOT IMPLEMENTED'
		STOP
	END SELECT
	CALL PDEPrim2Cons(Q,up)
END SUBROUTINE InitialData

RECURSIVE SUBROUTINE PDElimitervalue(limiter_value,xx,numberOfObservables, observablesMin, observablesMax)
	USE, INTRINSIC :: ISO_C_BINDING
	USE MainVariables, ONLY : nVar, nDim
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: xx(3)        ! 
	INTEGER, INTENT(IN)					:: numberOfObservables
	INTEGER, INTENT(OUT)              :: limiter_value        !
	REAL, INTENT(IN)					:: observablesMin(numberOfObservables), observablesMax(numberOfObservables)
	real	:: rr	

	limiter_value=0

   IF((observablesMax(1)<0.999 .and. observablesMin(1)>0.001) .or. observablesMax(1)>1.001) THEN  ! was -0.001   .or. levelmin<0.001
		limiter_value = 1 
   ENDIF 
   
   IF(observablesMax(2)>1.e-5  .and. observablesMax(1)>1.e-4) THEN 
       limiter_value = 1
   ENDIF 
   
   IF(observablesMax(3)>0.5 ) THEN 
       limiter_value = 1
   ENDIF 
   
   !limiter_value = 1
END SUBROUTINE PDElimitervalue

#endif
