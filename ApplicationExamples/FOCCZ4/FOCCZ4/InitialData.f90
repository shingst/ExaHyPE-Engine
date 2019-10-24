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
	ICType   = 'CCZ4MinkowskiSrc'
	select case(ICType)
		case('CCZ4MinkowskiSrc')
			EQN%CCZ4GLMc0   = 1.5   ! 0.1      
			EQN%CCZ4GLMc    = 1.2   ! 2.0    
			EQN%CCZ4GLMd    = 2.0   ! 1.0     
			EQN%CCZ4GLMepsA = 1.0   ! 5. 
			EQN%CCZ4GLMepsP = 1.0   ! 5.  
			EQN%CCZ4GLMepsD = 1.0   ! 0.1 
			!
			EQN%CCZ4k1  = 0.0  
			EQN%CCZ4k2  = 0.0 
			EQN%CCZ4k3  = 0.0 
			EQN%CCZ4eta = 0.0 
			EQN%CCZ4f   = 0.0 
			EQN%CCZ4g   = 0.0 
			EQN%CCZ4xi  = 0.0 
			EQN%CCZ4e   = 2.0 
			EQN%CCZ4c   = 0.0 
			EQN%CCZ4mu  = 0.0 
			EQN%CCZ4ds  = 1.0 
			EQN%CCZ4sk  = 0.0
			EQN%CCZ4bs   = 0.0      ! set bs=1 if you want to activate the shift convection for beta, b and B (standard CCZ4 formulation). set it to bs=0 to switch off shift convection for those quantities 
			EQN%CCZ4LapseType   = 0 ! harmonic lapse 
			EQN%EinsteinAutoAux = 0 
	end select
END SUBROUTINE InitParameters

RECURSIVE SUBROUTINE InitialData(xGP, tGp, Q)
	USE, INTRINSIC :: ISO_C_BINDING
	USE MainVariables, ONLY : nVar, nDim, EQN, ICType
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: xGP(3), tGp        ! 
	REAL, INTENT(OUT)              :: Q(nVar)        ! 
	REAL :: up(nVar),V0(nVar),r

	
	select case(ICType)
		case('CCZ4MinkowskiSrc')
 ! 
       V0(:) = 0.0 
       ! 
       V0(1)  = 1.0 
       V0(4)  = 1.0 
       V0(6)  = 1.0 
       V0(17) = 1.0 
       V0(55) = 1.0 
       !
       !V0(54) = 0.2*EXP(-0.5*( xGP(1)**2+xGP(2)**2)/1.0**2 ) 
       !
       r = SQRT( xGP(1)**2 + xGP(2)**2 ) 
       !
       V0(60) =     0.001*EXP(-0.5*( (xGP(1)-2.0)**2+xGP(2)**2+0*xGP(3)**2)/1.0**2 ) +     0.001*EXP(-0.5*( (xGP(1)+2.0)**2+xGP(2)**2+0*xGP(3)**2)/1.0**2 )
       IF( r>5.0 .AND. r<10. ) THEN 
           V0(61) = -0.2*xGP(2)*( 10.0 - r )/5.0  
           V0(62) = +0.2*xGP(1)*( 10.0 - r )/5.0  
       ELSEIF( r<=5.0 ) THEN
           V0(61) = -0.2*xGP(2)
           V0(62) = +0.2*xGP(1)           
       ELSE 
           V0(61:62) = 0.0  
       ENDIF 
       V0(63) = 0.0 
       V0(64) = 0.5*0.001*EXP(-0.5*( (xGP(1)-2.0)**2+xGP(2)**2+0*xGP(3)**2)/1.0**2 ) + 0.5*0.001*EXP(-0.5*( (xGP(1)+2.0)**2+xGP(2)**2+0*xGP(3)**2)/1.0**2 )
       ! 
	CASE DEFAULT
		PRINT *, 'NOT IMPLEMENTED'
		STOP
	END SELECT
	CALL PDEPrim2Cons(Q,V0)
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

   !limiter_value = 1
END SUBROUTINE PDElimitervalue

#endif
