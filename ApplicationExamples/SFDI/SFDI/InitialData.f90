! DIM Initial Data


RECURSIVE SUBROUTINE InitParameters() 
	USE MainVariables, ONLY : nVar, nDim, EQN, ICType, MHDRotomega, EPRadius,ICdelta, ICRot,ICuL,ICA,EPCenter
	IMPLICIT NONE     
	EQN%Pi=acos(-1.0)
	select case(ICType)
		case('NBodies')
		EQN%gamma        =  1.4 
		EQN%nu           =  0.0 
		EQN%mu           =  0.0 
		EQN%cv           =  2.5 
		EQN%Pr           =  0.75
		EQN%R     = EQN%cv*(EQN%gamma-1.)
        EQN%kappa = EQN%gamma*EQN%cv*EQN%mu/EQN%Pr
		
		ICuL(1:5)        =   (/ 1.4 , 0.0 , 0.0 , 0.0 , 1.0 /)
		ICrot            =   3
		EPCenter(1:nDim) =   (/ 0.0 , 0.0 /)
		EPRadius         =   1.0   
		ICA              =   0.2   
		ICdelta          =   0.001 
		MHDRotomega      =   3.0   
	END SELECT	
	
END SUBROUTINE InitParameters

RECURSIVE SUBROUTINE InitialData(xGP, tGp, Q)
	USE, INTRINSIC :: ISO_C_BINDING
	USE MainVariables, ONLY : nVar, nDim, EQN, ICType, MHDRotomega, EPRadius,ICdelta, ICRot, ICuL,ICA,EPCenter
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: xGP(3), tGp        ! 
	REAL, INTENT(OUT)              :: Q(nVar)        ! 
	REAL :: up(nVar), r
	INTEGER :: i
	
	select case(ICType)
		case('NBodies')
			up(1:5) = ICuL(1:5)
			up(6)   = 1.0 - ICdelta
			DO i = 1, ICrot     
				r = SQRT( (xGP(1)-EPRadius*COS(i*2*EQN%Pi/ICrot))**2+(xGP(2)-EPRadius*SIN(i*2*EQN%Pi/ICRot))**2) 
				IF(r<=ICA) THEN 
					up(6) = ICdelta
				ENDIF 
			ENDDO                           
			up(7)   =  MHDRotomega/EPRadius*xGP(2) 
			up(8)   = -MHDRotomega/EPRadius*xGP(1) 
			up(9)   = 0.0 
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
   IF((observablesMin(1)<1.e-5) then
      limiter_value = 1
   end if

	
   IF((observablesMax(6)<0.999 .and. observablesMin(6)>0.001) .or. observablesMax(6)>1.001) THEN  ! was -0.001   .or. levelmin<0.001
		limiter_value = 1 
   ENDIF 
   !
	!limiter_value = 1
END SUBROUTINE PDElimitervalue

  


