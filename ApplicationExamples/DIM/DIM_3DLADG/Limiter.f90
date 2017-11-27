! DIM Limiter functions


RECURSIVE SUBROUTINE PDEAssurePositivity(status, Qmin, Qmax)
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar, nDim
	IMPLICIT NONE 
	! Argument list 
    REAL, INTENT(IN)               :: Qmin(nVar), Qmax(nVar)
    INTEGER, INTENT(OUT)           :: status

    status=+1   ! By default the solution is Physically admissible 
 
    if(Qmin(13) .lt.  -1.e-12) then
        status=-1
        return    
    end if
    
    if(Qmin(13) .lt.  1.0+1.e-12) then
        status=-1
        return    
    end if
    
    if(Qmin(13) .lt. 1.0-1.e-7 .and. Qmax(13) .gt. 1.e-7) then
        status=-1
        return
    else
        status=-1
        return        
    end if
END SUBROUTINE PDEAssurePositivity

