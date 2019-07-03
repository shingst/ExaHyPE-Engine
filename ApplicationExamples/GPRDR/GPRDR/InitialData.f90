! DIM Initial Data


RECURSIVE SUBROUTINE InitParameters() 
	USE MainVariables, ONLY: nVar , nDim, EQN 
	IMPLICIT NONE     

END SUBROUTINE InitParameters

RECURSIVE SUBROUTINE InitialData(xGP, tGp, Q)
	USE, INTRINSIC :: ISO_C_BINDING
	USE MainVariables, ONLY : nVar, nDim, EQN, ICType
    USE SpecificVarEqn99
    USE GPRmaterials
#ifdef EQNTYPED99
    use expintegrator_ode       , only : convert_parameters
    use expintegrator_type      , only : t_ode_parameters
#endif
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: xGP(nDim), tGp        ! 
	REAL, INTENT(OUT)              :: Q(nVar)        ! 
	REAL :: up(nVar)
	REAL :: rr,xi,LamCoeff(3),LEsigma(6), LL_GPR, MM_GPR
	type(t_ode_parameters) :: ODE
	
	select case(ICType)
		case('NLOPRUPTURE')
#if defined(EQNTYPEC99) || defined(EQNTYPED99)

#ifdef EQNTYPED99
        rr=xGP(2)-2.0*xGP(1)-4000.0
        xi = 0.5+0.5*ERF(rr/100.0)
        up(22)=1-xi
        LamCoeff=getLameCoefficients(EQN%MATERIALS(1))
            up(1)     = LamCoeff(1)
            up(20)    = LamCoeff(2)
            up(19)    = LamCoeff(3)
#else
        call AssignMaterialProperties(up,ICtype2)  ! Assign the material properties as from the parfile
        rr=xGP(2)-2.0*xGP(1)-4000.0
        xi = 0.5+0.5*ERF(rr/300.0)
        call AssignMaterialPropertiesMix(up,(/'ROCK1' ,  'ROCK3'  /),(/1-xi, xi/),2,.TRUE.) 
#endif
        ! Initial velocity vector ----- !
        up(2) = 0.0                     !
        up(3) = 0.0                     ! 
        up(4) = 0.0                     !
        ! ----------------------------- !
        ! ------------ Limiter properties ----------------------
        EQN%Y0lim=1.0  ! 99% of the illness stress Y0
        ! ------------------------------------------------------
        ! Thermal properties of the material ------------------------------------------------------------
        EQN%T0    = 0.0
        EQN%ch    = 0.!0.1!1.0e-4
        EQN%tau2  = 1.e+22!1.e+5
        EQN%cv    = 1.0
        EQN%gamma = 2.0
        EQN%p0    = 0.0
        EQN%tau0  = 1.e-6
        EQN%Yeq_mode =1
        ! -----------------------------------------------------------------------------------------------
        ! Compute the initial stress --------------------------------------------------------------------------------
         LEsigma=0.
        LEsigma(2)=-120.0*1.e+6
        if(abs(xGP(1)) .le. 5000.0 .and. abs(xGP(2)) .le. 300.0/3.0 .and. abs(xGP(3)) .le. 1500.0) then 
            LEsigma(4)=70.0*1e+6
        else
            LEsigma(4)=70.0*1e+6
        end if
        !LEsigma=0.
#ifdef EQNTYPED99
        call convert_parameters(up, ODE, EQN)
        call compute_l_m_mix(LL_GPR, MM_GPR,up(19),up(20),ODE%lam1lam2ratio,ODE%mu1mu2ratio,up(21))
#else
        call compute_l_m_mix(LL_GPR, MM_GPR,up(19),up(20),up(23),up(22),up(21))
#endif
        up(9:17)= GPRsigma2ASGEOS(LEsigma,0.0,LL_GPR,MM_GPR,up(1),EQN%gamma,up(5),EQN%cv,EQN%p0,1.e-3,EQN%EOS,EQN%T0)
        !up(9:17)= (/0.999623180012770, 0.0, 0.0, -2.177843294032983D-003, 1.00149401820627,0.0,0.0,0.0,0.999627936634102/)
        ! -----------------------------------------------------------------------------------------------------------
        ! Rupture properties ----------------------------------------------------------------------------------------
        up(21)=0.0

        !if(nDim .eq. 3) then
        !    if(abs(xGP(1)) .le. 100.0/3.0 .and. abs(xGP(2)) .le. 100.0/3.0 .and. abs(xGP(3)) .le. 100.0/3.0) then 
        !        up(21)=0.5    
        !    end if            
        !else
        !    if(abs(xGP(1)) .le. 1000.0 .and. abs(xGP(2)) .le. 100.0/3.0 .and. abs(xGP(3)) .le. 100.0/3.0) then 
        !        up(21)=1.0    
        !    end if            
        !end if
        !if(abs(xGP(1)) .le. 100.0/3.0 .and. abs(xGP(2)) .le. 100.0/3.0 .and. abs(xGP(3)) .le. 100.0/3.0) then 
        !    up(21)=0.1    
        !end if
        ! -----------------------------------------------------------------------------------------------------------
        ! Complex geometry ------------------------------------------------------------------------------------------
        up(18)=1     ! alpha
        up(5)=0     ! alpha
        ! -----------------------------------------------------------------------------------------------------------
        

#endif		

	END SELECT
	CALL PDEPrim2Cons(Q,up)
END SUBROUTINE InitialData

RECURSIVE SUBROUTINE PDElimitervalue(limiter_value,xx,numberOfObservables, observablesMin, observablesMax)
	USE, INTRINSIC :: ISO_C_BINDING
	USE MainVariables, ONLY : nVar, nDim
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: xx(nDim)        ! 
	INTEGER, INTENT(IN)					:: numberOfObservables
	INTEGER, INTENT(OUT)              :: limiter_value        !
	REAL, INTENT(IN)					:: observablesMin(numberOfObservables), observablesMax(numberOfObservables)
	real	:: rr	

	limiter_value=0
	!if(observablesMin(1).ge.1.1 .and. observablesMax(1).le. 2.6) then
	!	limiter_value=1
	!end if
END SUBROUTINE PDElimitervalue

RECURSIVE SUBROUTINE ShuVortex2D(x, t, Q)
    USE, INTRINSIC :: ISO_C_BINDING
    USE MainVariables, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 

	INTEGER	:: iErr
    REAL    :: up(nVar), Pi = ACOS(-1.0), epsilon, r, du, dv, dT, drho, dp

       CALL PDEPrim2Cons(Q,up)
    !Q=up
END SUBROUTINE ShuVortex2D
  


