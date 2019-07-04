
RECURSIVE SUBROUTINE UpdateSolutionODE(Qnew,Qold,loc_dt)
  USE MainVariables, ONLY : nVar, nDim, EQN
#ifdef ODESOLVER
    use expintegrator_type, only : ODESETTINGS, nVarODE                                     , t_ode_parameters
    use expintegrator, only : expintegrator_adaptive
    USE expintegrator_ode, only : export_expintegrator_state,extract_expintegrator_state     , convert_parameters
#endif
  USE iso_c_binding
  IMPLICIT NONE
  REAL :: Qold(nVar), Qnew(nVar),loc_dt
  INTENT(IN)  :: Qold,loc_dt
  INTENT(OUT) :: Qnew
  ! Local varialbes
#ifdef ODESOLVER
    real :: AC(nVarODE),AC0(nVarODE)
    INTEGER :: LOCiter, substep_count
    
    real ::V0(nVar),Q0(nVar)
    type(t_ode_parameters) :: ODE
#endif
	REAL :: x00(3),time
	INTEGER :: iErr
	
	x00=0.
	time=0.
	Q0=Qold
    !Q0 = subuh1(:,ii,jj,kk)    
    call extract_expintegrator_state(AC0, Q0)
    
    CALL PDECons2Prim(V0,Q0,x00,time,iErr)
    CALL convert_parameters(V0, ODE, EQN)

    CALL expintegrator_adaptive(AC0, EQN, loc_dt, ODESETTINGS,AC, substep_count, LOCiter, ODE)
	!AC=AC0
    call export_expintegrator_state(Q0, AC,(/ substep_count, LOCiter /))
    !subuh1(:,ii,jj,kk) = Q0
    Qnew=Q0	
  END SUBROUTINE UpdateSolutionODE

