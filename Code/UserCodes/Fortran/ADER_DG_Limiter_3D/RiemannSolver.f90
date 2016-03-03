SUBROUTINE ADERRiemannSolver(lQbndL,lFbndL,lparbndL,lQbndR,lFbndR,lparbndR,nv)
    USE typesDef
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)     :: lQbndL(nVar,nDOF(2),nDOF(3))              ! left space-time degrees of freedom 
    REAL, INTENT(IN)     :: lQbndR(nVar,nDOF(2),nDOF(3))              ! right space-time degrees of freedom 
    REAL, INTENT(IN)     :: lparbndL(nParam,nDOF(2),nDOF(3))          ! left material parameters 
    REAL, INTENT(IN)     :: lparbndR(nParam,nDOF(2),nDOF(3))          ! right material parameters 
    REAL, INTENT(IN)     :: nv(d)                                     ! normal vector 
    REAL, INTENT(INOUT)  :: lFbndL(nVar,nDOF(2),nDOF(3))              ! left flux 
    REAL, INTENT(INOUT)  :: lFbndR(nVar,nDOF(2),nDOF(3))              ! right flux 
    ! Local variables 
    INTEGER           :: i,j,k,l
    REAL              :: aux(d), QavL(nVar), QavR(nVar), smax
    REAL              :: parL(nParam), parR(nParam) 
    REAL              :: LL(nVar), LR(nVar) 
    !
    ! Compute the average states from the left and the right, which we need to compute the numerical dissipation 
    QavL = 0. 
    QavR = 0. 
    DO k = 1, nDOF(3)
      DO j = 1, nDOF(2)
            aux = (/ 1., wGPN(j), wGPN(k) /) 
            QavL = QavL + PRODUCT(aux(1:nDim))*lQbndL(:,j,k) 
            QavR = QavR + PRODUCT(aux(1:nDim))*lQbndR(:,j,k) 
        ENDDO
    ENDDO
    !
    ! If necessary, compute the left and right averages of the material parameters 
    ! 
    IF(nParam.GT.0) THEN
        parL = 0. 
        parR = 0. 
        DO k = 1, nDOF(3)
          DO j = 1, nDOF(2)
                aux = (/ 1., wGPN(j), wGPN(k) /) 
                parL = parL + PRODUCT(aux(1:nDim))*lparbndL(:,j,k) 
                parR = parR + PRODUCT(aux(1:nDim))*lparbndR(:,j,k) 
            ENDDO
        ENDDO        
    ENDIF    
    !
    ! Here, we implement a very simple Rusanov scheme with scalar dissipation (smax*Id). 
    ! We can change this into a more sophisticated Osher or HLLEM Riemann solver whenever needed! 
    !
    CALL PDEEigenvalues(LL,QavL,parL,nv) 
    CALL PDEEigenvalues(LR,QavR,parR,nv) 
    smax = MAX( MAXVAL(ABS(LL)), MAXVAL(ABS(LR)) ) 
    !
    ! We now compute the numerical flux. Note that the scheme is at the moment written in 
    ! CONSERVATION FORM => no fluctuations, but real fluxes. 
    ! Later, this will be converted into the left and right fluctuations. 
    !
    DO k = 1, nDOF(3)
      DO j = 1, nDOF(2)
            lFbndL(:,j,k) = 0.5*( lFbndR(:,j,k) + lFbndL(:,j,k) ) - 0.5*smax*( lQbndR(:,j,k) - lQbndL(:,j,k) ) 
            lFbndR(:,j,k) = lFbndL(:,j,k) 
        ENDDO
    ENDDO
    !
END SUBROUTINE ADERRiemannSolver 
    
    