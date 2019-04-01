
!    
MODULE GREOS
    IMPLICIT NONE
    !
    REAL, PARAMETER :: gamma = 2.0
    private
    ! Functions needed to generate the EOS table  
    !
    INTERFACE InitGREOS
        MODULE PROCEDURE InitGREOS
    END INTERFACE
    !
    INTERFACE ReadEOSCoefficients
        MODULE PROCEDURE ReadEOSCoefficients
    END INTERFACE
    !
    INTERFACE EvaluateEOS
        MODULE PROCEDURE EvaluateEOS
    END INTERFACE
    !
    INTERFACE DP2TE
        MODULE PROCEDURE DP2TE
    END INTERFACE
    !
    INTERFACE DP2TE_scalar
        MODULE PROCEDURE DP2TE_scalar
    END INTERFACE
    !
    INTERFACE EvaluateEOS_scalar
        MODULE PROCEDURE EvaluateEOS_scalar
    END INTERFACE
    !
    PUBLIC ::  InitGREOS,ReadEOSCoefficients, EvaluateEOS, DP2TE,DP2TE_scalar, EvaluateEOS_scalar
    !
    CONTAINS


RECURSIVE SUBROUTINE InitGREOS
    USE EOSTypes
    IMPLICIT NONE
    !
    InputFile  = 'eos.table'                ! input file 
    EOSFile    = 'eos.dg'                   ! output file in DG format  
    !
END SUBROUTINE

! Pressure on temperature and density (POTD) 

RECURSIVE SUBROUTINE POTD(lp,lT,lrho) 
    USE EOSTypes
    IMPLICIT NONE 
    ! Argument list 
    REAL(8) :: lp(VECTORLENGTH), lT(VECTORLENGTH), lrho(VECTORLENGTH) 
    ! Local variables 
    INTEGER(8) :: i(VECTORLENGTH),   j(VECTORLENGTH), iv, ii, jj, kk, ll   
    REAL(8)    :: Lagrange(VECTORLENGTH) 
    !
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64:: i,j,Lagrange
    !DIR$ ASSUME_ALIGNED lp         : 64
    !DIR$ ASSUME_ALIGNED lt         : 64
    !DIR$ ASSUME_ALIGNED lrho       : 64
#else
    !DIR$ ATTRIBUTES ALIGN: 32:: i,j,Lagrange
    !DIR$ ASSUME_ALIGNED lp         : 32
    !DIR$ ASSUME_ALIGNED lt         : 32
    !DIR$ ASSUME_ALIGNED lrho       : 32
#endif 
    ! debug 
    !lp(:) = (4./3.-1.0)*lrho(:)*0.001*lT(:)  
    !RETURN 
    !
    ! Here we implement a simple centered sliding Lagrange interpolation of piecewise polynomials of degree 3 
    i = MIN( MAX(2, CEILING( (lT(:)  -EOS%Tmin  )/EOS%dTin   ) ), EOS%NT-2   ) 
    j = MIN( MAX(2, CEILING( (lrho(:)-EOS%rhomin)/EOS%drhoin ) ), EOS%NRHO-2 ) 
    !
    lp(:) = 0.0
    DO jj = -1, 2
     DO ii = -1, 2 
         Lagrange = 1.0 
         ! Lagrange polynomial in the first direction (T) 
         DO kk = -1, 2
           IF( kk==ii) CYCLE
           Lagrange(:) = Lagrange(:) * ( lT(:) - EOS%T(i(:)+kk) )/( EOS%T(i(:)+ii) - EOS%T(i(:)+kk) ) 
         ENDDO
         ! Lagrange polynomial in the second direction (rho) 
         DO ll = -1, 2
           IF( ll==jj) CYCLE 
           Lagrange(:) = Lagrange(:) * ( lrho(:) - EOS%rho(j(:)+ll) )/( EOS%rho(j(:)+jj) - EOS%rho(j(:)+ll) )  
         ENDDO
         ! 
         DO iv = 1, VECTORLENGTH 
            lp(iv) = lp(iv) + EOS%p(i(iv)+ii,j(iv)+jj)*Lagrange(iv) 
        ENDDO 
     ENDDO
    ENDDO    
    !
END SUBROUTINE POTD 
    
! Pressure on temperature and density (POTD) 

RECURSIVE SUBROUTINE POTD_scalar(lp,lT,lrho) 
    USE EOSTypes
    IMPLICIT NONE 
    ! Argument list 
    REAL :: lp, lT, lrho 
    ! Local variables 
    INTEGER :: i,   j, iv, ii, jj, kk, ll   
    REAL    :: Lagrange 
    !
    ! debug 
    !lp(:) = (4./3.-1.0)*lrho(:)*0.001*lT(:)  
    !RETURN 
    !
    ! Here we implement a simple centered sliding Lagrange interpolation of piecewise polynomials of degree 3 
    i = MIN( MAX(2, CEILING( (lT  -EOS%Tmin  )/EOS%dTin   ) ), EOS%NT-2   ) 
    j = MIN( MAX(2, CEILING( (lrho-EOS%rhomin)/EOS%drhoin ) ), EOS%NRHO-2 ) 
    !
    lp = 0.0
    DO jj = -1, 2
     DO ii = -1, 2 
         Lagrange = 1.0 
         ! Lagrange polynomial in the first direction (T) 
         DO kk = -1, 2
           IF( kk==ii) CYCLE
           Lagrange = Lagrange * ( lT - EOS%T(i+kk) )/( EOS%T(i+ii) - EOS%T(i+kk) ) 
         ENDDO
         ! Lagrange polynomial in the second direction (rho) 
         DO ll = -1, 2
           IF( ll==jj) CYCLE 
           Lagrange = Lagrange * ( lrho - EOS%rho(j+ll) )/( EOS%rho(j+jj) - EOS%rho(j+ll) )  
         ENDDO
         ! 
         lp = lp + EOS%p(i+ii,j+jj)*Lagrange 
     ENDDO
    ENDDO    
    !
END SUBROUTINE POTD_scalar
    
! Internal energy on temperature and density (UOTD)      

RECURSIVE SUBROUTINE UOTD(le,lT,lrho) 
    USE EOSTypes
    IMPLICIT NONE 
    ! Argument list 
    REAL(8) :: le(VECTORLENGTH), lT(VECTORLENGTH), lrho(VECTORLENGTH) 
    ! Local variables 
    INTEGER(8) :: i(VECTORLENGTH),   j(VECTORLENGTH), iv, ii, jj, kk, ll   
    REAL(8)    :: Lagrange(VECTORLENGTH) 
    !
#ifdef AVX512
    !DIR$ attributes align: 64:: i,j,Lagrange
    !DIR$ ASSUME_ALIGNED le         : 64
    !DIR$ ASSUME_ALIGNED lt         : 64
    !DIR$ ASSUME_ALIGNED lrho       : 64
#else
    !DIR$ ATTRIBUTES ALIGN: 32:: i,j,Lagrange
    !DIR$ ASSUME_ALIGNED le         : 32
    !DIR$ ASSUME_ALIGNED lt         : 32
    !DIR$ ASSUME_ALIGNED lrho       : 32
#endif 
    ! debug 
    !le(:) = 0.001*lT(:)  
    !RETURN    
    !
    ! Here we implement a simple centered sliding Lagrange interpolation of piecewise polynomials of degree 3 
    i = MIN( MAX(2, CEILING( (lT(:)  -EOS%Tmin  )/EOS%dTin   ) ), EOS%NT-2   ) 
    j = MIN( MAX(2, CEILING( (lrho(:)-EOS%rhomin)/EOS%drhoin ) ), EOS%NRHO-2 ) 
    !
    le(:) = 0.0
    DO jj = -1, 2
     DO ii = -1, 2 
         Lagrange = 1.0 
         ! Lagrange polynomial in the first direction (T) 
         DO kk = -1, 2
           IF( kk==ii) CYCLE
           Lagrange(:) = Lagrange(:) * ( lT(:) - EOS%T(i(:)+kk) )/( EOS%T(i(:)+ii) - EOS%T(i(:)+kk) ) 
         ENDDO
         ! Lagrange polynomial in the second direction (rho) 
         DO ll = -1, 2
           IF( ll==jj) CYCLE 
           Lagrange(:) = Lagrange(:) * ( lrho(:) - EOS%rho(j(:)+ll) )/( EOS%rho(j(:)+jj) - EOS%rho(j(:)+ll) )  
         ENDDO
         ! 
         DO iv = 1, VECTORLENGTH 
            le(iv) = le(iv) + EOS%e(i(iv)+ii,j(iv)+jj)*Lagrange(iv) 
        ENDDO 
     ENDDO
    ENDDO    
    !
END SUBROUTINE UOTD 

    
! Internal energy on temperature and density (UOTD)      

RECURSIVE SUBROUTINE UOTD_scalar(le,lT,lrho) 
    USE EOSTypes
    IMPLICIT NONE 
    ! Argument list 
    REAL :: le, lT, lrho 
    ! Local variables 
    INTEGER :: i,   j, iv, ii, jj, kk, ll   
    REAL    :: Lagrange 
    !
    ! debug 
    !le = 0.001*lT  
    !RETURN    
    !
    ! Here we implement a simple centered sliding Lagrange interpolation of piecewise polynomials of degree 3 
    i = MIN( MAX(2, CEILING( (lT  -EOS%Tmin  )/EOS%dTin   ) ), EOS%NT-2   ) 
    j = MIN( MAX(2, CEILING( (lrho-EOS%rhomin)/EOS%drhoin ) ), EOS%NRHO-2 ) 
    !
    le = 0.0
    DO jj = -1, 2
     DO ii = -1, 2 
         Lagrange = 1.0 
         ! Lagrange polynomial in the first direction (T) 
         DO kk = -1, 2
           IF( kk==ii) CYCLE
           Lagrange = Lagrange * ( lT - EOS%T(i+kk) )/( EOS%T(i+ii) - EOS%T(i+kk) ) 
         ENDDO
         ! Lagrange polynomial in the second direction (rho) 
         DO ll = -1, 2
           IF( ll==jj) CYCLE 
           Lagrange = Lagrange * ( lrho - EOS%rho(j+ll) )/( EOS%rho(j+jj) - EOS%rho(j+ll) )  
         ENDDO
         ! 
         le = le + EOS%e(i+ii,j+jj)*Lagrange
     ENDDO
    ENDDO    
    !
END SUBROUTINE UOTD_scalar 

    
    
!
! This subroutine converts the conservative quantitites density (rho) and internal energy (e) to the primitive variables
! temperature (T) and pressure (p). This is the EOS that we typically need in a Godunov-type scheme. 
!

RECURSIVE SUBROUTINE DE2TP(T,p,rho,e) 
    USE EOSTypes 
    !-------------------------------------       
    IMPLICIT NONE 
    !-------------------------------------       
    ! Argument list declaration
    REAL(8) :: T(VECTORLENGTH),p(VECTORLENGTH),rho(VECTORLENGTH),e(VECTORLENGTH) 
    ! Local variable declaration 
    INTEGER :: iNewton, inner 
    REAL(8)    :: u(VECTORLENGTH), u1(VECTORLENGTH), u2(VECTORLENGTH), f(VECTORLENGTH), df(VECTORLENGTH), dT(VECTORLENGTH), res(VECTORLENGTH), res2(VECTORLENGTH), dudT(VECTORLENGTH), delta(VECTORLENGTH) 
    REAL(8)    :: Ta(VECTORLENGTH),Tb(VECTORLENGTH),Tc(VECTORLENGTH),fa(VECTORLENGTH),fb(VECTORLENGTH),fc(VECTORLENGTH) 
    ! Temperature range of the EOS 
    REAL                :: TMIN, TMAX 
    ! Numerical parameters 
    REAL, PARAMETER     :: epsilon   = 1e-7 
    REAL, PARAMETER     :: tol       = 1e-10
    INTEGER, PARAMETER  :: MAXNEWTON = 50 
    !-------------------------------------       
    ! Input and output variables 
    INTENT(IN)  :: rho,e 
    INTENT(OUT) :: T,p         
    !-------------------------------------       
    !
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64:: u, u1, u2, f, df, dT, res, res2, dudT, delta,Ta,Tb,Tc,fa,fb,fc 
    !DIR$ ASSUME_ALIGNED t         : 64
    !DIR$ ASSUME_ALIGNED p         : 64
    !DIR$ ASSUME_ALIGNED rho       : 64
    !DIR$ ASSUME_ALIGNED e         : 64
#else
    !DIR$ ATTRIBUTES ALIGN: 32:: u, u1, u2, f, df, dT, res, res2, dudT, delta,Ta,Tb,Tc,fa,fb,fc 
    !DIR$ ASSUME_ALIGNED t         : 32
    !DIR$ ASSUME_ALIGNED p         : 32
    !DIR$ ASSUME_ALIGNED rho       : 32
    !DIR$ ASSUME_ALIGNED e         : 32
#endif 
!#ifdef DEBUG    
!    ! debug 
!    p(:) = (gamma-1.0)*rho(:)*e(:) 
!    T = 0.001*e(:) 
!    return 
!#endif    
    !
    TMIN = EOS%Tmin
    TMAX = EOS%Tmax 
    !
    ! Part 1. Compute the temperature using Newton's method, 
    ! starting from an arbitrary initial temperature 
    !
    T(:) = 0.5*( TMIN + TMAX ) 
    dT = T 
    DO iNewton = 1, MAXNEWTON
       CALL UOTD(u,T,rho)
       f = u-e
       res  =  ABS(f)/(ABS(e)+1.)  
       res2 =  ABS(dT)  
       IF( MAXVAL(res).LT.tol ) THEN
            EXIT
       ENDIF  
       CALL UOTD(u2,T*(1+epsilon),rho)
       CALL UOTD(u1,T*(1-epsilon),rho)
       dudT = (u2-u1)/(2.0*T*epsilon)
       df   = dudT 
       dT   = -f/dudT
       delta = 1.
       ! Globalization strategy: linesearch   
       DO inner = 1, 10
           CALL UOTD(u,T+delta*dT,rho) 
           WHERE( ABS(u-e)/(ABS(e)+1.) .LT. res ) 
                delta = delta 
           ELSEWHERE 
                delta = 0.5*delta
           ENDWHERE 
       ENDDO
       T    = T + delta*dT            
    ENDDO
    !
    CONTINUE
    !
    ! If Newton does not work, we may want to use a more robust bisection algorithm... 
    ! 
    IF( ANY(ISNAN(T)).OR.(iNewton.GT.MAXNEWTON-1)) THEN
        ! Try the more robust Bisection instead.... 
        !    PRINT *, ' newton warning ', iNewton, res, T
        !T = -1.
        !p = -1.
        Ta = EOS%Tmin 
        Tb = EOS%Tmax 
        CALL UOTD(u,Ta,rho)
        fa = u-e 
        CALL UOTD(u,Tb,rho)
        fb = u-e
        IF( ANY(fa*fb.GT.0.) ) THEN
            ! PRINT *, ' Bisection error. Signs are equal. ', fa, fb   
            T = -1. 
            p = -1. 
            RETURN            
        ENDIF
        DO iNewton = 1, MAXNEWTON
            Tc = 0.5*(Ta+Tb)          ! Bisection 
            !Tc = (fa*Tb-Ta*fb)/(fa-fb)  ! Regula falsi 
            CALL UOTD(u,Tc,rho)
            fc = u-e 
            res  =  ABS(fc)/(ABS(e)+1.)  
            IF(ALL(res.LT.tol)) THEN  
                T = Tc 
                EXIT
            ENDIF 
            CALL UOTD(u,Tc,rho) 
            WHERE( fa*fc.GT.0. ) 
                Ta = Tc 
                fa = u-e
            ELSEWHERE(fb*fc.GT.0.)  
                Tb = Tc 
                fb = u-e
            ENDWHERE 
        ENDDO 
    ENDIF    
    !
    WHERE( ISNAN(T) ) 
        T = -1.
        p = -1. 
    ENDWHERE 
    !
    ! Part 2. Compute the pressure from the (now) known temperature and density 
    !
    CALL POTD(p,T,rho)   
    ! 
END SUBROUTINE DE2TP 
    

!
! This subroutine converts the quantitites density (rho) and pressure (p) to the variables 
! temperature (T) and internal energy (e). Needed for initial condition setup.  
!

RECURSIVE SUBROUTINE DP2TE(T,e,rho,p) 
    USE EOSTypes 
    !-------------------------------------       
    IMPLICIT NONE 
    !-------------------------------------       
    ! Argument list declaration
    REAL(8) :: T(VECTORLENGTH),p(VECTORLENGTH),rho(VECTORLENGTH),e(VECTORLENGTH) 
    ! Local variable declaration 
    INTEGER :: iNewton, inner 
    REAL(8)    :: u(VECTORLENGTH), u1(VECTORLENGTH), u2(VECTORLENGTH), f(VECTORLENGTH), df(VECTORLENGTH), dT(VECTORLENGTH), res(VECTORLENGTH), res2(VECTORLENGTH), dudT(VECTORLENGTH), delta(VECTORLENGTH) 
    REAL(8)    :: Ta(VECTORLENGTH),Tb(VECTORLENGTH),Tc(VECTORLENGTH),fa(VECTORLENGTH),fb(VECTORLENGTH),fc(VECTORLENGTH) 
    ! Temperature range of the EOS 
    REAL                :: TMIN, TMAX 
    ! Numerical parameters 
    REAL, PARAMETER     :: epsilon   = 1e-7 
    REAL, PARAMETER     :: tol       = 1e-12   
    INTEGER, PARAMETER  :: MAXNEWTON = 50 
    !-------------------------------------       
    ! Input and output variables 
    INTENT(IN)  :: rho,p 
    INTENT(OUT) :: T,e          
    !-------------------------------------       
    ! 
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64:: u, u1, u2, f, df, dT, res, res2, dudT, delta,Ta,Tb,Tc,fa,fb,fc
    !DIR$ ASSUME_ALIGNED t         : 64
    !DIR$ ASSUME_ALIGNED p         : 64
    !DIR$ ASSUME_ALIGNED rho       : 64
    !DIR$ ASSUME_ALIGNED e         : 64
#else
    !DIR$ ATTRIBUTES ALIGN: 32:: u, u1, u2, f, df, dT, res, res2, dudT, delta,Ta,Tb,Tc,fa,fb,fc
    !DIR$ ASSUME_ALIGNED t         : 32
    !DIR$ ASSUME_ALIGNED p         : 32
    !DIR$ ASSUME_ALIGNED rho       : 32
    !DIR$ ASSUME_ALIGNED e         : 32
#endif 
!
#ifdef GEOSDEBUG    
    e = (gamma-1.0)*p/rho
    T = 0.001*e
    RETURN
#endif    
    !
    TMIN = EOS%Tmin
    TMAX = EOS%Tmax 
    !
    ! Part 1. Compute the temperature using Newton's method, 
    ! starting from an arbitrary initial temperature 
    !
    T(:) = 0.5*( TMIN + TMAX ) 
    dT = T 
    DO iNewton = 1, MAXNEWTON
       CALL POTD(u,T,rho) 
       f = u-p 
       res  =  ABS(f)/(ABS(p)+1.)  
       res2 =  ABS(dT)  
       IF( MAXVAL(res).LT.tol ) THEN
            EXIT
       ENDIF  
       CALL POTD(u2,T*(1+epsilon),rho)
       CALL POTD(u1,T*(1-epsilon),rho)
       dudT = (u2-u1)/(2.0*T*epsilon)
       df   = dudT 
       dT   = -f/dudT
       delta = 1.
       ! Globalization strategy: linesearch   
       !DO inner = 1, 10
       !    CALL POTD(u,T+delta*dT,rho) 
       !    WHERE( ABS(u-p)/(ABS(p)+1.) .LT. res ) 
       !         delta = delta 
       !    ELSEWHERE 
       !         delta = 0.5*delta
       !    ENDWHERE 
       !ENDDO
       T    = T + delta*dT            
    ENDDO
    !
    CONTINUE
    !
    ! Part 2. Compute the energy from the (now) known temperature and density 
    !
    CALL UOTD(e,T,rho)   
    ! 
END SUBROUTINE DP2TE      
    


RECURSIVE SUBROUTINE DP2TE_scalar(T,e,rho,p) 
    USE EOSTypes 
    !-------------------------------------       
    IMPLICIT NONE 
    !-------------------------------------       
    ! Argument list declaration
    REAL :: T,p,rho,e 
    ! Local variable declaration 
    INTEGER :: iNewton, inner 
    REAL    :: u, u1, u2, f, df, dT, res, res2,res3, dudT, delta 
    REAL    :: Ta,Tb,Tc,fa,fb,fc 
    ! Temperature range of the EOS 
    REAL                :: TMIN, TMAX 
    ! Numerical parameters 
    REAL, PARAMETER     :: epsilon   = 1e-7 
    REAL, PARAMETER     :: tol       = 1e-12   
    INTEGER, PARAMETER  :: MAXNEWTON = 50 
    !-------------------------------------       
    ! Input and output variables 
    INTENT(IN)  :: rho,p 
    INTENT(OUT) :: T,e          
    !-------------------------------------       
    !
#ifdef GEOSDEBUG    
    e = (gamma-1.0)*p/rho
    T = 0.001*e
    RETURN
#endif    
    !
    TMIN = EOS%Tmin
    TMAX = EOS%Tmax 
    !
    ! Part 1. Compute the temperature using Newton's method, 
    ! starting from an arbitrary initial temperature 
    !
    T = 0.5*( TMIN + TMAX ) 
    dT = T 
    DO iNewton = 1, MAXNEWTON
       CALL POTD_scalar(u,T,rho) 
       f = u-p 
       res  =  ABS(f)/(ABS(p)+1.)  
       res2 =  ABS(dT)  
       IF(res.LT.tol) THEN
           continue
            EXIT
       ENDIF  
       CALL POTD_scalar(u2,T*(1+epsilon),rho)
       CALL POTD_scalar(u1,T*(1-epsilon),rho)
       dudT = (u2-u1)/(2.0*T*epsilon)
       df   = dudT 
       dT   = -f/dudT
       delta = 1.
       ! Globalization strategy: linesearch   
       DO inner = 1, 50
           CALL POTD_scalar(u,T+delta*dT,rho) 
           res3 = ABS(u-p)/(ABS(p)+1.)
           IF( res3  .LT. res )  THEN
                delta = delta 
                exit
           ELSE 
                delta = 0.5*delta
                continue
           ENDIF 
       ENDDO
       T    = T + delta*dT            
    ENDDO
    !
    CONTINUE
    !
    ! Part 2. Compute the energy from the (now) known temperature and density 
    !
    CALL UOTD_scalar(e,T,rho)   
    ! 
END SUBROUTINE DP2TE_scalar    


RECURSIVE SUBROUTINE GenerateTable
    USE EOSTypes
    IMPLICIT NONE 
    INTEGER :: i, j, k, iGP, jGP, ii, jj, count
    INTEGER(8) :: ix(VECTORLENGTH), jx(VECTORLENGTH)   
    REAL(8)    :: xiGP(VECTORLENGTH), etaGP(VECTORLENGTH)
    REAL(8)    :: rhoGP(VECTORLENGTH), eGP(VECTORLENGTH) 
    REAL(8)    :: TGP(VECTORLENGTH), pGP(VECTORLENGTH) 
    !-----------------------------------       
    ! 
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64:: xiGP, etaGP,rhoGP, eGP,TGP, pGP ,ix,jx
#else
    !DIR$ ATTRIBUTES ALIGN: 32:: xiGP, etaGP,rhoGP, eGP,TGP, pGP ,ix,jx
#endif 
!

    PRINT *, ' Generating DG table... ' 
    !
    EOS%drho = (EOS%rhomax-EOS%rhomin) / REAL(EOS%IMAX) 
    EOS%de   = (EOS%emax  -EOS%emin  ) / REAL(EOS%JMAX) 
    !
    ALLOCATE( EOS%DGp(EOS%N+1,EOS%N+1,EOS%IMAX,EOS%JMAX) ) 
    ALLOCATE( EOS%DGT(EOS%N+1,EOS%N+1,EOS%IMAX,EOS%JMAX) ) 
    !
    DO j = 1, EOS%JMAX
    if(j==20) then
        continue
    endif
    
     DO i = 1, EOS%IMAX
         count = 0 
         DO jGP = 1, EOS%N+1
          DO iGP = 1, EOS%N+1
              count = count + 1 
              xiGP(count)  = EOS%xiGP(iGP)
              etaGP(count) = EOS%xiGP(jGP) 
              ix(count)    = iGP 
              jx(count)    = jGP 
              rhoGP(count) = EOS%rhomin + (REAL(i-1)+xiGP(count) )*EOS%drho 
              eGP(count)   = EOS%emin   + (REAL(j-1)+etaGP(count))*EOS%de 
              IF( MOD(count,VECTORLENGTH)==0 ) THEN
                  CALL DE2TP(TGP,pGP,rhoGP,eGP)
                  count = 0 
                  DO k = 1, VECTORLENGTH
                      ii = ix(k)
                      jj = jx(k) 
                      EOS%DGp(ii,jj,i,j) = pGP(k) 
                      EOS%DGT(ii,jj,i,j) = TGP(k) 
                  ENDDO                  
                  CONTINUE
              ENDIF                            
          ENDDO
         ENDDO         
     ENDDO
    ENDDO    
    !
    PRINT *, ' Done. '    
    !
    CONTINUE 
    !
END SUBROUTINE GenerateTable 


RECURSIVE SUBROUTINE EvaluateEOS(p,dpdrho,dpde,T,rho,e) 
    USE EOSTypes
    IMPLICIT NONE 
    REAL(8)    :: p(VECTORLENGTH), dpdrho(VECTORLENGTH), dpde(VECTORLENGTH), T(VECTORLENGTH), rho(VECTORLENGTH), e(VECTORLENGTH) 
    INTEGER(8) :: i(VECTORLENGTH), j(VECTORLENGTH)
    REAL(8)    :: xiGP(VECTORLENGTH), etaGP(VECTORLENGTH)
    REAL(8)    :: phi1(VECTORLENGTH,EOS%N+1), phi1_xi(VECTORLENGTH,EOS%N+1) 
    REAL(8)    :: phi2(VECTORLENGTH,EOS%N+1), phi2_xi(VECTORLENGTH,EOS%N+1) 
    INTEGER    :: ii,jj,iv
    !
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64 :: i,j,xiGP,etaGP,phi1,phi1_xi,phi2,phi2_xi
    !DIR$ ASSUME_ALIGNED p      : 64
    !DIR$ ASSUME_ALIGNED dpdrho : 64
    !DIR$ ASSUME_ALIGNED dpde   : 64
    !DIR$ ASSUME_ALIGNED T      : 64
    !DIR$ ASSUME_ALIGNED rho    : 64
    !DIR$ ASSUME_ALIGNED e      : 64
#else
    !DIR$ ATTRIBUTES ALIGN: 32 :: i,j,xiGP,etaGP,phi1,phi1_xi,phi2,phi2_xi
    !DIR$ ASSUME_ALIGNED p      : 32
    !DIR$ ASSUME_ALIGNED dpdrho : 32
    !DIR$ ASSUME_ALIGNED dpde   : 32
    !DIR$ ASSUME_ALIGNED T      : 32
    !DIR$ ASSUME_ALIGNED rho    : 32
    !DIR$ ASSUME_ALIGNED e      : 32
#endif 
    !
    i = CEILING( (rho(:)-EOS%rhomin)/EOS%drho   ) 
    j = CEILING( (e(:)  -EOS%emin  )/EOS%de     ) 
    xiGP(:)  = ( rho(:) - (EOS%rhomin+(i-1)*EOS%drho) )/EOS%drho 
    etaGP(:) = ( e(:)   - (EOS%emin  +(j-1)*EOS%de)   )/EOS%de  
    !CALL EOSBaseFunc1D(phi1,phi1_xi,xiGP, EOS%xiGP,EOS%N) 
    !CALL EOSBaseFunc1D(phi2,phi2_xi,etaGP,EOS%xiGP,EOS%N) 
    !CALL EOSBaseFunc1D_bis(phi1,phi1_xi,xiGP,phi2,phi2_xi,etaGP,EOS%xiGP,EOS%N) 
    CALL EOSBaseFunc1D_bis7(phi1,phi1_xi,xiGP,phi2,phi2_xi,etaGP,EOSxiGP) 
    !
    p = 0.0
    dpdrho = 0.0
    dpde = 0.0
    T = 0.0 
    !
    DO jj = 1, EOS%N+1
     DO ii = 1, EOS%N+1
         DO iv = 1, VECTORLENGTH 
             p(iv) = p(iv) + phi1(iv,ii)*phi2(iv,jj)*EOSDGp(ii,jj,i(iv),j(iv)) 
             dpdrho(iv) = dpdrho(iv) + phi1_xi(iv,ii)*phi2(iv,jj)*EOSDGp(ii,jj,i(iv),j(iv))/EOSdrho
             dpde(iv) = dpde(iv) + phi1(iv,ii)*phi2_xi(iv,jj)*EOSDGp(ii,jj,i(iv),j(iv))/EOSde 
             T(iv) = T(iv) + phi1(iv,ii)*phi2(iv,jj)*EOSDGT(ii,jj,i(iv),j(iv))  
         ENDDO         
     ENDDO
    ENDDO 
    !
    CONTINUE 
    !
END SUBROUTINE EvaluateEOS 



RECURSIVE SUBROUTINE EvaluateEOS_scalar(p,dpdrho,dpde,T,rho,e) 
    USE EOSTypes
    IMPLICIT NONE 
    REAL    :: p, dpdrho, dpde, T, rho, e 
    INTEGER :: i, j, ii, jj, iv  
    REAL    :: xiGP, etaGP
    REAL    :: phi1(EOS%N+1), phi1_xi(EOS%N+1) 
    REAL    :: phi2(EOS%N+1), phi2_xi(EOS%N+1) 
    !
    i = CEILING( (rho-EOS%rhomin)/EOS%drho   ) 
    j = CEILING( (e  -EOS%emin  )/EOS%de     ) 
    xiGP  = ( rho - (EOS%rhomin+(i-1)*EOS%drho) )/EOS%drho 
    etaGP = ( e   - (EOS%emin  +(j-1)*EOS%de)   )/EOS%de  
    CALL EOSBaseFunc1D_scalar(phi1,phi1_xi,xiGP, EOS%xiGP,EOS%N) 
    CALL EOSBaseFunc1D_scalar(phi2,phi2_xi,etaGP,EOS%xiGP,EOS%N) 
    !
    p = 0.0 
    dpdrho = 0.
    dpde = 0.
    T = 0.0 
    DO jj = 1, EOS%N+1
     DO ii = 1, EOS%N+1
        p = p + phi1(ii)*phi2(jj)*EOS%DGp(ii,jj,i,j) 
        dpdrho = dpdrho + phi1_xi(ii)*phi2(jj)*EOS%DGp(ii,jj,i,j)/EOS%drho
        dpde = dpde + phi1(ii)*phi2_xi(jj)*EOS%DGp(ii,jj,i,j)/EOS%de 
        T = T + phi1(ii)*phi2(jj)*EOS%DGT(ii,jj,i,j) 
     ENDDO
    ENDDO 
    !
    CONTINUE 
    !
    !p = 0.0 
    !dpdrho = 0.
    !dpde = 0.
    !T = 0.0
    !DO jj = 1, nP
    !    phiarray_p(:,jj)    =    phi1(:)   *phi2(jj)  !*EOSDGp(:,jj,i,j)
    !    phiarray_prho(:,jj) = phi1_xi(:)   *phi2(jj)  !*EOSDGp(:,jj,i,j)
    !    phiarray_pe(:,jj)   =    phi1(:)*phi2_xi(jj)  !*EOSDGp(:,jj,i,j)
    !    phiarray_T(:,jj)    =    phi1(:)   *phi2(jj)  !*EOSDGT(:,jj,i,j)
    !ENDDO  
    !!
    !DO jj = 1, nP  
    !    p       = p         + DOT_PRODUCT(phiarray_p(:,jj),EOSDGp(:,jj,i,j))
    !    dpdrho  = dpdrho    + DOT_PRODUCT(phiarray_prho(:,jj),EOSDGp(:,jj,i,j))
    !    dpde    = dpde      + DOT_PRODUCT(phiarray_pe(:,jj),EOSDGp(:,jj,i,j))
    !    T       = T         + DOT_PRODUCT(phiarray_T(:,jj),EOSDGT(:,jj,i,j))
    !ENDDO 
    

END SUBROUTINE EvaluateEOS_scalar


RECURSIVE SUBROUTINE WriteEOSCoefficients
    USE EOSTypes
    IMPLICIT NONE
    INTEGER :: i, j, ii, jj
    !
    PRINT *, ' Writing EOS table to file... ', TRIM(EOSFile) 
    OPEN(UNIT=333, FILE=TRIM(EOSFile), ACCESS='STREAM' )
    WRITE(333) EOS%N, EOS%IMAX, EOS%JMAX 
    WRITE(333) EOS%rhomin, EOS%rhomax, EOS%emin, EOS%emax, EOS%Tmin, EOS%Tmax, EOS%pmin, EOS%pmax 
    DO j = 1, EOS%JMAX
     DO i = 1, EOS%IMAX
         WRITE(333) EOS%DGp(:,:,i,j) 
         WRITE(333) EOS%DGT(:,:,i,j) 
     ENDDO
    ENDDO 
    ! If original table is also requested, then write it 
    IF(OriginalTable==1) THEN
        WRITE(333) EOS%NRHO, EOS%NT 
        WRITE(333) EOS%T
        WRITE(333) EOS%rho 
        WRITE(333) EOS%p
        WRITE(333) EOS%e         
    ENDIF         
    CLOSE(333) 
    PRINT *, ' Done. ' 
    !
END SUBROUTINE WriteEOSCoefficients 



RECURSIVE SUBROUTINE ReadEOSCoefficients
    USE EOSTypes
    IMPLICIT NONE
    INTEGER :: i, j, ii, jj
    !
    PRINT *, ' Reading EOS table to file... ', TRIM(EOSFile) 
    OPEN(UNIT=333, FILE=TRIM(EOSFile), ACCESS='STREAM' )
    READ(333) EOS%N, EOS%IMAX, EOS%JMAX 
    ALLOCATE( EOS%xiGP(EOS%N+1), EOS%wGP(EOS%N+1) ) 
    CALL gauleg(0.0,1.0,EOS%xiGP,EOS%wGP,EOS%N+1) 
    READ(333) EOS%rhomin, EOS%rhomax, EOS%emin, EOS%emax, EOS%Tmin, EOS%Tmax, EOS%pmin, EOS%pmax 
    ALLOCATE( EOS%DGp(EOS%N+1,EOS%N+1,EOS%IMAX,EOS%JMAX) ) 
    ALLOCATE( EOS%DGT(EOS%N+1,EOS%N+1,EOS%IMAX,EOS%JMAX) ) 
    EOS%drho = (EOS%rhomax-EOS%rhomin) / REAL(EOS%IMAX) 
    EOS%de   = (EOS%emax  -EOS%emin  ) / REAL(EOS%JMAX) 
    !
    EOSN=EOS%N
    EOSIMAX=EOS%IMAX
    EOSJMAX=EOS%JMAX
    ALLOCATE( EOSxiGP(EOS%N+1), EOSwGP(EOS%N+1) ) 
    CALL gauleg(0.0,1.0,EOSxiGP,EOSwGP,EOS%N+1) 
    ALLOCATE( EOSDGp(EOSN+1,EOSN+1,EOSIMAX,EOSJMAX) ) 
    ALLOCATE( EOSDGT(EOSN+1,EOSN+1,EOSIMAX,EOSJMAX) )  
#ifdef AVX512
    !DIR$ ASSUME_ALIGNED  EOSDGT(1,1,1,1):64, EOSDGT(1,1,1,1):64
#else                                         
    !DIR$ ASSUME_ALIGNED  EOSDGT(1,1,1,1):32, EOSDGT(1,1,1,1):32
#endif
    EOSdrho = EOS%drho 
    EOSde   = EOS%de  
    !
    DO j = 1, EOS%JMAX
     DO i = 1, EOS%IMAX
         READ(333) EOS%DGp(:,:,i,j) 
         READ(333) EOS%DGT(:,:,i,j) 
         EOSDGp(:,:,i,j)=EOS%DGp(:,:,i,j) 
         EOSDGT(:,:,i,j)=EOS%DGT(:,:,i,j)
     ENDDO
    ENDDO 
    ! If original table is also requested, then read it 
    IF(OriginalTable==1) THEN
        READ(333) EOS%NRHO, EOS%NT 
        ALLOCATE( EOS%T(EOS%NT), EOS%rho(EOS%NRHO) ) 
        ALLOCATE( EOS%p(EOS%NRHO, EOS%NT), EOS%e(EOS%NRHO, EOS%NT) ) 
        EOS%drhoin = (EOS%rhomax-EOS%rhomin)/REAL(EOS%NRHO-1) 
        EOS%dTin   = (EOS%Tmax-EOS%Tmin)/REAL(EOS%NT-1) 
        READ(333) EOS%T
        READ(333) EOS%rho 
        READ(333) EOS%p
        READ(333) EOS%e         
    ENDIF         
    CLOSE(333) 
    PRINT *, ' Done. ' 
    !
    Lag_U_c = 0.
    Lag_L_c = 0.
    !
    Lag_U_c(1,1) = 1.0
    Lag_U_c(2,1) = -EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6) - EOSxiGP(5) - EOSxiGP(4) - EOSxiGP(3) - EOSxiGP(2)
    Lag_U_c(3,1) = EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6)) * EOSxiGP(5) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6) - EOSxiGP(5)) * EOSxiGP(4) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6) - EOSxiGP(5) - EOSxiGP(4)) * EOSxiGP(3) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6) - EOSxiGP(5) - EOSxiGP(4) - EOSxiGP(3)) * EOSxiGP(2)
    Lag_U_c(4,1) = -EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6)) * EOSxiGP(5) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6)) * EOSxiGP(5)) * EOSxiGP(4) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6)) * EOSxiGP(5) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6) - EOSxiGP(5)) * EOSxiGP(4)) * EOSxiGP(3) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6)) * EOSxiGP(5) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6) - EOSxiGP(5)) * EOSxiGP(4) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6) - EOSxiGP(5) - EOSxiGP(4)) * EOSxiGP(3)) * EOSxiGP(2)
    Lag_U_c(5,1) = EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) * EOSxiGP(5) - (-EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6)) * EOSxiGP(5)) * EOSxiGP(4) - (-EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6)) * EOSxiGP(5) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6)) * EOSxiGP(5)) * EOSxiGP(4)) * EOSxiGP(3) - (-EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6)) * EOSxiGP(5) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6)) * EOSxiGP(5)) * EOSxiGP(4) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6)) * EOSxiGP(5) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6) - EOSxiGP(5)) * EOSxiGP(4)) * EOSxiGP(3)) * EOSxiGP(2)
    Lag_U_c(6,1) = -EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) * EOSxiGP(5) * EOSxiGP(4) - (EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) * EOSxiGP(5) - (-EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6)) * EOSxiGP(5)) * EOSxiGP(4)) * EOSxiGP(3) - (EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) * EOSxiGP(5) - (-EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6)) * EOSxiGP(5)) * EOSxiGP(4) - (-EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6)) * EOSxiGP(5) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6)) * EOSxiGP(5)) * EOSxiGP(4)) * EOSxiGP(3)) * EOSxiGP(2)
    Lag_U_c(7,1) = EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) * EOSxiGP(5) * EOSxiGP(4) * EOSxiGP(3) - (-EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) * EOSxiGP(5) * EOSxiGP(4) - (EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) * EOSxiGP(5) - (-EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6)) * EOSxiGP(5)) * EOSxiGP(4)) * EOSxiGP(3)) * EOSxiGP(2)
    Lag_U_c(8,1) = -EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) * EOSxiGP(5) * EOSxiGP(4) * EOSxiGP(3) * EOSxiGP(2)
    Lag_U_c(2,2) = 1.0
    Lag_U_c(3,2) = -EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6) - EOSxiGP(5) - EOSxiGP(4) - EOSxiGP(3)
    Lag_U_c(4,2) = EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6)) * EOSxiGP(5) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6) - EOSxiGP(5)) * EOSxiGP(4) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6) - EOSxiGP(5) - EOSxiGP(4)) * EOSxiGP(3)
    Lag_U_c(5,2) = -EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6)) * EOSxiGP(5) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6)) * EOSxiGP(5)) * EOSxiGP(4) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6)) * EOSxiGP(5) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6) - EOSxiGP(5)) * EOSxiGP(4)) * EOSxiGP(3)
    Lag_U_c(6,2) = EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) * EOSxiGP(5) - (-EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6)) * EOSxiGP(5)) * EOSxiGP(4) - (-EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6)) * EOSxiGP(5) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6)) * EOSxiGP(5)) * EOSxiGP(4)) * EOSxiGP(3)
    Lag_U_c(7,2) = -EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) * EOSxiGP(5) * EOSxiGP(4) - (EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) * EOSxiGP(5) - (-EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6)) * EOSxiGP(5)) * EOSxiGP(4)) * EOSxiGP(3)
    Lag_U_c(8,2) = EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) * EOSxiGP(5) * EOSxiGP(4) * EOSxiGP(3)
    Lag_U_c(3,3) = 1.0
    Lag_U_c(4,3) = -EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6) - EOSxiGP(5) - EOSxiGP(4)
    Lag_U_c(5,3) = EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6)) * EOSxiGP(5) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6) - EOSxiGP(5)) * EOSxiGP(4)
    Lag_U_c(6,3) = -EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6)) * EOSxiGP(5) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6)) * EOSxiGP(5)) * EOSxiGP(4)
    Lag_U_c(7,3) = EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) * EOSxiGP(5) - (-EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6)) * EOSxiGP(5)) * EOSxiGP(4)
    Lag_U_c(8,3) = -EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) * EOSxiGP(5) * EOSxiGP(4)
    Lag_U_c(4,4) = 1.0
    Lag_U_c(5,4) = -EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6) - EOSxiGP(5)
    Lag_U_c(6,4) = EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6) - (-EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6)) * EOSxiGP(5)
    Lag_U_c(7,4) = -EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) - (EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6)) * EOSxiGP(5)
    Lag_U_c(8,4) = EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6) * EOSxiGP(5)
    Lag_U_c(5,5) = 1.0
    Lag_U_c(6,5) = -EOSxiGP(8) - EOSxiGP(7) - EOSxiGP(6)
    Lag_U_c(7,5) = EOSxiGP(8) * EOSxiGP(7) - (-EOSxiGP(8) - EOSxiGP(7)) * EOSxiGP(6)
    Lag_U_c(8,5) = -EOSxiGP(8) * EOSxiGP(7) * EOSxiGP(6)
    Lag_U_c(6,6) = 1.0
    Lag_U_c(7,6) = -EOSxiGP(8) - EOSxiGP(7)
    Lag_U_c(8,6) = EOSxiGP(8) * EOSxiGP(7)
    Lag_U_c(7,7) = 1.0
    Lag_U_c(8,7) = -EOSxiGP(8)
    Lag_U_c(8,8) = 1.0
    !
    Lag_L_c(8,1) = 1.0
    Lag_L_c(7,2) = 1.0
    Lag_L_c(8,2) = -EOSxiGP(1)
    Lag_L_c(6,3) = 1.0
    Lag_L_c(7,3) = -EOSxiGP(1) - EOSxiGP(2)
    Lag_L_c(8,3) = EOSxiGP(1) * EOSxiGP(2)
    Lag_L_c(5,4) = 1.0
    Lag_L_c(6,4) = -EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3)
    Lag_L_c(7,4) = EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3)
    Lag_L_c(8,4) = -EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3)
    Lag_L_c(4,5) = 1.0
    Lag_L_c(5,5) = -EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3) - EOSxiGP(4)
    Lag_L_c(6,5) = EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3)) * EOSxiGP(4)
    Lag_L_c(7,5) = -EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3)) * EOSxiGP(4)
    Lag_L_c(8,5) = EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) * EOSxiGP(4)
    Lag_L_c(3,6) = 1.0
    Lag_L_c(4,6) = -EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3) - EOSxiGP(4) - EOSxiGP(5)
    Lag_L_c(5,6) = EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3)) * EOSxiGP(4) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3) - EOSxiGP(4)) * EOSxiGP(5)
    Lag_L_c(6,6) = -EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3)) * EOSxiGP(4) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3)) * EOSxiGP(4)) * EOSxiGP(5)
    Lag_L_c(7,6) = EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) * EOSxiGP(4) - (-EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3)) * EOSxiGP(4)) * EOSxiGP(5)
    Lag_L_c(8,6) = -EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) * EOSxiGP(4) * EOSxiGP(5)
    Lag_L_c(2,7) = 1.0
    Lag_L_c(3,7) = -EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3) - EOSxiGP(4) - EOSxiGP(5) - EOSxiGP(6)
    Lag_L_c(4,7) = EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3)) * EOSxiGP(4) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3) - EOSxiGP(4)) * EOSxiGP(5) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3) - EOSxiGP(4) - EOSxiGP(5)) * EOSxiGP(6)
    Lag_L_c(5,7) = -EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3)) * EOSxiGP(4) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3)) * EOSxiGP(4)) * EOSxiGP(5) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3)) * EOSxiGP(4) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3) - EOSxiGP(4)) * EOSxiGP(5)) * EOSxiGP(6)
    Lag_L_c(6,7) = EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) * EOSxiGP(4) - (-EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3)) * EOSxiGP(4)) * EOSxiGP(5) - (-EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3)) * EOSxiGP(4) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3)) * EOSxiGP(4)) * EOSxiGP(5)) * EOSxiGP(6)
    Lag_L_c(7,7) = -EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) * EOSxiGP(4) * EOSxiGP(5) - (EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) * EOSxiGP(4) - (-EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3)) * EOSxiGP(4)) * EOSxiGP(5)) * EOSxiGP(6)
    Lag_L_c(8,7) = EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) * EOSxiGP(4) * EOSxiGP(5) * EOSxiGP(6)
    Lag_L_c(1,8) = 1.0
    Lag_L_c(2,8) = -EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3) - EOSxiGP(4) - EOSxiGP(5) - EOSxiGP(6) - EOSxiGP(7)
    Lag_L_c(3,8) = EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3)) * EOSxiGP(4) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3) - EOSxiGP(4)) * EOSxiGP(5) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3) - EOSxiGP(4) - EOSxiGP(5)) * EOSxiGP(6) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3) - EOSxiGP(4) - EOSxiGP(5) - EOSxiGP(6)) * EOSxiGP(7)
    Lag_L_c(4,8) = -EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3)) * EOSxiGP(4) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3)) * EOSxiGP(4)) * EOSxiGP(5) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3)) * EOSxiGP(4) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3) - EOSxiGP(4)) * EOSxiGP(5)) * EOSxiGP(6) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3)) * EOSxiGP(4) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3) - EOSxiGP(4)) * EOSxiGP(5) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3) - EOSxiGP(4) - EOSxiGP(5)) * EOSxiGP(6)) * EOSxiGP(7)
    Lag_L_c(5,8) = EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) * EOSxiGP(4) - (-EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3)) * EOSxiGP(4)) * EOSxiGP(5) - (-EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3)) * EOSxiGP(4) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3)) * EOSxiGP(4)) * EOSxiGP(5)) * EOSxiGP(6) - (-EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3)) * EOSxiGP(4) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3)) * EOSxiGP(4)) * EOSxiGP(5) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3)) * EOSxiGP(4) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3) - EOSxiGP(4)) * EOSxiGP(5)) * EOSxiGP(6)) * EOSxiGP(7)
    Lag_L_c(6,8) = -EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) * EOSxiGP(4) * EOSxiGP(5) - (EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) * EOSxiGP(4) - (-EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3)) * EOSxiGP(4)) * EOSxiGP(5)) * EOSxiGP(6) - (EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) * EOSxiGP(4) - (-EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3)) * EOSxiGP(4)) * EOSxiGP(5) - (-EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3)) * EOSxiGP(4) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3) - (-EOSxiGP(1) - EOSxiGP(2) - EOSxiGP(3)) * EOSxiGP(4)) * EOSxiGP(5)) * EOSxiGP(6)) * EOSxiGP(7)
    Lag_L_c(7,8) = EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) * EOSxiGP(4) * EOSxiGP(5) * EOSxiGP(6) - (-EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) * EOSxiGP(4) * EOSxiGP(5) - (EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) * EOSxiGP(4) - (-EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) - (EOSxiGP(1) * EOSxiGP(2) - (-EOSxiGP(1) - EOSxiGP(2)) * EOSxiGP(3)) * EOSxiGP(4)) * EOSxiGP(5)) * EOSxiGP(6)) * EOSxiGP(7)
    Lag_L_c(8,8) = -EOSxiGP(1) * EOSxiGP(2) * EOSxiGP(3) * EOSxiGP(4) * EOSxiGP(5) * EOSxiGP(6) * EOSxiGP(7)
    !
    Lag_A(1) = 0.1D1 / (EOSxiGP(1) - EOSxiGP(2)) / (EOSxiGP(1) - EOSxiGP(3)) / (EOSxiGP(1) - EOSxiGP(4)) / (EOSxiGP(1) - EOSxiGP(5)) / (EOSxiGP(1) - EOSxiGP(6)) / (EOSxiGP(1) - EOSxiGP(7)) / (EOSxiGP(1) - EOSxiGP(8))
    Lag_A(2) = 0.1D1 / (EOSxiGP(2) - EOSxiGP(1)) / (EOSxiGP(2) - EOSxiGP(3)) / (EOSxiGP(2) - EOSxiGP(4)) / (EOSxiGP(2) - EOSxiGP(5)) / (EOSxiGP(2) - EOSxiGP(6)) / (EOSxiGP(2) - EOSxiGP(7)) / (EOSxiGP(2) - EOSxiGP(8))
    Lag_A(3) = 0.1D1 / (EOSxiGP(3) - EOSxiGP(1)) / (EOSxiGP(3) - EOSxiGP(2)) / (EOSxiGP(3) - EOSxiGP(4)) / (EOSxiGP(3) - EOSxiGP(5)) / (EOSxiGP(3) - EOSxiGP(6)) / (EOSxiGP(3) - EOSxiGP(7)) / (EOSxiGP(3) - EOSxiGP(8))
    Lag_A(4) = 0.1D1 / (EOSxiGP(4) - EOSxiGP(1)) / (EOSxiGP(4) - EOSxiGP(2)) / (EOSxiGP(4) - EOSxiGP(3)) / (EOSxiGP(4) - EOSxiGP(5)) / (EOSxiGP(4) - EOSxiGP(6)) / (EOSxiGP(4) - EOSxiGP(7)) / (EOSxiGP(4) - EOSxiGP(8))
    Lag_A(5) = 0.1D1 / (EOSxiGP(5) - EOSxiGP(1)) / (EOSxiGP(5) - EOSxiGP(2)) / (EOSxiGP(5) - EOSxiGP(3)) / (EOSxiGP(5) - EOSxiGP(4)) / (EOSxiGP(5) - EOSxiGP(6)) / (EOSxiGP(5) - EOSxiGP(7)) / (EOSxiGP(5) - EOSxiGP(8))
    Lag_A(6) = 0.1D1 / (EOSxiGP(6) - EOSxiGP(1)) / (EOSxiGP(6) - EOSxiGP(2)) / (EOSxiGP(6) - EOSxiGP(3)) / (EOSxiGP(6) - EOSxiGP(4)) / (EOSxiGP(6) - EOSxiGP(5)) / (EOSxiGP(6) - EOSxiGP(7)) / (EOSxiGP(6) - EOSxiGP(8))
    Lag_A(7) = 0.1D1 / (EOSxiGP(7) - EOSxiGP(1)) / (EOSxiGP(7) - EOSxiGP(2)) / (EOSxiGP(7) - EOSxiGP(3)) / (EOSxiGP(7) - EOSxiGP(4)) / (EOSxiGP(7) - EOSxiGP(5)) / (EOSxiGP(7) - EOSxiGP(6)) / (EOSxiGP(7) - EOSxiGP(8))
    Lag_A(8) = 0.1D1 / (EOSxiGP(8) - EOSxiGP(1)) / (EOSxiGP(8) - EOSxiGP(2)) / (EOSxiGP(8) - EOSxiGP(3)) / (EOSxiGP(8) - EOSxiGP(4)) / (EOSxiGP(8) - EOSxiGP(5)) / (EOSxiGP(8) - EOSxiGP(6)) / (EOSxiGP(8) - EOSxiGP(7))
    !
END SUBROUTINE ReadEOSCoefficients 
        
    
    

RECURSIVE SUBROUTINE PlotEOS 
    USE EOSTypes 
    IMPLICIT NONE 
    INTEGER :: nRealNodes, nSubPlotElem, i, j, k, c, nc, ii, jj, count, iNode, idxn(EOS%N+1,EOS%N+1), subtri(4,EOS%N**2)  
    INTEGER :: pDim = 2, pVar = 5 
    REAL(8) :: eGP(VECTORLENGTH), rhoGP(VECTORLENGTH), pGP(VECTORLENGTH), dpdrhoGP(VECTORLENGTH), dpdeGP(VECTORLENGTH), TGP(VECTORLENGTH), pe(VECTORLENGTH), Te(VECTORLENGTH)  
    REAL    :: xin(EOS%N+1) 
    ! 
    INTEGER*4, POINTER :: NData(:,:)  
    REAL(4), POINTER   :: DataArray(:,:) 
    CHARACTER(LEN=200) :: PlotFile 
    !
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64:: eGP, rhoGP, pGP, dpdrhoGP, dpdeGP, TGP, pe, Te
#else
    !DIR$ ATTRIBUTES ALIGN: 32:: eGP, rhoGP, pGP, dpdrhoGP, dpdeGP, TGP, pe, Te
#endif 
    !
    PRINT *, ' Plotting DG table and original EOS... '     
    !
    DO i = 1, EOS%N+1
        xin(i) = REAL(i-1)/REAL(EOS%N) 
    ENDDO
    xin(1) = 1e-7 
    xin(EOS%N+1) = 1.0-1e-7 
    !
    idxn = 0 
    c = 0 
    DO j = 1, EOS%N+1 
     DO i = 1, EOS%N+1 
       c = c + 1 
       idxn(i,j) = c
     ENDDO
    ENDDO
    c = 0 
    DO j = 1, EOS%N
     DO i = 1, EOS%N 
       c = c + 1 
       subtri(:,c) = (/ idxn(i,j), idxn(i+1,j), idxn(i+1,j+1), idxn(i,j+1)  /)         
     ENDDO 
    ENDDO
    !
    nRealNodes   = EOS%IMAX*EOS%JMAX*(EOS%N+1)**2 
    nSubPlotElem = EOS%IMAX*EOS%JMAX*EOS%N**2 
    !
    ALLOCATE(NData(4,nSubPlotElem))  
    ALLOCATE(DataArray(nRealNodes,pDim+pVar))     
    !
    PlotFile = 'EOS_plot.dat' 
    OPEN(UNIT=333,FILE=TRIM(PlotFile),RECL=1000) 
    ! 
    WRITE(333,*) ' VARIABLES = "rho" "e" "p" "pe" "T" "Te" "p_error"  ' 
    WRITE(333,*) ' ZONE N= ', nRealNodes, ' E= ', nSubPlotElem, 'F=FEPOINT  ET=QUADRILATERAL '   
    !
    c  = 0 
    nc = 0 
    DO j = 1, EOS%JMAX
     DO i = 1, EOS%IMAX 
        DO k = 1, EOS%N**2
            c = c + 1  
            NData(:,c) = nc + subtri(1:4,k)
        ENDDO 
        nc = nc + (EOS%N+1)**2 
     ENDDO 
    ENDDO   
    !
    iNode = 0 
    DO j = 1, EOS%JMAX
     DO i = 1, EOS%IMAX
         !
         count = 0 
         DO jj = 1, EOS%N+1
          DO ii = 1, EOS%N+1
              count = count + 1 
              rhoGP(count) = EOS%rhomin + (i-1)*EOS%drho + xin(ii)*EOS%drho 
              eGP(count)   = EOS%emin   + (j-1)*EOS%de   + xin(jj)*EOS%de 
              IF(MOD(count,VECTORLENGTH)==0) THEN
                  count = 0
                  CALL EvaluateEOS(pGP,dpdrhoGP,dpdeGP,TGP,rhoGP,eGP)
                  CALL DE2TP(Te,pe,rhoGP,eGP)
                  DO k = 1, VECTORLENGTH
                      iNode = iNode + 1 
                      DataArray(iNode,:) = (/ rhoGP(k), eGP(k), pGP(k), pe(k), TGP(k), Te(k), pe(k)-pGP(k) /) 
                  ENDDO                  
              ENDIF              
          ENDDO          
         ENDDO         
         !
     ENDDO
    ENDDO 
    ! 
    DO i = 1, nRealNodes
       WRITE(333, *) DataArray(i,:)
    ENDDO
    DO i = 1, nSubPlotElem
        WRITE(333,*) NData(:,i) 
    ENDDO    
    !
    CLOSE(333) 
    !
    DEALLOCATE(NData, DataArray) 
    !
    PRINT *, ' Done. ' 
    !
END SUBROUTINE PlotEOS 
    
    

RECURSIVE SUBROUTINE EOSBaseFunc1D(phi,phi_xi,xi,xin,N)
   USE EOSTypes  
   IMPLICIT NONE
   ! Argument list 
   INTEGER, INTENT(IN) :: N 
   REAL(8), INTENT(IN )   :: xi(VECTORLENGTH), xin(N+1)                            ! coordinate in [0,1] where to evaluate the basis and nodes that define the basis  
   REAL(8), INTENT(OUT)   :: phi(VECTORLENGTH,N+1), phi_xi(VECTORLENGTH,N+1)       ! the basis and its derivative w.r.t. xi 
   ! Local variables 
   INTEGER     :: i,j,m
   REAL(8)        :: tmp(VECTORLENGTH)   
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64:: tmp 
    !DIR$ ASSUME_ALIGNED phi   : 64
    !DIR$ ASSUME_ALIGNED phi_xi: 64
    !DIR$ ASSUME_ALIGNED xi    : 64
#else
    !DIR$ ATTRIBUTES ALIGN: 32:: tmp
    !DIR$ ASSUME_ALIGNED phi   : 32
    !DIR$ ASSUME_ALIGNED phi_xi: 32
    !DIR$ ASSUME_ALIGNED xi    : 32
#endif 
    ! 
   ! 
   ! Initialize variables 
   phi      = 1. 
   phi_xi   = 0. 
   ! Lagrange polynomial and its derivative 
   DO m = 1, N+1
      DO j = 1, N+1
         IF(j.EQ.m) CYCLE 
         phi(:,m) = phi(:,m)*(xi(:)-xin(j))/(xin(m)-xin(j)) 
      ENDDO 
      DO i = 1, N+1
         IF(i.EQ.m) CYCLE
         tmp = 1. 
         DO j = 1, N+1
            IF(j.EQ.i) CYCLE 
            IF(j.EQ.m) CYCLE 
            tmp(:) = tmp(:)*(xi(:)-xin(j))/(xin(m)-xin(j))    
         ENDDO 
         phi_xi(:,m) = phi_xi(:,m) + tmp(:)/(xin(m)-xin(i)) 
      ENDDO 
   ENDDO 
   !
END SUBROUTINE EOSBaseFunc1D




RECURSIVE SUBROUTINE EOSBaseFunc1D_bis(phi,phi_xi,xi,phi1,phi1_xi,xi1,xin,N)
   USE EOSTypes  
   IMPLICIT NONE
   ! Argument list 
   INTEGER, INTENT(IN) :: N 
   REAL(8), INTENT(IN )   :: xi(VECTORLENGTH),xin(N+1),xi1(VECTORLENGTH)                            ! coordinate in [0,1] where to evaluate the basis and nodes that define the basis  
   REAL(8), INTENT(OUT)   :: phi(VECTORLENGTH,N+1),phi_xi(VECTORLENGTH,N+1),phi1(VECTORLENGTH,N+1),phi1_xi(VECTORLENGTH,N+1)       ! the basis and its derivative w.r.t. xi 
   ! Local variables 
   INTEGER :: i,j,m
   REAL(8) :: tmp(VECTORLENGTH) ,tmp1(VECTORLENGTH)   
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64   :: tmp,tmp1 
    !DIR$ ASSUME_ALIGNED phi     : 64
    !DIR$ ASSUME_ALIGNED phi_xi  : 64
    !DIR$ ASSUME_ALIGNED xi      : 64
    !DIR$ ASSUME_ALIGNED phi1    : 64
    !DIR$ ASSUME_ALIGNED phi1_xi : 64
    !DIR$ ASSUME_ALIGNED xi1     : 64
#else
    !DIR$ ATTRIBUTES ALIGN: 32   :: tmp,tmp1
    !DIR$ ASSUME_ALIGNED phi     : 32
    !DIR$ ASSUME_ALIGNED phi_xi  : 32
    !DIR$ ASSUME_ALIGNED xi      : 32
    !DIR$ ASSUME_ALIGNED phi1    : 32
    !DIR$ ASSUME_ALIGNED phi1_xi : 32
    !DIR$ ASSUME_ALIGNED xi1     : 32
#endif 
    ! 
   ! 
   ! Initialize variables 
   phi      = 1. 
   phi_xi   = 0. 
   phi1      = 1. 
   phi1_xi   = 0. 
   ! Lagrange polynomial and its derivative 
   DO m = 1, N+1
      DO j = 1, N+1
         IF(j.EQ.m) CYCLE 
         phi(:,m) = phi(:,m)*(xi(:)-xin(j))/(xin(m)-xin(j))
         phi1(:,m) = phi1(:,m)*(xi1(:)-xin(j))/(xin(m)-xin(j)) 
      ENDDO 
      DO i = 1, N+1
         IF(i.EQ.m) CYCLE
         tmp = 1. 
         tmp1 = 1. 
         DO j = 1, N+1
            IF(j.EQ.i) CYCLE 
            IF(j.EQ.m) CYCLE 
            tmp(:) = tmp(:)*(xi(:)-xin(j))/(xin(m)-xin(j)) 
            tmp1(:) = tmp1(:)*(xi1(:)-xin(j))/(xin(m)-xin(j))    
         ENDDO 
         phi_xi(:,m) = phi_xi(:,m) + tmp(:)/(xin(m)-xin(i)) 
         phi1_xi(:,m) = phi1_xi(:,m) + tmp1(:)/(xin(m)-xin(i)) 
      ENDDO 
   ENDDO 
   !
END SUBROUTINE EOSBaseFunc1D_bis


RECURSIVE SUBROUTINE EOSBaseFunc1D_bis7(phi,phi_xi,xi,phi1,phi1_xi,xi1,xin)
   USE EOSTypes, ONLY: VECTORLENGTH,Lag_U_c,Lag_L_c,Lag_A
   IMPLICIT NONE
   ! Argument list 
   INTEGER, PARAMETER ::N=7
   REAL(8), INTENT(IN )   :: xi(VECTORLENGTH),xin(N+1),xi1(VECTORLENGTH)                            ! coordinate in [0,1] where to evaluate the basis and nodes that define the basis  
   REAL(8), INTENT(OUT)   :: phi(VECTORLENGTH,N+1),phi_xi(VECTORLENGTH,N+1),phi1(VECTORLENGTH,N+1),phi1_xi(VECTORLENGTH,N+1)       ! the basis and its derivative w.r.t. xi 
   ! Local variables 
   INTEGER :: i,j,m
   REAL(8) :: tmp(VECTORLENGTH) ,tmp1(VECTORLENGTH)
   REAL(8) :: Lag_U(VECTORLENGTH,N+1),Lag_L(VECTORLENGTH,N+1),Lag_U_xi(VECTORLENGTH,N+1),Lag_L_xi(VECTORLENGTH,N+1)
   REAL(8) :: Lag_U1(VECTORLENGTH,N+1),Lag_L1(VECTORLENGTH,N+1),Lag_U1_xi(VECTORLENGTH,N+1),Lag_L1_xi(VECTORLENGTH,N+1)
#ifdef AVX512                                         
    !DIR$ ATTRIBUTES ALIGN: 64   :: tmp,tmp1,Lag_U,Lag_L,Lag_U_xi,Lag_L_xi,Lag_U1,Lag_L1,Lag_U1_xi,Lag_L1_xi      
    !DIR$ ASSUME_ALIGNED phi     : 64               
    !DIR$ ASSUME_ALIGNED phi_xi  : 64
    !DIR$ ASSUME_ALIGNED xi      : 64
    !DIR$ ASSUME_ALIGNED phi1    : 64
    !DIR$ ASSUME_ALIGNED phi1_xi : 64
    !DIR$ ASSUME_ALIGNED xi1     : 64
#else
    !DIR$ ATTRIBUTES ALIGN: 32   :: tmp,tmp1,Lag_U,Lag_L,Lag_U_xi,Lag_L_xi,Lag_U1,Lag_L1,Lag_U1_xi,Lag_L1_xi  
    !DIR$ ASSUME_ALIGNED phi     : 32
    !DIR$ ASSUME_ALIGNED phi_xi  : 32
    !DIR$ ASSUME_ALIGNED xi      : 32
    !DIR$ ASSUME_ALIGNED phi1    : 32
    !DIR$ ASSUME_ALIGNED phi1_xi : 32
    !DIR$ ASSUME_ALIGNED xi1     : 32
#endif
   !
   ! Initialize variables 
   !phi      = 1. 
   !phi_xi   = 0. 
   !phi1      = 1. 
   !phi1_xi   = 0. 
   !! Lagrange polynomial and its derivative 
   !DO m = 1, N+1
   !   DO j = 1, N+1
   !      IF(j.EQ.m) CYCLE 
   !      phi(:,m) = phi(:,m)*(xi(:)-xin(j))/(xin(m)-xin(j))
   !      phi1(:,m) = phi1(:,m)*(xi1(:)-xin(j))/(xin(m)-xin(j)) 
   !   ENDDO 
   !   DO i = 1, N+1
   !      IF(i.EQ.m) CYCLE
   !      tmp = 1. 
   !      tmp1 = 1. 
   !      DO j = 1, N+1
   !         IF(j.EQ.i) CYCLE 
   !         IF(j.EQ.m) CYCLE 
   !         tmp(:) = tmp(:)*(xi(:)-xin(j))/(xin(m)-xin(j)) 
   !         tmp1(:) = tmp1(:)*(xi1(:)-xin(j))/(xin(m)-xin(j))    
   !      ENDDO 
   !      phi_xi(:,m) = phi_xi(:,m) + tmp(:)/(xin(m)-xin(i)) 
   !      phi1_xi(:,m) = phi1_xi(:,m) + tmp1(:)/(xin(m)-xin(i)) 
   !   ENDDO 
   !ENDDO 
   !
   DO m = 1, N+1
      Lag_U(:,m)    = Lag_U_c(1,m)  ! this is the coefficient of xi^(N+1)
      Lag_L(:,m)    = Lag_L_c(1,m)  ! this is the coefficient of xi^(N+1)
      Lag_U_xi(:,m) = real(N)*Lag_U_c(1,m)  ! this is the coefficient of xi^(N+1)
      Lag_L_xi(:,m) = real(N)*Lag_L_c(1,m)  ! this is the coefficient of xi^(N+1)
      Lag_U1(:,m)   = Lag_U_c(1,m)  ! this is the coefficient of xi^(N+1)
      Lag_L1(:,m)   = Lag_L_c(1,m)  ! this is the coefficient of xi^(N+1)
      Lag_L1_xi(:,m)= real(N)*Lag_L_c(1,m)  ! this is the coefficient of xi^(N+1)
      Lag_U1_xi(:,m)= real(N)*Lag_U_c(1,m)  ! this is the coefficient of xi^(N+1)
      DO j = 2, N
            Lag_U(:,m)    = Lag_U_c(j,m) + xi(:)*Lag_U(:,m)
            Lag_L(:,m)    = Lag_L_c(j,m) + xi(:)*Lag_L(:,m)
            Lag_U_xi(:,m) = real(N+1-j)*Lag_U_c(j,m) + xi(:)*Lag_U_xi(:,m)                
            Lag_L_xi(:,m) = real(N+1-j)*Lag_L_c(j,m) + xi(:)*Lag_L_xi(:,m)    
            Lag_U1(:,m)   = Lag_U_c(j,m) + xi1(:)*Lag_U1(:,m)
            Lag_L1(:,m)   = Lag_L_c(j,m) + xi1(:)*Lag_L1(:,m)
            Lag_U1_xi(:,m)= real(N+1-j)*Lag_U_c(j,m) + xi1(:)*Lag_U1_xi(:,m)                
            Lag_L1_xi(:,m)= real(N+1-j)*Lag_L_c(j,m) + xi1(:)*Lag_L1_xi(:,m)                
      ENDDO
      j=N+1
      Lag_U(:,m) = Lag_U_c(j,m) + xi(:)*Lag_U(:,m)
      Lag_L(:,m) = Lag_L_c(j,m) + xi(:)*Lag_L(:,m)
      Lag_U1(:,m) = Lag_U_c(j,m) + xi1(:)*Lag_U1(:,m)
      Lag_L1(:,m) = Lag_L_c(j,m) + xi1(:)*Lag_L1(:,m)
      !
      phi(:,m)      = Lag_U(:,m)*Lag_L(:,m)*Lag_A(m)
      phi_xi(:,m)   = (Lag_U_xi(:,m)*Lag_L(:,m) + Lag_U(:,m)*Lag_L_xi(:,m)) *Lag_A(m) 
      phi1(:,m)     = Lag_U1(:,m)*Lag_L1(:,m)*Lag_A(m)
      phi1_xi(:,m)  = (Lag_U1_xi(:,m)*Lag_L1(:,m) + Lag_U1(:,m)*Lag_L1_xi(:,m)) *Lag_A(m) 
   ENDDO 
   !
END SUBROUTINE EOSBaseFunc1D_bis7





RECURSIVE SUBROUTINE EOSBaseFunc1D_scalar(phi,phi_xi,xi,xin,N)
   USE EOSTypes  
   IMPLICIT NONE
   ! Argument list 
   INTEGER, INTENT(IN) :: N 
   REAL, INTENT(IN )   :: xi, xin(N+1)                            ! coordinate in [0,1] where to evaluate the basis and nodes that define the basis  
   REAL, INTENT(OUT)   :: phi(N+1), phi_xi(N+1)       ! the basis and its derivative w.r.t. xi 
   ! Local variables
   INTEGER     :: i,j,m
   REAL        :: tmp
   ! 
   ! Initialize variables 
   phi      = 1. 
   phi_xi   = 0. 
   ! Lagrange polynomial and its derivative 
   DO m = 1, N+1
      DO j = 1, N+1
         IF(j.EQ.m) CYCLE 
         phi(m) = phi(m)*(xi-xin(j))/(xin(m)-xin(j)) 
      ENDDO 
      DO i = 1, N+1
         IF(i.EQ.m) CYCLE
         tmp = 1. 
         DO j = 1, N+1
            IF(j.EQ.i) CYCLE 
            IF(j.EQ.m) CYCLE 
            tmp = tmp*(xi-xin(j))/(xin(m)-xin(j))    
         ENDDO 
         phi_xi(m) = phi_xi(m) + tmp/(xin(m)-xin(i)) 
      ENDDO 
   ENDDO 
   !
END SUBROUTINE EOSBaseFunc1D_scalar




    
END MODULE GREOS