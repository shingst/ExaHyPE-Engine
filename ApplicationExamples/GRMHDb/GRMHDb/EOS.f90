
#define C2PFF    
    
    
MODULE EOS_mod
    USE Parameters, ONLY: VECTORLENGTH,NSTOV_kappa,EQN  
    USE EOS_par 
    IMPLICIT NONE 
    private
    !
    !INTERFACE EOS_init
    !    MODULE PROCEDURE EOS_init
    !END INTERFACE 
    !! 
    !INTERFACE BaseFunc1D_Neos
    !    MODULE PROCEDURE BaseFunc1D_Neos
    !END INTERFACE 
    ! 
    INTERFACE C2P_GEOS
        MODULE PROCEDURE C2P_GEOS
    END INTERFACE 
    !  
    INTERFACE FUNC_eps_rho_C2P_GEOS
        MODULE PROCEDURE FUNC_eps_rho_C2P_GEOS
    END INTERFACE 
    ! 
    INTERFACE C2P_GEOS_3D 
        MODULE PROCEDURE C2P_GEOS_3D
    END INTERFACE 
    ! 
    INTERFACE FUNC_C2P_GEOS_3D
        MODULE PROCEDURE FUNC_C2P_GEOS_3D
    END INTERFACE 
    ! 
    INTERFACE EOS_eps_depsdrho_depsdp
        MODULE PROCEDURE EOS_eps_depsdrho_depsdp
    END INTERFACE 
    !
    INTERFACE EOS_eps_depsdrho_depsdp_vector
        MODULE PROCEDURE EOS_eps_depsdrho_depsdp_vector
    END INTERFACE 
    !  
    INTERFACE FUNC_C2P_GEOS
        MODULE PROCEDURE FUNC_C2P_GEOS
    END INTERFACE 
    !  
    INTERFACE dFUNC_C2P_GEOS
        MODULE PROCEDURE dFUNC_C2P_GEOS
    END INTERFACE 
    !
    INTERFACE EOS_p_dprho_dpeps
        MODULE PROCEDURE EOS_p_dprho_dpeps
    END INTERFACE 
    ! 
    INTERFACE EOS_eps
        MODULE PROCEDURE EOS_eps
    END INTERFACE 
    !  
    INTERFACE EOS_EPS_P2C
        MODULE PROCEDURE EOS_EPS_P2C
    END INTERFACE 
    !   
    INTERFACE RTSAFE_C2P_RHD1
            MODULE PROCEDURE RTSAFE_C2P_RHD1
    END INTERFACE
    !
    INTERFACE FUNC_C2P_RHD1
        MODULE PROCEDURE FUNC_C2P_RHD1
    END INTERFACE
    !
    ! VECTORIZED VERSION:
    ! 
    INTERFACE C2P_GEOS_vector
        MODULE PROCEDURE C2P_GEOS_vector
    END INTERFACE 
    !  
    INTERFACE FUNC_eps_rho_C2P_GEOS_vector
        MODULE PROCEDURE FUNC_eps_rho_C2P_GEOS_vector
    END INTERFACE 
    ! 
    
    INTERFACE FUNC_C2P_GEOS_vector
        MODULE PROCEDURE FUNC_C2P_GEOS_vector
    END INTERFACE 
    !  
    INTERFACE dFUNC_C2P_GEOS_vector
        MODULE PROCEDURE dFUNC_C2P_GEOS_vector
    END INTERFACE 
    !
    INTERFACE EOS_p_dprho_dpeps_vector
        MODULE PROCEDURE EOS_p_dprho_dpeps_vector
    END INTERFACE 
    ! 
    INTERFACE EOS_eps_vector
        MODULE PROCEDURE EOS_eps_vector
    END INTERFACE 
    !  
    INTERFACE EOS_EPS_P2C_vector
        MODULE PROCEDURE EOS_EPS_P2C_vector
    END INTERFACE 
    !
    INTERFACE RTSAFE_C2P_RHD1_VECTOR
        MODULE PROCEDURE RTSAFE_C2P_RHD1_VECTOR
    END INTERFACE
    !
    INTERFACE FUNC_C2P_RHD1_VECTORF
        MODULE PROCEDURE FUNC_C2P_RHD1_VECTORF
    END INTERFACE
    !
    !
    PUBLIC :: C2P_GEOS          ,EOS_p_dprho_dpeps          ,EOS_eps            ,dFUNC_C2P_GEOS,         EOS_EPS_P2C,                   &
              C2P_GEOS_vector   ,EOS_p_dprho_dpeps_vector   ,EOS_eps_vector     ,dFUNC_C2P_GEOS_vector,  EOS_EPS_P2C_vector,            &
              C2P_GEOS_3D, FUNC_C2P_GEOS_3D,  EOS_eps_depsdrho_depsdp,          &
                                            EOS_eps_depsdrho_depsdp_vector, &
            RTSAFE_C2P_RHD1,FUNC_C2P_RHD1, &
            RTSAFE_C2P_RHD1_VECTOR,FUNC_C2P_RHD1_VECTORF  
    !
    CONTAINS



 

RECURSIVE SUBROUTINE C2P_GEOS_3D(w,p,v2,gam,d,e,s2,b2,sb2,FAILED)
  !***********************************************************
  ! given w=rho*h*W^2, solve the 2d nonlinear system for x=(w,p) 
  !*****
  USE Astromod
  IMPLICIT NONE
  ! arguments
  REAL :: w,p,v2,gam,d,e,s2,b2,sb2
  LOGICAL :: FAILED
  INTENT (INOUT) :: v2,w,p
  !INTENT (IN)    :: gam,d,e,s2,b2,sb2
  INTENT (INOUT)    :: gam,d,e,s2,b2,sb2
  ! local variables
  INTEGER:: i,iloc
  REAL :: alpha,xnew(3),fnew(3),fmod,fmodnew
  REAL :: dx(3),x(3),xtmp(3),f(3),df(3,3),idf(3,3),det,drhodx,dydx,depsdrho,depsdp,eps,v20,p0,w0
  REAL :: dFtmp,dx1tmp,idet,df2(2,2),idf2(2,2),dx2(3)
  INTEGER :: inner
  INTEGER, PARAMETER :: innerMaxNewton = 500
  !  
  REAL, PARAMETER :: third=1./3.
  INTEGER, PARAMETER :: MaxNewton=500
  REAL, PARAMETER :: tolerance = 1e-12 
  REAL, PARAMETER :: tolthreshold = 1e-10
  FAILED = .FALSE.
  !v2=x0
  ! check the value fo 
  x(1) = v2
  x(2) = w
  x(3) = p
  v20 = v2
  w0 = w
  p0 = p
  !
  DO i=1,MaxNewton
        ! 
        CALL FUNC_C2P_GEOS_3D(x,f,df,gam,d,e,s2,b2,sb2)
        !
        CALL MatrixInverse3x3(df,idf,det) 
        !
        dx(1) = -idf(1,1)*f(1)-idf(1,2)*f(2)-idf(1,3)*f(3)     ! = dv2
        dx(2) = -idf(2,1)*f(1)-idf(2,2)*f(2)-idf(2,3)*f(3)     ! = dp
        dx(3) = -idf(3,1)*f(1)-idf(3,2)*f(2)-idf(3,3)*f(3)     ! = dp
        ! 
        !IF(SUM(ABS(dx)).LT.tolerance) THEN
        !    continue
        !    RETURN
        !ENDIF
        !
        fmod=SUM(ABS(f))
        ! 
        alpha = 1.0 
        !
        DO iloc=1,20 
            xnew = x + alpha*dx
            IF(ANY(xnew.LT.0).OR.xnew(1).GT.1.0) THEN
                alpha = alpha*0.5
            ELSE
                exit
            ENDIF
            !
        END DO
        !
        !xnew=abs(xnew)
        !
        DO iloc=1,20  
            xnew = ABS(x + alpha*dx)
            CALL FUNC_C2P_GEOS_3D(xnew,fnew,df,gam,d,e,s2,b2,sb2)
            fmodnew=SUM(ABS(fnew))
            IF(fmodnew.GT.fmod) THEN 
                alpha = alpha*0.5 
            ELSE
                exit
            ENDIF 
        END DO
        !
        f=fnew
        x=xnew
        fmod=fmodnew
        !
        IF(fmodnew.GT.tolthreshold.AND.alpha.LT.1e-5) THEN
            ! 2d search: 
              DO inner=1,40 !innerMaxNewton
                    ! 
                    CALL FUNC_C2P_GEOS_3D(x,f,df,gam,d,e,s2,b2,sb2)
                    !
                    det = df(1,1)*df(2,2)-df(1,2)*df(2,1)
                    idet = 1.0/det
                    !
                    idf2(1,1) = df(2,2)
                    idf2(1,2) = -df(1,2)
                    idf2(2,1) = -df(2,1)
                    idf2(2,2) = df(1,1)  
                    !
                    dx2(1) = -idf2(1,1)*f(1)-idf2(1,2)*f(2) !-idf(1,3)*f(3)     ! = dv2
                    dx2(2) = -idf2(2,1)*f(1)-idf2(2,2)*f(2) !-idf(2,3)*f(3)     ! = dw 
                    dx2(3) = 0.
                    !  
                    fmod=SUM(ABS(f))
                    !
                    alpha = 1.0 
                    !
                    DO iloc=1,20 
                        xnew = x + alpha*dx2
                        IF(ANY(xnew(1:2).LT.0).OR.xnew(1).GT.1.0) THEN
                            alpha = alpha*0.5
                        ELSE
                            exit
                        ENDIF
                        !
                    END DO
                    !
                    !xnew=abs(xnew)
                    !
                    DO iloc=1,20  
                        xnew = ABS(x + alpha*dx2)
                        CALL FUNC_C2P_GEOS_3D(xnew,fnew,df,gam,d,e,s2,b2,sb2)
                        fmodnew=SUM(ABS(fnew))
                        IF(fmodnew.GT.fmod) THEN 
                            alpha = alpha*0.5 
                        ELSE
                            exit
                        ENDIF 
                    END DO
                    !
                    f=fnew
                    x=xnew
                    fmod=fmodnew
                    !
              ENDDO
              !
        ENDIF
        !
        IF(ANY(ISNAN(dx))) THEN
            continue
        ENDIF
        !
        IF(ANY(ISNAN(f))) THEN
            continue
        ENDIF
        !
        CALL MatrixInverse3x3(df,idf,det) 
        !
        IF(ANY(ISNAN(df))) THEN
            continue
        ENDIF
        !
        IF(ANY(ISNAN(idf))) THEN
            continue
        ENDIF
        !
        dx(1) = -idf(1,1)*f(1)-idf(1,2)*f(2)-idf(1,3)*f(3)     ! = dv2
        dx(2) = -idf(2,1)*f(1)-idf(2,2)*f(2)-idf(2,3)*f(3)     ! = dp
        dx(3) = -idf(3,1)*f(1)-idf(3,2)*f(2)-idf(3,3)*f(3)     ! = dp
        !
        IF(SUM(ABS(dx)).LT.tolerance) THEN
            continue
            IF(SUM(ABS(f)).LT.tolerance) THEN
                v2= x(1)
                w = x(2)
                p = x(3)
                continue
                RETURN
            ELSE
                continue
            ENDIF
            !
        ENDIF
  ENDDO
  !
  IF(SUM(ABS(f)).GT.tolthreshold) THEN
      continue
      FAILED = .TRUE.
  ENDIF
  !
  v2= x(1)
  w = x(2)
  p = x(3)
  !
  !
  continue
  !
END SUBROUTINE C2P_GEOS_3D

                           

RECURSIVE SUBROUTINE FUNC_C2P_GEOS_3D(x,f,df,gam,d,e,s2,b2,sb2)
  IMPLICIT NONE
  ! arguments
  REAL :: x(3),f(3),df(3,3)
  REAL :: gam,d,e,s2,b2,sb2 
  INTENT (INOUT) :: x,f,df
  INTENT (IN)    :: gam,d,e,s2,b2,sb2
  ! local variables
  REAL, PARAMETER :: third=1./3. 
  REAL, PARAMETER ::  i27=1./27. 
  INTEGER    :: i
  REAL :: v2,w,p,det,drhodx,dydx,depsdrho,depsdp,eps 
  REAL :: rho,c0,c2,c3,dw,dc2,dc3,dlogw,wb,vb2,wb2 
  REAL :: sqrtdet,q,w0,wtest,c2ic3,c2ic33,Lor
  v2=x(1)
  w=x(2)
  p=x(3)
  rho=d*sqrt(1.-v2) 
  wb    = w + b2
  wb2   = wb**2
  vb2   = sb2 / w**2
  Lor = 1.0/SQRT(1.-v2)
  ! 
  f(1) = wb2*v2-vb2*(wb+w)-s2
  !
  c3=1. !-gam*(1.-v2)
  c2=.5*b2*(1.+v2)-e-p
  c0=-0.5*sb2
  !
  f(2) = c3*w**3 + c2*w**2 + c0
  !
  CALL EOS_eps_depsdrho_depsdp(eps,depsdrho,depsdp,rho,p) 
  !
  f(3) = p + eps*rho + rho - (1.-v2)*w 
  !
  ! compute the jacobian: (d/dw, d/dp)[ (/ f(1), f(2) /) ]
  !
  df(1,1) = wb2
  df(1,2) = 2.0*( wb*v2 + vb2*(wb+w)/w -vb2)
  df(1,3) = 0.0
  !
  df(2,1) = 0.5*b2*w**2
  df(2,2) = (3*c3*w + 2*c2)*w
  df(2,3) = -w**2
  !
  df(3,1) = w-0.5*D*Lor*(1.0 + eps + rho*depsdrho )
  df(3,2) = v2-1.0
  df(3,3) = 1.0+rho*depsdp
  !
        !IF(ANY(ISNAN(f))) THEN
        !    continue
        !ENDIF
  continue
  !
END SUBROUTINE FUNC_C2P_GEOS_3D
    


RECURSIVE SUBROUTINE dFUNC_C2P_GEOS(dF1,dp1rho,dp1eps,p,rho,eps,tau,D,s2)
    USE EOS_par
    !USE EOS_mod
    IMPLICIT NONE
    ! input variables
    REAL, INTENT(IN) :: dp1rho,dp1eps,p,rho,eps,tau,D,s2
    REAL, INTENT(out) :: dF1
    REAL :: drhop,depsp,sqrtX,den,den2
    !
    !CALL EOS_p_dprho_dpeps(p,rho,eps) 
    !
    den = tau+p+D
    den2 = den**2
    sqrtX = SQRT(den2-s2)
    !
    drhop = D*s2/den2/sqrtX
    depsp = p*s2/(den2-s2)/den/rho
    !
    SELECT CASE(EOSfunc)
    CASE(0) ! p-p[rho,eps] 
        dF1 = 1.0 - dp1rho*drhop - dp1eps*depsp
    CASE(1) ! p^2-F[rho,eps] 
        dF1 = 2*p - dp1rho*drhop - dp1eps*depsp
     CASE DEFAULT
    	dF1 = 1.0 - dp1rho*drhop - dp1eps*depsp
    END SELECt
    
    !
END SUBROUTINE dFUNC_C2P_GEOS


RECURSIVE SUBROUTINE C2P_GEOS(p,rho,eps,p0,tau,D,s2,FAILED)
    !
    ! this subroutine  (tau,D,s2) => (p,rho,eps)
    ! given an equation of state (analytical or tabulated)
    ! by means of a Newton method, looking for the root of the functional F1=F1[p,rho(p),eps(p)]:= p - EOS(p,rho(p),eps(p))
    !
    IMPLICIT NONE
    ! input variables
    REAL, INTENT(IN) :: tau,D,s2
    REAL, INTENT(INOUT) :: p0
    !REAL, INTENT(OUT) :: p0
    ! output variables
    REAL :: p,rho,eps
    ! local variables
    REAL :: F1,dF1,dp,dp1rho,dp1eps,eps0,rho0,alpha,res,Newres
    LOGICAL :: FAILED
    INTEGER :: Nt,inNt
    INTEGER, PARAMETER :: maxNt = 500
    REAL :: tol,tolthreshold
    FAILED = .FALSE.
    ! here is the initial guess:
    p=p0
    CALL FUNC_eps_rho_C2P_GEOS(eps0,rho0,p,tau,D,s2)  
    !
    rho=rho0
    eps=eps0
    tol = 1e-18
    tolthreshold = 1e-10
    !
    ! Newton algorithm:
    DO Nt = 1, maxNt
        CALL FUNC_C2P_GEOS(F1,dp1rho,dp1eps,p,rho,eps)        
        IF(ABS(F1).LT.tol) THEN
            continue
            RETURN
        ENDIF
        !
        CALL dFUNC_C2P_GEOS(dF1,dp1rho,dp1eps,p,rho,eps,tau,D,s2)
        !
        dp=-F1/dF1
        p=p+dp
        CALL FUNC_eps_rho_C2P_GEOS(eps,rho,p,tau,D,s2)  
        !
		CYCLE
		!
        res = ABS(F1)
        IF(res.LT.tol) THEN
            continue
            RETURN
        ENDIF
        !
        CALL dFUNC_C2P_GEOS(dF1,dp1rho,dp1eps,p,rho,eps,tau,D,s2)
        IF(ISNAN(dF1)) THEN
            EXIT
        ENDIF
        dp=-F1/dF1
        alpha=1.0
        DO inNt=1,20
            IF(p+alpha*dp<1e-10) THEN
                alpha = alpha / 2.0 
            ELSE
                CALL FUNC_C2P_GEOS(F1,dp1rho,dp1eps,p+alpha*dp,rho,eps)
                newres = ABS(F1) 
                IF(newres.GE.res) THEN
                    alpha = alpha / 2.0
                ELSE
                    EXIT 
                ENDIF 
            ENDIF          
        ENDDO      
        p=ABS(p+alpha*dp)
        CALL FUNC_eps_rho_C2P_GEOS(eps,rho,p,tau,D,s2)  
        !
    ENDDO
    ! 
    IF(ABS(F1).GT.tolthreshold.OR.ISNAN(F1)) THEN
        continue
    	FAILED = .TRUE.
    ENDIF
    !
    !
    continue
    !
END SUBROUTINE C2P_GEOS




RECURSIVE SUBROUTINE FUNC_eps_rho_C2P_GEOS(eps,rho,p,tau,D,s2) 
    ! this function evaluate rho =rho(p), eps = eps(p)
    USE EOS_par
    !USE EOS_mod
    IMPLICIT NONE
    ! input variables
    REAL, INTENT(IN) :: p,tau,D,s2
    REAL, INTENT(out) :: eps,rho
    real :: den,tmp
    ! 
    ! this is true only for pure GR-Hydro, independently on the EOS.
    den = tau+p+D 
    tmp = SQRT(den**2-s2)
    !
    rho =D*tmp/den
    eps =tmp-p*den/tmp-D
    eps=eps/D
    !
    continue
    !
END SUBROUTINE FUNC_eps_rho_C2P_GEOS




RECURSIVE SUBROUTINE FUNC_C2P_GEOS(F1,dprho,dpeps,p,rho,eps)
    USE EOS_par
    !USE EOS_mod
    IMPLICIT NONE
    ! input variables
    REAL, INTENT(IN) :: p, rho, eps
    !REAL, INTENT(INout) :: p
    REAL, INTENT(out) :: F1,dprho,dpeps
    REAL :: p1
    !
    ! define the functional and its partial derivatives.
    !
    CALL EOS_p_dprho_dpeps(p1,dprho,dpeps,rho,eps) 
    !
    SELECT CASE(EOSfunc)
    CASE(0) ! p-p[rho,eps] 
    	F1 = p - p1
    CASE(1)
        F1 = p**2 - p1
    CASE DEFAULT     ! p-p[rho,eps] 
        F1 = p - p1
    END SELECT
    !
    !
END SUBROUTINE FUNC_C2P_GEOS
 



RECURSIVE SUBROUTINE EOS_p_dprho_dpeps(p,dprho,dpeps,rho,eps) 
    USE Parameters, ONLY: EQN 
    USE EOS_par
    USE GREOS
    !USE EOS_mod
    IMPLICIT NONE
    ! input variables
    REAL, INTENT(IN) :: rho, eps
    REAL, INTENT(out) :: p,dprho,dpeps
    ! local variables::
    INTEGER :: ii,jj,kk !,N,nDOF
    REAL :: T
    !REAL :: psi_rho(EOS%N+1),psi_rho_xi
    !REAL :: psi_eps(EOS%N+1),psi_eps_xi
    !
    SELECT CASE(EOStype)
    CASE(-1)
        ! polytopic GASES:
        p = NSTOV_kappa*rho**EQN%gamma  
        dprho = EQN%gamma*NSTOV_kappa*rho**(EQN%gamma-1.0)
        dpeps = 0.0
        !
    CASE(0)
        ! IDEAL GASES:
        p = (EQN%gamma -1.0)*eps*rho
        dprho = (EQN%gamma -1.0)*eps
        dpeps = (EQN%gamma -1.0)*rho
        !
    CASE(1)
        ! NON-IDEAL GASES:  p^2 = (gamma -1)*eps*rho^1.5
        p = SQRT((EQN%gamma -1.0)*eps*rho*SQRT(rho))
        dprho = SQRT((EQN%gamma -1.0)*eps)*0.75/rho**0.25
        dpeps = SQRT((EQN%gamma -1.0)*rho*SQRT(rho))*0.5/SQRT(eps)
    CASE(2)
#ifdef GEOSDEBUG
        p = (EQN%gamma -1.0)*eps*rho
        dprho = (EQN%gamma -1.0)*eps
        dpeps = (EQN%gamma -1.0)*rho
        RETURN
#endif
        ! tabulated-DG
        !PRINT *,'***ERROR:  DG-mode EOS is available only in vector-mode!'
        !STOP
        ! tabulated DG (only vectorized)
        CALL EvaluateEOS_scalar(p,dprho,dpeps,T,rho,eps)
        !
        !
    CASE DEFAULT
        ! IDEAL GASES:
        p = (EQN%gamma -1.0)*eps*rho
        dprho = (EQN%gamma -1.0)*eps
        dpeps = (EQN%gamma -1.0)*rho
        !
    END SELECT
    !
    SELECT CASE(EOSfunc)
    CASE(0) ! p-p[rho,eps] 
        continue
        !
    CASE(1) ! p^2-p^2[rho,eps]  
        dprho=2*p*dprho
        dpeps=2*p*dpeps
        p=p**2
        !
    CASE DEFAULT ! p-p[rho,eps]  
        continue
        !
    END SELECT
    
    
    !
    ! TABULATED :
    !
    !CALL DGEOS_p_dprho_dpeps(p,dprho,dpeps,rho,eps)  
    !
    continue
    !
    !!N=EOS%N
    !!nDOF=EOS%nDOF
    !! locate point in the table
    !ii = MIN( IMAX, MAX( 1, CEILING( (rho-EOS%rhoL)/EOS%drho ) ) ) 
    !jj = MIN( IMAX, MAX( 1, CEILING( (eps-EOS%epsL)/EOS%deps ) ) )
    !!
    !iElem = EOS%idxe(ii,jj)
    !! interpolate p = p(rho,eps) from the EOS-DG
    !!
    !xi_rho = (rho - EOS%rho(ii))/EOS%drho
    !xi_eps = (rho - EOS%eps(jj))/EOS%deps
    !!
    !CALL BaseFunc1D_Neos(psi_rho,psi_rho_xi,xi_rho)  
    !CALL BaseFunc1D_Neos(psi_eps,psi_eps_xi,xi_eps)   
    !!
    !! evaluate p and its partial derivative with DG direct interpolation/evaluation
    !p=0
    !dprho=0
    !dpeps=0
    !DO jj=1,EOS%N
    !    DO ii=1,EOS%N
    !        i=EOS%rse(ii,jj)
    !        p = p + EOS%p(i,iElem)*psi_rho(ii)*psi_eps(jj)
    !        dprho = dprho + EOS%p(i,iElem)*psi_rho_xi(ii)*psi_eps(jj)
    !        dpeps = dpeps + EOS%p(i,iElem)*psi_rho(ii)*psi_eps_xi(jj)
    !    ENDDO
    !ENDDO
    !
    !continue
    !
END SUBROUTINE EOS_p_dprho_dpeps
 


RECURSIVE SUBROUTINE EOS_eps_P2C(rhoeps,rho,p)
	USE Parameters, ONLY: EQN,nVar
    USE GREOS
	IMPLICIT NONE
	! input variables:
	REAL, INTENT(IN) :: rho,p
	! output variables:
	REAL, INTENT(OUT) :: rhoeps
	! local var:
    REAL :: T,eps
    !
    SELECT CASE(EOStype)
    CASE(-1)
        ! polytopic GASES:
        ! rho*eps = p/(EQN%gamma-1)
        !p = NSTOV_kappa*rho**EQN%gamma  
	    rhoeps = NSTOV_kappa*rho**EQN%gamma/(EQN%gamma-1)     !p
        !dprho = EQN%gamma*NSTOV_kappa*rho**(EQN%gamma-1.0)
        !dpeps = 0.0
        !
    CASE(0)
        ! IDEAL GASES: p = (EQN%gamma-1)*rho*eps 
	    rhoeps = p/(EQN%gamma-1)
	    !
    CASE(1)
        ! NON-IDEAL GASES:  p^2 = (gamma -1)*eps*rho^1.5
        !                   eps*rho^1.5 =  p^2/(gamma -1)
        !                   eps*rho =  p^2/(gamma -1)/SQRT(rho) 
        rhoeps =  p**2/(EQN%gamma -1)/SQRT(rho) 
    CASE(2)
#ifdef GEOSDEBUG
	    rhoeps = p/(EQN%gamma-1)
        RETURN
#endif
        ! EOS-DG tabulated:
        ! first: evaluate the temperature; second: evaluate the energy
        ! PRINT *,'EOS-DG implemented only in vector mode!'
        ! STOP
        !
        CALL DP2TE_scalar(T,eps,rho,p)
	    rhoeps = rho*eps
        !
    CASE DEFAULT
        ! IDEAL GASES:
	    rhoeps = p/(EQN%gamma-1)
        !
    END SELECT
        
	!
    !
    ! TABULATED :
    !
    !CALL DGEOS_eps(eps,rho,p)  
    !rhoeps = rho*eps
    !
    continue
    !
END SUBROUTINE EOS_eps_P2C
	




RECURSIVE SUBROUTINE EOS_eps(rhoeps,rho,p,dd,tau,s2)
	USE Parameters, ONLY: EQN,nVar
	IMPLICIT NONE
	! input variables:
	REAL, INTENT(IN) :: rho,p,dd,tau,s2
	! output variables:
	REAL, INTENT(OUT) :: rhoeps
	! auxiliary variables
    REAL :: tmp,eps
    SELECT CASE(EOStype)
    CASE(0)
        ! IDEAL GASES: p = (EQN%gamma-1)*rho*eps 
	    rhoeps = p/(EQN%gamma-1)
	    !
    CASE(1)
        ! NON-IDEAL GASES:  p^2 = (gamma -1)*eps*rho^1.5
        !                   eps*rho^1.5 =  p^2/(gamma -1)
        !                   eps*rho =  p^2/(gamma -1)/SQRT(rho) 
        rhoeps =  p**2/(EQN%gamma -1)/SQRT(rho) 
    CASE(2)
#ifdef GEOSDEBUG
        rhoeps =  p**2/(EQN%gamma -1)/SQRT(rho) 
        RETURN
#endif
        ! EOS-DG tabulated:
        ! pure hydrodynamics: use a direct equation for epsilon
        !PRINT *,'EOS-DG implemented only in vector mode!'
        !STOP 
        ! EOS-DG tabulated:
        ! first: evaluate the temperature; second: evaluate the energy  
        tmp = SQRT((tau+p+dd)**2-s2)
        eps = tmp-p*(tau+p+dd)/tmp-dd
	    rhoeps = rho*eps/dd
        !
    CASE DEFAULT
        ! IDEAL GASES:
	    rhoeps = p/(EQN%gamma-1)
        !
    END SELECT
        
	!
    !
    ! TABULATED :
    !
    !CALL DGEOS_eps(eps,rho,p)  
    !rhoeps = rho*eps
    !
    continue
    !
END SUBROUTINE EOS_eps
	






RECURSIVE SUBROUTINE dFUNC_C2P_GEOS_vector(dF1,dp1rho,dp1eps,p,rho,eps,tau,D,s2)
    USE EOS_par
    !USE EOS_mod
    IMPLICIT NONE
    ! input variables
    REAL(8), INTENT(IN) :: dp1rho(VECTORLENGTH),dp1eps(VECTORLENGTH),p(VECTORLENGTH),rho(VECTORLENGTH),eps(VECTORLENGTH),tau(VECTORLENGTH),D(VECTORLENGTH),s2(VECTORLENGTH)
    REAL(8), INTENT(out) :: dF1(VECTORLENGTH)
    REAL(8) :: drhop(VECTORLENGTH),depsp(VECTORLENGTH),sqrtX(VECTORLENGTH),den(VECTORLENGTH),den2(VECTORLENGTH)
    !   
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64:: drhop,depsp,sqrtX,den,den2 
    !DIR$ ASSUME_ALIGNED dF1   : 64
    !DIR$ ASSUME_ALIGNED dp1rho: 64
    !DIR$ ASSUME_ALIGNED dp1eps: 64
    !DIR$ ASSUME_ALIGNED p     : 64
    !DIR$ ASSUME_ALIGNED rho   : 64
    !DIR$ ASSUME_ALIGNED eps   : 64
    !DIR$ ASSUME_ALIGNED tau   : 64
    !DIR$ ASSUME_ALIGNED D     : 64
    !DIR$ ASSUME_ALIGNED s2    : 64 
#else
    !DIR$ ATTRIBUTES ALIGN: 32:: drhop,depsp,sqrtX,den,den2 
    !DIR$ ASSUME_ALIGNED dF1   : 32
    !DIR$ ASSUME_ALIGNED dp1rho: 32
    !DIR$ ASSUME_ALIGNED dp1eps: 32
    !DIR$ ASSUME_ALIGNED p     : 32
    !DIR$ ASSUME_ALIGNED rho   : 32
    !DIR$ ASSUME_ALIGNED eps   : 32
    !DIR$ ASSUME_ALIGNED tau   : 32
    !DIR$ ASSUME_ALIGNED D     : 32
    !DIR$ ASSUME_ALIGNED s2    : 32 
#endif 
    !
    !CALL EOS_p_dprho_dpeps(p,rho,eps) 
    !
    den = tau+p+D
    den2 = den**2
    sqrtX = SQRT(den2-s2)
    !
    drhop = D*s2/den2/sqrtX
    depsp = p*s2/(den2-s2)/den/rho
    !
    SELECT CASE(EOSfunc)
    CASE(0) ! p-p[rho,eps] 
        dF1 = 1.0 - dp1rho*drhop - dp1eps*depsp
    CASE(1) ! p^2-F[rho,eps] 
        dF1 = 2.0*p - dp1rho*drhop - dp1eps*depsp
     CASE DEFAULT
    	dF1 = 1.0 - dp1rho*drhop - dp1eps*depsp
    END SELECt
    !
END SUBROUTINE dFUNC_C2P_GEOS_vector


RECURSIVE SUBROUTINE C2P_GEOS_vector(p,rho,eps,p0,tau,D,s2,FAILED)
    !
    ! this subroutine  (tau,D,s2) => (p,rho,eps)
    ! given an equation of state (analytical or tabulated)
    ! by means of a Newton method, looking for the root of the functional F1=F1[p,rho(p),eps(p)]:= p - EOS(p,rho(p),eps(p))
    !
    IMPLICIT NONE
    ! input variables
    REAL(8), INTENT(IN) :: tau(VECTORLENGTH),D(VECTORLENGTH),s2(VECTORLENGTH)
    REAL(8), INTENT(INOUT) :: p0(VECTORLENGTH) 
    ! output variables
    REAL(8), INTENT(OUT)  :: p(VECTORLENGTH),rho(VECTORLENGTH),eps(VECTORLENGTH)
    ! local variables
    REAL(8) :: F1(VECTORLENGTH),dF1(VECTORLENGTH),dp(VECTORLENGTH),dp1rho(VECTORLENGTH),dp1eps(VECTORLENGTH),eps0(VECTORLENGTH),rho0(VECTORLENGTH)
    REAL(8) :: alpha(VECTORLENGTH),res(VECTORLENGTH),Newres(VECTORLENGTH),ptmp(VECTORLENGTH)
    LOGICAL(8) :: FAILED(VECTORLENGTH),DONE(VECTORLENGTH),innerDONE(VECTORLENGTH)
    INTEGER :: Nt,inNt
    INTEGER, PARAMETER :: maxNt = 200 
    REAL :: tol,tolthreshold 
    !
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64:: F1,dF1,dp,dp1rho,dp1eps,eps0,rho0,alpha,res,Newres,DONE,ptmp,innerDONE
    !DIR$ ASSUME_ALIGNED p      : 64
    !DIR$ ASSUME_ALIGNED rho    : 64
    !DIR$ ASSUME_ALIGNED eps    : 64
    !DIR$ ASSUME_ALIGNED p0     : 64
    !DIR$ ASSUME_ALIGNED tau    : 64
    !DIR$ ASSUME_ALIGNED D      : 64
    !DIR$ ASSUME_ALIGNED s2     : 64
    !DIR$ ASSUME_ALIGNED FAILED : 64
#else                         
    !DIR$ ATTRIBUTES ALIGN: 32:: F1,dF1,dp,dp1rho,dp1eps,eps0,rho0,alpha,res,Newres,DONE,ptmp,innerDONE
    !DIR$ ASSUME_ALIGNED p      : 32
    !DIR$ ASSUME_ALIGNED rho    : 32
    !DIR$ ASSUME_ALIGNED eps    : 32
    !DIR$ ASSUME_ALIGNED p0     : 32
    !DIR$ ASSUME_ALIGNED tau    : 32
    !DIR$ ASSUME_ALIGNED D      : 32
    !DIR$ ASSUME_ALIGNED s2     : 32
    !DIR$ ASSUME_ALIGNED FAILED : 32
#endif 
    FAILED = .FALSE.
    DONE = .FALSE.
    ! here is the initial guess:
    p=p0                                                         
    CALL FUNC_eps_rho_C2P_GEOS_vector(eps0,rho0,p,tau,D,s2)      
    !
    rho=rho0
    eps=eps0
    tol = 1e-14
    tolthreshold = 1e-10
    !
    ! Newton algorithm:
    DO Nt = 1, maxNt
  !      CALL FUNC_C2P_GEOS_vector(F1,dp1rho,dp1eps,p,rho,eps)
  !      IF(SQRT(SUM(F1**2)).LT.tol*VECTORLENGTH) THEN
  !          continue
  !          RETURN
  !      ENDIF
  !      !
  !      CALL dFUNC_C2P_GEOS_vector(dF1,dp1rho,dp1eps,p,rho,eps,tau,D,s2)
  !      !
  !      dp=-F1/dF1
  !      WHERE(ABS(F1).GT.tol)  
  !          p=p+dp
  !      ENDWHERE
  !      !
  !      CALL FUNC_eps_rho_C2P_GEOS_vector(eps,rho,p,tau,D,s2)  
  !      !
		!CYCLE
		!
        CALL FUNC_C2P_GEOS_vector(F1,dp1rho,dp1eps,p,rho,eps)
        res = ABS(F1)
        !WHERE(res.LT.tol*VECTORLENGTH)
        !    DONE =.TRUE.
        !END WHERE
        !
        CALL dFUNC_C2P_GEOS_vector(dF1,dp1rho,dp1eps,p,rho,eps,tau,D,s2)
        !WHERE(ISNAN(dF1)) 
        !    FAILED = .TRUE.
        !    DONE =.TRUE.
        !ENDWHERE
        !
        IF(ALL(res.LT.tol*VECTORLENGTH)) THEN !all done!
            continue
            RETURN
        ENDIF
        !
        dp=-F1/dF1
        !WHERE(ABS(F1).GT.tol)  
        !    p=p+dp
        !ENDWHERE
        !
        alpha=1.0
        innerDONE = DONE
        DO inNt=1,20
            ptmp =p+alpha*dp
            IF(ANY(ptmp<0.0)) THEN
                alpha = 0.5*alpha
                CYCLE
            ENDIF
            ptmp = ABS(ptmp)
            CALL FUNC_C2P_GEOS_vector(F1,dp1rho,dp1eps,ptmp,rho,eps)
            newres = ABS(F1) 
            IF(ANY(newres.GE.res)) THEN
                alpha = 0.5*alpha
            ELSE 
                EXIT 
            ENDIF
            !         
        ENDDO      
        !WHERE(.NOT.DONE)
        !    p=ABS(p+alpha*dp)
        !ENDWHERE
        !
        CALL FUNC_eps_rho_C2P_GEOS_vector(eps,rho,p,tau,D,s2)  
        !
    ENDDO
    ! 
    WHERE(ABS(F1).GT.tolthreshold.OR.ISNAN(F1)) 
    	FAILED = .TRUE.
    ENDWHERE
    !
    continue
    if(any(failed)) then
        CONTINUE
    endif
    !
    !
END SUBROUTINE C2P_GEOS_vector




RECURSIVE SUBROUTINE FUNC_eps_rho_C2P_GEOS_vector(eps,rho,p,tau,D,s2) 
    ! this function evaluate rho =rho(p), eps = eps(p)
    USE EOS_par
    !USE EOS_mod
    IMPLICIT NONE
    ! input variables
    REAL(8), INTENT(IN) :: p(VECTORLENGTH),tau(VECTORLENGTH),D(VECTORLENGTH),s2(VECTORLENGTH)
    REAL(8), INTENT(out) :: eps(VECTORLENGTH),rho(VECTORLENGTH)
    real(8) :: den(VECTORLENGTH),tmp(VECTORLENGTH)
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64:: den,tmp 
    !DIR$ ASSUME_ALIGNED eps     : 64 
    !DIR$ ASSUME_ALIGNED rho     : 64 
    !DIR$ ASSUME_ALIGNED p       : 64 
    !DIR$ ASSUME_ALIGNED tau     : 64 
    !DIR$ ASSUME_ALIGNED D       : 64 
    !DIR$ ASSUME_ALIGNED s2      : 64 
#else                         
    !DIR$ ATTRIBUTES ALIGN: 32:: den,tmp
    !DIR$ ASSUME_ALIGNED eps     : 32 
    !DIR$ ASSUME_ALIGNED rho     : 32 
    !DIR$ ASSUME_ALIGNED p       : 32 
    !DIR$ ASSUME_ALIGNED tau     : 32 
    !DIR$ ASSUME_ALIGNED D       : 32 
    !DIR$ ASSUME_ALIGNED s2      : 32 
#endif 
    ! 
    ! this is true only for pure GR-Hydro, independently on the EOS.
    den = tau+p+D  
    tmp = SQRT(den**2-s2) 
    !
    rho =D*tmp/den
    eps =tmp-p*den/tmp-D
    eps=eps/D
    !
    continue
    !
END SUBROUTINE FUNC_eps_rho_C2P_GEOS_vector




RECURSIVE SUBROUTINE FUNC_C2P_GEOS_vector(F1,dprho,dpeps,p,rho,eps)
    USE EOS_par
    !USE EOS_mod
    IMPLICIT NONE
    ! input variables
    REAL(8), INTENT(IN) :: rho(VECTORLENGTH), eps(VECTORLENGTH)
    REAL(8), INTENT(INout) :: p(VECTORLENGTH)
    REAL(8), INTENT(out) :: F1(VECTORLENGTH),dprho(VECTORLENGTH),dpeps(VECTORLENGTH)
    REAL(8) :: p1(VECTORLENGTH)
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64:: p1
    !DIR$ ASSUME_ALIGNED F1    : 64 
    !DIR$ ASSUME_ALIGNED dprho : 64 
    !DIR$ ASSUME_ALIGNED dpeps : 64 
    !DIR$ ASSUME_ALIGNED p     : 64 
    !DIR$ ASSUME_ALIGNED rho   : 64 
    !DIR$ ASSUME_ALIGNED eps   : 64 
#else                         
    !DIR$ ATTRIBUTES ALIGN: 32:: p1
    !DIR$ ASSUME_ALIGNED F1    : 32 
    !DIR$ ASSUME_ALIGNED dprho : 32 
    !DIR$ ASSUME_ALIGNED dpeps : 32 
    !DIR$ ASSUME_ALIGNED p     : 32 
    !DIR$ ASSUME_ALIGNED rho   : 32 
    !DIR$ ASSUME_ALIGNED eps   : 32 
#endif 
    !
    ! define the functional and its partial derivatives.
    !
    CALL EOS_p_dprho_dpeps_vector(p1,dprho,dpeps,rho,eps) 
    !
    SELECT CASE(EOSfunc)
    CASE(0) ! p-p[rho,eps] 
    F1 = p - p1
    CASE(1)
        F1 = p**2 - p1
    CASE DEFAULT     ! p-p[rho,eps] 
        F1 = p - p1
    END SELECT
    !
    !
END SUBROUTINE FUNC_C2P_GEOS_vector
 



RECURSIVE SUBROUTINE EOS_p_dprho_dpeps_vector(p,dprho,dpeps,rho,eps) 
    ! ****
    ! given rho and eps, this subroutine evaluates p=p(rho,eps); and its partial derivatives
    ! ****
    USE Parameters, ONLY: EQN
    USE GREOS
    USE EOS_par
    !USE EOS_mod
    IMPLICIT NONE
    ! input variables
    REAL(8), INTENT(IN) :: rho(VECTORLENGTH), eps(VECTORLENGTH)
    REAL(8), INTENT(out) :: p(VECTORLENGTH),dprho(VECTORLENGTH),dpeps(VECTORLENGTH)
    ! local variables::
    INTEGER :: ii,jj,kk,iv
    REAL(8), DIMENSION(VECTORLENGTH) :: T,p1,dprho1,dpeps1
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64:: T,p1,dprho1,dpeps1
    !DIR$ ASSUME_ALIGNED p     : 64
    !DIR$ ASSUME_ALIGNED dprho : 64
    !DIR$ ASSUME_ALIGNED dpeps : 64
    !DIR$ ASSUME_ALIGNED rho   : 64
    !DIR$ ASSUME_ALIGNED eps   : 64 
#else                         
    !DIR$ ATTRIBUTES ALIGN: 32:: T,p1,dprho1,dpeps1
    !DIR$ ASSUME_ALIGNED p     : 32
    !DIR$ ASSUME_ALIGNED dprho : 32
    !DIR$ ASSUME_ALIGNED dpeps : 32
    !DIR$ ASSUME_ALIGNED rho   : 32
    !DIR$ ASSUME_ALIGNED eps   : 32 
#endif 
    !
    SELECT CASE(EOStype)
    CASE(0)
        ! IDEAL GASES:
        p = (EQN%gamma -1.0)*eps*rho
        dprho = (EQN%gamma -1.0)*eps
        dpeps = (EQN%gamma -1.0)*rho
        !
    CASE(1)
        ! NON-IDEAL GASES:  p^2 = (gamma -1)*eps*rho^1.5
        p = SQRT((EQN%gamma -1.0)*eps*rho*SQRT(rho))
        dprho = SQRT((EQN%gamma -1.0)*eps)*0.75/rho**0.25
        dpeps = SQRT((EQN%gamma -1.0)*rho*SQRT(rho))*0.5/SQRT(eps)
    CASE(2)
#ifdef GEOSDEBUG
        p = (EQN%gamma -1.0)*eps*rho
        dprho = (EQN%gamma -1.0)*eps
        dpeps = (EQN%gamma -1.0)*rho
        RETURN
#endif
        ! tabulated DG (only vectorized) 
        !p1 = (EQN%gamma -1.0)*eps*rho+1e-14
        !dprho1 = (EQN%gamma -1.0)*eps+1e-13
        !dpeps1 = (EQN%gamma -1.0)*rho+1e-13
        !p=p1
        !dprho=dprho1
        !dpeps=dpeps1
        !p = (EQN%gamma -1.0)*eps*rho
        !dprho = (EQN%gamma -1.0)*eps
        !dpeps = (EQN%gamma -1.0)*rho
        CALL EvaluateEOS(p,dprho,dpeps,T,rho,eps) 
        !    IF(MAXVAL(p1-p).gt.1e-13) THEN
        !        continue
        !    ENDIF
        !    IF(MAXVAL(dprho1-dprho).gt.1e-13) THEN
        !        continue
        !    ENDIF
        !    IF(MAXVAL(dpeps1-dpeps).gt.1e-13) THEN
        !        continue
        !    ENDIF
        !    
        !
    CASE DEFAULT
        ! IDEAL GASES:
        p = (EQN%gamma -1.0)*eps*rho
        dprho = (EQN%gamma -1.0)*eps
        dpeps = (EQN%gamma -1.0)*rho
        !
    END SELECT
    !
    SELECT CASE(EOSfunc)
    CASE(0) ! p-p[rho,eps] 
        continue
        !
    CASE(1) ! p^2-p^2[rho,eps]  
        dprho=2*p*dprho
        dpeps=2*p*dpeps
        p=p**2
        !
    CASE DEFAULT ! p-p[rho,eps]  
        continue
        !
    END SELECT
    !
    ! TABULATED :
    !
    !CALL DGEOS_p_dprho_dpeps(p,dprho,dpeps,rho,eps)  
    !
    continue
    !
    !!N=EOS%N
    !!nDOF=EOS%nDOF
    !! locate point in the table
    !ii = MIN( IMAX, MAX( 1, CEILING( (rho-EOS%rhoL)/EOS%drho ) ) ) 
    !jj = MIN( IMAX, MAX( 1, CEILING( (eps-EOS%epsL)/EOS%deps ) ) )
    !!
    !iElem = EOS%idxe(ii,jj)
    !! interpolate p = p(rho,eps) from the EOS-DG
    !!
    !xi_rho = (rho - EOS%rho(ii))/EOS%drho
    !xi_eps = (rho - EOS%eps(jj))/EOS%deps
    !!
    !CALL BaseFunc1D_Neos(psi_rho,psi_rho_xi,xi_rho)  
    !CALL BaseFunc1D_Neos(psi_eps,psi_eps_xi,xi_eps)   
    !!
    !! evaluate p and its partial derivative with DG direct interpolation/evaluation
    !p=0
    !dprho=0
    !dpeps=0
    !DO jj=1,EOS%N
    !    DO ii=1,EOS%N
    !        i=EOS%rse(ii,jj)
    !        p = p + EOS%p(i,iElem)*psi_rho(ii)*psi_eps(jj)
    !        dprho = dprho + EOS%p(i,iElem)*psi_rho_xi(ii)*psi_eps(jj)
    !        dpeps = dpeps + EOS%p(i,iElem)*psi_rho(ii)*psi_eps_xi(jj)
    !    ENDDO
    !ENDDO
    !
    !continue
    !
END SUBROUTINE EOS_p_dprho_dpeps_vector



RECURSIVE SUBROUTINE EOS_eps_depsdrho_depsdp(eps,depsdrho,depsdp,rho,p) 
    ! ****
    ! given rho and eps, this subroutine evaluates p=p(rho,eps); and its partial derivatives
    ! ****
    USE Parameters, ONLY: EQN
    USE EOS_par
    !USE EOS_mod
    IMPLICIT NONE
    ! input variables
    REAL, INTENT(IN) :: p,rho
    REAL, INTENT(out) ::  eps,depsdp,depsdrho
    ! local variables::
    INTEGER :: ii,jj,kk !,N,nDOF
    !REAL :: xi_rho,xi_eps
    !REAL :: psi_rho(EOS%N+1),psi_rho_xi
    !REAL :: psi_eps(EOS%N+1),psi_eps_xi
    ! 
    !
    ! IDEAL GASES:
    eps         =p/(EQN%gamma -1.0)/rho 
    depsdrho    =-p/(EQN%gamma -1.0)/rho**2 
    depsdp      =1.0/(EQN%gamma -1.0)/rho 
    !
    ! TABULATED :
    !
    !CALL DGEOS_p_dprho_dpeps(p,dprho,dpeps,rho,eps)  
    !
    continue
    ! 
    !continue
    !
END SUBROUTINE EOS_eps_depsdrho_depsdp




RECURSIVE SUBROUTINE EOS_eps_depsdrho_depsdp_vector(eps,depsdrho,depsdp,rho,p) 
    ! ****
    ! given rho and eps, this subroutine evaluates p=p(rho,eps); and its partial derivatives
    ! ****
    USE Parameters, ONLY: EQN
    USE EOS_par
    !USE EOS_mod
    IMPLICIT NONE
    ! input variables
    REAL(8), INTENT(IN) :: p(VECTORLENGTH),rho(VECTORLENGTH)
    REAL(8), INTENT(out) ::  eps(VECTORLENGTH),depsdp(VECTORLENGTH),depsdrho(VECTORLENGTH)
    ! local variables::
    INTEGER :: ii,jj,kk !,N,nDOF
    !REAL :: xi_rho,xi_eps
    !REAL :: psi_rho(EOS%N+1),psi_rho_xi
    !REAL :: psi_eps(EOS%N+1),psi_eps_xi
#ifdef AVX512
    !DIR$ ASSUME_ALIGNED p          : 64
    !DIR$ ASSUME_ALIGNED depsdp     : 64
    !DIR$ ASSUME_ALIGNED depsdrho   : 64
    !DIR$ ASSUME_ALIGNED rho        : 64
    !DIR$ ASSUME_ALIGNED eps        : 64 
#else                         
    !DIR$ ASSUME_ALIGNED p          : 32
    !DIR$ ASSUME_ALIGNED depsdp     : 32
    !DIR$ ASSUME_ALIGNED depsdrho   : 32
    !DIR$ ASSUME_ALIGNED rho        : 32
    !DIR$ ASSUME_ALIGNED eps        : 32
#endif 
    ! 
    !
    ! IDEAL GASES:
    eps         =p/(EQN%gamma -1.0)/rho 
    depsdrho    =-p/(EQN%gamma -1.0)/rho**2 
    depsdp      =1.0/(EQN%gamma -1.0)/rho 
    !
    ! TABULATED :
    !
    !CALL DGEOS_p_dprho_dpeps(p,dprho,dpeps,rho,eps)  
    !
    continue
    ! 
    !continue
    !
END SUBROUTINE EOS_eps_depsdrho_depsdp_vector





RECURSIVE SUBROUTINE EOS_eps_vector(rhoeps,rho,p,dd,tau,s2)
    ! ****
    ! given rho and p, this subroutine evaluates rho*eps:   where eps = eps(rho,p)
    ! ****
	USE Parameters, ONLY: EQN,nVar
    USE GREOS
	IMPLICIT NONE
	! input variables:
	REAL(8), INTENT(IN) :: rho(VECTORLENGTH),p(VECTORLENGTH),dd(VECTORLENGTH),tau(VECTORLENGTH),s2(VECTORLENGTH)
	! output variables:
	REAL(8), INTENT(OUT) :: rhoeps(VECTORLENGTH) 
    ! local variables
    REAL(8)                :: tmp(VECTORLENGTH),eps(VECTORLENGTH)
	!
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64:: tmp
    !DIR$ ASSUME_ALIGNED rho    : 64 
    !DIR$ ASSUME_ALIGNED p      : 64 
    !DIR$ ASSUME_ALIGNED rhoeps : 64  
    !DIR$ ASSUME_ALIGNED dd     : 64  
    !DIR$ ASSUME_ALIGNED tau    : 64  
    !DIR$ ASSUME_ALIGNED s2     : 64  
#else
    !DIR$ ATTRIBUTES ALIGN: 32:: tmp
    !DIR$ ASSUME_ALIGNED rho    : 32 
    !DIR$ ASSUME_ALIGNED p      : 32 
    !DIR$ ASSUME_ALIGNED rhoeps : 32  
    !DIR$ ASSUME_ALIGNED dd     : 32  
    !DIR$ ASSUME_ALIGNED tau    : 32  
    !DIR$ ASSUME_ALIGNED s2     : 32  
#endif
	!
    SELECT CASE(EOStype)
    CASE(0)
        ! IDEAL GASES: p = (EQN%gamma-1)*rho*eps 
	    rhoeps(:) = p(:)/(EQN%gamma-1.0)
	    !
    CASE(1)
        ! NON-IDEAL GASES:  p^2 = (gamma -1)*eps*rho^1.5
        !                   eps*rho^1.5 =  p^2/(gamma -1)
        !                   eps*rho =  p^2/(gamma -1)/SQRT(rho) 
        rhoeps(:) = p(:)**2/(EQN%gamma -1.0)/SQRT(rho) 
    CASE(2)
#ifdef GEOSDEBUG
        rhoeps =  p**2/(EQN%gamma -1.0)/SQRT(rho) 
        RETURN
#endif
        ! EOS-DG tabulated:
        ! first: evaluate the temperature; second: evaluate the energy  
        tmp(:) = SQRT((tau(:)+p(:)+dd(:))**2-s2(:))
        eps(:) = tmp(:)-p(:)*(tau(:)+p(:)+dd(:))/tmp(:)-dd(:)
	    rhoeps(:) = rho(:)*eps(:)/dd(:)
        !
    CASE DEFAULT
        ! IDEAL GASES:
	    rhoeps(:) = p(:)/(EQN%gamma-1.0)
        !
    END SELECT
    !
    ! TABULATED :
    !
    !CALL DGEOS_eps(eps,rho,p)  
    !rhoeps = rho*eps
    !
    continue
    !
END SUBROUTINE EOS_eps_vector
	
   
   
   

RECURSIVE SUBROUTINE EOS_eps_P2C_vector(rhoeps,rho,p)
	USE Parameters, ONLY: EQN,nVar
    USE GREOS
	IMPLICIT NONE
	! input variables:
	REAL(8), INTENT(IN) :: rho(VECTORLENGTH),p(VECTORLENGTH)
	! output variables:
	REAL(8), INTENT(OUT) :: rhoeps(VECTORLENGTH)
	! local var:
    REAL(8) :: T(VECTORLENGTH),eps(VECTORLENGTH)
    !
	!
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64:: T,eps
    !DIR$ ASSUME_ALIGNED rho    : 64 
    !DIR$ ASSUME_ALIGNED p      : 64 
    !DIR$ ASSUME_ALIGNED rhoeps : 64  
#else
    !DIR$ ATTRIBUTES ALIGN: 32:: T,eps
    !DIR$ ASSUME_ALIGNED rho    : 32 
    !DIR$ ASSUME_ALIGNED p      : 32 
    !DIR$ ASSUME_ALIGNED rhoeps : 32  
#endif
	!
    SELECT CASE(EOStype)
    CASE(0)
        ! IDEAL GASES: p = (EQN%gamma-1)*rho*eps 
	    rhoeps = p/(EQN%gamma-1.0)
	    !
    CASE(1)
        ! NON-IDEAL GASES:  p^2 = (gamma -1)*eps*rho^1.5
        !                   eps*rho^1.5 =  p^2/(gamma -1)
        !                   eps*rho =  p^2/(gamma -1)/SQRT(rho) 
        rhoeps =  p**2/(EQN%gamma -1.0)/SQRT(rho) 
    CASE(2)
#ifdef GEOSDEBUG
	    rhoeps = p/(EQN%gamma-1.0)
        RETURN
#endif
        ! EOS-DG tabulated:
        ! first: evaluate the temperature; second: evaluate the energy  
        CALL DP2TE(T,eps,rho,p)
	    rhoeps = rho*eps
        !
    CASE DEFAULT
        ! IDEAL GASES:
	    rhoeps = p/(EQN%gamma-1.0)
        !
    END SELECT
        
	!
    !
    ! TABULATED :
    !
    !CALL DGEOS_eps(eps,rho,p)  
    !rhoeps = rho*eps
    !
    continue
    !
END SUBROUTINE EOS_eps_P2C_vector
	
   
   
   
   




!    
!SUBROUTINE BaseFunc1D_Neos(phi,phi_xi,xi)    !compute the vector of the basis at a given point x=xi
!   !USE Parameters, ONLY : tin
!   IMPLICIT NONE
!   ! Argument list 
!   REAL        :: phi(EOS%N+1), phi_xi(EOS%N+1), xi  
!   INTEGER     :: N
!   ! Local variables 
!   INTEGER     :: i,j,m
!   REAL        :: tmp   
!   INTENT(IN)  :: xi
!   INTENT(OUT) :: phi, phi_xi 
!   ! 
!   N=EOS%N
!   ! Initialize variables 
!   phi      = 1. 
!   phi_xi   = 0. 
!   ! Lagrange polynomial and its derivative 
!   DO m = 1,N+1
!      DO j = 1,N+1
!         IF(j.EQ.m) CYCLE 
!         phi(m) = phi(m)*(xi-tin(j))/(tin(m)-tin(j))    ! phi_m(x) is zero in every x=tin(j), but tin(m)
!      ENDDO                                             ! note: the integration is done along more points xi, nGP!=N+1=nDOF
!      DO i = 1,N+1
!         IF(i.EQ.m) CYCLE
!         tmp = 1. 
!         DO j = 1,N+1
!            IF(j.EQ.i) CYCLE 
!            IF(j.EQ.m) CYCLE 
!            tmp = tmp*(xi-tin(j))/(tin(m)-tin(j))    
!         ENDDO 
!         phi_xi(m) = phi_xi(m) + tmp/(tin(m)-tin(i)) 
!      ENDDO 
!   ENDDO 
!   !
!END SUBROUTINE BaseFunc1D_Neos
             



RECURSIVE SUBROUTINE FUNC_C2P_RHD1(p,f1,df1,d,tau,s2,FAILED)
  !
  ! This is the CONS2PRIM strategy adopted in Whisky and originally described in Baiotti's PhD thesis. 
  ! See also Appendix D of the book "Relativistic hydrodynamics" by Rezzolla & Zanotti (2013)
  !
  IMPLICIT NONE
  INTEGER    :: iter
  REAL       :: p,f1,df1
  REAL       :: d,tau,s2
  
  LOGICAL    :: FAILED
  REAL       :: sqrtX,rho,epsilon,drhodp
  REAL       :: den,eps,den2,drhop,depsp,p1,dprho,dpeps,T
  
  INTENT(IN) :: p,d,tau,s2
  INTENT(OUT):: f1,df1,FAILED
  
  ! x = pressure, the unknown
  FAILED   = .FALSE.
  ! 
  !CALL FUNC_eps_rho_C2P_GEOS(eps0,rho0,p,tau,D,s2)  
  !
  den = tau+p+D 
  sqrtX = den**2 - s2
  IF(sqrtX.GT.0.0) THEN
     sqrtX    = SQRT(sqrtX)
  ELSE
     FAILED = .TRUE.
     RETURN
  ENDIF 
  !
  rho =D*sqrtX/den
  eps =sqrtX-p*den/sqrtX-D
  eps=eps/D
  !
  den2 = den**2 
  !
  drhop = D*s2/den2/sqrtX
  depsp = p*s2/(den2-s2)/den/rho
  !
  ! define the functional and its partial derivatives.
  !
  !CALL EOS_p_dprho_dpeps(p1,dprho,dpeps,rho,eps)  
  !
  SELECT CASE(EOStype)
  CASE(-1)
      ! polytopic GASES:
      p1 = NSTOV_kappa*rho**EQN%gamma  
      dprho = EQN%gamma*NSTOV_kappa*rho**(EQN%gamma-1.0)
      dpeps = 0.0
      !
  CASE(0)
      ! IDEAL GASES:
      p1 = (EQN%gamma -1.0)*eps*rho
      dprho = (EQN%gamma -1.0)*eps
      dpeps = (EQN%gamma -1.0)*rho
      !
  CASE(1)
      ! NON-IDEAL GASES:  p^2 = (gamma -1)*eps*rho^1.5
      p1 = SQRT((EQN%gamma -1.0)*eps*rho*SQRT(rho))
      dprho = SQRT((EQN%gamma -1.0)*eps)*0.75/rho**0.25
      dpeps = SQRT((EQN%gamma -1.0)*rho*SQRT(rho))*0.5/SQRT(eps)
  CASE(2)
#ifdef GEOSDEBUG
      p1 = (EQN%gamma -1.0)*eps*rho
      dprho = (EQN%gamma -1.0)*eps
      dpeps = (EQN%gamma -1.0)*rho
      RETURN
#endif
      ! tabulated-DG
      !PRINT *,'***ERROR:  DG-mode EOS is available only in vector-mode!'
      !STOP
      ! tabulated DG (only vectorized)
      CALL EvaluateEOS_scalar(p1,dprho,dpeps,T,rho,eps)
      !
      !
  CASE DEFAULT
      ! IDEAL GASES:
      p1 = (EQN%gamma -1.0)*eps*rho
      dprho = (EQN%gamma -1.0)*eps
      dpeps = (EQN%gamma -1.0)*rho
      !
  END SELECT
  !
  SELECT CASE(EOSfunc)
  CASE(0) ! p-p[rho,eps] 
      continue
      !
  	  F1 = p - p1
      dF1 = 1.0 - dprho*drhop - dpeps*depsp
	  !
  CASE(1) ! p^2-p^2[rho,eps]  
      dprho=2*p1*dprho
      dpeps=2*p1*dpeps
      p1=p1**2
      !
      F1 = p**2 - p1
      dF1 = 2*p - dprho*drhop - dpeps*depsp
	  !
  CASE DEFAULT ! p-p[rho,eps]  
      continue
      !
      F1 = p - p1
      dF1 = 1.0 - dprho*drhop - dpeps*depsp
  END SELECT 
  !
  continue
  !
END SUBROUTINE FUNC_C2P_RHD1


            

RECURSIVE SUBROUTINE RTSAFE_C2P_RHD1(p,X1,X2,XACC,d,tau,s2,FAILED)
  USE Parameters, ONLY : nVar
  IMPLICIT NONE
  ! argument variables
  REAL, INTENT(INOUT)       :: p,X1,X2,XACC
  REAL, INTENT(IN)          :: d,tau,s2
  LOGICAL, INTENT(INOUT)    :: FAILED
  ! auxiliary variables
  INTEGER, PARAMETER   :: MAXIT=400
  INTEGER              :: J,i
  REAL                 :: FL,FH,DF,XH,XL,SWAP,DXOLD,XACC2,DX,F,TEMP,Xtmp,Vtmp(7),Ftmp,dFtmp,dXtmp
  ! 
  !IF(DebugPrint) THEN 
  !    WRITE(DebugFile,'(a)') 'DebugFile.dat' 
  !    OPEN(UNIT=4321,FILE=TRIM(DebugFile), RECL=800) 
  !    WRITE(4321,*) '  VARIABLES = "i"   "Xtmp"   "Ftmp"   "dFtmp"   "gam" "d" "e" "s2" "b2" "sb2" "w" '
  !    WRITE(4321,*)
  !    Vtmp(1:7) = (/ gam,d,e,s2,b2,sb2,w /)
  !    dXtmp = (1.0-1.0e-12)/300
  !    DO i =1, 350
  !        Xtmp = 0.5*dXtmp + real(i)*dXtmp
  !      CALL FUNC_C2P_RHD1(Xtmp,Ftmp,dFtmp,Vtmp(1),Vtmp(2),Vtmp(3),Vtmp(4),Vtmp(5),Vtmp(6),Vtmp(7))
  !      WRITE(4321,'(i9.4,10(E16.6))') i,Xtmp,Ftmp,dFtmp,Vtmp(1:7)
  !    ENDDO
  !    CLOSE(4321)   
  !ENDIF
  !
  !PRINT *,e
	!
  XACC2 = 1e-60 !XACC(:)**2
  FAILED = .FALSE. 
  CALL FUNC_C2P_RHD1(X1,FL,DF,d,tau,s2,FAILED)
  IF( ISNAN(FL) ) THEN
         FL = 1e20
  ENDIF
#ifdef C2PFF  
  IF(ABS(FL).LT.XACC2) THEN
      continue
#else
  IF(FL.EQ.0.) THEN
#endif
     p=X1
     RETURN
  ENDIF
  CALL FUNC_C2P_RHD1(X2,FH,DF,d,tau,s2,FAILED)
  IF( ISNAN(FH) ) THEN
         FH = 1e20
  ENDIF
#ifdef C2PFF  
  IF(ABS(FH).LT.XACC2) THEN
      continue
#else
  IF(FH.EQ.0.) THEN
#endif
     p=X2
     RETURN
  ENDIF
  IF(FL*FH.GT.0.) THEN
     FAILED = .TRUE.
     p = 0. ! assign value even if it fails
     RETURN
  ENDIF
  IF(FL.LT.0.)THEN
     XL=X1
     XH=X2
  ELSE
     XH=X1
     XL=X2
     SWAP=FL
     FL=FH
     FH=SWAP
  ENDIF
  p=.5*(X1+X2)
  DXOLD=ABS(X2-X1)
  DX=DXOLD
  CALL FUNC_C2P_RHD1(p,F,DF,d,tau,s2,FAILED)
  IF( ISNAN(F) ) THEN
      F = 1e20
  ENDIF 
  DO 11 J=1,MAXIT
     IF(((p-XH)*DF-F)*((p-XL)*DF-F).GE.0. &
          .OR. ABS(2.*F).GT.ABS(DXOLD*DF) ) THEN
        DXOLD=DX
        DX=0.5*(XH-XL)
        p=XL+DX
        IF(XL.EQ.p) THEN
            continue
            RETURN
        ENDIF
     ELSE
        DXOLD=DX
        DX=F/DF
        TEMP=p
        p=p-DX
        IF(TEMP.EQ.p) THEN
            continue
            RETURN
        ENDIF
     ENDIF
     IF(ABS(DX).LT.XACC) THEN
            !PRINT *,j,DX
            continue
            RETURN
     ENDIF
     CALL FUNC_C2P_RHD1(p,F,DF,d,tau,s2,FAILED)
     IF( ISNAN(F) ) THEN
            F = 1e20
     ENDIF
     IF(F.LT.0.) THEN
        XL=p
        FL=F
     ELSE
        XH=p
        FH=F
     ENDIF
#ifdef C2PFF  
     IF(ABS(F).LT.XACC2) THEN
         continue
#else
     IF(F.EQ.0.) THEN
#endif 
        RETURN
     ENDIF
11   CONTINUE
     FAILED = .TRUE.
     p = 0. ! assign value even if it fails
     RETURN
  END SUBROUTINE RTSAFE_C2P_RHD1
   !
  
  
  
  
  

RECURSIVE SUBROUTINE FUNC_C2P_RHD1_VECTORF(p,f1,df1,d,tau,s2,FAILED)
  !
  ! This is the CONS2PRIM strategy adopted in Whisky and originally described in Baiotti's PhD thesis. 
  ! See also Appendix D of the book "Relativistic hydrodynamics" by Rezzolla & Zanotti (2013)
  !
  IMPLICIT NONE 
  INTEGER    :: iter
  REAL(8), DIMENSION(VECTORLENGTH)  :: p,f1,df1
  REAL(8), DIMENSION(VECTORLENGTH)  :: d,tau,s2
  LOGICAL(8), DIMENSION(VECTORLENGTH)  :: FAILED
  REAL(8), DIMENSION(VECTORLENGTH)  :: sqrtX,rho,epsilon,drhodp
  REAL(8), DIMENSION(VECTORLENGTH)  :: den,eps,den2,drhop,depsp,p1,dprho,dpeps,T
  !
  INTENT(IN) :: p,d,tau,s2
  INTENT(OUT):: f1,df1,FAILED
  
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64:: sqrtX,rho,epsilon,drhodp,den,eps,den2,drhop,depsp,p1,dprho,dpeps,T
    !DIR$ ASSUME_ALIGNED p : 64 
    !DIR$ ASSUME_ALIGNED f1 : 64 
    !DIR$ ASSUME_ALIGNED df1 : 64  
    !DIR$ ASSUME_ALIGNED d    : 64 
    !DIR$ ASSUME_ALIGNED tau    : 64 
    !DIR$ ASSUME_ALIGNED s2   : 64 
    !DIR$ ASSUME_ALIGNED FAILED    : 64 
#else
    !DIR$ ATTRIBUTES ALIGN: 32 :: sqrtX,rho,epsilon,drhodp,den,eps,den2,drhop,depsp,p1,dprho,dpeps,T
    !DIR$ ASSUME_ALIGNED p   : 32 
    !DIR$ ASSUME_ALIGNED f1   : 32 
    !DIR$ ASSUME_ALIGNED df1   : 32  
    !DIR$ ASSUME_ALIGNED d    : 32 
    !DIR$ ASSUME_ALIGNED tau    : 32 
    !DIR$ ASSUME_ALIGNED s2   : 32 
    !DIR$ ASSUME_ALIGNED FAILED    : 32  
#endif 
  !
  ! x = pressure, the unknown
  FAILED   = .FALSE.
  ! 
  !CALL FUNC_eps_rho_C2P_GEOS(eps0,rho0,p,tau,D,s2)  
  !
  den = tau+p+D 
  sqrtX = den**2 - s2
  IF(ANY(sqrtX.LT.0.0)) THEN
     FAILED = .TRUE.
     RETURN
  ELSE
     sqrtX    = SQRT(sqrtX)
  ENDIF 
  !
  rho =D*sqrtX/den
  eps =sqrtX-p*den/sqrtX-D
  eps=eps/D
  !
  den2 = den**2 
  !
  drhop = D*s2/den2/sqrtX
  depsp = p*s2/(den2-s2)/den/rho
  !
  ! define the functional and its partial derivatives.
  !
  !CALL EOS_p_dprho_dpeps(p1,dprho,dpeps,rho,eps)  
  !
  SELECT CASE(EOStype)
  CASE(-1)
      ! polytopic GASES:
      p1 = NSTOV_kappa*rho**EQN%gamma  
      dprho = EQN%gamma*NSTOV_kappa*rho**(EQN%gamma-1.0)
      dpeps = 0.0
      !
  CASE(0)
      ! IDEAL GASES:
      p1 = (EQN%gamma -1.0)*eps*rho
      dprho = (EQN%gamma -1.0)*eps
      dpeps = (EQN%gamma -1.0)*rho
      !
  CASE(1)
      ! NON-IDEAL GASES:  p^2 = (gamma -1)*eps*rho^1.5
      p1 = SQRT((EQN%gamma -1.0)*eps*rho*SQRT(rho))
      dprho = SQRT((EQN%gamma -1.0)*eps)*0.75/rho**0.25
      dpeps = SQRT((EQN%gamma -1.0)*rho*SQRT(rho))*0.5/SQRT(eps)
  CASE(2)
#ifdef GEOSDEBUG
      p1 = (EQN%gamma -1.0)*eps*rho
      dprho = (EQN%gamma -1.0)*eps
      dpeps = (EQN%gamma -1.0)*rho
      RETURN
#endif
      ! tabulated-DG
      !PRINT *,'***ERROR:  DG-mode EOS is available only in vector-mode!'
      !STOP
      ! tabulated DG (only vectorized)
      CALL EvaluateEOS_scalar(p1,dprho,dpeps,T,rho,eps)
      !
      !
  CASE DEFAULT
      ! IDEAL GASES:
      p1 = (EQN%gamma -1.0)*eps*rho
      dprho = (EQN%gamma -1.0)*eps
      dpeps = (EQN%gamma -1.0)*rho
      !
  END SELECT
  !
  SELECT CASE(EOSfunc)
  CASE(0) ! p-p[rho,eps] 
      continue
      !
  	  F1 = p - p1
      dF1 = 1.0 - dprho*drhop - dpeps*depsp
	  !
  CASE(1) ! p^2-p^2[rho,eps]  
      dprho=2*p1*dprho
      dpeps=2*p1*dpeps
      p1=p1**2
      !
      F1 = p**2 - p1
      dF1 = 2*p - dprho*drhop - dpeps*depsp
	  !
  CASE DEFAULT ! p-p[rho,eps]  
      continue
      !
      F1 = p - p1
      dF1 = 1.0 - dprho*drhop - dpeps*depsp
  END SELECT 
  !
  continue
  !
END SUBROUTINE FUNC_C2P_RHD1_VECTORF


 

RECURSIVE SUBROUTINE RTSAFE_C2P_RHD1_VECTOR(p,X1,X2,XACC,d,tau,s2,FAILED)
  USE VSMM_mod 
  IMPLICIT NONE
  INTEGER, PARAMETER    :: MAXIT=400
  INTEGER               :: i,j 
  REAL(8)               :: p(VECTORLENGTH), XG(VECTORLENGTH), FG(VECTORLENGTH), Vtmp1(VECTORLENGTH), Vtmp2(VECTORLENGTH) 
  REAL(8)               :: X1(VECTORLENGTH), X2(VECTORLENGTH), XACC(VECTORLENGTH), XACC2(VECTORLENGTH), gam(VECTORLENGTH)
  REAL(8)               :: d(VECTORLENGTH),tau(VECTORLENGTH),s2(VECTORLENGTH)
  REAL(8)               :: FL(VECTORLENGTH),FH(VECTORLENGTH),DF(VECTORLENGTH),XH(VECTORLENGTH),XL(VECTORLENGTH),SWAP(VECTORLENGTH)
  REAL(8)               :: DXOLD(VECTORLENGTH),DX(VECTORLENGTH),F(VECTORLENGTH),TEMP(VECTORLENGTH),Xtmp(VECTORLENGTH),Ftmp(VECTORLENGTH),dFtmp(VECTORLENGTH),dXtmp(VECTORLENGTH) 
  LOGICAL(8)            :: FAILED(VECTORLENGTH) 
  !
#ifdef AVX512
    !DIR$ ATTRIBUTES ALIGN: 64:: FL,FH,DF,XH,XL,SWAP,DXOLD,DX,F,TEMP,Xtmp,Ftmp,dFtmp,dXtmp,XACC2
    !DIR$ ASSUME_ALIGNED p : 64 
    !DIR$ ASSUME_ALIGNED X1 : 64 
    !DIR$ ASSUME_ALIGNED X2 : 64 
    !DIR$ ASSUME_ALIGNED XACC : 64 
    !DIR$ ASSUME_ALIGNED d    : 64 
    !DIR$ ASSUME_ALIGNED tau    : 64 
    !DIR$ ASSUME_ALIGNED s2   : 64 
    !DIR$ ASSUME_ALIGNED FAILED    : 64 
#else
    !DIR$ ATTRIBUTES ALIGN: 32 :: FL,FH,DF,XH,XL,SWAP,DXOLD,DX,F,TEMP,Xtmp,Ftmp,dFtmp,dXtmp,XACC2
    !DIR$ ASSUME_ALIGNED p   : 32 
    !DIR$ ASSUME_ALIGNED X1   : 32 
    !DIR$ ASSUME_ALIGNED X2   : 32 
    !DIR$ ASSUME_ALIGNED XACC : 32 
    !DIR$ ASSUME_ALIGNED d    : 32 
    !DIR$ ASSUME_ALIGNED tau    : 32 
    !DIR$ ASSUME_ALIGNED s2   : 32 
    !DIR$ ASSUME_ALIGNED FAILED    : 32  
#endif 
  !
  XACC2(:) = 1e-60 !XACC(:)**2
  FAILED(:) = .FALSE.
  !!
  XG = p ! initial guess
  ! do bracketing with the initial guess: 
  CALL FUNC_C2P_RHD1_VECTORF(XG,FG,DF,d,tau,s2,FAILED)
  WHERE( ISNAN(FG) )
      FG = 1e20
  ENDWHERE 
  IF(ANY(ABS(FG(:)).GT.XACC2(:))) THEN
      continue
  ELSE 
     p(:)=XG(:)
     RETURN
  ENDIF
  !
  
  !! First try four Newton iterations. 
  !!   
  !DO j = 1 , 5 
  !    CALL FUNC_C2P_RHD1_VECTORF(p,F,DF,d,tau,s2,FAILED)
  !    WHERE( ISNAN(F) )
  !        F = 1e20
  !    ENDWHERE 
  !    IF( ALL(ABS(F).LT.XACC) ) THEN
  !        RETURN
  !    ENDIF 
  !    p(:) = p(:) - F(:)/DF(:)  
  !ENDDO 
  !!
  CALL FUNC_C2P_RHD1_VECTORF(X1,FL,DF,d,tau,s2,FAILED)
  WHERE( ISNAN(FL) )
      FL = 1e20
  ENDWHERE 
        !IF(ANY(ABS(FL(:)).GT.XACC2(:))) THEN
        !  WHERE(ABS(FL(:)).LT.XACC2(:))  
        !     p(:) = X1(:) 
        !     X2(:) = X1(:) 
        !  ENDWHERE
        !ELSE 
        !   p(:)=X1(:)
        !   RETURN
        !ENDIF
        ! 
  CALL FUNC_C2P_RHD1_VECTORF(X2,FH,DF,d,tau,s2,FAILED)
  WHERE( ISNAN(FH) )
      FH = 1e20
  ENDWHERE 
        !IF(ANY(ABS(FH(:)).GT.XACC2(:))) THEN
        !      WHERE(ABS(FH(:)).LT.XACC2(:))  
        !         p(:) = X2(:) 
        !         X1(:) = X2(:) 
        !      ENDWHERE 
        !ELSE 
        !   p(:)=X2(:)
        !   RETURN
        !ENDIF
  WHERE(FL(:)*FH(:).GT.0.)  
     FAILED(:) = .TRUE.
     p(:)     = 0. ! assign value even if it fails     
  ENDWHERE 
  WHERE(FL(:).LT.0.0)
        WHERE(FG(:).LT.0.0)  
            XL(:) = XG(:) 
            XH(:) = X2(:)
        ELSEWHERE   
            XL(:) = X1(:) 
            XH(:) = XG(:)
        ENDWHERE
  ELSEWHERE 
     XH(:) = X1(:)
     XL(:) = X2(:) 
     SWAP(:) = FL(:) 
     FL(:) = FH(:) 
     FH(:) = SWAP(:)
        WHERE(FG(:).LT.0.0)  
            XL(:) = XG(:)  
        ELSEWHERE    
            XH(:) = XG(:)
        ENDWHERE
  ENDWHERE 
  !
  ! First try four Newton iterations. 
  !   
  Vtmp1 = ABS(FL)-ABS(FH)
  WHERE(Vtmp1.GT.0.) 
      p = XH
  ELSEWHERE
      p = XL
  ENDWHERE
  ! Do 5 Newton iterations
  DO j = 1,5 
      CALL FUNC_C2P_RHD1_VECTORF(p,F,DF,d,tau,s2,FAILED)
      WHERE( ISNAN(F) )
          F = 1e20
      ENDWHERE 
      !IF(ALL(ABS(F).LT.XACC2)) THEN
      !    RETURN
      !ENDIF 
      DX(:) = F(:)/DF(:)
      IF( ALL(ABS(DX).LT.XACC) ) THEN
        CONTINUE 
        RETURN
      ENDIF
      p(:) = MIN(p(:) - DX(:),XH)
  ENDDO
  !
  !p(:)=0.5*(X1(:)+X2(:))
  !
  DXOLD(:)=DX(:) !ABS(X2(:)-X1(:))
  DX(:)=DXOLD(:)
  CALL FUNC_C2P_RHD1_VECTORF(p,F,DF,d,tau,s2,FAILED)  
  WHERE( ISNAN(F) )
      F = 1e20
  ENDWHERE 
  DO j = 1, MAXIT 
      Vtmp1 = ((p-XH)*DF-F)*((p-XL)*DF-F)
      Vtmp2 = ABS(2.*F)-ABS(DXOLD*DF) 
     WHERE(Vtmp1.GE.0. .OR. Vtmp2.GT.0. )  
        DXOLD(:) = DX(:)
        DX(:) = 0.5*(XH(:)-XL(:)) 
        p(:)=XL(:)+DX(:) 
     ELSEWHERE 
        DXOLD(:) = DX(:)
        DX(:)=F(:)/DF(:) 
        TEMP(:)=p(:)
        p(:)=p(:)-DX(:) 
     ENDWHERE 
     IF( ALL(ABS(DX).LT.XACC) ) THEN
        CONTINUE 
        RETURN
     ENDIF
     CALL FUNC_C2P_RHD1_VECTORF(p,F,DF,d,tau,s2,FAILED)
     WHERE( ISNAN(F) )
         F = 1e20
     ENDWHERE 
     WHERE(F.LT.0.)  
        XL(:) = p(:)
        FL(:) = F(:) 
     ELSEWHERE 
        XH(:)=p(:)
        FH(:)=F(:)
     ENDWHERE 
#ifdef C2PFF  
     WHERE(ABS(F).LT.XACC2)
#else
     WHERE(F.EQ.0.)
#endif 
         DF = 0. 
     ENDWHERE
    ! 
  ENDDO 
  !
  WHERE(ABS(DX).GT.1e-10)
      FAILED = .TRUE.
      p = 0. ! assign value even if it fails
  ELSEWHERE
      FAILED = .FALSE.
  ENDWHERE
  !
  !IF(ANY(FAILED)) THEN
  !    continue
  !ENDIF 
  RETURN
  ! 
END SUBROUTINE   RTSAFE_C2P_RHD1_VECTOR
   !
  
  
  


END MODULE EOS_mod