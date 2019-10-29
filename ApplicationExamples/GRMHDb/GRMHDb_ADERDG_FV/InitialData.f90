! GRMHDb Initial Data
 

RECURSIVE SUBROUTINE PDESetup(myrank)
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters 
#ifdef RNSTOV
    USE NSTOV_mod
#endif  
#ifdef GEOS
    USE EOS_par
    USE GREOS
#endif
    IMPLICIT NONE
	INTEGER, INTENT(IN)            :: myrank
!#ifndef RNSTOV   
    SELECT CASE(TRIM(ICTYPE))    ! 
    CASE('GRMHDTOV','CCZ4TOV','GRMHDTOV_perturbed') 
        !
        continue
        !
        CASE DEFAULT 
            NSTOV_rho_atmo = 1e-10
            NSTOV_p_atmo = 1e-10 
            !
    END SELECT
    !
!#endif  
    !
#ifdef GEOS
    !
    SELECT CASE(EOStype)
    CASE(2)
        CALL InitGREOS
        CALL ReadEOSCoefficients  
    CASE DEFAULT
        continue
    END SELECT
    !
#endif
    ! 
	myrank_F90 = myrank
    !
    EQN%Pi    = ACOS(-1.0)                                          ! Pi 
    !EQN%gamma = 1.4                                                 ! ratio of specific heats for compressible Euler
	!EQN%ch = 1.0
	EQN%DivCleaning_a = 1.0
    !            
    !ExcisionRadius = -1.0     
    ! 
    !GRMHD-BondiAccretion
    EQN%gamma = 4./3. 
    EQN%cs = 1.0                                 ! cs  
    EQN%tau = 0.2                                 ! tau 
    !GRMHD-TOV
   ! EQN%gamma = 2.0
    !EQN%cs = 0.2                                 ! cs  
    !EQN%tau = 0.2                             ! tau  
    !
    !!AlfvenWave SRMHD limit
    !EQN%gamma = 2.0  
	! 
	!
	SELECT CASE(TRIM(ICTYPE))
    CASE('Sod')
    	EQN%gamma = 5.0/3.0
    CASE('Riemann0')
        EQN%gamma = 5.0/3.0
	CASE('GRMHDAccretion')
		EQN%gamma = 4./3. 
	CASE('GRMHDTOV')
		EQN%gamma = 2.0
		!ICTYPE   = 'GRMHDTOV'
		!ICTYPE   = 'GRMHDTOV_perturbed'
	CASE('GRMHDTorus')
		EQN%gamma = 4./3.
		!    
		CASE DEFAULT
			continue
    END SELECT
	!
	!IF(myrank.eq.0) THEN
    PRINT *, "<<<<<<<<<<<<<<<<<<<<<<<<<<<------------------------------------"
		PRINT *, 'EQN%gamma=',EQN%GAMMA, myrank, myrank_f90 
    PRINT *, "<<<<<<<<<<<<<<<<<<<<<<<<<<<------------------------------------"
	!ENDIF
    !
    ! 
    SELECT CASE(TRIM(ICTYPE))    ! 
    CASE('GRMHDTOV','CCZ4TOV','GRMHDTOV_perturbed') 
#ifdef RNSTOV
        IF(NSTOVVar%Computed.NE.12345) THEN
            NSTOVVar%Computed=12345   
            CALL NSTOV_Main
        ENDIF
#endif  
	CASE DEFAULT
		continue
    ! 
    END SELECT 
    !
    !
    continue
    !
END SUBROUTINE PDESetup



RECURSIVE SUBROUTINE InitialData(x, t, Q)
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar, nDim,myrank_f90
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: x(nDim), t        ! 
	REAL, INTENT(OUT)              :: Q(nVar)        ! 
	! local variables
	INTEGER :: i
	!
	! DO i=1,nVar 
		! Q(i) = 0.
	! ENDDO
	! Q(1) = 1.0
	! Q(5) = 1.0
	! Q(10) = 1.0
	! Q(14) = 1.0
	! Q(17) = 1.0
	! Q(19) =  1.0
	
	! return
	
	! STOP
	
	call InitialField(x, t, Q);
	!!!call MichelAccretionSpherical(x, t, Q);
	!IF(myrank_f90.eq.0) THEN
	!	PRINT *,'myrank_f90=', myrank_f90
	!	PRINT *, 'x=', x
	!	PRINT *, 'Q =', Q
	!ENDIF
	!
	!IF(MAXval(abs(Q(:))).lt.1e-8) then 
	!		print *, 'fatal error: InitialData' 
	!		write(*,*) Q(1:8)
	!		stop 
	!endif
	!IF(ANY(ISNAN(Q(:)))) then 
	!		print *, 'fatal error:ISNAN InitialData' 
	!		write(*,*) x
	!		write(*,*) t
	!		write(*,*) Q(1:8)
	!		write(*,*) Q(9:19)
	!		stop 
	!endif
	!WRITE(*,*) t,x,Q(1)
	!
END SUBROUTINE InitialData


RECURSIVE SUBROUTINE PDEdefineobservables(numberOfObservables,observables,Q)
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar, nDim
	IMPLICIT NONE 
	! Argument list 
	INTEGER, INTENT(IN)	:: numberOfObservables
	REAL, INTENT(INOUT) :: observables(numberOfObservables)        !
	REAL, INTENT(IN)  :: Q(nVar)   
	
	observables = 1.0
	! compute observables from Q
	
	
END SUBROUTINE PDEdefineobservables


RECURSIVE SUBROUTINE PDEassurepositivity(limiter_value,xx,numberOfObservables, observablesMin, observablesMax)
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar, nDim
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: xx(nDim)        ! 
	INTEGER, INTENT(IN)					:: numberOfObservables
	INTEGER, INTENT(OUT)              :: limiter_value        !
	REAL, INTENT(IN)					:: observablesMin(numberOfObservables), observablesMax(numberOfObservables)
	real	:: rr	

	limiter_value=0
	return
	
END SUBROUTINE PDEassurepositivity


RECURSIVE SUBROUTINE InitialField(xGP,tGP,u0) 
    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim, d, ICType, EQN, ICType2,ExcisionRadius,NSTOVVar,NSTOVVar_bar,NSTOV_kappa  , aom, Mbh,NSTOV_ATMO,NSTOV_t_atm,myrank_f90
    USE AstroMod
	!USE Bessel_mod 
#ifdef TWOPUNCTURES  
	USE TwoPunctures_mod, ONLY : TwoPunctures_Interpolate
#endif 
#ifdef RNSID 
    USE RNSID_mod, ONLY : RNSID_Interpolate
#endif
#ifdef RNSTOV
  USE NSTOV_mod
#endif
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN ) :: xGP(nDim), tGP        ! spatial position vector and time 
    REAL, INTENT(OUT) :: u0(nVar)           ! initial data vector in terms of conserved variables 
    ! Local variables 
    INTEGER :: i,j,k,l,nm,iNewton,Maxnewton
    REAL :: VBase(nVar), ampl(nVar), sigma(d) ,xGP3D(d)
    REAL :: V0(nVar),r,VLL(nVar),VRR(nVar), VZ4(54),V70(70)  
    REAL :: de,u0_test(nVar),v0_tesT(nVar)
    REAL :: du,dv,drho,dTemp,dp,epsilon,xc(d),rbar
    REAL :: omega,igamma,gamma2,gamma1,rho_atm,h,t,p_atm
    REAL :: r1(9), lambda, lambdat, lambdax, mu, rho, cs, cp, ICA, HH, dxH, Kxx,A11
    REAL :: theta, a,  M, zz, r_hz, delta, g_cov(3,3), g_contr(3,3), g_tmp(3,3), Kex(3,3), Aex(3,3) 
    REAL :: traceK, alpha, fa, ggg, fff, k0, sig, AA(3), BB(3,3), BMat(3,3), DD(3,3,3), PP(3)  
    REAL :: rho0, p0, eta, b0, tempaa, tempab, tempac, va2, vax, bv(3), BV_contr(3), vv(3), gm, gp, ddet(3) 
    REAL :: beta(3), betadown(3), dbetadown(3,3), nablabeta(3,3), test(3,3), Gtilde(3)   
    REAL :: Christoffel(3,3,3), dxphi, bbb, dxb, dtb, kyy, kzz, Aref, x0, re0, ms, uu0  
    REAL :: vv_cov(3), shift(3), p, lapse, gammaij(6), ng, up, phi,phi2, dphi(3), rc, vc2, tc, pc, rhoc, lf, vtheta, vphi, vx, vy, vz 
    REAL :: urr, tt, c1, c2, df, vr, ut, f, urc, g_tt, betaru, dtt, detg, psi , xloc(d), minR
    REAL :: BJ(0:1),DJ(0:1),BY(0:1),DY(0:1) 
    REAL :: IJ(0:1),DIJ(0:1),IY(0:1),DIY(0:1) 
    REAL, PARAMETER :: Pi = ACOS(-1.0) 
    REAL :: xGP_sph(3),xGP_loc(3),A_coord(3,3), detA,iA_coord(3,3),A_coord_contr(3,3),iA_coord_contr(3,3),v_cov(3),v_contr(3),v2
    INTEGER :: iErr
    !
    Maxnewton = 50 
    u0  = 0.0 
    !
    xGP3D = 0.
    xGP3D(1:nDim) = xGP(1:nDim)
    ! no local material parameters for Euler equations 
    ! Gaussian perturbation 
    SELECT CASE(TRIM(ICType)) 
#if defined(GRMHD) 
    CASE('Sod')
        !CALL Curved2Cartesian(xGP_loc,xGP3D)          ! no effects if you are in Cartesian coordiantes
        !CALL Cartesian2Spherical(xGP_sph,xGP_loc)   ! here we convert Cartesian to Spherical, independently on the chosen coordiantes
        !PRINT *,"SOD"
        r=xGP3D(1)
        if(r.LT.0.5) THEN
            v0 = 0.
            v0(1) = 1.0
            v0(5) = 1.0
            v0(10) = 1.0
            v0(14) = 1.0
            v0(17) = 1.0
            v0(19) = 1.0  
        ELSE
            v0 = 0.
            v0(1) = 0.125
            v0(5) = 0.1
            v0(10) = 1.0
            v0(14) = 1.0
            v0(17) = 1.0
            v0(19) = 1.0   
        ENDIF
        !
        !v0(2) = 0.05        
        
       ! PRINT *,"r, rho", r,V0(1)
        !
    CASE('Riemann0')
        !CALL Curved2Cartesian(xGP_loc,xGP3D)          ! no effects if you are in Cartesian coordiantes
        !CALL Cartesian2Spherical(xGP_sph,xGP_loc)   ! here we convert Cartesian to Spherical, independently on the chosen coordiantes
        r=xGP3D(1)
        if(r.LT.0.5) THEN
            v0 = 0.
            v0(1) = 1.0
            v0(5) = 1.0
            v0(10) = 1.0
            v0(14) = 1.0
            v0(17) = 1.0
            v0(19) = 1.0 
        ELSE
            v0 = 0.
            v0(1) = 0.125
            v0(5) = 1.0
            v0(10) = 1.0
            v0(14) = 1.0
            v0(17) = 1.0
            v0(19) = 1.0 
        ENDIF
        v0(2)=0.2
        !
    CASE('GRMHDAlfvenWave') 
       rho0 = 1.
       p0   = 1.
       eta  = 1. 
       B0   = 1.0  
       !
       hh = 1.0 + EQN%gamma / ( EQN%gamma - 1.0) * p0 / rho0
       tempaa = rho0 * hh + B0**2 * ( 1.0 + eta**2)
       tempab = 2.0 * eta * B0**2 / tempaa
       tempac = 0.5 * ( 1.0 + sqrt ( 1.0 - tempab**2))
       va2 = b0**2 / ( tempaa * tempac)
       vax = sqrt ( va2) 
       !
       BV(1) = B0
       BV(2) = eta * B0 * COS( xGP3D(1) - vax*tGP)
       BV(3) = eta * B0 * SIN( xGP3D(1) - vax*tGP)
       !
       VV(1)   = 0.0
       VV(2:3) = - vax * BV(2:3) / B0
       !
       alpha=1.0
       gammaij = 0.0 
       gammaij(1) = 1.0 
       gammaij(4) = 1.0 
       gammaij(6) = 1.0  
       ! Now convert to conservative variables 
       !
       V0(1) = rho0
       V0(2) = VV(1)
       V0(3) = VV(2)
       V0(4) = VV(3)
       V0(5) = p0
       V0(6) = BV(1)
       V0(7) = BV(2)
       V0(8) = BV(3)
       V0(9) = 0.       
       !
        V0(10)  = alpha 
        !
        DO i=1,3 
            V0(10+i) = 0.
        ENDDO
        !
        DO i=1,6
            V0(13+i) = gammaij(i)
        ENDDO
        !
       CONTINUE
       !
    
       !
    CASE('GRMHDTOV')
#ifdef RNSTOV
        CALL Curved2Cartesian(xGP_loc,xGP3D)          ! no effects if you are in Cartesian coordiantes
        CALL Cartesian2Spherical(xGP_sph,xGP_loc)   ! here we convert Cartesian to Spherical, independently on the chosen coordiantes
        r=XGP_sph(1)
        ! 
#ifdef SPHERICAL  ! then we are in the original coordinate system        
        CALL NSTOV_x(r,NSTOVVar%qloc)
#elif CYLINDRICAL
        PRINT *, 'CYLINDRICAL COORDINATES NOT TESTED FOR RNSTOV'
#else   
        IF(r.LT.0.) THEN
            PRINT *, 'FATAL ERROR: negative radius, R=',r
        ENDIF
        !
        CALL NSTOV_rbar(r,NSTOVVar%qloc)
#endif
        !
        !p = NSTOVVar%qloc(3)
        p = MAX( NSTOVVar%p_R, NSTOVVar%qloc(3) ) 
        igamma = 1.0/EQN%gamma 
        rho=(p/NSTOV_kappa)**igamma 
        VV_cov(1:3) = 0.
        BV_contr(1:3) = 0.
        !
        CALL METRIC(  xGP3D, alpha,  gp, gm, shift, Kex, g_cov, g_contr, phi )
        CALL DMETRIC( xGP3D, AA, BB, DD, ddet )
        V0(20)  = AA(1) 
        !! 
        gammaij(1) = g_cov(1,1)
        gammaij(2) = g_cov(1,2)
        gammaij(3) = g_cov(1,3)
        gammaij(4) = g_cov(2,2)
        gammaij(5) = g_cov(2,3)
        gammaij(6) = g_cov(3,3) 
        !
        ! BUILD THE ATMOSPHERE:
        IF(r.GT.NSTOVVar_bar%radius) THEN
            SELECT CASE(NSTOV_ATMO)
            CASE(1) ! hydrostatic (Landau) atmo.
                ! -------------------------------
                ! Static (Landau) atmosphere 
                rho_atm= 0.01*NSTOVVar%rho_R   !jump_l * rho_c
                gamma1 = EQN%gamma / ( EQN%gamma - 1.0)
                h      = (1.0 + gamma1*NSTOV_t_atm) * NSTOVVar%lapse_C/alpha
                t      = (h - 1.0) / gamma1
                gamma2 = 1.0 / ( EQN%gamma - 1.0)
                rho    = rho_atm * ( t / NSTOV_t_atm)**gamma2
                p      = rho*t
                !vr     = 0.0
                !vtheta = 0.0
                !vphi   = 0.0
                !br     = 0.0
                !btheta = 0.0
                !bphi   = 0.0 
                continue
                !
            CASE(2) ! hydrostatic (Landau) atmo. B
                ! -------------------------------
                ! Static (Landau) atmosphere 
                gamma1 = EQN%gamma / ( EQN%gamma - 1.0)
                h      = (1.0 + gamma1*NSTOV_t_atm) * NSTOVVar%lapse_C/alpha
                t      = (h - 1.0) / gamma1
                gamma2 = 1.0 / ( EQN%gamma - 1.0) 
                
                p_atm= 0.01*NSTOVVar%p_R   !jump_l * rho_c
                p    = p_atm * ( t / NSTOV_t_atm)**gamma2 
                
                !p      = rho*t
                rho = p/t
                !vr     = 0.0
                !vtheta = 0.0
                !vphi   = 0.0
                !br     = 0.0
                !btheta = 0.0
                !bphi   = 0.0 
                continue
                !
                CASE DEFAULT
                    ! const. floor. atmo 
                    !
            END SELECt
        ENDIF
        !
        !
        !
        !
        BV(1:3) = 0.
        !
        !V0 = (/ rho, VV_cov(1:3), p ,BV(1:3), 0.0, alpha, shift(1:3), gammaij(1:6)/)
        !
        V0(1) = rho
        V0(5) = p
        V0(9) = 0.
        V0(10)  = alpha 
        !
        DO i=1,3
            V0(1+i) = VV_cov(i) 
            V0(5+i) = BV(i) 
            V0(10+i) = 0.
        ENDDO
        !
        DO i=1,6
            V0(13+i) = gammaij(i)
        ENDDO
        !
        !V0(9)=exp(-(xGP3D(1)**2 + xGP3D(2)**2 + xGP3D(3)**2) / 8.0)*0.002
        !PRINT *,'V0=',V0
#else
        PRINT *,' GRMHDTOV not implemented for your PDE'
        STOP
#endif
       continue
       ! 
       !
    CASE('GRMHDTOV_perturbed')
#ifdef RNSTOV
        CALL Curved2Cartesian(xGP_loc,xGP3D)          ! no effects if you are in Cartesian coordiantes
        CALL Cartesian2Spherical(xGP_sph,xGP_loc)   ! here we convert Cartesian to Spherical, independently on the chosen coordiantes
        r=XGP_sph(1)
        ! 
#ifdef SPHERICAL  ! then we are in the original coordinate system        
        CALL NSTOV_x(r,NSTOVVar%qloc)
#elif CYLINDRICAL
        PRINT *, 'CYLINDRICAL COORDINATES NOT TESTED FOR RNSTOV'
#else  
        CALL NSTOV_rbar(r,NSTOVVar%qloc)
#endif
        !
        !p = NSTOVVar%qloc(3)
        p = MAX( NSTOVVar%p_R, NSTOVVar%qloc(3) ) 
        igamma = 1.0/EQN%gamma 
        rho=(p/NSTOV_kappa)**igamma 
        ! then add a density perturbation and recompute the pressure.
        rho = rho*(1 + 0.1* cos(pi/2.0 * r/NSTOVVar_bar%radius))
        p=rho**EQN%gamma*NSTOV_kappa
        p = MAX( NSTOVVar%p_R,p ) 
        
        VV_cov(1:3) = 0.
        BV_contr(1:3) = 0.
        !
        CALL METRIC(  xGP3D, alpha,  gp, gm, shift, Kex, g_cov, g_contr, phi )
        CALL DMETRIC( xGP3D, AA, BB, DD, ddet )
        V0(20)  = AA(1) 
        !! 
        gammaij(1) = g_cov(1,1)
        gammaij(2) = g_cov(1,2)
        gammaij(3) = g_cov(1,3)
        gammaij(4) = g_cov(2,2)
        gammaij(5) = g_cov(2,3)
        gammaij(6) = g_cov(3,3) 
        !
        ! BUILD THE ATMOSPHERE:
        IF(r.GT.NSTOVVar_bar%radius) THEN
            SELECT CASE(NSTOV_ATMO)
            CASE(1) ! hydrostatic (Landau) atmo.
                ! -------------------------------
                ! Static (Landau) atmosphere 
                rho_atm= 0.01*NSTOVVar%rho_R   !jump_l * rho_c
                gamma1 = EQN%gamma / ( EQN%gamma - 1.0)
                h      = (1.0 + gamma1*NSTOV_t_atm) * NSTOVVar%lapse_C/alpha
                t      = (h - 1.0) / gamma1
                gamma2 = 1.0 / ( EQN%gamma - 1.0)
                rho    = rho_atm * ( t / NSTOV_t_atm)**gamma2
                p      = rho*t
                !vr     = 0.0
                !vtheta = 0.0
                !vphi   = 0.0
                !br     = 0.0
                !btheta = 0.0
                !bphi   = 0.0 
                continue
                !
            CASE(2) ! hydrostatic (Landau) atmo. B
                ! -------------------------------
                ! Static (Landau) atmosphere 
                gamma1 = EQN%gamma / ( EQN%gamma - 1.0)
                h      = (1.0 + gamma1*NSTOV_t_atm) * NSTOVVar%lapse_C/alpha
                t      = (h - 1.0) / gamma1
                gamma2 = 1.0 / ( EQN%gamma - 1.0) 
                
                p_atm= 0.01*NSTOVVar%p_R   !jump_l * rho_c
                p    = p_atm * ( t / NSTOV_t_atm)**gamma2 
                
                !p      = rho*t
                rho = p/t
                !vr     = 0.0
                !vtheta = 0.0
                !vphi   = 0.0
                !br     = 0.0
                !btheta = 0.0
                !bphi   = 0.0 
                continue
                !
                CASE DEFAULT
                    ! const. floor. atmo 
                    !
            END SELECt
        ENDIF
        !
        !
        !
        !
        BV(1:3) = 0.
        !
        !V0 = (/ rho, VV_cov(1:3), p ,BV(1:3), 0.0, alpha, shift(1:3), gammaij(1:6)/)
        !
        V0(1) = rho
        V0(5) = p
        V0(9) = 0.
        V0(10)  = alpha 
        !
        DO i=1,3
            V0(1+i) = VV_cov(i) 
            V0(5+i) = BV(i) 
            V0(10+i) = 0.
        ENDDO
        !
        DO i=1,6
            V0(13+i) = gammaij(i)
        ENDDO
        !
#else
        PRINT *,' GRMHDTOV not implemented for your PDE'
        STOP
#endif
       continue
       ! 
       !
#ifdef TORUS        
    CASE('GRMHDTorus')
       xloc = xGP3D 
       !
       minR = 0.8*ExcisionRadius
       !
#ifdef SPHERICAL  
       !
       r=xloc(1)
       IF ( r .LT. minR) THEN ! RETURN
           r = minR
           xloc(1) = r
       ENDIF
       !
#else
       !
       r=SQRT(SUM(xloc(1:nDim)**2))
       !
       IF ( r .LT. minR) THEN ! RETURN
           r=minR
           xloc(1) = r*DSIN(THETA)*DCOS(PHI2)
           xloc(2) = r*DSIN(THETA)*DSIN(PHI2)
           xloc(3) = r*DCOS(THETA) 
       ENDIF
       !
#endif 
       !
       CALL METRIC ( xloc , lapse, gp, gm, shift, Kex, g_cov, g_contr, phi )
       ng     = 1.0/(EQN%gamma - 1.0)

       gammaij(1) = g_cov(1,1)
       gammaij(2) = g_cov(1,2)
       gammaij(3) = g_cov(1,3)
       gammaij(4) = g_cov(2,2)
       gammaij(5) = g_cov(2,3)
       gammaij(6) = g_cov(3,3) 
       !
       CALL MatrixInverse3x3(g_cov,g_contr,gp)
       !
#ifdef SPHERICAL  
           !
           CALL init_disk_isentropic(xloc, V0(1:8),METRIC_KSS) 
           V0(9) = 0.
           V0(10:19) = (/  lapse, shift(1:3), gammaij(1:6)/)
           ! 
#else  
       ! The Following is for Kerr-Schild Cartesian coordinates 
       r      = SQRT( xloc(1)**2 + xloc(2)**2 + xloc(3)**2) 
       theta  = ACOS( xloc(3)/r)
       phi    = ATAN2( xloc(2), xloc(1))
       ! 
       minR = 0.8*ExcisionRadius
       IF(r.LT.minR) THEN  ! RETURN
           r=ExcisionRadius
           xloc(3) = r*DCOS(theta)
           xloc(2) = r*DSIN(theta)*DSIN(phi)
           xloc(1) = r*DSIN(theta)*DCOS(phi)
       ENDIF 
       !
       xGP_sph = (/r,theta, phi   /) 
       !
       CALL init_disk_isentropic(xGP_sph, V0(1:8),METRIC_KSS)
       !    
        CALL Cart2SphMatrix_cov(A_coord,xGP_sph)
        ! 
        CALL MatrixInverse3x3(A_coord,iA_coord,detA)      ! iA = Sph2CartMatrix_cov    = Cart2SphMatrix_contr
        A_coord_contr(:,:) = TRANSPOSE(iA_coord(:,:))     ! Cart2SphMatrix_contr = TRANSPOSE( Sph2CartMatrix_cov )
        iA_coord_contr(:,:) = TRANSPOSE(A_coord(:,:))    ! Sph2CartMatrix_contr = TRANSPOSE( Cart2SphMatrix_cov )
        !
        VV_cov(1:3) = MATMUL(iA_coord,V0(2:4))  ! SPHERICAL to Cartesian: Velocity is covariant
        BV_contr(1:3) = MATMUL(iA_coord_contr,V0(6:8))  !  SPHERICAL to Cartesian:  Magnetic field is contravariant  
        !
        V0(2:4)=VV_cov(1:3)
        V0(6:8)=BV_contr(1:3)
        !   
        V0(10:19) = (/  lapse, shift(1:3), gammaij(1:6)/)
        !
        IF(MAXVAL(ABS(V0(2:4))).GT.1e-10) THEN
            continue
        ENDIF
        IF(MAXVAL(ABS(V0(10:19))).LT.1e-10) THEN
            continue
        ENDIF
        ! 
#endif        
#endif 

    CASE('GRMHDRNSID')
        ! 
        V0 = 0.0
        ! Compute the metric and its derivatives 
        CALL METRIC(  xGP3D, alpha,  gp, gm, beta, Kex, g_cov, g_contr, phi )
        CALL DMETRIC( xGP3D, AA, BB, DD, PP )
        test = matmul(g_cov,g_contr)
        detg = g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)*g_cov(3,2)-g_cov(2,1)*g_cov(1,2)*g_cov(3,3)+g_cov(2,1)*g_cov(1,3)*g_cov(3,2)+g_cov(3,1)*g_cov(1,2)*g_cov(2,3)-g_cov(3,1)*g_cov(1,3)*g_cov(2,2)
        !phi = detg**(-1./6.) 
        !
#ifdef RNSID   
        CALL RNSID_Interpolate(xGP3D, V70) 
        rho         = MAX( 1e-12, V70(60) ) 
        VV_cov(1:3) = V70(61:63) 
        p           = MAX( 1e-12, V70(64) ) 
        BV          = 0.0 
#else
        PRINT *, ' RNSID not available. Please compile with -DRNSID flag. '
        STOP 
#endif      
        
       !
       V0(1:9) = (/ rho, VV_cov(1:3), p, BV(1:3), 0.0 /)
       V0(10)  = alpha 
       V0(11:13) = 0.0 
       V0(14:19) = (/ g_cov(1,1), g_cov(1,2), g_cov(1,3), g_cov(2,2), g_cov(2,3), g_cov(3,3) /)  
       !
       CONTINUE
       !

    CASE('GRMHDAlfven')
       rho = 1.
       p   = 1.
       eta  = 1. 
       B0   = 1.0  
       !
       hh = 1.0 + EQN%gamma / ( EQN%gamma - 1.0) * p / rho
       tempaa = rho * hh + B0**2 * ( 1.0 + eta**2)
       tempab = 2.0 * eta * B0**2 / tempaa
       tempac = 0.5 * ( 1.0 + sqrt ( 1.0 - tempab**2))
       va2 = b0**2 / ( tempaa * tempac)
       vax = sqrt ( va2)       
       !
       ! flat metric: contr and cov are the same in Cartesian
       BV_contr(1) = B0
       BV_contr(2) = eta * B0 * COS( xGP3D(1) - vax*tGP)
       BV_contr(3) = eta * B0 * SIN( xGP3D(1) - vax*tGP)
       !
       VV_cov(1)   = 0.0
       VV_cov(2:3) = - vax * BV_contr(2:3) / B0
       !
       ! Now convert to conservative variables 
       !
       lapse = 1.0
       shift(1:3) = 0.
       gammaij = 0.
       gammaij(1) = 1.
       gammaij(4) = 1.
       gammaij(6) = 1. 
       !       
       V0(1:19) = (/ rho, VV_cov(1:3), p ,BV_contr(1:3), 0.0, lapse, shift(1:3), gammaij(1:6)/) 
       !V0 = (/ 1.0, 0.0, 0.0, 0.0, 1.0 ,0.0, 0.0, 0.0, 0.0, lapse, shift(1:3), gammaij(1:6)/)  
       ! 
    CASE('GRMHDAccretion')       
       !aom = 0.0
       !Mbh = 1.0  
       rc   = 8.0 
       rhoc = 1./16.
       B0 = 0. 
       xloc = xGP3D
       IF ( aom > 0.0) THEN 
           WRITE(*,*)'SPHERICAL Accretion solution is only for a = 0!'
           STOP
       ENDIF 
       !
       ! The Following is for Kerr-Schild Cartesian coordinates
       r      = SQRT( xGP3D(1)**2 + xGP3D(2)**2 + xGP3D(3)**2) 
       theta  = ACOS( xGP3D(3)/r)
       phi2    = ATAN2( xGP3D(2), xGP3D(1))
       !
       !IF ( r .LT. 0.0001) THEN ! RETURN
       !    lapse = 1.0
       !    shift(1:3) = 0.
       !    !
       !    rho = 1.0
       !    p = 1.0
       !    VV_cov = 0.
       !    BV(1:3) = 0.
       !    !
       !    !shift_contr = MATMUL(g_contr,shift)  !shift is controvariant. See fluxes....
       !    !
       !    gammaij = 0.
       !    gammaij(1) = 1.0
       !    gammaij(4) = 1.0
       !    gammaij(6) = 1.0
       !    !  
       !    V0 = (/ rho, VV_cov(1:3), p ,BV(1:3), 0.0, lapse, shift(1:3), gammaij(1:6)/)
       !    !
       !    CALL PDEPrim2Cons(u0,V0)
       !    RETURN
       !ELSE
       !minR = MAX(MAXVAL(dx),0.5*ExcisionRadius) 
       minR = 0.8*ExcisionRadius
       IF ( r .LT. minR) THEN ! RETURN
           r = minR
           xloc(1) = r*DSIN(THETA)*DCOS(PHI2)
           xloc(2) = r*DSIN(THETA)*DSIN(PHI2)
           xloc(3) = r*DCOS(THETA) 
       ENDIF
       !
       CALL METRIC ( xloc , lapse, gp, gm, shift, Kex, g_cov, g_contr, phi )
       ng     = 1.0/(EQN%gamma - 1.0)  
       !
       !
       !
       !phi2    = ACOS( xGP3D(1) / (r*SIN(theta)))
       zz     = 2.0/r   ! we are computing the solution at theta=pi/2
       betaru = zz/(1.0 + zz)
       g_tt   = zz - 1.0
       !
       urc = sqrt(1.0 / (2.0*rc))
       vc2 = urc**2 / (1.0 - 3.0*urc**2)
       tc  = ng*vc2 / ((1.0 + ng)*(1.0 - ng*vc2))
       pc  = rhoc*tc
       
       c1 = urc*tc**ng*rc**2
       c2 = (1.0 + ( 1.0 + ng)*tc)**2*(1.0 - 2.0/rc+urc**2)
       !
       tt = tc
       DO iNewton = 1, MAXNEWTON   
          urr = c1 / (r**2*tt**ng)
          f   = (1.0 + (1.0 + ng)*tt)**2*(1.0 - 2.0/r + urr**2) - c2
          df  = 2.0 * (1.0 + ng)*(1.0 + (1.0 + ng)*tt)*(1.0 - 2.0/r + urr**2) - 2.0*ng*urr**2/tt*(1.0 + (1.0 + ng)*tt)**2
          dtt  = -f/df
          IF (abs(dtt) < 1.e-14) EXIT
          tt = tt + dtt
       ENDDO
       IF(iNewton.gE.MAXNEWTON) THEN
           continue
       ENDIF
       !
       ut     = (-zz*urr + sqrt(urr**2 - zz + 1.0))/(zz - 1.0)
       LF     = lapse*ut
       vr     = ( urr / LF + betaru / lapse)
       vtheta = 0.0
       vphi   = 0.0
       !
       vx  = DSIN(theta)*DCOS(phi2)*vr
       vy  = DSIN(theta)*DSIN(phi2)*vr
       vz  = DCOS(theta)*vr
       !
       VV(1:3) = (/ vx, vy, vz /)
       ! Convert to covariant velocities
       VV_cov = MATMUL(g_cov, VV)
       !
       rho = rhoc*(tt/tc)**ng
       p   = rho*tt
       !
       BV(1:3) = 0.
       BV_contr(1:3) = 0.
       !
       lapse = lapse
       !
       !shift_contr = MATMUL(g_contr,shift)  !shift is controvariant. See fluxes....
       !
       gammaij(1) = g_cov(1,1)
       gammaij(2) = g_cov(1,2)
       gammaij(3) = g_cov(1,3)
       gammaij(4) = g_cov(2,2)
       gammaij(5) = g_cov(2,3)
       gammaij(6) = g_cov(3,3) 
       !
       SELECT CASE(TRIM(ICType2))
       CASE('GRMHD-BondiAccretion')
           !
#ifdef SPHERICAL
           PRINT *, 'GRMHD-BondiAccretio only for Cartesian coordinates!'
           continue
           !
           STOP
#endif           
           !
           BV_contr(1:3) = 0.
           !
           IF(B0.GT.0) THEN
                BV_contr(3) = 2.2688/Mbh*SQRT(B0)
           ENDIF
           !
           BV_contr(1) = gm*BV_contr(3)*Mbh**2/r**3*xloc(1)
           BV_contr(2) = gm*BV_contr(3)*Mbh**2/r**3*xloc(2)
           BV_contr(3) = gm*BV_contr(3)*Mbh**2/r**3*xloc(3)
           !
           !BV(1:3) = MATMUL(g_cov,BV_contr(1:3))
           ! 
           !
           ! 
       CASE('GRMHDCBS') 
           !
#ifndef SPHERICAL
           PRINT *, 'GRMHDCBS only for spherical coordinates!'
           continue
           !
           STOP
#endif           
           !
           BV_contr(1:3) = 0.
           BV_contr(1) = lapse*B0*DCOS(xGP3D(2))
           BV_contr(2) = -  lapse*B0*DSIN(xGP3D(2))/xGP3D(1) 
           !
           !BV(1:3) = MATMUL(g_cov,BV(1:3))
           ! 
           !
           CASE DEFAULT
           !
           BV_contr(1:3) = 0.
           continue
           !
       END SELECT
       !   
       V0(1) = rho
       V0(2:4) = VV_cov(1:3)
       V0(5) = p
       V0(6:8) = BV_contr(1:3)
       V0(9) = 0.0
       V0(10) = lapse
       V0(11:13) = shift(1:3)
       V0(14:19) = gammaij(1:6) 
       ! 
       CONTINUE
#endif
    CASE DEFAULT
       !
       V0(1) = 0.2
       V0(2:4) =0.
       V0(5) = 0.1
       V0(6:8) = 0.
       V0(9) = 0.0
       lapse = 1.0
       shift(1:3) = 0.
       gammaij = 0.
       gammaij(1) = 1.
       gammaij(4) = 1.
       gammaij(6) = 1. 
       V0(10) = lapse
       V0(11:13) = shift(1:3)
       V0(14:19) = gammaij(1:6) 
       !
    END SELECT 
    !
	IF(ANY(ABS(v0(6:8)).GT.1e-14)) THEN
        PRINT *,v0(6:8)
        PRINT *,v0(6:8)
        PRINT *,v0(6:8)
        PRINT *, "I feel magnetized :-( InitialField  v0"
        ERROR STOP
    ENDIF
    !
    !
    CALL PDEPrim2Cons(u0,V0) 
    !!
	IF(ANY(ABS(v0(6:8)).GT.1e-14)) THEN
        PRINT *,v0(6:8)
        PRINT *,v0(6:8)
        PRINT *,v0(6:8)
        PRINT *, "I feel magnetized :-( InitialField  v02"
        ERROR STOP
    ENDIF
    !
	IF(ANY(ABS(u0(6:8)).GT.1e-14)) THEN
        PRINT *,u0(6:8)
        PRINT *,u0(6:8)
        PRINT *,u0(6:8)
        PRINT *, "I feel magnetized :-( InitialField u0 "
        ERROR STOP
    ENDIF
    !
        !PRINT *,"xGP(1), u0", xGP(1), u0
    !V0_test=v0
    !CALL PDECons2Prim(V0_test,u0,iErr) 
    !!
    !CALL PDEPrim2Cons(u0,V0) 
    !!
    !CALL PDEPrim2Cons(u0_test,V0_test) 
    !!
    !de=SQRT(SUM((u0_test-u0)**2))
    !IF(de.GT.1e-10) THEN
    !    continue
    !ENDIF
    !de=SQRT(SUM((v0_test-v0)**2))
    !IF(de.GT.1e-10) THEN
    !    continue
    !ENDIF
    
    CONTINUE
    !
END SUBROUTINE InitialField
    
