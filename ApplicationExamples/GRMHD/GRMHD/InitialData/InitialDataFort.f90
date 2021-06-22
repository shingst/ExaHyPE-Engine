! GRMHD Initial Data
#define GRMHD



RECURSIVE SUBROUTINE InitialData(x, t, Q)
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar, nDim
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: x(nDim), t        ! 
	REAL, INTENT(OUT)              :: Q(nVar)        ! 

	! Call here one of
	! CALL InitialBlast(x, 0.0,  Q)
  !	Call AlfenWave(x, t, Q)
	! Call InitialRotor(x, 0.0, Q)
	! Call InitialBlast(x, 0.0, Q)
	! Call InitialOrsagTang(x, 0.0 , Q)
	
	! CALL InitialAccretionDisc(x, 0.0,  Q)
	!CALL InitialAccretionDisc3D(x, 0.0, Q)
	CALL InitialDataTN(x, t,  Q)
    
END SUBROUTINE InitialData

RECURSIVE SUBROUTINE InitialDataTN(x, t, Q)
	USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim, d, ICType,  ICType2,ExcisionRadius,NSTOVVar,NSTOVVar_bar,NSTOV_kappa,  aom, Mbh,NSTOV_ATMO,NSTOV_t_atm, gamma
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
    REAL, INTENT(IN ) :: x(nDim), t        ! spatial position vector and time 
    REAL, INTENT(OUT) :: Q(nVar)        ! 

    INTEGER , PARAMETER :: nParam = 1
    REAL :: xGP(d),tGP,u0(nVar)
    ! Local variables 
    INTEGER :: i,j,k,l,nm,iNewton, Maxnewton = 50    
    REAL :: VBase(nVar), ampl(nVar), sigma(d) 
    REAL :: V0(nVar),r,VLL(nVar),VRR(nVar), VZ4(54),V70(70)  
    REAL :: du,dv,drho,dTemp,dp,epsilon,xc(d),rbar
    REAL :: omega, parL(nParam), parR(nParam), igamma,gamma2,gamma1,rho_atm,h,p_atm
    REAL :: r1(9), lambda, lambdat, lambdax, mu, rho, cs, cp, ICA, HH, dxH, Kxx
#ifndef GRMHD
    REAL :: aom,Mbh
#endif
    REAL :: theta, a,  M, zz, r_hz, delta, g_cov(3,3), g_contr(3,3), g_tmp(3,3), Kex(3,3), Aex(3,3) 
    REAL :: traceK, alpha, fa, ggg, fff, k0, sig, AA(3), BB(3,3), BMat(3,3), DD(3,3,3), PP(3)  
    REAL :: rho0, p0, eta, b0, tempaa, tempab, tempac, va2, vax, bv(3), BV_contr(3), vv(3), gm, gp, ddet(3) 
    REAL :: beta(3), betadown(3), dbetadown(3,3), nablabeta(3,3), test(3,3), Gtilde(3)   
    REAL :: Christoffel(3,3,3), dxphi, bbb, dxb, dtb, kyy, kzz, Aref, x0, re0, ms, uu0  
    REAL :: vv_cov(3), shift(3), p,  gammaij(6), ng, up, phi,phi2, dphi(3), rc, vc2, tc, pc, rhoc, lf, vtheta, vphi, vx, vy, vz 
    REAL :: urr, tt, c1, c2, df, vr, ut, f, urc, g_tt, betaru, dtt, detg, psi , xloc(d), minR
    REAL :: BJ(0:1),DJ(0:1),BY(0:1),DY(0:1) 
    REAL :: IJ(0:1),DIJ(0:1),IY(0:1),DIY(0:1) 
    REAL, PARAMETER :: Pi = ACOS(-1.0) 
    REAL :: xGP_sph(3),xGP_loc(3),A_coord(3,3), detA,iA_coord(3,3),A_coord_contr(3,3),iA_coord_contr(3,3)
    !
    IF(nDIM.eq.3) THEN
		xGP=x
	ELSE
		xGP(1:2) = x(1:2)
		xGP(3) = 0.
	ENDIF
	tGP = t

    u0  = 0.0
    !
#ifdef GRMHD
    SELECT CASE(ICType)
    CASE('GRMHDAlfvenWave') 
       rho0 = 1.
       p0   = 1.
       eta  = 1. 
       B0   = 1.0  
       !
       hh = 1.0 + gamma / ( gamma - 1.0) * p0 / rho0
       tempaa = rho0 * hh + B0**2 * ( 1.0 + eta**2)
       tempab = 2.0 * eta * B0**2 / tempaa
       tempac = 0.5 * ( 1.0 + sqrt ( 1.0 - tempab**2))
       va2 = b0**2 / ( tempaa * tempac)
       vax = sqrt ( va2) 
       !
       BV(1) = B0
       BV(2) = eta * B0 * COS( xGP(1) - vax*tGP)
       BV(3) = eta * B0 * SIN( xGP(1) - vax*tGP)
       !
       VV_cov(1)   = 0.0
       VV_cov(2:3) = - vax * BV(2:3) / B0
       !
       ! Now convert to conservative variables 
       !
       !V0 = (/ rho0, VV_cov(1:3), p0, BV(1:3), 0.0 /)  
       gammaij = 0.
		gammaij(1)=1.0
		gammaij(4)=1.0
		gammaij(6)=1.0
       !
        BV(1:3) = 0.
        !
        !V0 = (/ rho, VV_cov(1:3), p ,BV(1:3), 0.0, alpha, shift(1:3), gammaij(1:6)/)
        !
        V0(1) = rho
        V0(5) = p
        V0(9) = 0.
        V0(10) = alpha 
        !
        DO i=1,3
            V0(1+i) = VV_cov(i) 
            V0(5+i) = BV(i) 
            V0(10+i) = shift(i)
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
        CALL Curved2Cartesian(xGP_loc,xGP)          ! no effects if you are in Cartesian coordiantes
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
        igamma = 1.0/gamma 
        rho=(p/NSTOV_kappa)**igamma 
        VV_cov(1:3) = 0.
        BV_contr(1:3) = 0.
        !
        CALL METRIC(  xGP, alpha,  gp, gm, shift, Kex, g_cov, g_contr, phi )
        CALL DMETRIC( xGP, AA, BB, DD, ddet )
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
                gamma1 = gamma / ( gamma - 1.0)
                h      = (1.0 + gamma1*NSTOV_t_atm) * NSTOVVar%lapse_C/alpha
                t      = (h - 1.0) / gamma1
                gamma2 = 1.0 / ( gamma - 1.0)
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
                gamma1 = gamma / ( gamma - 1.0)
                h      = (1.0 + gamma1*NSTOV_t_atm) * NSTOVVar%lapse_C/alpha
                t      = (h - 1.0) / gamma1
                gamma2 = 1.0 / ( gamma - 1.0) 
                
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
            V0(10+i) = shift(i)
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
    CASE('GRMHDTOV_perturbed')
#ifdef RNSTOV
        CALL Curved2Cartesian(xGP_loc,xGP)          ! no effects if you are in Cartesian coordiantes
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
        igamma = 1.0/gamma 
        rho=(p/NSTOV_kappa)**igamma 
        ! then add a density perturbation and recompute the pressure.
        rho = rho*(1 + 0.1* cos(pi/2.0 * r/NSTOVVar_bar%radius))
        p=rho**gamma*NSTOV_kappa
        p = MAX( NSTOVVar%p_R,p ) 
        
        VV_cov(1:3) = 0.
        BV_contr(1:3) = 0.
        !
        CALL METRIC(  xGP, alpha,  gp, gm, shift, Kex, g_cov, g_contr, phi )
        CALL DMETRIC( xGP, AA, BB, DD, ddet )
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
                gamma1 = gamma / ( gamma - 1.0)
                h      = (1.0 + gamma1*NSTOV_t_atm) * NSTOVVar%lapse_C/alpha
                t      = (h - 1.0) / gamma1
                gamma2 = 1.0 / ( gamma - 1.0)
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
                gamma1 = gamma / ( gamma - 1.0)
                h      = (1.0 + gamma1*NSTOV_t_atm) * NSTOVVar%lapse_C/alpha
                t      = (h - 1.0) / gamma1
                gamma2 = 1.0 / ( gamma - 1.0) 
                
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
            V0(10+i) = shift(i)
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
       xloc = xGP 
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
       CALL METRIC ( xloc , alpha, gp, gm, shift, Kex, g_cov, g_contr, phi )
       ng     = 1.0/(gamma - 1.0)

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
        BV(1:3) = MATMUL(iA_coord_contr,V0(6:8))  !  SPHERICAL to Cartesian:  Magnetic field is contravariant  
        !
		
		
        V0(1) = rho
        V0(5) = p
        V0(9) = 0.
        V0(10)  = alpha 
        !
        DO i=1,3
            V0(1+i) = VV_cov(i) 
            V0(5+i) = BV(i) 
            V0(10+i) = shift(i)
        ENDDO
        !
        DO i=1,6
            V0(13+i) = gammaij(i)
        ENDDO
        !  
#endif        
#endif 

    CASE('GRMHDRNSID')
        ! 
        V0 = 0.0
        ! Compute the metric and its derivatives 
        CALL METRIC(  xGP, alpha,  gp, gm, beta, Kex, g_cov, g_contr, phi )
        CALL DMETRIC( xGP, AA, BB, DD, PP )
        test = matmul(g_cov,g_contr)
        detg = g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)*g_cov(3,2)-g_cov(2,1)*g_cov(1,2)*g_cov(3,3)+g_cov(2,1)*g_cov(1,3)*g_cov(3,2)+g_cov(3,1)*g_cov(1,2)*g_cov(2,3)-g_cov(3,1)*g_cov(1,3)*g_cov(2,2)
        !phi = detg**(-1./6.) 
        !
#ifdef RNSID   
        CALL RNSID_Interpolate(xGP, V70) 
        rho         = MAX( 1e-12, V70(60) ) 
        VV_cov(1:3) = V70(61:63) 
        p           = MAX( 1e-12, V70(64) ) 
        BV          = 0.0 
#else
        PRINT *, ' RNSID not available. Please compile with -DRNSID flag. '
        STOP 
#endif      
        
       ! 
       !
       gammaij(1) = g_cov(1,1)
       gammaij(2) = g_cov(1,2)
       gammaij(3) = g_cov(1,3)
       gammaij(4) = g_cov(2,2)
       gammaij(5) = g_cov(2,3)
       gammaij(6) = g_cov(3,3) 
	   
	   
        V0(1) = rho
        V0(5) = p
        V0(9) = 0.
        V0(10)  = alpha 
        !
        DO i=1,3
            V0(1+i) = VV_cov(i) 
            V0(5+i) = BV(i) 
            V0(10+i) = beta(i)
        ENDDO
        !
        DO i=1,6
            V0(13+i) = gammaij(i)
        ENDDO
        !
       CONTINUE
       !

    CASE('GRMHDAlfven')
       rho = 1.
       p   = 1.
       eta  = 1. 
       B0   = 1.0  
       !
       hh = 1.0 + gamma / ( gamma - 1.0) * p / rho
       tempaa = rho * hh + B0**2 * ( 1.0 + eta**2)
       tempab = 2.0 * eta * B0**2 / tempaa
       tempac = 0.5 * ( 1.0 + sqrt ( 1.0 - tempab**2))
       va2 = b0**2 / ( tempaa * tempac)
       vax = sqrt ( va2)       
       !
       ! flat metric: contr and cov are the same in Cartesian
       BV_contr(1) = B0
       BV_contr(2) = eta * B0 * COS( xGP(1) - vax*tGP)
       BV_contr(3) = eta * B0 * SIN( xGP(1) - vax*tGP)
       !
       VV_cov(1)   = 0.0
       VV_cov(2:3) = - vax * BV_contr(2:3) / B0
       !
       ! Now convert to conservative variables 
       !
       alpha = 1.0
       shift(1:3) = 0.
       gammaij = 0.
       gammaij(1) = 1.
       gammaij(4) = 1.
       gammaij(6) = 1. 
       !        
        V0(1) = rho
        V0(5) = p
        V0(9) = 0.
        V0(10)  = alpha 
        !
        DO i=1,3
            V0(1+i) = VV_cov(i) 
            V0(5+i) = BV_contr(i) 
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
    CASE('GRMHDAccretion','GRMHD-SphericalAccretion')       
       !aom = 0.0
       !Mbh = 1.0  
       rc   = 8.0 
       rhoc = 1./16.
       B0 = 4. 
       xloc = xGP
       IF ( aom > 0.0) THEN 
           WRITE(*,*)'SPHERICAL Accretion solution is only for a = 0!'
           STOP
       ENDIF 
       !
       ! The Following is for Kerr-Schild Cartesian coordinates
       r      = SQRT( xGP(1)**2 + xGP(2)**2 + xGP(3)**2) 
       theta  = ACOS( xGP(3)/r)
       phi2    = ATAN2( xGP(2), xGP(1))
       !
       !IF ( r .LT. 0.0001) THEN ! RETURN
       !    alpha = 1.0
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
       !    V0 = (/ rho, VV_cov(1:3), p ,BV(1:3), 0.0, alpha, shift(1:3), gammaij(1:6)/)
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
       CALL METRIC ( xloc , alpha, gp, gm, shift, Kex, g_cov, g_contr, phi )
       ng     = 1.0/(gamma - 1.0)  
       !
       !
       !
       !phi2    = ACOS( xGP(1) / (r*SIN(theta)))
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
       LF     = alpha*ut
       vr     = ( urr / LF + betaru / alpha)
       vtheta = 0.0
       vphi   = 0.0
       !
       vx  = DSIN(theta)*DCOS(phi2)*vr
       vy  = DSIN(theta)*DSIN(phi2)*vr
       vz  = DCOS(theta)*vr
       !
        VV(1) = vx
        VV(2) = vy
        VV(3) = vz

       ! Convert to covariant velocities
       VV_cov = MATMUL(g_cov, VV)
       !
       rho = rhoc*(tt/tc)**ng
       p   = rho*tt
       !
       BV(1:3) = 0.
       BV_contr(1:3) = 0.
       ! 
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
           BV_contr(1) = alpha*B0*DCOS(xGP(2))
           BV_contr(2) = -  alpha*B0*DSIN(xGP(2))/xGP(1) 
           !
           !BV(1:3) = MATMUL(g_cov,BV(1:3))
           ! 
           !
           CASE DEFAULT
           !
           continue
           !
       END SELECT
       !  
       gammaij(1) = g_cov(1,1)
       gammaij(2) = g_cov(1,2)
       gammaij(3) = g_cov(1,3)
       gammaij(4) = g_cov(2,2)
       gammaij(5) = g_cov(2,3)
       gammaij(6) = g_cov(3,3) 
	   
	   !
        V0(1) = rho
        V0(5) = p
        V0(9) = 0.
        V0(10)  = alpha 
        !
        DO i=1,3
            V0(1+i) = VV_cov(i) 
            V0(5+i) = BV_contr(i) 
            V0(10+i) = shift(i)
        ENDDO
        !
        DO i=1,6
            V0(13+i) = gammaij(i)
        ENDDO
        !
       CONTINUE
       !
       CONTINUE
       !
    END SELECT 
#endif
    !
    !
    CALL PDEPrim2Cons(Q,V0) 
    !
    CONTINUE
    !
END SUBROUTINE InitialDataTN


RECURSIVE SUBROUTINE AlfenWave(x, t, Q)
    ! Computes the AlfenWave conserved variables (Q) at a given time t.
    ! Use it ie. with t=0 for initial data
    ! Use it for any other time ie. for comparison
    
    ! GRID FOR ALFENWAVE:
    !     dimension const                = 2
    !     width                          = 1.0, 0.3
    !     offset                         = 0.0, 0.0
    !     end-time                       = 2.1
    !
    !  maximum-mesh-size              = 0.04

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim, gamma
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 
    REAL, PARAMETER                :: t_offset = 1.0 ! time offset

    REAL :: rho0, p0, eta, B0, hh, tempaa, tempab, tempac, va2, vax
    REAL :: V(nVar), BV(3), VV(3), Pi = ACOS(-1.0)

    rho0 = 1.
    p0   = 1.
    eta  = 1.
    B0   = 1.0 
    !
    hh = 1.0 + gamma / ( gamma - 1.0) * p0 / rho0
    tempaa = rho0 * hh + B0**2 * ( 1.0 + eta**2)
    tempab = 2.0 * eta * B0**2 / tempaa
    tempac = 0.5 * ( 1.0 + sqrt ( 1.0 - tempab**2))
    va2 = b0**2 / ( tempaa * tempac)
    vax = sqrt ( va2)
    !
    BV(1) = B0
!    BV(2) = eta * B0 * COS(2*Pi*( x(1) - vax*(t-t_offset)))
!    BV(3) = eta * B0 * SIN(2*Pi*( x(1) - vax*(t-t_offset)))
    BV(2) = eta * B0 * COS(x(1) - vax*t)
    BV(3) = eta * B0 * SIN(x(1) - vax*t)
    !
    VV(1)   = 0.0
    VV(2:3) = - vax * BV(2:3) / B0
    !
    ! Now convert to conservative variables
    !
    V(1:9) = (/ rho0, VV(1:3), p0, BV(1:3), 0.0 /) ! psi
    !            lapse, shift(3),    gamma11, ...  gamma22, -, gamma33
    V(10:19) = (/ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0 /)
    
    ! create primitive ID
    !Q = V
    CALL PDEPrim2Cons(Q,V)
END SUBROUTINE AlfenWave

RECURSIVE SUBROUTINE InitialBlast(x, t, Q)
    ! Blast wave initial data in conserved variables (Q):
    ! Simulation domain:  -6 .. +6

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 
    
    REAL :: rho0, p0, p1, rho1, r, r0, r1, taper
    REAL :: V(nVAR)


    p0   = 5.0e-4
    p1   = 1.0
    rho0 = 1.0e-4 
    rho1 = 1.0e-2
    ! Here, we set the primitive variables 
    V(1) = 0.0   ! rho
    V(2) = 0.0   ! vx
    V(3) = 0.0   ! vy
    V(4) = 0.0   ! vz
    V(5) = 0.0   ! p
    V(6) = 0.1   ! bx
    V(7) = 0.0   
    V(8) = 0.0   
    V(9) = 0.0   ! preserving
    r = SQRT(x(1)**2 + x(2)**2)
    r0 = 0.8
    r1 = 1.0
    !
    IF(r .LE. 0.8) THEN
        V(1) = rho1
        V(5) = p1  
    ELSEIF ( r.GT.0.8 .AND. r.LE.1.0) THEN
        taper = (1.-(r-r0)/(r1-r0))
        V(1) = rho0 + taper*(rho1-rho0)
        taper = (1.-(r-r0)/(r1-r0))
        V(5) = p0 + taper*(p1-p0)
    ELSE  
        V(1) = rho0
        V(5) = p0 
    ENDIF
    V(10:19) = (/ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0 /)

    CALL PDEPrim2Cons(Q,V)
END SUBROUTINE

RECURSIVE SUBROUTINE InitialOrsagTang(x, t, Q)
    ! Orsang-Tang initial data in conserved variables (Q)
    ! Simulation Domain: 0 .. 2*pi = 6.283185307179586

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 

    REAL :: rho0, p0, vel0, B0
    REAL :: V(nVAR), VV(3), BV(3)

    ! RMHDOrszagTang
    rho0 = 1.0
    p0   = 1.0
    vel0 = 0.75
    B0   = 1.0
    !
    VV(1) = -1./SQRT(2.0)*vel0*SIN(x(2))
    VV(2) =  1./SQRT(2.0)*vel0*SIN(x(1))
    VV(3) = 0. 
    !
    BV(1) = -B0*SIN(x(2))
    BV(2) =  B0*SIN(2.0*x(1))
    BV(3) =  0.
    !
    V(1:9) = (/ rho0, VV(1:3), p0, BV(1:3), 0.0 /)
    V(10:19) = (/ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0 /)

    Call PDEPrim2Cons(Q,V)
END SUBROUTINE InitialOrsagTang

RECURSIVE SUBROUTINE InitialRotor(x,t,Q)
    ! GRMHD Rotor initial data in conserved variables (Q)
    ! Needs Limiting, otherwise crash.
    ! Domain: 0.0 .. 1.0 square

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 
    
    REAL ::  rho, p, r
    REAL :: rho0, p0, rvel0, rho1, p1, B0, v0
    REAL :: V(nVAR), VV(3), BV(3)
    
    REAL, PARAMETER :: MHDRotomega = 0.95
    REAL :: EPCenter(2), EPRadius
    
    EPCenter = (/ 0.5, 0.5 /)
    EPRadius = 0.1
    
    rho0 = 1.0
    p0   = 1.0
    rho1 = 10.0
    p1   = 1.0
    B0   = 1.0
    v0   = MHDRotomega  ! this is omega
    r = sqrt ( (x(1)-EPCenter(1))**2 + (x(2)-EPCenter(2))**2)   ! cylindrical
    IF  (r < EPRadius) THEN
        rho   =  rho1                   
        VV(1) = -v0 * x(2)
        VV(2) =  v0 * x(1) 
        VV(3) =  0.
        p     =  p1
    ELSE
        rho   = rho0
        VV    = 0.
        p     = p0
    ENDIF
    !
    BV(1) = B0
    BV(2) = 0.
    BV(3) = 0.
    !
    V(1:9) = (/ rho, VV(1:3), p, BV(1:3), 0.0 /)
    V(10:19) = (/ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0 /)
    CALL PDEPrim2Cons(Q, V)
END SUBROUTINE InitialRotor



RECURSIVE SUBROUTINE InitialAccretionDisc(x,t,Q)
    ! 2D Accretion disk, with simulation domain:
    !    dimension const                = 2
    !    width                          = 8.5, 2.0
    !    offset                         = 1.5, 0.5
    !    end-time                       = 2.1
    ! mesh:     maximum-mesh-size              = 0.1
    ! Timestep size is around repeat = 0.000785674, for plotter.

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim, gamma
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 
    ! Local variables
    
    REAL :: rho0, p0, eta, B0, hh, tempaa, tempab, tempac, va2, vax
    REAL :: V(nVar), BV(3), VV(3), Pi = ACOS(-1.0)
    REAL :: r, zz, urc, vc2, tc, pc,tt, c1, c2, urr, f
    REAL :: df, dt, ut, LF, vr, vtheta, vphi, rho, p, VV_cov(3), Kex(3,3), g_cov(3,3), g_contr(3,3)
    REAL :: phi_c
    REAL :: gp, gm, shift(3), lapse
    
    ! PARAMETERS:
    REAL :: rhoc = 0.0625  ! Critical radius
    REAL :: rc = 8.0
    INTEGER :: MAXNEWTON = 50, iNewton
    REAL :: ng = 1.0 / (gamma-1.0)

    CALL METRIC ( x, lapse, gp, gm, shift, Kex, g_cov, g_contr,phi_c)

    
     ! The Following is for Kerr-Schild spherical coordinates
       r      = x(1)
       !
       zz     = 2.0/r               ! we are computing the solution at theta=pi/2
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
          dt  = -f/df
          IF (abs(dt) < 1.e-10) EXIT
          tt = tt + dt
       ENDDO
       ut     = (-zz*urr + sqrt(urr**2 - zz + 1.0))/(zz - 1.0)
       LF     = lapse*ut
       vr     = ( urr / LF + shift(1) / lapse)
       vtheta = 0.0
       vphi   = 0.0
       !
       VV(1:3) = (/ vr, vtheta, vphi /)
       ! Convert to covariant velocities
       VV_cov = MATMUL(g_cov, VV)
       !
       rho = rhoc*(tt/tc)**ng
       p   = rho*tt      

       V(1:9) = (/ rho, VV_cov(1:3), p, 0., 0., 0., 0. /)
       V(10:19) = (/ 1., 0., 0., 0., 1., 0., 0., 1., 0., 1. /)
       CALL PDEPrim2Cons(Q,V)
END SUBROUTINE InitialAccretionDisc


RECURSIVE SUBROUTINE InitialAccretionDisc3D(x,t,Q)
    ! 3D Accretion disk, with simulation domain
    !     width                          = 2.0, 2.0, 2.0
    !     offset                         = 0.0, 0.0, 0.0
    ! 

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim, gamma
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 
    ! Local variables
    
    REAL :: rho0, p0, eta, B0, hh, tempaa, tempab, tempac, va2, vax
    REAL :: V(nVar), BV(3), BV_contr(3), VV(3), Pi = ACOS(-1.0)
    REAL :: r, zz, urc, vc2, tc, pc,tt, c1, c2, urr, f
    REAL :: df, dt, ut, LF, vr, vtheta, vphi, rho, p, VV_cov(3), Kex(3,3),g_cov(3,3), g_contr(3,3)
    REAL :: gp, gm, shift(3), lapse, gammaij(6), betaru, g_tt, phi, theta, vx, vy, vz
    REAL :: phi_c
    REAL :: detgamma    
    ! PARAMETERS:
  !  REAL :: rhoc = 0.625  ! Critical radius
    REAL :: rhoc = 0.0625  ! Critical radius
    
    REAL :: rc = 8.0
    INTEGER :: MAXNEWTON = 50, iNewton
    INTEGER :: i 
    REAL :: ng = 1.0 / (gamma-1.0)

    
     CALL METRIC ( x, lapse, gp, gm, shift, Kex , g_cov, g_contr,phi_c)

     ng     = 1.0/(gamma - 1.0)

       ! The Following is for Kerr-Schild Cartesian coordinates
       r      = SQRT( x(1)**2 + x(2)**2 + x(3)**2)
       !
       !IF ( x(1) .LT. 0.1) RETURN
       !IF ( x(2) .LT. 0.1) RETURN
       !IF ( x(3) .LT. 0.1) RETURN
       IF ( r .LT. 0.8) THEN
          ! To avoid division by zero, never used for evolution or BC
          rho = 1.0 !1.0
          VV_cov(1:3) = 0.0
          p = 1.0
          BV = 0.0 
          V(1:9) = (/ rho, VV_cov(1:3), p, BV(1:3), 0. /)
          V(10:19) = (/ 1.0, 0.0,0.0,0.0, 1.0,0.0,0.0,1.0,0.0,1.0 /)
          CALL PDEPrim2Cons(Q,V)       
          RETURN
       ENDIF
       !
       theta  = ACOS( x(3)/r)
       phi    = ATAN2( x(2), x(1))
       !phi    = ACOS( x(1) / (r*SIN(theta)))
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
          dt  = -f/df
          IF (abs(dt) < 1.e-10) EXIT
          tt = tt + dt
       ENDDO
       ut     = (-zz*urr + sqrt(urr**2 - zz + 1.0))/(zz - 1.0)
       LF     = lapse*ut
       vr     = ( urr / LF + betaru / lapse)
       vtheta = 0.0
       vphi   = 0.0
       !
       vx  = SIN(theta)*COS(phi)*vr
       vy  = SIN(theta)*SIN(phi)*vr
       vz  = COS(theta)*         vr
       !
!       VV(1:3) = (/ vx, vy, vz /)
       VV(1) = vx
       VV(2) = vy
       VV(3) = vz
       ! Convert to covariant velocities
       VV_cov = MATMUL(g_cov, VV)
       !
       rho = rhoc*(tt/tc)**ng
       p   = rho*tt
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

       !If hydro (Clasical) Michel accretion                                                    

       !BV(1:3) = 0.                                                                            

       ! MHD  mICHEl accretion                                                                  

!       detgamma= gammaij(1)*( gammaij(4)*gammaij(6)-gammaij(5)*gammaij(5)) &
!               - gammaij(2)*( gammaij(2)*gammaij(6)-gammaij(5)*gammaij(3)) &
!               + gammaij(3)*( gammaij(2)*gammaij(5)-gammaij(3)*gammaij(4))
!       detgamma = sqrt(detgamma)




       ! B0 = (2.2688/M)*(sqrt(b^2/rho0)), here sqrt(b^2/rho0)=4.0                              

       B0 = 0

       !---------------------------------------------------------------------- 
       BV_contr = B0*x(1:3)/(sqrt(gp) * r*r*r)

!       BV = BV_contr
       
       BV = MATMUL(g_cov,BV_contr)

        !

    V(1) = rho
    V(2) = VV_cov(1)
    V(3) = VV_cov(2)
    V(4) = VV_cov(3)
    V(5) = p  
    V(6) = BV(1)
    V(7) = BV(2)
    V(8) = BV(3)
    V(9) = 0
    V(10) = lapse
    V(11) = shift(1)
    V(12) = shift(2)
    V(13) = shift(3)
        DO i=1,6
            V(13+i) = gammaij(i)
        ENDDO
   
  ! DO i=1,3
   !     V(1+i) = VV_cov(i)
   !     V(5+i) = BV(i)
   !     V(10+i) = shift(i)
   ! ENDDO
    !
    !
       !V(1:9) = (/ rho, VV_cov(1:3), p, BV(1:3), 0. /)
       !V(10:19) = (/ lapse, shift(1:3), gammaij(1:6) /)
       CALL PDEPrim2Cons(Q,V)


 END SUBROUTINE InitialAccretionDisc3D

 
 
 
 
    
! SUBROUTINE PDESetup
!     USE Parameters 
! #ifdef TWOPUNCTURES 
! 	USE TwoPunctures_mod, ONLY : TwoPunctures_Setup, TwoPunctures_Run, TwoPunctures_AmendParameters
! #endif         
! #ifdef RNSID
!     USE RNSID_mod, ONLY : RNSID_Setup, RNSID_Run, RNSID_AmendParameters 
! #endif          
! #ifdef RNSTOV
!     USE NSTOV_mod
! #endif  
!     IMPLICIT NONE
! #ifdef PARALLEL
!     INCLUDE 'mpif.h'
! #endif    
!     ! Local variables 
!     INTEGER          :: i, j, k, ii, jj, kk, l, c, iGP, iElem, iVar, VMAX(d), count, cnt, iErr, iDim, idxs(d), iper(d) 
!     INTEGER          :: Left, Right, cInnerFaces, cOuterFaces, itemp(N+1), CPUN(-1:1,-1:1,-1:1)    
!     REAL             :: phi(N+1), phi0(N+1), phi1(N+1), phi_xi(N+1), phi_xixi(N+1)
!     REAL             :: phi_i(N+1), phi_j(N+1), phi_k(N+1) 
!     REAL             :: ux(nVar,N+1,N+1,N+1), uy(nVar,N+1,N+1,N+1), uz(nVar,N+1,N+1,N+1) 
!     REAL             :: u0(nVar), v0(nVar), xGP(d), x0(d), subxi(N+1), aux(d), nv(d), par0(nParam)  
!     REAL             :: lparambnd(nParam,6,nDOF(2),nDOF(3))    
!     REAL             :: xi, dxi, xi1, xi2, rtemp(N+1), ImLambda(N+1)  
!     REAL             :: TestMatrix(N+1,N+1), TestMatrix2(nSubLim,nSubLim), test(nSubLim)
!     REAL             :: r1(9), wGPNMat(N+1,N+1) 
!     REAL, POINTER    :: LSQM(:,:), iLSQM(:,:), LSQrhs(:,:) 
!     LOGICAL          :: dmpresult 
!     REAL, PARAMETER  :: Pi = ACOS(-1.0) 
!     ! 
!     SELECT CASE(ICType)
!     CASE('CCZ4TwoPunctures','Z4TwoPunctures') 
! #ifdef TWOPUNCTURES  
!         CALL TwoPunctures_Setup()	
! 	    ! After Setup, you can change the parameters of TwoPunctures, such as Black Hole
! 	    ! mass, etc. They are stored as C++ class members, so there is no way of changing
! 	    ! them from Fortran. However, you can do so if you call the following function
! 	    ! and changes it's implementation
! 	    !CALL TwoPunctures_AmendParameters()	
! 	    PRINT *, ' Starting the TwoPunctures algorithm... '
! 	    CALL TwoPunctures_Run() 
! 	    PRINT *, ' Done. '
! #ifdef PARALLEL
!         CALL MPI_BARRIER(MPI_COMM_WORLD,mpiErr) 
! #endif
! #else
!         PRINT *, ' Twopunctures not available. Please compile with -DTWOPUNCTURES. ' 
!         STOP 
! #endif         
!     CASE('CCZ4RNSID','GRMHDRNSID') 
! #ifdef RNSID 
!         CALL RNSID_Setup()	 
!         CALL RNSID_AmendParameters()  
!         PRINT *, ' Starting the RNSID algorithm... '
!         CALL RNSID_Run() 
! 	    PRINT *, ' Done. '        
! #ifdef PARALLEL
!         CALL MPI_BARRIER(MPI_COMM_WORLD,mpiErr) 
! #endif
!         ! 
! #else
!         !PRINT *, ' RNSID not available. Please compile with -DRNSID. ' 
!         STOP 
! #endif 
!     CASE('GRMHDTorus')
! #ifdef TORUS
!         CALL init_disk_isentropic_par(METRIC_KSS)
! #endif
!     ! 
!     CASE('GRMHDTOV','CCZ4TOV','GRMHDTOV_perturbed') 
! #ifdef RNSTOV
!         !IF(myrank.eq.0) PRINT *, ' Starting the RNSTOV algorithm... '
!             CALL NSTOV_Main
! 	    !IF(myrank.eq.0) PRINT *, ' Done. '  
! #ifdef PARALLEL 
!         CALL MPI_BARRIER(MPI_COMM_WORLD,mpiErr)  
! #endif     
! #else
!         PRINT *, ' RNSTOV not available. Please compile with -DRNSTOV. ' 
!         STOP 
! #endif
!     ! 
!     END SELECT 
!     ! 
!     !
!     continue
!     !
! END SUBROUTINE PDESetUp
