! Kerr, Twopunctures, etc. from METRIC and DMETRIC derived ID

RECURSIVE SUBROUTINE InitialField(u0,xGP,tGP) 
    USE typesDef, ONLY : nVar, d, ICType, IC, EQN  
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN ) :: xGP(d), tGP        ! spatial position vector and time 
    REAL, INTENT(OUT) :: u0(nVar)           ! initial data vector in terms of conserved variables 
    ! Local variables 
    INTEGER :: i,j,k,l, iNewton, Maxnewton = 50    
    REAL :: VBase(nVar), ampl(nVar), sigma(d) 
    REAL :: V0(nVar),r,VLL(nVar),VRR(nVar), VZ4(54) 
    REAL :: du,dv,drho,dTemp,dp,epsilon,xc(d) 
    REAL :: omega
    REAL :: r1(9), lambda, mu, rho, cs, cp, ICA, HH, dxH, Kxx
    REAL :: theta, a, aom, M, Mbh, zz, r_hz, delta, g_cov(3,3), g_contr(3,3), Kex(3,3), Aex(3,3) 
    REAL :: traceK, alpha, fa, ggg, fff, k0, sig, AA(3), BB(3,3), BMat(3,3), DD(3,3,3), PP(3)  
    REAL :: rho0, p0, eta, b0, tempaa, tempab, tempac, va2, vax, bv(3), vv(3), gm, gp, ddet(3) 
    REAL :: beta(3), betadown(3), dbetadown(3,3), nablabeta(3,3), test(3,3), Gtilde(3)   
    REAL :: Christoffel(3,3,3), dxphi   
    REAL :: vv_cov(3), shift(3), p, lapse, gammaij(6), ng, up, phi, dphi(3), rc, vc2, tc, pc, rhoc, lf, vtheta, vphi, vx, vy, vz 
    REAL :: urr, tt, c1, c2, df, vr, ut, f, urc, g_tt, betaru, dtt, detg, psi     
    REAL, PARAMETER :: Pi = ACOS(-1.0) 
    !
    u0  = 0.0
    V0 = 0.0
    IF(IC.EQ.ICType%CCZ4GaugeWave) THEN
        V0 = 0.0
        !
        ICA = 0.1 
        ! Gauge wave propagating along x 
        HH     = 1.0 - ICA*SIN( 2.0*Pi*( xGP(1) - tGP)) 
        dxH    = -2.0*Pi*ICA*COS( 2.0 * Pi*(xGP(1) - tGP)) 
        dxphi  = - HH**(-7.0/6.0)*dxH/6.0
        phi    = ( 1.0 / HH)**(1.0/6.0)
        Kxx    = - Pi*ICA*COS( 2.0 * Pi*(xGP(1) - tGP))/SQRT( 1.0 - ICA*SIN( 2.0*Pi*( xGP(1) - tGP))  )  ! extrinsic curvature  
        traceK = Kxx/HH
        !
!        HH = 1.0 
!        phi = 1.0 
!        traceK = 0.0
!        Kxx = 0.0 
!        traceK = 0.0 
!        dxH = 0.0
!        dxphi = 0.0 
        !
        V0(:)  = 0.0
        V0(1)  = phi**2*HH                          ! \tilde\gamma_xx
        V0(4)  = phi**2                             ! \tilde\gamma_yy
        V0(6)  = phi**2                             ! \tilde\gamma_zz
        V0(7)  = phi**2*(Kxx - 1.0/3.0*traceK*HH )  ! \tilde A_{xx}
        V0(10) = phi**2*(0.0 - 1.0/3.0*traceK*1.0)  ! \tilde A_{yy}
        V0(12) = phi**2*(0.0 - 1.0/3.0*traceK*1.0)  ! \tilde A_[zz}
        V0(17) = SQRT(HH)
        V0(14) = 2.0/(3.0*HH**(5.0/3.0))*dxH        ! Gtilde 
        !
        ! Auxiliary variables
        V0(24) = 1.0/(2.0*HH)*dxH                ! A_x  
        V0(36) = HH**(-1.0/3.0)*dxH/3.0          ! D_xxx
        V0(39) = phi*dxphi                       ! D_xyy
        V0(41) = phi*dxphi                       ! D_xzz
        ! 
        V0(54) = traceK
        V0(55) = phi
        V0(56) = dxphi/phi                       ! P_x
        
    ELSE IF(IC.EQ.ICType%CCZ4Puncture) THEN
        V0 = 0.0
        !
        ! Compute the metric and its derivatives 
        CALL METRIC(  xGP, alpha,  gp, gm, beta, Kex, g_cov, g_contr, phi )
        CALL DMETRIC( xGP, AA, BB, DD, PP )
        test = matmul(g_cov,g_contr)
        detg = g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)*g_cov(3,2)-g_cov(2,1)*g_cov(1,2)*g_cov(3,3)+g_cov(2,1)*g_cov(1,3)*g_cov(3,2)+g_cov(3,1)*g_cov(1,2)*g_cov(2,3)-g_cov(3,1)*g_cov(1,3)*g_cov(2,2)
        !phi = detg**(-1./6.) 
        AA = AA/alpha   ! convert pure derivatives to aux. variables 
        DD = 0.5*DD     ! convert pure derivatives to aux. variables 
        PP = PP/phi     ! convert pure derivatives to aux. variables  
        ! Christoffel symbols 
        Christoffel = 0.0
        DO k = 1, 3
         DO j = 1, 3
          DO i = 1, 3
              DO l = 1, 3 
                Christoffel(i,j,k) = Christoffel(i,j,k) + g_contr(k,l)*( DD(i,j,l) + DD(j,i,l) - DD(l,i,j) ) 
              ENDDO 
          ENDDO
         ENDDO
        ENDDO 
        ! Extrinsic curvature 
        Kex = 0.0
        ! Trace of K 
        traceK = SUM( Kex*g_contr ) 
        ! The conformal traceless part of the extrinsic curvature
        Aex = phi**2*( Kex - 1./3*traceK*g_cov ) 
        ! The contracted connection coefficients 
        Gtilde = 0.0
        DO i = 1, 3
         DO j = 1, 3
          DO k = 1, 3
           DO l = 1, 3
               Gtilde(i) = Gtilde(i) + 1./phi**2*( g_contr(i,j)*g_contr(k,l)*( 2*DD(l,j,k) + 2*PP(l)*g_cov(j,k) )  ) 
           ENDDO
          ENDDO 
         ENDDO
        ENDDO 
        SELECT CASE(EQN%CCZ4LapseType) 
        CASE(0)  ! harmonic 
           fa = 1.0 
        CASE DEFAULT  ! 1 + log 
           fa = 2.0/alpha
        END SELECT 
        ! K0 to make the PDE for alpha stationary 
        K0 = 0.0 
        !
        ! The metric tensor 
        V0(1) = phi**2*g_cov(1,1)    ! \gamma_11 
        V0(2) = phi**2*g_cov(1,2)    ! \gamma_12
        V0(3) = phi**2*g_cov(1,3)    ! \gamma_13
        V0(4) = phi**2*g_cov(2,2)    ! \gamma_22
        V0(5) = phi**2*g_cov(2,3)    ! \gamma_23
        V0(6) = phi**2*g_cov(3,3)    ! \gamma_33 
        ! The extrinsic curvature         
        V0(7)  = Aex(1,1)  
        V0(8)  = Aex(1,2) 
        V0(9)  = Aex(1,3) 
        V0(10) = Aex(2,2) 
        V0(11) = Aex(2,3) 
        V0(12) = Aex(3,3) 
        ! The cleaning variables Z and Theta 
        V0(13) = 0.0          ! Theta 
        V0(14) = Gtilde(1)    ! G1 
        V0(15) = Gtilde(2)    ! G2 
        V0(16) = Gtilde(3)    ! G3  
        ! The logarithm of the lapse 
        V0(17) = alpha 
        ! The shift 
        V0(18) = beta(1)  
        V0(19) = beta(2) 
        V0(20) = beta(3) 
        !
        ! Auxiliary variables
        ! 
        !! vector A_i = \partial_i \alpha / \alpha 
        !!  
        V0(24) = AA(1)     ! A1 = alpha_1/alpha  
        V0(25) = AA(2)     ! A2 = alpha_2/alpha    
        V0(26) = AA(3)     ! A3 = alpha_3/alpha   
        !
        ! Matrix B_ik = \partial_i \beta_k 
        !
        V0(27) = BB(1,1) 
        V0(28) = BB(2,1)
        V0(29) = BB(3,1)
        V0(30) = BB(1,2)
        V0(31) = BB(2,2)
        V0(32) = BB(3,2)
        V0(33) = BB(1,3)
        V0(34) = BB(2,3)
        V0(35) = BB(3,3) 
        !
        ! tensor D_ijk = 0.5 \partial i \tilde \gamma_jk   
        !
        DO j = 1, 3
         DO i = 1, 3 
          DO k = 1, 3 
            DD(k,i,j) = phi**2*( DD(k,i,j) + PP(k)*g_cov(i,j) ) 
          ENDDO
         ENDDO
        ENDDO 
        !
        V0(36) = DD(1,1,1) 
        V0(37) = DD(1,1,2) 
        V0(38) = DD(1,1,3) 
        V0(39) = DD(1,2,2) 
        V0(40) = DD(1,2,3) 
        V0(41) = DD(1,3,3) 
        V0(42) = DD(2,1,1) 
        V0(43) = DD(2,1,2) 
        V0(44) = DD(2,1,3) 
        V0(45) = DD(2,2,2) 
        V0(46) = DD(2,2,3) 
        V0(47) = DD(2,3,3) 
        V0(48) = DD(3,1,1) 
        V0(49) = DD(3,1,2) 
        V0(50) = DD(3,1,3) 
        V0(51) = DD(3,2,2) 
        V0(52) = DD(3,2,3) 
        V0(53) = DD(3,3,3) 
        ! trace of K 
        V0(54) = traceK 
        ! logarithm of the conformal factor phi 
        V0(55) = phi 
        ! derivative of phi 
        V0(56) = PP(1) 
        V0(57) = PP(2) 
        V0(58) = PP(3) 
        ! The trace of the extrinsic curvature at the initial time 
        V0(59) = K0 
        ! 
        
    ELSE IF(IC.EQ.ICType%CCZ4TwoPunctures) THEN
        !
        V0 = 0.0
        ! Compute the metric and its derivatives 
        CALL METRIC(  xGP, alpha,  gp, gm, beta, Kex, g_cov, g_contr, phi )
        CALL DMETRIC( xGP, AA, BB, DD, PP )
        test = matmul(g_cov,g_contr)
        detg = g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)*g_cov(3,2)-g_cov(2,1)*g_cov(1,2)*g_cov(3,3)+g_cov(2,1)*g_cov(1,3)*g_cov(3,2)+g_cov(3,1)*g_cov(1,2)*g_cov(2,3)-g_cov(3,1)*g_cov(1,3)*g_cov(2,2)
        !phi = detg**(-1./6.) 
        AA = AA/alpha   ! convert pure derivatives to aux. variables 
        DD = 0.5*DD     ! convert pure derivatives to aux. variables 
        PP = PP/phi     ! convert pure derivatives to aux. variables  
        ! Trace of K 
        traceK = SUM( Kex*g_contr ) 
        ! The conformal traceless part of the extrinsic curvature 
        Aex = phi**2*( Kex - 1./3*traceK*g_cov ) 
        ! The contracted connection coefficients 
        Gtilde = 0.0
        DO i = 1, 3
         DO j = 1, 3
          DO k = 1, 3
           DO l = 1, 3
               Gtilde(i) = Gtilde(i) + 1./phi**2*( g_contr(i,j)*g_contr(k,l)*( 2*DD(l,j,k) + 2*PP(l)*g_cov(j,k) )  ) 
           ENDDO
          ENDDO 
         ENDDO
        ENDDO 
        !
        ! The metric tensor 
        V0(1) = phi**2*g_cov(1,1)    ! \gamma_11 
        V0(2) = phi**2*g_cov(1,2)    ! \gamma_12
        V0(3) = phi**2*g_cov(1,3)    ! \gamma_13
        V0(4) = phi**2*g_cov(2,2)    ! \gamma_22
        V0(5) = phi**2*g_cov(2,3)    ! \gamma_23
        V0(6) = phi**2*g_cov(3,3)    ! \gamma_33 
        ! The extrinsic curvature         
        V0(7)  = Aex(1,1)  
        V0(8)  = Aex(1,2) 
        V0(9)  = Aex(1,3) 
        V0(10) = Aex(2,2) 
        V0(11) = Aex(2,3) 
        V0(12) = Aex(3,3) 
        ! The cleaning variables Z and Theta 
        V0(13) = 0.0          ! Theta 
        V0(14) = Gtilde(1)    ! G1 
        V0(15) = Gtilde(2)    ! G2 
        V0(16) = Gtilde(3)    ! G3  
        ! The logarithm of the lapse 
        V0(17) = alpha 
        ! The shift 
        V0(18) = beta(1)  
        V0(19) = beta(2) 
        V0(20) = beta(3) 
        !
        ! Auxiliary variables
        ! 
        !! vector A_i = \partial_i \alpha / \alpha 
        !!  
        V0(24) = AA(1)     ! A1 = alpha_1/alpha  
        V0(25) = AA(2)     ! A2 = alpha_2/alpha    
        V0(26) = AA(3)     ! A3 = alpha_3/alpha   
        !
        ! Matrix B_ik = \partial_i \beta_k 
        !
        V0(27) = BB(1,1) 
        V0(28) = BB(2,1)
        V0(29) = BB(3,1)
        V0(30) = BB(1,2)
        V0(31) = BB(2,2)
        V0(32) = BB(3,2)
        V0(33) = BB(1,3)
        V0(34) = BB(2,3)
        V0(35) = BB(3,3) 
        !
        ! tensor D_ijk = 0.5 \partial i \tilde \gamma_jk   
        !
        DO j = 1, 3
         DO i = 1, 3 
          DO k = 1, 3 
            DD(k,i,j) = phi**2*( DD(k,i,j) + PP(k)*g_cov(i,j) ) 
          ENDDO
         ENDDO
        ENDDO 
        !
        V0(36) = DD(1,1,1) 
        V0(37) = DD(1,1,2) 
        V0(38) = DD(1,1,3) 
        V0(39) = DD(1,2,2) 
        V0(40) = DD(1,2,3) 
        V0(41) = DD(1,3,3) 
        V0(42) = DD(2,1,1) 
        V0(43) = DD(2,1,2) 
        V0(44) = DD(2,1,3) 
        V0(45) = DD(2,2,2) 
        V0(46) = DD(2,2,3) 
        V0(47) = DD(2,3,3) 
        V0(48) = DD(3,1,1) 
        V0(49) = DD(3,1,2) 
        V0(50) = DD(3,1,3) 
        V0(51) = DD(3,2,2) 
        V0(52) = DD(3,2,3) 
        V0(53) = DD(3,3,3) 
        ! trace of K 
        V0(54) = traceK 
        ! logarithm of the conformal factor phi 
        V0(55) = phi 
        ! derivative of phi 
        V0(56) = PP(1) 
        V0(57) = PP(2) 
        V0(58) = PP(3)         
        !
        CONTINUE 
        !         
     ELSE IF(IC.EQ.ICType%CCZ4GRMHDAccretion) THEN
        V0 = 0.0
        ! Compute the metric and its derivatives
        CALL METRIC(  xGP, alpha,  gp, gm, beta, Kex, g_cov, g_contr, phi )
        CALL DMETRIC( xGP, AA, BB, DD, PP )
        test = matmul(g_cov,g_contr)
        detg = g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)*g_cov(3,2)-g_cov(2,1)*g_cov(1,2)*g_cov(3,3)+g_cov(2,1)*g_cov(1,3)*g_cov(3,2)+g_cov(3,1)*g_cov(1,2)*g_cov(2,3)-g_cov(3,1)*g_cov(1,3)*g_cov(2,2)
        !phi = detg**(-1./6.)
        AA = AA/alpha   ! convert pure derivatives to aux. variables
        DD = 0.5*DD     ! convert pure derivatives to aux. variables
        PP = PP/phi     ! convert pure derivatives to aux. variables 
        betadown = MATMUL(g_cov, beta)
        ! Derivative of beta down
        dbetadown = 0.0
        DO j = 1, 3
         DO i = 1, 3
          DO k = 1, 3
            dbetadown(k,i) = dbetadown(k,i) + 2.0*DD(k,i,j)*beta(j) + g_cov(i,j)*BB(k,j)
          ENDDO
         ENDDO
        ENDDO
        ! Christoffel symbols
        Christoffel = 0.0
        DO k = 1, 3
         DO j = 1, 3
          DO i = 1, 3
              DO l = 1, 3
                Christoffel(i,j,k) = Christoffel(i,j,k) + g_contr(k,l)*( DD(i,j,l) + DD(j,i,l) - DD(l,i,j) )
              ENDDO
          ENDDO
         ENDDO
        ENDDO
        ! covariant derivative of the shift beta
        DO j = 1, 3
         DO i = 1, 3
           nablabeta(i,j) = dbetadown(i,j)
           DO k = 1, 3
              nablabeta(i,j) = nablabeta(i,j) - Christoffel(i,j,k)*betadown(k)
           ENDDO
         ENDDO
        ENDDO
        ! Extrinsic curvature
        Kex = 0.0
        DO j = 1, 3
         DO i = 1, 3
             Kex(i,j) = 1.0/(2*alpha)*( nablabeta(i,j) + nablabeta(j,i) )
         ENDDO
        ENDDO
        ! Trace of K
        traceK = SUM( Kex*g_contr )
        ! The conformal traceless part of the extrinsic curvature
        Aex = phi**2*( Kex - 1./3*traceK*g_cov )
        ! The contracted connection coefficients
        Gtilde = 0.0
        DO i = 1, 3
         DO j = 1, 3
          DO k = 1, 3
           DO l = 1, 3
               Gtilde(i) = Gtilde(i) + 1./phi**2*( g_contr(i,j)*g_contr(k,l)*( 2*DD(l,j,k) + 2*PP(l)*g_cov(j,k) )  )
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        SELECT CASE(EQN%CCZ4LapseType)
        CASE(0)  ! harmonic
           fa = 1.0
        CASE DEFAULT  ! 1 + log
           fa = 2.0/alpha
        END SELECT
        ! K0 to make the PDE for alpha stationary
        K0 = traceK - 1.0/(alpha*fa)*SUM(beta*AA)  
        !
        ! The metric tensor
        V0(1) = phi**2*g_cov(1,1)    ! \gamma_11
        V0(2) = phi**2*g_cov(1,2)    ! \gamma_12
        V0(3) = phi**2*g_cov(1,3)    ! \gamma_13
        V0(4) = phi**2*g_cov(2,2)    ! \gamma_22
        V0(5) = phi**2*g_cov(2,3)    ! \gamma_23
        V0(6) = phi**2*g_cov(3,3)    ! \gamma_33
        ! The extrinsic curvature        
        V0(7)  = Aex(1,1) 
        V0(8)  = Aex(1,2)
        V0(9)  = Aex(1,3)
        V0(10) = Aex(2,2)
        V0(11) = Aex(2,3)
        V0(12) = Aex(3,3)
        ! The cleaning variables Z and Theta
        V0(13) = 0.0          ! Theta
        V0(14) = Gtilde(1)    ! G1
        V0(15) = Gtilde(2)    ! G2
        V0(16) = Gtilde(3)    ! G3 
        ! The logarithm of the lapse
        V0(17) = alpha
        ! The shift
        V0(18) = beta(1) 
        V0(19) = beta(2)
        V0(20) = beta(3)
        !
        ! Auxiliary variables
        !
        !! vector A_i = \partial_i \alpha / \alpha
        !! 
        V0(24) = AA(1)     ! A1 = alpha_1/alpha 
        V0(25) = AA(2)     ! A2 = alpha_2/alpha   
        V0(26) = AA(3)     ! A3 = alpha_3/alpha  
        !
        ! Matrix B_ik = \partial_i \beta_k
        !
        V0(27) = BB(1,1)
        V0(28) = BB(2,1)
        V0(29) = BB(3,1)
        V0(30) = BB(1,2)
        V0(31) = BB(2,2)
        V0(32) = BB(3,2)
        V0(33) = BB(1,3)
        V0(34) = BB(2,3)
        V0(35) = BB(3,3)
        !
        ! tensor D_ijk = 0.5 \partial i \tilde \gamma_jk  
        !
        DO j = 1, 3
         DO i = 1, 3
          DO k = 1, 3
            DD(k,i,j) = phi**2*( DD(k,i,j) + PP(k)*g_cov(i,j) )
          ENDDO
         ENDDO
        ENDDO
        !
        V0(36) = DD(1,1,1)
        V0(37) = DD(1,1,2)
        V0(38) = DD(1,1,3)
        V0(39) = DD(1,2,2)
        V0(40) = DD(1,2,3)
        V0(41) = DD(1,3,3)
        V0(42) = DD(2,1,1)
        V0(43) = DD(2,1,2)
        V0(44) = DD(2,1,3)
        V0(45) = DD(2,2,2)
        V0(46) = DD(2,2,3)
        V0(47) = DD(2,3,3)
        V0(48) = DD(3,1,1)
        V0(49) = DD(3,1,2)
        V0(50) = DD(3,1,3)
        V0(51) = DD(3,2,2)
        V0(52) = DD(3,2,3)
        V0(53) = DD(3,3,3)
        ! trace of K
        V0(54) = traceK
        ! logarithm of the conformal factor phi
        V0(55) = phi
        ! derivative of phi
        V0(56) = PP(1)
        V0(57) = PP(2)
        V0(58) = PP(3)
        ! The trace of the extrinsic curvature at the initial time
        V0(59) = K0 

    END IF
    
    CALL PDEPrim2Cons(u0,V0) 
    !PRINT *, "U0(1)=", U0(1)
END SUBROUTINE InitialField
    
