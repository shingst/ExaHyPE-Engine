! DIM Parameters
  
  MODULE MainVariables 
    IMPLICIT NONE  
    PUBLIC  

    ! This typesDef.f90 is a relict still from the old Fortran interface in Exahype.
    ! However, the following two variables are needed in the Fortran code. They
    ! should provided by the glue code generator in an appropriate way.

    ! If you modify the SRHD.exahype, please do the following mapping by hand:
    !
    ! solver ADER-DG SRHDSolver / unknowns   ->  goes to ->  nVar
    ! computational-domain / dimension       ->  goes to ->  nDim
    !

    ! Here, we obtain DIMENSIONS from the C++ preprocessor
#if defined(Dim3)
    INTEGER, PARAMETER             	:: nDim = 3                   ! The number of space dimensions
#elif defined(Dim2)
    INTEGER, PARAMETER             	:: nDim = 2                   ! The number of space dimensions
#endif
	INTEGER, PARAMETER             	:: nAux = 0
    INTEGER, PARAMETER             	:: nVar = 9                           ! The number of variables of the PDE system 
    INTEGER, PARAMETER 				:: nLin = 7
    CHARACTER(LEN=20), PARAMETER	:: ICType='NBodies'
	
TYPE tEquation 
     REAL    :: gamma, gamma2, g, alpha, beta, nu, sigma, pi, rho0, rho1, c0, p0, p2, k0, gz, cs, cv, cv1, cv2, tau, tau2, ch, cp
     REAL    :: lambda, mu, chi, D12, m1, m2, Kbeta, phi, Cf, Cr, theta   
     REAL    :: gamma1, pi1, pi2, rho2, epsilon, g1, g2
     REAL    :: lambda1, lambda2, mu1, mu2, kappa1, kappa2, pr1, pr2, K1, K2, KK, Y0, aexp, epsilon1, cl, cl2        
     REAL    :: Su, T0, Pr, Sc, kappa, R, eta, taul, taus, cleaning_a, cleaning_eps
  
     REAL    :: rhoref(2) 
     REAL    :: CCZ4k1, CCZ4k2, CCZ4k3, CCZ4eta, CCZ4f, CCZ4g, CCZ4xi, CCZ4e, CCZ4c, CCZ4mu, CCZ4bs, CCZ4itau, CCZ4ds
     REAL    :: switch_threshold, switch_exponent
     INTEGER :: CCZ4sk = 1 
     INTEGER :: SubType                          ! Equation subtype
     INTEGER :: CCZ4LapseType = 1                ! default is 1+log. Set this value to 0 for harmonic slicing. 
  END TYPE tEquation
	
	TYPE(tEquation) :: EQN
	REAL :: ICuL(nVar), EPCenter(1:nDim), EPRadius, ICA, ICdelta, MHDRotomega
	INTEGER :: ICrot
  ! 3-point Gaussian quadrature 
  REAL, PARAMETER     :: sGP3(3) = (/ 0.5-sqrt(15.)/10.,0.5,0.5+sqrt(15.)/10. /) 
  REAL, PARAMETER     :: wGP3(3) = (/ 5./18., 8./18., 5./18. /) 
  INTEGER, PARAMETER  :: nGP3 = 3 

  END MODULE MainVariables  
