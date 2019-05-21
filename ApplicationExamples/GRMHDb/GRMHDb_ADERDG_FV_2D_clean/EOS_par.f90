
MODULE EOS_par 
      !INTEGER, PARAMETER :: EOSDG_N = 3
      !INTEGER, PARAMETER :: EOSDG_IMAX_rho = 200
      !INTEGER, PARAMETER :: EOSDG_JMAX_eps = 200
      !TYPE EOS_type 
      !   INTEGER                :: N,nDOF,IMAX_rho,JMAX_eps
      !   REAL                   :: rhoL,rhoR,epsL,epsR,pL,pR
      !   INTEGER,   ALLOCATABLE :: rsa(:,:) ca
      !END TYPE EOS_type
      !TYPE(EOS_type) :: EOS 
    INTEGER, PARAMETER :: EOSfunc = 0   ! functional type:  0 = pressure eos  p - p[rho(p),eps(p)]; 1 =  p^2 - p[rho(p),eps(p)]^2
    INTEGER, PARAMETER :: EOStype = 0   ! 0 = ideal gases; 1 = non-ideal gas ; 2 = EOS-DG tabulated
    
    !
END MODULE EOS_par
