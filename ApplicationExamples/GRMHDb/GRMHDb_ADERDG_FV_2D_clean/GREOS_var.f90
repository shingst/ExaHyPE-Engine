MODULE EOSTypes
    IMPLICIT NONE 
    PUBLIC 
#ifdef AVX512    
    INTEGER, PARAMETER :: VECTORLENGTH = 8 
#else
    INTEGER, PARAMETER :: VECTORLENGTH = 4 
#endif
    INTEGER            :: OriginalTable = 1 
    INTEGER,  PARAMETER:: EOSNpoly=7 
    CHARACTER(LEN=200) :: InputFile, EOSFile    
    REAL(8)            :: Lag_U(EOSNpoly+1),Lag_L(EOSNpoly+1),Lag_A(EOSNpoly+1)
    REAL(8)            :: Lag_U_xi(EOSNpoly+1),Lag_L_xi(EOSNpoly+1) 
    REAL(8)            :: Lag_U_c(EOSNpoly+1,EOSNpoly+1),Lag_L_c(EOSNpoly+1,EOSNpoly+1) 
    !
    ! EOS Table 
    TYPE tEOSTable 
        INTEGER           :: N                                                                          ! polynomial approximation degree in the table 
        INTEGER           :: NRHO, NT                                                                   ! number of points in the input table 
        INTEGER           :: IMAX, JMAX                                                                 ! number of elements in the output table         
        REAL              :: Tmin, Tmax, rhomin, rhomax, emin, emax, pmin, pmax, drho, de
        REAL              :: drhoin, dTin     
        REAL, ALLOCATABLE :: xiGP(:), wGP(:)    
        REAL, ALLOCATABLE :: T(:), rho(:) 
        REAL, ALLOCATABLE :: p(:,:), e(:,:) 
        REAL, ALLOCATABLE :: DGp(:,:,:,:), DGT(:,:,:,:), DGe(:,:,:,:) 
    END TYPE tEOSTable 
    
    TYPE(tEOSTable)       :: EOS 
    
    INTEGER :: EOSN,EOSIMAX,EOSJMAX
    REAL(8), ALLOCATABLE :: EOSDGp(:,:,:,:),EOSDGT(:,:,:,:),EOSxiGP(:),EOSwGP(:)
    REAL(8) ::EOSdrho,EOSde 
#ifdef AVX512 
	  !DIR$ attributes align:64 :: EOSDGp, EOSDGT 
	  !!DIR$ ASSUME_ALIGNED EOSDGp:64, EOSDGT:64 
#else 
	  !DIR$ attributes align:32 :: EOSDGp, EOSDGT 
	  !!DIR$ ASSUME_ALIGNED  EOSDGp:32, EOSDGT:32 
#endif 
END MODULE EOSTypes
!
!    
!