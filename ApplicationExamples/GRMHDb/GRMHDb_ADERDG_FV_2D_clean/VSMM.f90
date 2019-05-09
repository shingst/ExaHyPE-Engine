#define INTRINSICS 
#define NPOLY3     
    
MODULE VSMM_mod

    IMPLICIT NONE 

INTERFACE PDEFluxVector
#if defined(EULER) && defined(INTRINSICS) 
    
RECURSIVE SUBROUTINE pdefluxvectoravx(FVec,VVec,QVec) bind ( c ) 
        use Parameters, ONLY : nVar, d, VECTORLENGTH  
        use iso_c_binding
        real ( c_double ) :: FVec(VECTORLENGTH,nVar,d)  
        real ( c_double ) :: VVec(VECTORLENGTH,nVar) 
        real ( c_double ) :: QVec(VECTORLENGTH,nVar) 
    END SUBROUTINE pdefluxvectoravx 
#else
    
RECURSIVE SUBROUTINE PDEFluxVectorF(FVec,VVec,QVec,QxVec,QyVec,QzVec) 
        use Parameters, ONLY : nVar, d, VECTORLENGTH  
        REAL :: FVec(VECTORLENGTH,nVar,d)  
        REAL :: VVec(VECTORLENGTH,nVar) 
        REAL :: QVec(VECTORLENGTH,nVar) 
        REAL :: QxVec(VECTORLENGTH,nVar) 
        REAL :: QyVec(VECTORLENGTH,nVar) 
        REAL :: QzVec(VECTORLENGTH,nVar) 
    END SUBROUTINE PDEFluxVectorF  
#endif 
END INTERFACE PDEFluxVector 

! INTERFACE FUNC_C2P_RMHD1_VECTOR
! #ifdef CINTRINSICS 
    ! SUBROUTINE FUNC_C2P_RMHD1_AVX(x,f,df,gam,dd,e,s2,b2,sb2,w) bind ( c ) 
        ! use Parameters, ONLY : VECTORLENGTH  
        ! use iso_c_binding
        ! real ( c_double ) :: x(VECTORLENGTH)  
        ! real ( c_double ) :: f(VECTORLENGTH)  
        ! real ( c_double ) :: df(VECTORLENGTH)  
        ! real ( c_double ) :: gam(VECTORLENGTH)  
        ! real ( c_double ) :: dd(VECTORLENGTH)  
        ! real ( c_double ) :: e(VECTORLENGTH)  
        ! real ( c_double ) :: s2(VECTORLENGTH)  
        ! real ( c_double ) :: b2(VECTORLENGTH)  
        ! real ( c_double ) :: sb2(VECTORLENGTH)  
        ! real ( c_double ) :: w(VECTORLENGTH)           
    ! END SUBROUTINE FUNC_C2P_RMHD1_AVX 
! #else
    ! SUBROUTINE FUNC_C2P_RMHD1_VECTORF(x,f,df,gam,dd,e,s2,b2,sb2,w,VECTORLENGTH) 
        ! INTEGER :: VECTORLENGTH  
        ! REAL(8) :: x(VECTORLENGTH)  
        ! REAL(8) :: f(VECTORLENGTH)  
        ! REAL(8) :: df(VECTORLENGTH)  
        ! REAL(8) :: gam(VECTORLENGTH)  
        ! REAL(8) :: dd(VECTORLENGTH)  
        ! REAL(8) :: e(VECTORLENGTH)  
        ! REAL(8) :: s2(VECTORLENGTH)  
        ! REAL(8) :: b2(VECTORLENGTH)  
        ! REAL(8) :: sb2(VECTORLENGTH)  
        ! REAL(8) :: w(VECTORLENGTH)    
! #ifdef AVX512 
            ! !DIR$ ASSUME_ALIGNED x      : 64 
            ! !DIR$ ASSUME_ALIGNED f      : 64 
            ! !DIR$ ASSUME_ALIGNED df     : 64 
            ! !DIR$ ASSUME_ALIGNED XACC   : 64 
            ! !DIR$ ASSUME_ALIGNED gam    : 64 
            ! !DIR$ ASSUME_ALIGNED dd     : 64 
            ! !DIR$ ASSUME_ALIGNED e      : 64 
            ! !DIR$ ASSUME_ALIGNED s2     : 64 
            ! !DIR$ ASSUME_ALIGNED b2     : 64 
            ! !DIR$ ASSUME_ALIGNED sb2    : 64 
            ! !DIR$ ASSUME_ALIGNED w      : 64 
! #else 
            ! !DIR$ ASSUME_ALIGNED x      : 32 
            ! !DIR$ ASSUME_ALIGNED f      : 32 
            ! !DIR$ ASSUME_ALIGNED df     : 32 
            ! !DIR$ ASSUME_ALIGNED XACC   : 32 
            ! !DIR$ ASSUME_ALIGNED gam    : 32 
            ! !DIR$ ASSUME_ALIGNED dd     : 32 
            ! !DIR$ ASSUME_ALIGNED e      : 32 
            ! !DIR$ ASSUME_ALIGNED s2     : 32 
            ! !DIR$ ASSUME_ALIGNED b2     : 32 
            ! !DIR$ ASSUME_ALIGNED sb2    : 32 
            ! !DIR$ ASSUME_ALIGNED w      : 32 
! #endif
    ! END SUBROUTINE FUNC_C2P_RMHD1_VECTORF  
! #endif 
! END INTERFACE FUNC_C2P_RMHD1_VECTOR  


INTERFACE ccz4test
#ifdef INTRINSICS 
    
RECURSIVE SUBROUTINE ccz4testc(dChristoffel_tilde,Z,Zup,Gtilde,Christoffel,Christoffel_tilde,dgup,QVec,QxV,QyV,QzV,ltime) bind ( c ) 
        use Parameters, ONLY : nVar, d, VECTORLENGTH  
        use iso_c_binding
        real ( c_double ) :: dChristoffel_tilde(VECTORLENGTH,3,3,3,3)
        real ( c_double ) :: Z(VECTORLENGTH,3)
        real ( c_double ) :: Zup(VECTORLENGTH,3)
        real ( c_double ) :: Gtilde(VECTORLENGTH,3)
        real ( c_double ) :: Christoffel(VECTORLENGTH,3,3,3)
        real ( c_double ) :: Christoffel_tilde(VECTORLENGTH,3,3,3)
        real ( c_double ) :: dgup(VECTORLENGTH,3,3,3)
        real ( c_double ) :: QVec(VECTORLENGTH,nVar)
        real ( c_double ) :: QxV(VECTORLENGTH,nVar)
        real ( c_double ) :: QyV(VECTORLENGTH,nVar)
        real ( c_double ) :: QzV(VECTORLENGTH,nVar)
        real ( c_double ) :: ltime 
    END SUBROUTINE ccz4testc 
#else
    
RECURSIVE SUBROUTINE ccz4testf(dChristoffel_tilde,Z,Zup,Gtilde,Christoffel,Christoffel_tilde,dgup,QVec,QxV,QyV,QzV,ltime) 
        use Parameters, ONLY : nVar, d, VECTORLENGTH  
        REAL :: dChristoffel_tilde(VECTORLENGTH,3,3,3,3) 
        REAL :: Z(VECTORLENGTH,3) 
        REAL :: Zup(VECTORLENGTH,3) 
        REAL :: Gtilde(VECTORLENGTH,3) 
        REAL :: Christoffel(VECTORLENGTH,3,3,3) 
        REAL :: Christoffel_tilde(VECTORLENGTH,3,3,3) 
        REAL :: dgup(VECTORLENGTH,3,3,3) 
        REAL :: QVec(VECTORLENGTH,nVar) 
        REAL :: QxV(VECTORLENGTH,nVar) 
        REAL :: QyV(VECTORLENGTH,nVar) 
        REAL :: QzV(VECTORLENGTH,nVar) 
        REAL :: ltime
    END SUBROUTINE ccz4testf   
#endif 
END INTERFACE ccz4test  

 INTERFACE VSMM0 
#ifdef INTRINSICS  
    
RECURSIVE SUBROUTINE vsmm4x4c( Qx, Q, dudx ) bind ( c )
      use Parameters, ONLY : nVar 
      use iso_c_binding
      real ( c_double ) :: Qx(nVar,4)  
      real ( c_double ) :: Q(nVar,4) 
      real ( c_double ) :: dudx(4,4) 
    END SUBROUTINE vsmm4x4c
#else
    MODULE PROCEDURE VSMMf 
#endif 
  END INTERFACE VSMM0 
    
 INTERFACE VSMMS 
#ifdef INTRINSICS  
    
RECURSIVE SUBROUTINE vsmm4x4cs( Qx, Q, dudx, s ) bind ( c )
      use Parameters, ONLY : nVar 
      use iso_c_binding
      real ( c_double ) :: Qx(nVar,4)  
      real ( c_double ) :: Q(nVar,4) 
      real ( c_double ) :: dudx(4,4) 
      real ( c_double ) :: s  
    END SUBROUTINE vsmm4x4cs
#else
    MODULE PROCEDURE VSMMsf   
#endif 
 END INTERFACE VSMMS      

INTERFACE cmatmultest
    
RECURSIVE SUBROUTINE cmatmultest( Qx, Q, dudx, s, nElem ) bind ( c )
      use Parameters, ONLY : nVar 
      use iso_c_binding
      real ( c_double ) :: Qx(nVar,4,4,4,4)  
      real ( c_double ) :: Q(nVar,4,4,4,4) 
      real ( c_double ) :: dudx(4,4)       
      real ( c_double ) :: s  
      integer (c_int )  :: nElem 
    END SUBROUTINE cmatmultest 
END INTERFACE cmatmultest  
 
INTERFACE VSMMs4dall_x
    
RECURSIVE SUBROUTINE vsmm4x4cs4d_x( Qx, Q, dudx, s ) bind ( c )
      use Parameters, ONLY : nVar 
      use iso_c_binding
      real ( c_double ) :: Qx(nVar,4,4,4,4)  
      real ( c_double ) :: Q(nVar,4,4,4,4) 
      real ( c_double ) :: dudx(4,4) 
      real ( c_double ) :: s  
#ifdef AVX512       
!DIR$ ASSUME_ALIGNED Qx         : 64 
!DIR$ ASSUME_ALIGNED Q          : 64 
!DIR$ ASSUME_ALIGNED dudx       : 64 
#else
!DIR$ ASSUME_ALIGNED Qx         : 32 
!DIR$ ASSUME_ALIGNED Q          : 32 
!DIR$ ASSUME_ALIGNED dudx       : 32 
#endif
    END SUBROUTINE vsmm4x4cs4d_x 
END INTERFACE VSMMs4dall_x 

INTERFACE VSMMs4dall_y
    
RECURSIVE SUBROUTINE vsmm4x4cs4d_y( Qy, Q, dudx, s ) bind ( c )
      use Parameters, ONLY : nVar 
      use iso_c_binding
      real ( c_double ) :: Qy(nVar,4,4,4,4)  
      real ( c_double ) :: Q(nVar,4,4,4,4) 
      real ( c_double ) :: dudx(4,4) 
      real ( c_double ) :: s  
    END SUBROUTINE vsmm4x4cs4d_y 
END INTERFACE VSMMs4dall_y 

INTERFACE VSMMs4dall_z
    
RECURSIVE SUBROUTINE vsmm4x4cs4d_z( Qz, Q, dudx, s ) bind ( c )
      use Parameters, ONLY : nVar 
      use iso_c_binding
      real ( c_double ) :: Qz(nVar,4,4,4,4)  
      real ( c_double ) :: Q(nVar,4,4,4,4) 
      real ( c_double ) :: dudx(4,4) 
      real ( c_double ) :: s  
    END SUBROUTINE vsmm4x4cs4d_z 
END INTERFACE VSMMs4dall_z 

INTERFACE VSMMs3dall_x
    
RECURSIVE SUBROUTINE vsmm4x4cs3d_x( Qx, Q, dudx, s ) bind ( c )
      use Parameters, ONLY : nVar 
      use iso_c_binding
      real ( c_double ) :: Qx(nVar,4,4,4)  
      real ( c_double ) :: Q(nVar,4,4,4) 
      real ( c_double ) :: dudx(4,4) 
      real ( c_double ) :: s  
    END SUBROUTINE vsmm4x4cs3d_x 
END INTERFACE VSMMs3dall_x 

INTERFACE VSMMs3dall_y
    
RECURSIVE SUBROUTINE vsmm4x4cs3d_y( Qy, Q, dudx, s ) bind ( c )
      use Parameters, ONLY : nVar 
      use iso_c_binding
      real ( c_double ) :: Qy(nVar,4,4,4)  
      real ( c_double ) :: Q(nVar,4,4,4) 
      real ( c_double ) :: dudx(4,4) 
      real ( c_double ) :: s  
    END SUBROUTINE vsmm4x4cs3d_y 
END INTERFACE VSMMs3dall_y 

INTERFACE VSMMs3dall_z
    
RECURSIVE SUBROUTINE vsmm4x4cs3d_z( Qz, Q, dudx, s ) bind ( c )
      use Parameters, ONLY : nVar 
      use iso_c_binding
      real ( c_double ) :: Qz(nVar,4,4,4)  
      real ( c_double ) :: Q(nVar,4,4,4) 
      real ( c_double ) :: dudx(4,4) 
      real ( c_double ) :: s  
    END SUBROUTINE vsmm4x4cs3d_z 
END INTERFACE VSMMs3dall_z 

INTERFACE VSMMs4dall_kx
    SUBROUTINE vsmm4x4c4d_kx( Qt, Q, dudx, s ) bind ( c )
      use Parameters, ONLY : nVar 
      use iso_c_binding
      real ( c_double ) :: Qt(nVar,4,4,4,4)  
      real ( c_double ) :: Q(nVar,4,4,4,4) 
      real ( c_double ) :: dudx(4,4) 
      real ( c_double ) :: s  
    END SUBROUTINE vsmm4x4c4d_kx 
END INTERFACE VSMMs4dall_kx 

INTERFACE VSMMs4dall_ky
    
RECURSIVE SUBROUTINE vsmm4x4c4d_ky( Qt, Q, dudx, s ) bind ( c )
      use Parameters, ONLY : nVar 
      use iso_c_binding
      real ( c_double ) :: Qt(nVar,4,4,4,4)  
      real ( c_double ) :: Q(nVar,4,4,4,4) 
      real ( c_double ) :: dudx(4,4) 
      real ( c_double ) :: s  
    END SUBROUTINE vsmm4x4c4d_ky 
END INTERFACE VSMMs4dall_ky 

INTERFACE VSMMs4dall_kz
    
RECURSIVE SUBROUTINE vsmm4x4c4d_kz( Qt, Q, dudx, s ) bind ( c )
      use Parameters, ONLY : nVar 
      use iso_c_binding
      real ( c_double ) :: Qt(nVar,4,4,4,4)  
      real ( c_double ) :: Q(nVar,4,4,4,4) 
      real ( c_double ) :: dudx(4,4) 
      real ( c_double ) :: s  
    END SUBROUTINE vsmm4x4c4d_kz 
END INTERFACE VSMMs4dall_kz 


INTERFACE VSMMs4dall_kt
    
RECURSIVE SUBROUTINE vsmm4x4c4d_kt( Qt, Q, dudx ) bind ( c )
      use Parameters, ONLY : nVar 
      use iso_c_binding
      real ( c_double ) :: Qt(nVar,4,4,4,4)  
      real ( c_double ) :: Q(nVar,4,4,4,4) 
      real ( c_double ) :: dudx(4,4) 
    END SUBROUTINE vsmm4x4c4d_kt 
END INTERFACE VSMMs4dall_kt 


INTERFACE VSMMreduceall_t
    
RECURSIVE SUBROUTINE vsmm4x4ra_t( Qi, Q, vec ) bind ( c )
      use Parameters, ONLY : nVar 
      use iso_c_binding
      real ( c_double ) :: Qi(nVar,4,4,4)  
      real ( c_double ) :: Q(nVar,4,4,4,4) 
      real ( c_double ) :: vec(4) 
    END SUBROUTINE vsmm4x4ra_t  
END INTERFACE VSMMreduceall_t  

INTERFACE VSMMreduceall_x
    
RECURSIVE SUBROUTINE vsmm4x4ra_x( Qi, Q, vec ) bind ( c )
      use Parameters, ONLY : nVar 
      use iso_c_binding
      real ( c_double ) :: Qi(nVar,4,4)  
      real ( c_double ) :: Q(nVar,4,4,4) 
      real ( c_double ) :: vec(4) 
    END SUBROUTINE vsmm4x4ra_x  
END INTERFACE VSMMreduceall_x  

INTERFACE VSMMreduceall_y
    
RECURSIVE SUBROUTINE vsmm4x4ra_y( Qi, Q, vec ) bind ( c )
      use Parameters, ONLY : nVar 
      use iso_c_binding
      real ( c_double ) :: Qi(nVar,4,4)  
      real ( c_double ) :: Q(nVar,4,4,4) 
      real ( c_double ) :: vec(4) 
    END SUBROUTINE vsmm4x4ra_y  
END INTERFACE VSMMreduceall_y

INTERFACE VSMMreduceall_z
    
RECURSIVE SUBROUTINE vsmm4x4ra_z( Qi, Q, vec ) bind ( c )
      use Parameters, ONLY : nVar 
      use iso_c_binding
      real ( c_double ) :: Qi(nVar,4,4)  
      real ( c_double ) :: Q(nVar,4,4,4) 
      real ( c_double ) :: vec(4) 
    END SUBROUTINE vsmm4x4ra_z  
END INTERFACE VSMMreduceall_z  

PUBLIC :: VSMM0, VSMMS, vsmm4x4cs4d_x, vsmm4x4cs4d_y, vsmm4x4cs4d_z, VSMMs4dall_Kt, VSMMreduceall_t, VSMMreduceall_x, VSMMreduceall_y, VSMMreduceall_z   
PUBLIC :: PDEFluxVector 
 
CONTAINS     
    

RECURSIVE SUBROUTINE VSMMsf(C,A,B,s) 
    USE Parameters, ONLY : nVar, N, nAdd, nMul  
    REAL(8)       :: C(nVar,N+1), A(nVar,N+1), B(N+1,N+1), s  
    REAL(8)       :: k1, k2, k3, k4, k5 
    INTEGER       :: j, k  
    INTENT(IN)    :: A, B, s 
    INTENT(INOUT) :: C 
#ifndef CRAYFTN 
#ifdef AVX512     
    !DIR$ ASSUME_ALIGNED A       : 64 
    !DIR$ ASSUME_ALIGNED B       : 64 
    !DIR$ ASSUME_ALIGNED C       : 64 
#else 
    !DIR$ ASSUME_ALIGNED A       : 32 
    !DIR$ ASSUME_ALIGNED B       : 32 
    !DIR$ ASSUME_ALIGNED C       : 32 
#endif
#endif
    !
    !C = 0.0 
    !RETURN 
    continue
    !
#if defined(NPOLY2) 
        C(:,1) = C(:,1) + A(:,1)*B(1,1) + A(:,2)*B(2,1) + A(:,3)*B(3,1)  
        C(:,2) = C(:,2) + A(:,2)*B(1,2) + A(:,2)*B(2,2) + A(:,3)*B(3,2)  
        C(:,3) = C(:,3) + A(:,3)*B(1,3) + A(:,2)*B(2,3) + A(:,3)*B(3,3)  
#elif defined(NPOLY3) 
        !k1 = B(1,1); k2=B(2,1); k3=B(3,1); k4=B(4,1)
        !C(:,1) = C(:,1) + A(:,1)*k1 + A(:,2)*k2 + A(:,3)*k3 + A(:,4)*k4 
        !k1 = B(1,2); k2=B(2,2); k3=B(3,2); k4=B(4,2)
        !C(:,2) = C(:,2) + A(:,1)*k1 + A(:,2)*k2 + A(:,3)*k3 + A(:,4)*k4 
        !k1 = B(1,3); k2=B(2,3); k3=B(3,3); k4=B(4,3)
        !C(:,3) = C(:,3) + A(:,1)*k1 + A(:,2)*k2 + A(:,3)*k3 + A(:,4)*k4 
        !k1 = B(1,4); k2=B(2,4); k3=B(3,4); k4=B(4,4)
        !C(:,4) = C(:,4) + A(:,1)*k1 + A(:,2)*k2 + A(:,3)*k3 + A(:,4)*k4 
        C(:,1) = A(:,1)*B(1,1) + A(:,2)*B(2,1) + A(:,3)*B(3,1) + A(:,4)*B(4,1)  
        C(:,2) = A(:,1)*B(1,2) + A(:,2)*B(2,2) + A(:,3)*B(3,2) + A(:,4)*B(4,2)  
        C(:,3) = A(:,1)*B(1,3) + A(:,2)*B(2,3) + A(:,3)*B(3,3) + A(:,4)*B(4,3)  
        C(:,4) = A(:,1)*B(1,4) + A(:,2)*B(2,4) + A(:,3)*B(3,4) + A(:,4)*B(4,4)  
        nAdd = nAdd + 3*4*nVar 
        nMul = nMul + 4*4*nVar + 4*4  ! the last 4*4 is the for scaling B with s 
#elif defined(NPOLY4) 
        C(:,1) = A(:,1)*B(1,1) + A(:,2)*B(2,1) + A(:,3)*B(3,1) + A(:,4)*B(4,1) + A(:,5)*B(5,1)
        C(:,2) = A(:,1)*B(1,2) + A(:,2)*B(2,2) + A(:,3)*B(3,2) + A(:,4)*B(4,2) + A(:,5)*B(5,2)
        C(:,3) = A(:,1)*B(1,3) + A(:,2)*B(2,3) + A(:,3)*B(3,3) + A(:,4)*B(4,3) + A(:,5)*B(5,3)
        C(:,4) = A(:,1)*B(1,4) + A(:,2)*B(2,4) + A(:,3)*B(3,4) + A(:,4)*B(4,4) + A(:,5)*B(5,4)
        C(:,5) = A(:,1)*B(1,5) + A(:,2)*B(2,5) + A(:,3)*B(3,5) + A(:,4)*B(4,5) + A(:,5)*B(5,5)
#elif defined(NPOLY5) 
        C(:,1) = A(:,1)*B(1,1) + A(:,2)*B(2,1) + A(:,3)*B(3,1) + A(:,4)*B(4,1) + A(:,5)*B(5,1) + A(:,6)*B(6,1)
        C(:,2) = A(:,1)*B(1,2) + A(:,2)*B(2,2) + A(:,3)*B(3,2) + A(:,4)*B(4,2) + A(:,5)*B(5,2) + A(:,6)*B(6,2)
        C(:,3) = A(:,1)*B(1,3) + A(:,2)*B(2,3) + A(:,3)*B(3,3) + A(:,4)*B(4,3) + A(:,5)*B(5,3) + A(:,6)*B(6,3)
        C(:,4) = A(:,1)*B(1,4) + A(:,2)*B(2,4) + A(:,3)*B(3,4) + A(:,4)*B(4,4) + A(:,5)*B(5,4) + A(:,6)*B(6,4)
        C(:,5) = A(:,1)*B(1,5) + A(:,2)*B(2,5) + A(:,3)*B(3,5) + A(:,4)*B(4,5) + A(:,5)*B(5,5) + A(:,6)*B(6,5)
        C(:,6) = A(:,1)*B(1,6) + A(:,2)*B(2,6) + A(:,3)*B(3,6) + A(:,4)*B(4,6) + A(:,5)*B(5,6) + A(:,6)*B(6,6)
#else 
        C(:,:) = 0.0 
        DO j = 1, N+1  
         DO k = 1, N+1 
            C(:,j) = C(:,j) + A(:,k)*B(k,j) 
         ENDDO
        ENDDO        
#endif 
    C(:,:) = s*C(:,:) 
    !
END SUBROUTINE VSMMsf      
    

RECURSIVE SUBROUTINE VSMMf(C,A,B) 
    USE Parameters, ONLY : nVar, N 
    REAL(8) :: C(nVar,N+1), A(nVar,N+1), B(N+1,N+1) 
    INTEGER :: j, k  
#ifndef CRAYFTN 
#ifdef AVX512     
    !DIR$ ASSUME_ALIGNED A       : 64 
    !DIR$ ASSUME_ALIGNED B       : 64 
    !DIR$ ASSUME_ALIGNED C       : 64 
#else 
    !DIR$ ASSUME_ALIGNED A       : 32 
    !DIR$ ASSUME_ALIGNED B       : 32 
    !DIR$ ASSUME_ALIGNED C       : 32 
#endif
#endif
    !
    !C = 0.0 
    !RETURN 
    ! 
#if defined(NPOLY2) 
        C(:,1) = A(:,1)*B(1,1) + A(:,2)*B(2,1) + A(:,3)*B(3,1)  
        C(:,2) = A(:,1)*B(1,2) + A(:,2)*B(2,2) + A(:,3)*B(3,2)  
        C(:,3) = A(:,1)*B(1,3) + A(:,2)*B(2,3) + A(:,3)*B(3,3)  
#elif defined(NPOLY3) 
        C(:,1) = A(:,1)*B(1,1) + A(:,2)*B(2,1) + A(:,3)*B(3,1) + A(:,4)*B(4,1) 
        C(:,2) = A(:,1)*B(1,2) + A(:,2)*B(2,2) + A(:,3)*B(3,2) + A(:,4)*B(4,2) 
        C(:,3) = A(:,1)*B(1,3) + A(:,2)*B(2,3) + A(:,3)*B(3,3) + A(:,4)*B(4,3) 
        C(:,4) = A(:,1)*B(1,4) + A(:,2)*B(2,4) + A(:,3)*B(3,4) + A(:,4)*B(4,4) 
#elif defined(NPOLY4) 
        C(:,1) = A(:,1)*B(1,1) + A(:,2)*B(2,1) + A(:,3)*B(3,1) + A(:,4)*B(4,1) + A(:,5)*B(5,1)
        C(:,2) = A(:,1)*B(1,2) + A(:,2)*B(2,2) + A(:,3)*B(3,2) + A(:,4)*B(4,2) + A(:,5)*B(5,2)
        C(:,3) = A(:,1)*B(1,3) + A(:,2)*B(2,3) + A(:,3)*B(3,3) + A(:,4)*B(4,3) + A(:,5)*B(5,3)
        C(:,4) = A(:,1)*B(1,4) + A(:,2)*B(2,4) + A(:,3)*B(3,4) + A(:,4)*B(4,4) + A(:,5)*B(5,4)
        C(:,5) = A(:,1)*B(1,5) + A(:,2)*B(2,5) + A(:,3)*B(3,5) + A(:,4)*B(4,5) + A(:,5)*B(5,5)
#elif defined(NPOLY5) 
        C(:,1) = A(:,1)*B(1,1) + A(:,2)*B(2,1) + A(:,3)*B(3,1) + A(:,4)*B(4,1) + A(:,5)*B(5,1) + A(:,6)*B(6,1)
        C(:,2) = A(:,1)*B(1,2) + A(:,2)*B(2,2) + A(:,3)*B(3,2) + A(:,4)*B(4,2) + A(:,5)*B(5,2) + A(:,6)*B(6,2)
        C(:,3) = A(:,1)*B(1,3) + A(:,2)*B(2,3) + A(:,3)*B(3,3) + A(:,4)*B(4,3) + A(:,5)*B(5,3) + A(:,6)*B(6,3)
        C(:,4) = A(:,1)*B(1,4) + A(:,2)*B(2,4) + A(:,3)*B(3,4) + A(:,4)*B(4,4) + A(:,5)*B(5,4) + A(:,6)*B(6,4)
        C(:,5) = A(:,1)*B(1,5) + A(:,2)*B(2,5) + A(:,3)*B(3,5) + A(:,4)*B(4,5) + A(:,5)*B(5,5) + A(:,6)*B(6,5)
        C(:,6) = A(:,1)*B(1,6) + A(:,2)*B(2,6) + A(:,3)*B(3,6) + A(:,4)*B(4,6) + A(:,5)*B(5,6) + A(:,6)*B(6,6)
#else 
        C(:,:) = 0.0 
        DO j = 1, N+1  
         DO k = 1, N+1 
            C(:,j) = C(:,j) + A(:,k)*B(k,j) 
         ENDDO
        ENDDO        
#endif 
    !
END SUBROUTINE VSMMf      
    


RECURSIVE SUBROUTINE VSMV(C,A,B) 
    USE Parameters, ONLY : nVar, N, nAdd, nMul  
    REAL(8) :: C(nVar), A(nVar,N+1), B(N+1) 
    INTEGER :: j, k  
#ifndef CRAYFTN  
#ifdef AVX512     
    !DIR$ ASSUME_ALIGNED A       : 64 
    !DIR$ ASSUME_ALIGNED B       : 64 
    !DIR$ ASSUME_ALIGNED C       : 64 
#else 
    !DIR$ ASSUME_ALIGNED A       : 32 
    !DIR$ ASSUME_ALIGNED B       : 32 
    !DIR$ ASSUME_ALIGNED C       : 32 
#endif
#endif
    !
#if defined(NPOLY2) 
        C(:) = A(:,1)*B(1) + A(:,2)*B(2) + A(:,3)*B(3)  
#elif defined(NPOLY3) 
        C(:) = A(:,1)*B(1) + A(:,2)*B(2) + A(:,3)*B(3) + A(:,4)*B(4)  
        nAdd = nAdd + 3*nVar
        nMul = nMul + 4*nVar 
#elif defined(NPOLY4) 
        C(:) = A(:,1)*B(1) + A(:,2)*B(2) + A(:,3)*B(3) + A(:,4)*B(4) + A(:,5)*B(5)  
#elif defined(NPOLY5) 
        C(:) = A(:,1)*B(1) + A(:,2)*B(2) + A(:,3)*B(3) + A(:,4)*B(4) + A(:,5)*B(5) + A(:,6)*B(6)
#else 
        C(:) = 0.0 
        DO k = 1, N+1 
            C(:) = C(:) + A(:,k)*B(k) 
        ENDDO        
#endif 
    !
END SUBROUTINE VSMV      
        

END MODULE VSMM_mod 