! GRMHD PDE.f90
! Trento (EQNTYPE4)


RECURSIVE SUBROUTINE PDEPrim2Cons(Q,V)
  USE Parameters, ONLY: nVar
  USE GRMHD_Mod
  IMPLICIT NONE
  ! Argument list declaration
  REAL :: Q(nVar), V(nVar)
  INTENT(IN)  :: V
  INTENT(OUT) :: Q 
  !
  CALL PDEPrim2ConsGRMHD(Q,V)
  !
  continue
  !
END SUBROUTINE PDEPrim2Cons

RECURSIVE SUBROUTINE PDECons2Prim(V,Q,iErr)
  USE Parameters, ONLY: nVar
  USE GRMHD_Mod
  IMPLICIT NONE
  !--------------------------------------------!
  ! Argument list declaration
  REAL :: Q(nVar), V(nVar)
  INTEGER :: iErr
  INTENT(IN)  :: Q 
  INTENT(OUT) :: V  
  INTENT(INOUT) :: iErr 
  !
  CALL PDECons2PrimGRMHD(V,Q,iErr)
  !
  continue
  !
END SUBROUTINE PDECons2Prim




!RECURSIVE SUBROUTINE PDEFlux(FF,Q)
RECURSIVE SUBROUTINE PDEFlux(f,g,h,Q)
  USE Parameters, ONLY : d,nVar,ndim
  USE GRMHD_Mod
  USE iso_c_binding
  IMPLICIT NONE
  REAL :: f(nVar), g(nVar), h(nVar), Q(nVar)
  !REAL :: FF(nVar,ndim) , Q(nVar)
  INTENT(IN)  :: Q
  INTENT(OUT) :: f, g, h
  !INTENT(OUT) :: FF
  ! Local varialbes
  INTEGER :: i
  REAL :: FF(nVar,d), V(nVar)
  !REAL :: V(nVar)
  REAL, PARAMETER :: epsilon = 1e-14 
  !
  !PRINT *, 'PDEFluxGRMHD'
  CALL PDEFluxGRMHD(FF,Q)
  !FF = 0.
  !
#if defined(Dim2)
  DO i=1,nVar
		f(i)=FF(i,1)
		g(i)=FF(i,2)
		h(i)=0.0
  ENDDO
#else
  DO i=1,nVar
		f(i)=FF(i,1)
		g(i)=FF(i,2)
		h(i)=FF(i,3)
  ENDDO
#endif
  
  ! DO i=1,nDim
	! IF(MAXVAL(ABS(FF(:,i))).LT.1e-8) THEN 
			! print *, 'FATAL ERROR: PDEFlux' 
			! WRITE(*,*) i,FF(1:8,i)
			! STOP 
	! ENDIF
  ! ENDDO
	!  
  !
	! IF( ANY( f.NE.0.0 )) THEN
		! PRINT *,' PDEFlux nonzero ' 
		! STOP
	! ENDIF  
	! IF( ANY( g.NE.0.0 )) THEN
		! PRINT *,' PDEFlux nonzero ' 
		! STOP
	! ENDIF  
!	IF( ANY( ISNAN(f) )) THEN
!		PRINT *,' PDEFlux NAN ' 
!		STOP
!	ENDIF 
!	IF( ANY( ISNAN(g) )) THEN
!		PRINT *,' PDEFlux NAN ' 
!		STOP
!	ENDIF 
!#if defined(Dim3)
!	IF( ANY( ISNAN(h) )) THEN
!		PRINT *,' PDEFlux NAN ' 
!		STOP
!	ENDIF
!#endif
!	IF(ANY(Q(6:8).NE.0.0)) THEN
!        PRINT *, Q(6:8)
!        PRINT *, "I feel magnetized :-("
!        ERROR STOP
!    ENDIF
    !
  END SUBROUTINE PDEFlux


RECURSIVE SUBROUTINE PDENCP(BgradQ,Q,gradQ) 
   USE Parameters, ONLY :  nVar, nDim
   USE GRMHD_Mod
   IMPLICIT NONE
   ! 11. Oct 21:40: This was a matrix BGradQ(nVar, nDim) but is a vector in spaceTimePredictorNonlinear
   REAL, INTENT(OUT) :: BgradQ(nVar)
   REAL, INTENT(IN)  :: gradQ(nVar, nDim)
   REAL, INTENT(IN)  :: Q(nVar)
   ! Local variables
   INTEGER :: i
   !
#ifdef DEBUGGONE
   BgradQ = 0.
  STOP
   RETURN
#endif
	!
   CALL PDENCPGRMHD(BgradQ,Q,gradQ)

!IF(ANY(BgradQ(6:8).NE.0)) THEN
!    PRINT *, "PDENCP BgradQ(6:8)",BgradQ(6:8) 
!    continue
!ENDIF
   !BgradQ = 0.
   !
	! IF(MAXVAL(ABS(BgradQ(:))).LT.1e-8) THEN 
			! print *, 'FATAL ERROR: PDENCP' 
			! WRITE(*,*) BgradQ(1:8)
			! STOP 
	! ENDIF
	! IF( ANY( BgradQ.NE.0.0 )) THEN
		! PRINT *,' PDENCP nonzero ' 
		! STOP
	! ENDIF   
!	IF( ANY( ISNAN(gradQ(:,1:2)) )) THEN
!		PRINT *,'gradQ, PDENCP NAN ' 
!		WRITE(*,*) Q(1:5)
!		WRITE(*,*) Q(6:10)
!		WRITE(*,*) Q(11:15)
!		WRITE(*,*) Q(16:19)
!		PRINT *,'---------------' 
!		WRITE(*,*) gradQ(1:5,1)
!		WRITE(*,*) gradQ(6:10,1)
!		WRITE(*,*) gradQ(11:15,1)
!		WRITE(*,*) gradQ(16:19,1)
!		PRINT *,'---------------' 
!		WRITE(*,*) gradQ(1:5,2)
!		WRITE(*,*) gradQ(6:10,2)
!		WRITE(*,*) gradQ(11:15,2)
!		WRITE(*,*) gradQ(16:19,2)
!		PRINT *,'---------------' 
!		STOP
!	ENDIF 
!#if defined(Dim3)
!	IF( ANY( ISNAN(gradQ(:,3)) )) THEN
!		PRINT *,'gradQ, PDENCP NAN ' 
!		STOP
!	ENDIF 
!#endif
!	IF( ANY( ISNAN(BgradQ) )) THEN
!		PRINT *,'BgradQ,  PDENCP NAN ' 
!		WRITE(*,*) Q(1:5)
!		WRITE(*,*) Q(6:10)
!		WRITE(*,*) Q(11:15)
!		WRITE(*,*) Q(16:19)
!		PRINT *,'---------------' 
!		WRITE(*,*) gradQ(1:5,1)
!		WRITE(*,*) gradQ(6:10,1)
!		WRITE(*,*) gradQ(11:15,1)
!		WRITE(*,*) gradQ(16:19,1)
!		PRINT *,'---------------' 
!		WRITE(*,*) gradQ(1:5,2)
!		WRITE(*,*) gradQ(6:10,2)
!		WRITE(*,*) gradQ(11:15,2)
!		WRITE(*,*) gradQ(16:19,2)
!		PRINT *,'---------------' 
!		WRITE(*,*) BgradQ(1:5)
!		WRITE(*,*) BgradQ(6:10)
!		WRITE(*,*) BgradQ(11:15)
!		WRITE(*,*) BgradQ(16:19)
!		PRINT *,'---------------' 
!		STOP
!		STOP
!	ENDIF 
	!IF( ANY( ISNAN(Q) )) THEN
	!	PRINT *,'Q,  PDENCP NAN ' 
	!	STOP
	!ENDIF
	!IF(ANY(Q(6:8).NE.0.0)) THEN
 !       PRINT *, Q(6:8)
 !       PRINT *, "I feel magnetized :-( NCP"
 !       ERROR STOP
 !   ENDIF
 !   !
	
   continue
   !
END SUBROUTINE PDENCP


RECURSIVE SUBROUTINE PDEEigenvalues(L,Q,n)
  USE Parameters, ONLY :  nVar, nDim 
  USE iso_c_binding
   USE GRMHD_Mod
  IMPLICIT NONE
  REAL :: L(nVar), n(nDim), Q(nVar), Vp(nVar)
  INTENT(IN)  :: Q,n
  INTENT(OUT) :: L 
  ! Local variables
  

  !
  L(:) = 0
  !
#ifdef DEBUGGONE
  L(:)= 1.0
  STOP
  return
#endif
  !
  ! L = 1.0
  CALL PDEEigenvaluesGRMHD(L,Q,n)
  !
	!IF( ANY( ISNAN(L) )) THEN
	!	PRINT *,' PDEEigenvalues NAN ' 
	!	STOP
	!ENDIF 
	!IF(ANY(Q(6:8).NE.0.0)) THEN
 !       PRINT *, Q(6:8)
 !       PRINT *, "I feel magnetized :-( Eigenvalues"
 !       ERROR STOP
 !   ENDIF
    !
  continue
  !
END SUBROUTINE PDEEigenvalues

RECURSIVE SUBROUTINE PDESource(S,Q) 
  USE Parameters, ONLY:  nVar
  USE iso_c_binding
   USE GRMHD_Mod
  IMPLICIT NONE
  ! --------------------------------------------
  ! Argument list declaration
  REAL :: S(nvar), Q(nvar)
  INTENT(IN)  :: Q 
  INTENT(OUT) :: S
  ! --------------------------------------------
  !
  CALL PDESourceGRMHD(S,Q) 
  !S = 0.
  !
	!IF( ANY( S.NE.0.0 )) THEN
	!	PRINT *,' PDESource nonzero ' 
	!	STOP
	!ENDIF   
	!IF( ANY( ISNAN(S) )) THEN
	!	PRINT *,' PDESource NAN ' 
	!	STOP
	!ENDIF 
	!IF(ANY(Q(6:8).NE.0.0)) THEN
 !       PRINT *, Q(6:8)
 !       PRINT *, "I feel magnetized :-( Source"
 !       ERROR STOP
 !   ENDIF
    !
  continue
  !
END SUBROUTINE PDESource

RECURSIVE SUBROUTINE PDEVarName(MyNameOUT,ind) 
  USE Parameters, ONLY: nVar  
  USE iso_c_binding
  USE GRMHD_Mod
  IMPLICIT NONE     
  CHARACTER(LEN=20):: MyName(nVar),MyNameOUT
  INTEGER			:: ind
	!
	CALL PDEVarNameGRMHD(MyName)
	!
	MyNameOUT=MyName(ind+1)
	!
END SUBROUTINE PDEVarName

RECURSIVE SUBROUTINE PDEMatrixB(An,Q,nv) 
  USE Parameters, ONLY : nVar, nDim
  USE iso_c_binding
  USE GRMHD_Mod
  IMPLICIT NONE
  ! Argument list 
  REAL :: An(nVar,nVar)
  REAL :: Q(nVar), nv(nDim) 
  INTENT(IN)  :: Q,nv
  INTENT(OUT) :: An  
  ! Local variables
  ! Linear elasticity variables
  REAL :: A(nVar,nVar), B(nVar,nVar), C(nVar,nVar), Vp(nVar)
  ! 
  An = 0.
#ifdef DEBUGGONE
  STOP
  RETURN 
#endif
  ! 
  ! we use this only for the Roe Matrix
	! --------------------------------------------
	!
	CALL PDEMatrixBGRMHD(An,Q,nv)  
	!
	continue
	!
END SUBROUTINE PDEMatrixB


    

RECURSIVE SUBROUTINE PDEAuxName(AuxName) 
  USE Parameters, ONLY: nAux
  USE iso_c_binding
  USE GRMHD_Mod  
  IMPLICIT NONE     
  CHARACTER(LEN=10):: AuxName(nAux)
  !
  IF(nAux.GT.0) THEN
		CALL PDEAuxNameGRMHD(AuxName)  
  ENDIF
  !
    END SUBROUTINE PDEAuxName     
  
      

RECURSIVE SUBROUTINE PDEAuxVar(aux,V,x) 
  USE Parameters, ONLY : nVar, nAux, d, nDim
  USE iso_c_binding
  USE GRMHD_Mod
  IMPLICIT NONE
  REAL :: aux(nAux), V(nVar) 
  REAL :: time, x(nDim)
  REAL, PARAMETER :: epsilon = 1e-14 
  ! Local Variables 
  !
  IF(nAux.GT.0) THEN
		CALL PDEAuxVarGRMHD(aux,V,x)
  ENDIF
  !
END SUBROUTINE PDEAuxVar     

!
!RECURSIVE SUBROUTINE PDEJacobian(An,Q,gradQ,nv) 
!  USE Parameters, ONLY : EQN, nVar, aom, d,nDim 
!  IMPLICIT NONE
!  ! Argument list 
!  REAL :: An(nVar,nVar)
!  REAL :: Q(nVar), gradQ(nVar,d), nv(d) 
!  INTENT(IN)  :: Q, gradQ, nv
!  INTENT(OUT) :: An  
!  !
!  A = 0. 
!  B = 0. 
!  C = 0.
!  An = 0. 
!  !
!  CALL PDEMatrixBGRMHD(An,Q,nv) 
!  !
!  continue
!  !
!END SUBROUTINE PDE
	

    
 
RECURSIVE SUBROUTINE RoeMatrix(ARoe,QL,QR,nv) 
  USE Parameters, ONLY : nVar, nDim, d,sGP3,wGP3,nGP3
  IMPLICIT NONE
  ! Argument list 
  REAL        :: ARoe(nVar,nVar), QL(nVar), QR(nVar), nv(nDim) 
  INTENT(IN)  :: QL, QR, nv  
  INTENT(OUT) :: ARoe 
  ! Local variables 
  INTEGER             :: iGP  
  REAL                :: psi(nVar) 
  REAL                :: A(nVar,nVar)
  REAL                :: gradQ(nVar,nDim) 
  ! Midpoint Rule 
  !REAL, PARAMETER     :: sGP(1) = (/ 0.5 /) 
  !REAL, PARAMETER     :: wGP(1) = (/ 1.0 /) 
  !INTEGER, PARAMETER :: nGP = 1 
  ! Trapezoidal Rule 
  !REAL, PARAMETER     :: sGP(2) = (/ 0.0, 1.0 /) 
  !REAL, PARAMETER     :: wGP(2) = (/ 0.5, 0.5 /) 
  !INTEGER, PARAMETER :: nGP = 2 
  ! 4-point Gaussian quadrature 
  !REAL, PARAMETER     :: sGP(4) = (/ 1.D0/2.D0-sqrt(525+70*sqrt(30.D0))/70,  1.D0/2.D0-sqrt(525-70*sqrt(30.D0))/70,  1.D0/2.D0+sqrt(525-70*sqrt(30.D0))/70, 1.D0/2.D0+sqrt(525+70*sqrt(30.D0))/70 /) 
  !REAL, PARAMETER     :: wGP(4) = (/ 1.D0/4.D0-sqrt(30.D0)/72, 1.D0/4.D0+sqrt(30.D0)/72,  1.D0/4.D0+sqrt(30.D0)/72,  1.D0/4.D0-sqrt(30.D0)/72 /) 
  !INTEGER, PARAMETER  :: nGP = 4  
  !
  gradQ = 0.0 
  !
  ARoe = 0. 
  DO iGP = 1, nGP3  
     psi = QL + sGP3(iGP)*(QR-QL)  
     CALL PDEMatrixB(A,psi,nv) 
     ARoe = ARoe + wGP3(iGP)*A   ! Numerical integration of the Roe-matrix  
  ENDDO
  !
END SUBROUTINE RoeMatrix 
    
 
RECURSIVE SUBROUTINE HLLEMFluxFV(FL,FR,QL,QR,QavL,QavR,NormalNonZero) 
  USE Parameters, ONLY : nVar, nDim, nLin
  USE iso_c_binding 
  IMPLICIT NONE
  ! Local variables
  INTEGER, INTENT(IN)   :: NormalNonZero
  REAL, INTENT(IN)     :: QL(nVar)
  REAL, INTENT(IN)     :: QR(nVar)
  REAL, INTENT(INOUT)  :: FL(nVar)
  REAL, INTENT(INOUT)  :: FR(nVar)
  REAL    :: QavL(nVar), QavR(nVar)  
    ! Local variables 
  INTEGER           :: i,j,k,l, ml(1)  
  REAL              :: smax, Qav(nVar),sL,sR
  REAL              ::  nv(nDim), flattener(nLin)
  REAL    :: absA(nVar,nVar), amax  
  REAL    :: QM(nVar),LL(nVar),LR(nVar),LM(nVar)
  REAL    :: deltaL(nLin,nLin),Lam(nLin,nLin),Lap(nLin,nLin) 
  REAL    :: RL(nVar,nLin),iRL(nLin,nVar),LLin(nLin,nLin) 
  REAL    :: Aroe(nVar,nVar),Aroep(nVar,nVar), Aroem(nVar,nVar), Dm(nVar), Dp(nVar), dQ(nVar)
  REAL :: f1R(nVar), g1R(nVar), h1R(nVar) 
  REAL :: f1L(nVar), g1L(nVar), h1L(nVar) 
  !  
  nv(:)=0.
  nv(NormalNonZero+1)=1.
  !
  flattener=1.
  !
  CALL PDEFlux(f1L,g1L,h1L,QL)
  CALL PDEFlux(f1R,g1R,h1R,QR)
  !
  fR = f1R*nv(1)+g1R*nv(2)+h1R*nv(3)
  fL = f1L*nv(1)+g1L*nv(2)+h1L*nv(3)
  !
IF(ANY(fR(6:8).NE.0)) THEN
     PRINT *,"f1R",f1R
     PRINT *,"g1R",g1R
     PRINT *,"h1R",h1R
     PRINT *,"f1L",f1L
     PRINT *,"g1L",g1L
     PRINT *,"h1L",h1L
     PRINT *,"dQ",dQ
     PRINT *,"QR",QR
     PRINT *,"QL",QL
     PRINT *,"dQ",dQ
     PRINT *,"BEFORE: fR(6:8).NE.0",fR(6:8)
    STOP
ENDIF
IF(ANY(fl(6:8).NE.0)) THEN
     PRINT *,"QR",QR
     PRINT *,"QL",QL
     PRINT *,"dQ",dQ
     PRINT *,"BEFORE: fl(6:8).NE.0",fl(6:8)
    STOP
ENDIF
  !USE Parameters, ONLY : d,nVar,ndim 
  QM = 0.5*(QL+QR) 
  CALL PDEEigenvalues(LL,QL,nv)  
  CALL PDEEigenvalues(LR,QR,nv)  
IF(ANY(QM(6:8).NE.0)) THEN
    PRINT *, "HLLEMFluxFV QM(6:8)",QM(6:8)
    STOP
ENDIF
  CALL PDEEigenvalues(LM,QM,nv)  
  sL = MIN( 0., MINVAL(LL(:)), MINVAL(LM(:)) ) 
  sR = MAX( 0., MAXVAL(LR(:)), MAXVAL(LM(:)) ) 
 ! PRINT *, "PDEIntermediateFields"
  !DO i=1,nVar
  !  WRITE(*,'(E16.6)'), QM(i)
  !ENDDO
  CALL PDEIntermediateFields(RL,LLin,iRL,QM,nv) 
  !PRINT *, "PDEIntermediateFields finished"
  Lam = 0.5*(LLin-ABS(LLin))
  Lap = 0.5*(LLin+ABS(LLin)) 
  deltaL = 0.0
  DO i = 1, nLin
      deltaL(i,i) = (1. - Lam(i,i)/(sL-1e-14) - Lap(i,i)/(sR+1e-14) )*flattener(i)  
  ENDDO    
#ifdef VISCOUS
  CALL PDEViscEigenvalues(LL,QL,nv)  
  CALL PDEViscEigenvalues(LR,QR,nv)
  amax = 2.0*MAX( MAXVAL(ABS(LL)), MAXVAL(ABS(LR)) )/dist 
#else
  amax = 0. 
#endif 
  absA = 0. 
  DO i = 1, nVar
      absA(i,i) = sR*sL/(sR-sL)  - 0.5*amax ! regular HLL diffusion, only on the diagonal 
  ENDDO  
  !
  IF(QR(1).LT.1e-9.OR.QL(1).LT.1e-9) THEN
      deltaL = 0.
  ENDIF
  !
  !deltaL = 0.
  absA = absA - sR*sL/(sR-sL)*MATMUL( RL, MATMUL(deltaL, iRL) )  ! HLLEM anti-diffusion  
  !    
  CALL RoeMatrix(ARoe,QL,QR,nv)
  !
  ARoem = -sL*ARoe/(sR-sL)
  ARoep = +sR*ARoe/(sR-sL)
  ! 
  !DR = ARoep
  !DL = ARoem
  !!!!FL(:) = 0.5*( FR(:) + FL(:) ) + MATMUL(absA, QR(:) - QL(:) )    ! purely conservative flux 
  !!!!FR(:) = FL(:) - 0.5*ncp(:)                                                              ! subtract the non-conservative product 
  !!!!FL(:) = FL(:) + 0.5*ncp(:)
  !
  dQ = QR - QL
  fL = (sR*fL - sL*fR)/(sR-sL) + MATMUL( absA, dQ ) 
  !
  !Dp = -MATMUL(Aroep,dQ)
  !Dm = -MATMUL(Aroem,dQ)        ! these are the path integral of the MatrixB from QL to QR. (the NCP as a first approximation)
  !
  Dp = MATMUL(Aroep,dQ)
  Dm = MATMUL(Aroem,dQ)        ! these are the path integral of the MatrixB from QL to QR. (the NCP as a first approximation)
  !
  fR = fL - Dp
  fL = fL + Dm
  ! 
IF(ANY(Dp(6:8).NE.0)) THEN
     PRINT *,"QR",QR
     PRINT *,"QL",QL
     PRINT *,"dQ",dQ
     PRINT *,"Dp(6:8).NE.0",Dp(6:8)
    STOP
ENDIF
IF(ANY(Dm(6:8).NE.0)) THEN
     PRINT *,"QR",QR
     PRINT *,"QL",QL
     PRINT *,"dQ",dQ
     PRINT *,"Dm(6:8).NE.0",Dm(6:8)
    STOP
ENDIF
IF(ANY(fR(6:8).NE.0)) THEN
     PRINT *,"QR",QR
     PRINT *,"QL",QL
     PRINT *,"dQ",dQ
     PRINT *,"fR(6:8).NE.0",fR(6:8)
    STOP
ENDIF
IF(ANY(fl(6:8).NE.0)) THEN
     PRINT *,"QR",QR
     PRINT *,"QL",QL
     PRINT *,"dQ",dQ
     PRINT *,"fl(6:8).NE.0",fl(6:8)
    STOP
ENDIF
  ! REMEMBER THE FOLLOWING: we are recursively updating qh as
  ! q_i^{n+1} = q_i^n - FL              .... i.e. F_=F_{i+1/2}_ right flux
  ! q_{i+1}^{n+1} = q_i^n + FR             .... i.e. FR=F_{i+1/2} left flux
  ! see musclhancock.cpph after "// 4. Solve Riemann problems"
  !
    END SUBROUTINE HLLEMFluxFV

    
RECURSIVE SUBROUTINE PDEEigenvectors(R,L,iR,Q,nv) 
  USE Parameters, ONLY: EQN, nVar, nDim, d 
  USE GRMHD_Mod
  IMPLICIT NONE
  ! Argument list 
  REAL :: R(nVar,nVar), L(nVar,nVar), iR(nVar,nVar)
  REAL :: Q(nVar), nv(nDim)
  ! 
IF(ANY(Q(6:8).NE.0)) THEN
    PRINT *, "PDEEigenvectors Q(6:8)",Q(6:8)
    STOP
ENDIF
  !
  CALL PDEEigenVectorsGRMHD(R,L,iR,Q,nv)
  !
  
  END SUBROUTINE PDEEigenvectors

RECURSIVE SUBROUTINE PDEIntermediateFields(RL,LL,iRL,Q,nv) 
  USE Parameters, ONLY: EQN, nVar, nLin, nDim
  IMPLICIT NONE
  ! Argument list 
  REAL :: RL(nVar,nLin), LL(nLin,nLin), iRL(nLin,nVar)
  REAL :: Q(nVar), nv(nDim)
  REAL :: x(nDim), time
  INTENT(IN)  :: Q,nv
  INTENT(OUT) :: RL,LL,iRL 
  ! Local variables
  INTEGER :: i,j,k, zero,iErr, minl(1), maxl(1) 
  REAL :: R(nVar,nVar), L(nVar,nVar), Lv(nVar), iR(nVar,nVar)
  ! 
 ! PRINT *, "HEY, I'm in the IntermediateFields subroutine!" 
  !
  !DO i=1,nVar
  !  WRITE(*,'(a,E16.6)'), "Q(i)",Q(i)
  !ENDDO
  CALL PDEEigenvectors(R,L,iR,Q,nv)
  DO i = 1, nVar
      Lv(i) = L(i,i)
  ENDDO
  minl = MINLOC(Lv) 
  maxl = MAXLOC(Lv)
  !
  i = 0 
  LL = 0. 
  DO j = 1, nVar
      IF( (j.NE.minl(1)).AND.(j.NE.maxl(1)) ) THEN 
          i = i + 1 
          RL(:,i) = R(:,j) 
          iRL(i,:) = iR(j,:) 
          LL(i,i)  = L(j,j) 
      ENDIF      
  ENDDO   
  !
END SUBROUTINE PDEIntermediateFields
    