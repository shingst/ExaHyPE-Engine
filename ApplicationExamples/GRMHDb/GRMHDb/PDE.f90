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





RECURSIVE SUBROUTINE PDEFlux(f,g,h,Q)
  USE Parameters, ONLY : d,nVar,ndim
  USE GRMHD_Mod
  USE iso_c_binding
  IMPLICIT NONE
  REAL :: f(nVar), g(nVar), h(nVar), Q(nVar)
  INTENT(IN)  :: Q
  INTENT(OUT) :: f, g, h
  ! Local varialbes
  INTEGER :: i
  REAL :: FF(nVar,d), V(nVar)
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
   REAL, INTENT(IN)  :: gradQ(nVar, 3)
   REAL, INTENT(IN)  :: Q(nVar)
   ! Local variables
   REAL :: TMPgradQ(nVar, 3)
   INTEGER :: i
   !REAL :: u(3),VP(nVar) 
   !REAL :: AQx(nVar), BQy(nVar), CQz(nVar) , Qx(nVar), Qy(nVar), Qz(nVar)
   !
  !PRINT *, 'PDENCPGRMHD'
#if defined(Dim2)
	DO i =1,nVar
		TMPgradQ(i,1) = gradQ(i,1)
		TMPgradQ(i,2) = gradQ(i,2)
		TMPgradQ(i,3) = 0.0
	ENDDO
#else
	DO i =1,nVar
		TMPgradQ(i,1) = gradQ(i,1)
		TMPgradQ(i,2) = gradQ(i,2)
		TMPgradQ(i,3) = gradQ(i,3)
	ENDDO
#endif
	!
   CALL PDENCPGRMHD(BgradQ,Q,TMPgradQ)
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
  REAL :: L(nVar), n(3), Q(nVar), Vp(nVar)
  INTENT(IN)  :: Q,n
  INTENT(OUT) :: L 
  ! Local variables
  REAL :: TMPn(3)
  !
  L(:) = 0
  !
#if defined(Dim2)
	TMPn(1) = n(1)
	TMPn(2) = n(2)
	TMPn(3) = 0.0
#else
	TMPn(1) = n(1)
	TMPn(2) = n(2)
	TMPn(3) = n(3)
#endif
  !PRINT *, 'PDEEigenvaluesGRMHD'
  ! WRITE(*,*) Q(1)
  ! WRITE(*,*) Q(2:4)
  ! WRITE(*,*) Q(5)
  ! WRITE(*,*) Q(6:8)
  ! WRITE(*,*) Q(9)
  ! WRITE(*,*) Q(10)
  ! WRITE(*,*) Q(11:13)
  ! WRITE(*,*) Q(14:19)
  ! L = 1.0
  CALL PDEEigenvaluesGRMHD(L,Q,TMPn)
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
  REAL :: Q(nVar), nv(3) 
  INTENT(IN)  :: Q,nv
  INTENT(OUT) :: An  
  ! Local variables
  ! Linear elasticity variables
  REAL :: A(nVar,nVar), B(nVar,nVar), C(nVar,nVar), Vp(nVar)
  !
  PRINT *, ' Impossible error! ' 
  STOP
	! --------------------------------------------
	!
	!CALL PDEMatrixBGRMHD(An,Q,nv)  
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
  USE Parameters, ONLY : nVar, nAux, d
  USE iso_c_binding
  USE GRMHD_Mod
  IMPLICIT NONE
  REAL :: aux(nAux), V(nVar) 
  REAL :: time, x(d)
  REAL, PARAMETER :: epsilon = 1e-14 
  ! Local Variables 
  !
  IF(nAux.GT.0) THEN
		CALL PDEAuxVarGRMHD(aux,V,x)
  ENDIF
  !
END SUBROUTINE PDEAuxVar     



