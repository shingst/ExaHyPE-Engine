MODULE TECPLOTPLOTTERmod
!USE Parameters, ONLY: nVar, nDim, N, M
IMPLICIT NONE
! Problem standard parameters
INTEGER :: MC, M,N
INTEGER :: nGPMC
INTEGER :: nDOFm
INTEGER :: nSubLim
INTEGER :: nSubLimV(3)
INTEGER :: nGPM, nGPVM(3)
INTEGER :: nFace, nVertex
REAL, ALLOCATABLE 	:: xiGPM(:), wGPM(:) 
REAL, ALLOCATABLE :: MPoly1D(:,:)
! 
REAL, ALLOCATABLE :: SubOutputMatrix(:,:), SubGradOutputMatrix(:,:,:)
!
! ------------------------------------------
! Variables needed for the tecplot plotter
!
INTEGER					:: nSubPlotElem,nRealNodes,nPlotElem
INTEGER					:: nElem_max,nRealNodes_max, nSubPlotElem_max
INTEGER, PARAMETER		:: td=4
INTEGER*4	, POINTER 	:: NData(:,:),NData_max(:,:)
REAL(td)	, POINTER  	:: DataArray(:,:), DataArray_max(:,:)
INTEGER, ALLOCATABLE    :: subtri(:,:)
REAL, ALLOCATABLE       :: allsubxi(:,:)
INTEGER :: Element_c, Element_nc 
INTEGER :: PlotIndex
REAL	:: PLOT_TIME
!
! ------------------------------------------
private :: MatrixInverse
public :: SetMainParameters
public :: ComputeOutputMatrices
private :: MakeMCPoly1D
private :: MBaseFunc1D
CONTAINS

SUBROUTINE InitTECPLOTPLOTTER(time_in)
    USE Parameters, ONLY: nDim,nVar,nAux
	IMPLICIT NONE
	REAL, INTENT(IN) :: time_in
	ALLOCATE(NData_max(nVertex,nSubPlotElem_max))  
	ALLOCATE(DataArray_max(nRealNodes_max,nDim+nVar+nAux+1+1))
	
	PlotIndex=PlotIndex+1
	PLOT_TIME=time_in
	nSubPlotElem=0
	nRealNodes=0
	nPlotElem=0
	Element_c=0 
	Element_nc=0
	
END SUBROUTINE InitTECPLOTPLOTTER

SUBROUTINE ElementTECPLOTPLOTTER(wh,lx0,ldx,limiter)
    USE Parameters, ONLY: nDim,nVar,nAux
	IMPLICIT NONE
	REAL, INTENT(IN) :: wh(nVar,nDOFm)
	INTEGER :: nSubNodes,J
	REAL    :: LocNode(nVar,(M+1)**nDim),xvec(3),lx0(3),ldx(3)
	REAL	:: VN(nVar),QN(nVar)
	integer :: limiter
	
	nPlotElem = nPlotElem + 1
	nSubPlotElem = nSubPlotElem + M**nDim  
	nSubNodes = (M+1)**nDim  
	nRealNodes = nRealNodes + nSubNodes
	
	Element_c = Element_c + 1

	!print *,'Element,',nPlotElem, '->', lx0(1:nDim),'dx=',ldx(1:nDim)
	
	DO j = 1, M**nDim
		NData_max(1:nVertex,(Element_c-1)*M**nDim+j) = (Element_c-1)*(M+1)**nDim + subtri(1:nVertex,j)
	END DO

	LocNode = MATMUL( wh(:,1:nDOFm), SubOutputMatrix(1:nDOFm,1:(M+1)**nDim) )
	DO j = 1, (M+1)**nDim  
	   QN(:) = LocNode(:,j)
	   xvec = lx0 + allsubxi(:,j)*ldx
	   CALL PDECons2Prim(VN,QN)
	   !print *, 'Node of elem=',Element_c,'=', xvec(1:nDim)
	    DataArray_max((Element_c-1)*(M+1)**nDim+j,:) = (/ xvec(1:nDim), VN, REAL(nPlotElem), REAL(limiter) /)  
	   !DataArray_max(Element_c,:) = (/ xvec(1:nDim), VN /) 	   
	END DO
		
	Element_nc = Element_nc + (M+1)**nDim
END SUBROUTINE ElementTECPLOTPLOTTER

SUBROUTINE FinalizeTECPLOTPLOTTER(Myrank)
    USE Parameters, ONLY: nDim,nVar,nAux
	USE teciomod
	USE ISO_C_BINDING
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: Myrank 	! CPU rank
	CHARACTER(LEN=200) :: Filename,ZoneTitle,Title,ScratchDir, BaseFile ! BaseFile is the folder	where place the plots
	CHARACTER(LEN=1000) :: VarString
	CHARACTER(LEN=10)  :: cmyrank,varname , AuxName
	INTEGER				:: ZoneType, iRet, StrandID,visdouble,i
	REAL				:: loctime
	REAL(td)           :: Test
	!POINTER(NullPtr,Null)
	!Integer*4 MyNull(*)
	integer*4,dimension(:), POINTER :: NullPtr => NULL ()
	
	SELECT CASE(KIND(Test))
	CASE(4)
		visdouble = 0
	CASE(8)
		visdouble = 1 
	CASE DEFAULT
		PRINT *, ' IO Kind error. ' 
		STOP 
	END SELECT 
	
	WRITE(cmyrank,'(I5.5)') myrank
	BaseFile="./output/TECPLOT"
	WRITE(FileName,'(a,a2,i1.1,a1,i1.1,a1,i8.8,a1,a,a)') TRIM(BaseFile),'-P',N,'P',M,'-',PlotIndex,'-',TRIM(cmyrank),'.plt'
!return
	!NullPtr = 0 
	
	! Now I know the number of elements, associate the proper Data
	ALLOCATE(NData(nVertex,nSubPlotElem))  
	ALLOCATE(DataArray(nRealNodes,nDim+nVar+nAux+1+1))
	
	! **********************************
	print *, "****************************************************"
	print *, "*********** TECPLOT PLOTTER INFO *******************"
	PRINT *, "Myrank		=",myrank
	print *, "NElem			=",nSubPlotElem/M**nDim
	print *, "nSubPlotElem	=",nSubPlotElem
	print *, "nRealNodes	=",nRealNodes
	print *, "****************************************************"
	

	DataArray=DataArray_max(1:nRealNodes,1:(nDim+nVar+nAux+1+1))
	
	 NData=NData_max(1:nVertex,1:nSubPlotElem)
	!do i=1,nSubPlotElem
	!	print *, "NData=",NData_max(:,i)
	!end do
	!do i=1,nRealNodes
	!	print *,"Vert:",i,DataArray_max(i,1:2)
	!end do
	DEALLOCATE(NData_max,DataArray_max)
	
	WRITE(Title,'(a,f9.4,a)') 'Time t = ', PLOT_TIME, ''//C_NULL_CHAR  
	WRITE(ScratchDir,'(a)') '.'//C_NULL_CHAR 
	SELECT CASE(nDim)
	CASE(2)
	WRITE(VarString,*) 'x y ' 
	ZoneType = 3 ! FEM Quad  
	CASE(3)
	WRITE(VarString,*) 'x y z ' 
	ZoneType = 5 ! FEM Brick 
	END SELECT 
	  
	DO i = 0, nVar-1
		CALL PDEVarName(varname,i)	
		WRITE(VarString,'(a,a,a,a)') TRIM(VarString), ' ', TRIM(varname) , ' '   
	ENDDO
	  
	!DO i = 1, nAux 
		!CALL PDEAuxName(AuxName,i)
		!WRITE(VarString,'(a,a,a,a)') TRIM(VarString), ' ', TRIM(AuxName()) , ' '   
	!ENDDO
	WRITE(VarString,'(a,a)') TRIM(VarString), ' iE lim ' 

	iret = TecIni112(TRIM(Title)//''//C_NULL_CHAR,TRIM(Varstring)//''//C_NULL_CHAR,TRIM(FileName)//''//C_NULL_CHAR,TRIM(ScratchDir)//''//C_NULL_CHAR,0,0,visdouble) 
	loctime = PLOT_TIME 
	ZoneTitle = Title
	StrandID=0
	iRet = TecZne112(TRIM(ZoneTitle)//C_NULL_CHAR, ZoneType, nRealNodes, nSubPlotElem, 0, 0, 0, 0, loctime, StrandID, 0, 1, 0, 0, 0, 0, 0, NullPtr, NullPtr, NullPtr, 0) 
	iRet = TecDat112( nRealNodes*(nDim+nVar+nAux+1+1), DataArray, visdouble )	
	iRet = TecNod112(NData)
	iRet = TecEnd112()
	!
	
	DEALLOCATE(NData,DataArray)
	!stop
END SUBROUTINE FinalizeTECPLOTPLOTTER

SUBROUTINE SetMainParameters(N_in,M_in)
	USE Parameters, ONLY: nDim
	USE recipies_mod
	IMPLICIT NONE
	INTEGER,INTENT(IN)	:: N_in,M_in
	
	M=M_in
	N=N_in
	MC = M+1
	nGPMC = MC + 1
	nDOFm = (M+1)**nDim
	nSubLim = 2*N+1
	nGPM  = M + 1
	if(nDim .eq. 3) then
		nSubLimV(1) = nSubLim
		nSubLimV(2) = nSubLim
		nSubLimV(3) = nSubLim
		nGPVM(1)=nGPM
		nGPVM(2)=nGPM
		nGPVM(3)=nGPM
	else
		nSubLimV(1) = nSubLim
		nSubLimV(2) = nSubLim
		nSubLimV(3) = 1
		nGPVM(1)=nGPM
		nGPVM(2)=nGPM
		nGPVM(3)=1
	end if
	nFace = 2*nDim
	nVertex = 2**nDim
	
	PlotIndex=0

	
	
	nElem_max=100000
	nSubPlotElem_max=nElem_max*M**nDim
	nRealNodes_max=nElem_max*(M+1)**nDim
	! Initialize matries and gauss points
	ALLOCATE(MPoly1D(0:MC,MC+1))
	
	ALLOCATE(xiGPM(nGPM), wGPM(nGPM) )
	
	
	CALL gauleg(0.,1.,xiGPM,wGPM,M+1)
	
	CALL MakeMCPoly1D()
	
	CALL ComputeOutputMatrices()
	
	print *, "***********************************************************"
	print *, "************ TECPLOT INITIALIZATION ***********************"
	PRINT *, "N,M=",N,M
	print *, "(nElem,nSubPlotElem,nRealNodes)_MAX=",nElem_max,nSubPlotElem_max,nRealNodes_max
	print *, "***********************************************************"
END SUBROUTINE SetMainParameters

SUBROUTINE ComputeOutputMatrices
	USE Parameters, ONLY: nDim
	IMPLICIT NONE
	! Local variables
	real :: psi_i(nGPM),psi_j(nGPM),psi_k(nGPM),psi_xi(nGPM),psi_xj(nGPM),psi_xk(nGPM),subxi(MC)
	integer :: cnt,i,j,k,kk, jj, ii, counter,c
	real	:: aux(3)
	integer, allocatable :: idxn(:,:,:)
	
	ALLOCATE(SubOutputMatrix((M+1)**3,(M+1)**3), SubGradOutputMatrix((M+1)**3,(M+1)**3,3)) ! First allocate the Outputmatrices
	
	DO i = 1, M+1 
		subxi(i) = REAL(i-1)/REAL(M) 
	ENDDO
	
	cnt = 0
     DO k = 1, M+1
        DO j = 1, M+1 
           DO i = 1, M+1  
              cnt = cnt + 1 
              CALL MBaseFunc1D(psi_i,psi_xi,subxi(i))
              CALL MBaseFunc1D(psi_j,psi_xj,subxi(j))
              CALL MBaseFunc1D(psi_k,psi_xk,subxi(k))
              counter = 0 
              DO kk = 1, nGPVM(3)  
                 DO jj = 1, nGPVM(2)  
                    DO ii = 1, nGPVM(1) 
                       counter = counter + 1 
					   aux(1)=psi_i(ii)
					   aux(2)=psi_j(jj)
					   aux(3)=psi_k(kk)
                       SubOutputMatrix(counter,cnt) = PRODUCT( aux(1:nDim) ) 
					   aux(1)=psi_xi(ii)
					   aux(2)=psi_j(jj)
					   aux(3)=psi_k(kk) 
                       SubGradOutputMatrix(counter,cnt,1) = PRODUCT( aux(1:nDim) )
					   aux(1)=psi_i(ii)
					   aux(2)=psi_xj(jj)
					   aux(3)=psi_k(kk)
                       SubGradOutputMatrix(counter,cnt,2) = PRODUCT( aux(1:nDim) )
					   aux(1)=psi_i(ii)
					   aux(2)=psi_j(jj)
					   aux(3)=psi_xk(kk) 
                       SubGradOutputMatrix(counter,cnt,3) = PRODUCT( aux(1:nDim) )
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
 			  
	! Compute subtri
	ALLOCATE( idxn(M+1,M+1,M+1),subtri(8,M**3),allsubxi(3,(M+1)**3) )
     idxn = 0 
     c = 0 
     DO k = 1, M+1
        DO j = 1, M+1 
           DO i = 1, M+1 
              c = c + 1 
              idxn(i,j,k) = c
			  allsubxi(1,c) = REAL(i-1)/REAL(M)
			  allsubxi(2,c) = REAL(j-1)/REAL(M)
			  allsubxi(3,c) = REAL(k-1)/REAL(M)			  
           ENDDO
        ENDDO
     ENDDO

     c = 0 
     DO k = 1, M 
        DO j = 1, M
           DO i = 1, M 
              c = c + 1 
			  subtri(1,c) =idxn(i,j,k)
			  subtri(2,c) =idxn(i+1,j,k)
			  subtri(3,c) =idxn(i+1,j+1,k)
			  subtri(4,c) =idxn(i,j+1,k)
			  subtri(5,c) =idxn(i,j,k+1)
			  subtri(6,c) =idxn(i+1,j,k+1)
			  subtri(7,c) =idxn(i+1,j+1,k+1)
			  subtri(8,c) =idxn(i,j+1,k+1)        
           ENDDO
        ENDDO
     ENDDO

     DEALLOCATE( idxn ) 
	 
END SUBROUTINE ComputeOutputMatrices

SUBROUTINE MakeMCPoly1D()
	IMPLICIT NONE
	REAL, ALLOCATABLE :: CoeffMat(:,:),iCoeffMat(:,:)
	INTEGER :: iDegFr, ii
  ALLOCATE( CoeffMat(nGPM,nGPM)  )  
  ALLOCATE( iCoeffMat(nGPM,nGPM) )
  CoeffMat  = 0.
  DO iDegFr = 1, nGPM 
     DO ii = 0, M  
        CoeffMat(iDegFr,ii+1) = xiGPM(iDegFr)**ii  
     ENDDO
  ENDDO
  CALL MatrixInverse(nGPM,CoeffMat,iCoeffMat)
  DO iDegFr = 1, nGPM  
     DO ii = 0, M    
        MPoly1D(ii,iDegFr) = iCoeffMat(ii+1,iDegFr)
     ENDDO
  ENDDO
  DEALLOCATE( CoeffMat  )  
  DEALLOCATE( iCoeffMat )
END SUBROUTINE MakeMCPoly1D

SUBROUTINE MBaseFunc1D(psi,psi_xi,xi) 
	IMPLICIT NONE     
	REAL             :: psi(M+1), psi_xi(M+1) 
	REAL             :: xi
	REAL             :: xipower(0:M)  
	INTEGER          :: i 
	!
	! 1D Lagrange basis function associated with the Gauss-Legendre points 
	psi = 0 
	xipower(0) = 1. 
	DO i = 1, M
		xipower(i) = xipower(i-1)*xi 
	ENDDO
	DO i = 0, M
		psi = psi + MPoly1D(i,:)*xipower(i)  
	ENDDO
	! Derivative of the basis functions 
	psi_xi = 0  
	DO i = 1, M
		psi_xi = psi_xi + i*MPoly1D(i,:)*xipower(i-1)  
	ENDDO
	!  
END SUBROUTINE MBaseFunc1D

SUBROUTINE MatrixInverse(NN,A,iA)
  IMPLICIT NONE
  INTEGER       :: NN
  REAL          :: A(NN,NN), iA(NN,NN)
  !
  INTEGER       :: i,j,flag,ml(1) 
  REAL          :: piv
  REAL          :: temp(2*NN)
  REAL          :: C(2*NN,NN)  
  !
  C(1:NN,:)     = TRANSPOSE(A(:,:))
  C(NN+1:2*NN,:) = 0. 
  DO i = 1, NN
     C(NN+i,i) = 1.
  ENDDO
  !    
  ! Forward elimination and row swapping (if necessary)
  ! 
  DO i = 1, NN
     ! If pivot element is zero, then swap rows 
     ml = MAXLOC(ABS(C(i,i:NN))) 
     j = i - 1 + ml(1) 
     temp   = C(:,j) 
     C(:,j) = C(:,i)
     C(:,i) = temp      
     IF(C(i,i).EQ.0.) THEN
        PRINT *, 'ERROR. Matrix is singular!'
        !DO j = 1, N
        !   PRINT *, A(j,:) 
        !ENDDO
        STOP
     ENDIF
     piv    = 1./C(i,i)
     C(:,i) = C(:,i)*piv 
     DO j = i+1, NN 
        C(:,j) = C(:,j) - C(i,j)*C(:,i)
     ENDDO
  ENDDO
  !
  ! Back substitution
  !
  DO i = NN,1,-1   
     DO j = i-1,1,-1
        C(:,j) = C(:,j) - C(i,j)*C(:,i)
     ENDDO
  ENDDO
  !
  iA = TRANSPOSE( C(NN+1:2*NN,:) ) 
  !
END SUBROUTINE MatrixInverse
 

END MODULE TECPLOTPLOTTERmod