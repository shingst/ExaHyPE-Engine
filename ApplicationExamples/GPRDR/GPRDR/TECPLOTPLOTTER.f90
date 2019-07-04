#ifndef HEADER_TECPLOTPLOTTERMOD
#define HEADER_TECPLOTPLOTTERMOD

#include "MainVariables.f90"
MODULE TECPLOTPLOTTERmod

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
private :: gauleg

INTERFACE 


SUBROUTINE tecforeign112 &
(OutputForeignByteOrder)
!MS$ATTRIBUTES STDCALL :: tecforeign112
!MS$ATTRIBUTES REFERENCE :: OutputForeignByteOrder
INTEGER(4) OutputForeignByteOrder
END SUBROUTINE tecforeign112

INTEGER(4) FUNCTION tecini112 &
(Title, &
Variables, &
FName, &
ScratchDir, &
FileType, &
Debug, &
VIsDouble)
!MS$ATTRIBUTES STDCALL :: tecini112
!MS$ATTRIBUTES REFERENCE :: Title,Variables,FName,ScratchDir
!MS$ATTRIBUTES REFERENCE :: FileType,Debug,VIsDouble
CHARACTER(LEN=*) Title
CHARACTER(LEN=*) Variables
CHARACTER(LEN=*) FName
CHARACTER(LEN=*) ScratchDir
INTEGER(4)       FileType
INTEGER(4)       Debug
INTEGER(4)       VIsDouble
END FUNCTION tecini112

INTEGER(4) FUNCTION teczne112 &
(ZoneTitle, &
ZoneType, &
IMxOrNumPts, &
JMxOrNumElements, &
KMxOrNumFaces, &
ICellMax, &
JCellMax, &
KCellMax, &
SolutionTime, &
StrandID, &
ParentZone, &
IsBlock, &
NumFaceConnections, &
FaceNeighborMode, &
TotalNumFaceNodes, &
NumConnectedBoundaryFaces, &
TotalNumBoundaryConnections, &
PassiveVarList, &
ValueLocation, &
ShareVarFromZone, &
ShareConnectivityFromZone)
!MS$ATTRIBUTES STDCALL :: teczne112
!MS$ATTRIBUTES REFERENCE :: ZoneTitle,ZoneType,IMxOrNumPts
!MS$ATTRIBUTES REFERENCE :: JMxOrNumElements,KMxOrNumFaces
!MS$ATTRIBUTES REFERENCE :: ICellMax,JCellMax,KCellMax
!MS$ATTRIBUTES REFERENCE :: SolutionTime,StrandID,ParentZone
!MS$ATTRIBUTES REFERENCE :: IsBlock,PassiveVarList
!MS$ATTRIBUTES REFERENCE :: NumFaceConnections,FaceNeighborMode
!MS$ATTRIBUTES REFERENCE :: TotalNumFaceNodes
!MS$ATTRIBUTES REFERENCE :: NumConnectedBoundaryFaces
!MS$ATTRIBUTES REFERENCE :: TotalNumBoundaryConnections
!MS$ATTRIBUTES REFERENCE :: ValueLocation,ShareVarFromZone
!MS$ATTRIBUTES REFERENCE :: ShareConnectivityFromZone
CHARACTER(LEN=*) ZoneTitle
INTEGER(4)       ZoneType
INTEGER(4)       IMxOrNumPts
INTEGER(4)       JMxOrNumElements
INTEGER(4)       KMxOrNumFaces
INTEGER(4)       ICellMax
INTEGER(4)       JCellMax
INTEGER(4)       KCellMax
REAL(8)          SolutionTime
INTEGER(4)       StrandID
INTEGER(4)       ParentZone
INTEGER(4)       IsBlock
INTEGER(4)       NumFaceConnections
INTEGER(4)       FaceNeighborMode
INTEGER(4)       TotalNumFaceNodes
INTEGER(4)       NumConnectedBoundaryFaces
INTEGER(4)       TotalNumBoundaryConnections
INTEGER(4)       PassiveVarList(*)
INTEGER(4)       ValueLocation(*)
INTEGER(4)       ShareVarFromZone(*)
INTEGER(4)       ShareConnectivityFromZone
END FUNCTION teczne112

#ifdef TECPLOTDOUBLE 
INTEGER(4) FUNCTION tecdat112 &
(N, &
FieldData, &
IsDouble)
!MS$ATTRIBUTES STDCALL :: tecdat112
!MS$ATTRIBUTES REFERENCE :: N,FieldData,IsDouble
INTEGER(4)  N
REAL(8)     FieldData(*)
INTEGER(4)  IsDouble
END FUNCTION tecdat112
#else
INTEGER(4) FUNCTION tecdat112 &
(N, &
FieldData, &
IsDouble)
!MS$ATTRIBUTES STDCALL :: tecdat112
!MS$ATTRIBUTES REFERENCE :: N,FieldData,IsDouble
INTEGER(4)  N
REAL(4)     FieldData(*)
INTEGER(4)  IsDouble
END FUNCTION tecdat112
#endif

INTEGER(4) FUNCTION tecnod112 &
(NData)
!MS$ATTRIBUTES STDCALL :: tecnod112
!MS$ATTRIBUTES REFERENCE :: NData
INTEGER(4)  NData(*)
END FUNCTION tecnod112

INTEGER(4) FUNCTION tecgeo112 &
(XPos, &
YPos, &
ZPos, &
PosCoordMode, &
AttachToZone, &
Zone, &
Color, &
FillColor, &
IsFilled, &
GeomType, &
LinePattern, &
PatternLength, &
LineThickness, &
NumEllipsePts, &
ArrowheadStyle, &
ArrowheadAttachment, &
ArrowheadSize, &
ArrowheadAngle, &
Scope, &
Clipping, &
NumSegments, &
NumSegPts, &
XGeomData, &
YGeomData, &
ZGeomData, &
mfc)
!MS$ATTRIBUTES STDCALL :: tecgeo112
!MS$ATTRIBUTES REFERENCE :: XPos,YPos,ZPos,PosCoordMode
!MS$ATTRIBUTES REFERENCE :: AttachToZone,Zone,Color,FillColor
!MS$ATTRIBUTES REFERENCE :: IsFilled,GeomType,LinePattern
!MS$ATTRIBUTES REFERENCE :: PatternLength,LineThickness
!MS$ATTRIBUTES REFERENCE :: NumEllipsePts,ArrowheadStyle
!MS$ATTRIBUTES REFERENCE :: ArrowheadAttachment,ArrowheadSize
!MS$ATTRIBUTES REFERENCE :: ArrowheadAngle,Scope,Clipping
!MS$ATTRIBUTES REFERENCE :: NumSegments,NumSegPts
!MS$ATTRIBUTES REFERENCE :: XGeomData,YGeomData
!MS$ATTRIBUTES REFERENCE :: ZGeomData,mfc
REAL(8)        XPos
REAL(8)        YPos
REAL(8)        ZPos
INTEGER(4)     PosCoordMode
INTEGER(4)     AttachToZone
INTEGER(4)     Zone
INTEGER(4)     Color
INTEGER(4)     FillColor
INTEGER(4)     IsFilled
INTEGER(4)     GeomType
INTEGER(4)     LinePattern
REAL(8)        PatternLength
REAL(8)        LineThickness
INTEGER(4)     NumEllipsePts
INTEGER(4)     ArrowheadStyle
INTEGER(4)     ArrowheadAttachment
REAL(8)        ArrowheadSize
REAL(8)        ArrowheadAngle
INTEGER(4)     Scope
INTEGER(4)     Clipping
INTEGER(4)     NumSegments
INTEGER(4)     NumSegPts(*)
REAL(4)        XGeomData(*)
REAL(4)        YGeomData(*)
REAL(4)        ZGeomData(*)
character(len=*) mfc
END FUNCTION tecgeo112

INTEGER(4) FUNCTION tectxt112 &
(XOrThetaPos, &
YOrRPos, &
ZOrUnusedPos, &
PosCoordMode, &
AttachToZone, &
Zone, &
Font, &
FontHeightUnits, &
FontHeight, &
BoxType, &
BoxMargin, &
BoxLineThickness, &
BoxColor, &
BoxFillColor, &
Angle, &
Anchor, &
LineSpacing, &
TextColor, &
Scope, &
Clipping, &
Text, &
mfc)
!MS$ATTRIBUTES STDCALL :: tectxt112
!MS$ATTRIBUTES REFERENCE :: XOrThetaPos,YOrRPos
!MS$ATTRIBUTES REFERENCE :: ZOrUnusedPos,PosCoordMode
!MS$ATTRIBUTES REFERENCE :: AttachToZone,Zone,Font
!MS$ATTRIBUTES REFERENCE :: FontHeightUnits
!MS$ATTRIBUTES REFERENCE :: FontHeight,BoxType,BoxMargin
!MS$ATTRIBUTES REFERENCE :: BoxLineThickness,BoxColor
!MS$ATTRIBUTES REFERENCE :: BoxFillColor,Angle,Anchor
!MS$ATTRIBUTES REFERENCE :: LineSpacing,TextColor,Scope,Clipping
!MS$ATTRIBUTES REFERENCE :: Text,mfc
REAL(8)          XOrThetaPos
REAL(8)          YOrRPos
REAL(8)          ZOrUnusedPos
INTEGER(4)       PosCoordMode
INTEGER(4)       AttachToZone
INTEGER(4)       Zone
INTEGER(4)       Font
INTEGER(4)       FontHeightUnits
REAL(8)          FontHeight
INTEGER(4)       BoxType
REAL(8)          BoxMargin
REAL(8)          BoxLineThickness
INTEGER(4)       BoxColor
INTEGER(4)       BoxFillColor
REAL(8)          Angle
INTEGER(4)       Anchor
REAL(8)          LineSpacing
INTEGER(4)       TextColor
INTEGER(4)       Scope
INTEGER(4)       Clipping
CHARACTER(LEN=*) Text
CHARACTER(LEN=*) mfc
END FUNCTION tectxt112

INTEGER(4) FUNCTION teclab112 &
(S)
!MS$ATTRIBUTES STDCALL :: teclab112
!MS$ATTRIBUTES REFERENCE :: S
character(len=*) S
END FUNCTION teclab112

INTEGER(4) FUNCTION tecfil112 &
(F)
!MS$ATTRIBUTES STDCALL :: tecfil112
!MS$ATTRIBUTES REFERENCE :: F
INTEGER(4)  F
END FUNCTION tecfil112

INTEGER(4) FUNCTION tecend112()
!MS$ATTRIBUTES STDCALL :: tecend112
END FUNCTION tecend112

INTEGER(4) FUNCTION tecusr112 &
(S)
!MS$ATTRIBUTES STDCALL :: tecusr112
!MS$ATTRIBUTES REFERENCE :: S
character(len=*) S
END FUNCTION tecusr112

INTEGER(4) FUNCTION tecauxstr112 &
(Name, &
Value)
!MS$ATTRIBUTES STDCALL :: tecauxstr112
!MS$ATTRIBUTES REFERENCE :: Name, Value
CHARACTER(LEN=*) Name 
CHARACTER(LEN=*) Value
END FUNCTION tecauxstr112

INTEGER(4) FUNCTION teczauxstr112 &
(Name, &
Value)
!MS$ATTRIBUTES STDCALL :: teczauxstr112
!MS$ATTRIBUTES REFERENCE :: Name, Value
CHARACTER(LEN=*) Name 
CHARACTER(LEN=*) Value
END FUNCTION teczauxstr112

INTEGER(4) FUNCTION tecvauxstr112 &
(Name, &
Value)
!MS$ATTRIBUTES STDCALL :: tecvauxstr112
!MS$ATTRIBUTES REFERENCE :: Name, Value
CHARACTER(LEN=*) Name 
CHARACTER(LEN=*) Value
END FUNCTION tecvauxstr112

INTEGER(4) FUNCTION tecface112 &
(FaceConnections)
!MS$ATTRIBUTES STDCALL :: tecface112
!MS$ATTRIBUTES REFERENCE :: FaceConnections
INTEGER(4) FaceConnections(*)
END FUNCTION tecface112

INTEGER(4) FUNCTION tecpoly112 &
(FaceNodeCounts, &
FaceNodes, &
FaceLeftElems, &
FaceRightElems, &
FaceBndryConnectionCounts, &
FaceBndryConnectionElems, &
FaceBndryConnectionZones)
!MS$ATTRIBUTES STDCALL   :: tecpoly112
!MS$ATTRIBUTES REFERENCE :: FaceNodes
!MS$ATTRIBUTES REFERENCE :: FaceLeftElems
!MS$ATTRIBUTES REFERENCE :: FaceRightElems
!MS$ATTRIBUTES REFERENCE :: FaceBndryConnectionCounts
!MS$ATTRIBUTES REFERENCE :: FaceBndryConnectionElems
!MS$ATTRIBUTES REFERENCE :: FaceBndryConnectionZones
INTEGER(4) FaceNodeCounts(*)
INTEGER(4) FaceNodes(*)
INTEGER(4) FaceLeftElems(*)
INTEGER(4) FaceRightElems(*)
INTEGER(4) FaceBndryConnectionCounts(*)
INTEGER(4) FaceBndryConnectionElems(*)
INTEGER(2) FaceBndryConnectionZones(*)
END FUNCTION tecpoly112


END INTERFACE

CONTAINS

SUBROUTINE InitTECPLOTPLOTTER(time_in)
    USE MainVariables, ONLY: nDim,nVar,nAux
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
    USE MainVariables, ONLY: nDim,nVar,nAux
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
    USE MainVariables, ONLY: nDim,nVar,nAux
	!USE teciomod
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
	USE MainVariables, ONLY: nDim
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
	USE MainVariables, ONLY: nDim
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
 
SUBROUTINE gauleg(x1,x2,x,w,n)
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  !--------------------------------------------------------------------------
  INTEGER     ::  n
  REAL        :: x1,x2,x(n),w(n)
  REAL        :: EPS
  INTEGER     :: i,j,m
  REAL        :: p1,p2,p3,pp,xl,xm,z,z1
  !--------------------------------------------------------------------------
  INTENT(IN)  :: x1,x2,n
  INTENT(OUT) :: x,w
  !--------------------------------------------------------------------------
  PARAMETER (EPS=3.E-14)
  !--------------------------------------------------------------------------

  m  = (n+1)/2
  xm = 0.5*(x2+x1)
  xl = 0.5*(x2-x1)

  DO i=1,m
     z = COS(3.141592654*(i-.25)/(n+.5))
1 CONTINUE
     p1 = 1.
     p2 = 0.
     DO j = 1,n
        p3 = p2
        p2 = p1
        p1 = ((2.*j-1.)*z*p2-(j-1.)*p3)/j
     END DO
     pp = n*(z*p1-p2)/(z*z-1.)
     z1 = z
     z  = z1-p1/pp
     IF(ABS(z-z1).GT.EPS)GOTO 1
     x(i)    = xm-xl*z
     x(n+1-i)= xm+xl*z
     w(i)    = 2.*xl/((1.-z*z)*pp*pp)
     w(n+1-i)= w(i)
  END DO
  RETURN
END SUBROUTINE gauleg
 

END MODULE TECPLOTPLOTTERmod
#endif
