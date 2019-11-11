! HERE every input variables refere to the "nDim" total space dimensions:
! at the moment, the externally defined TECPLOT routines are referred to 3D space dimensions.
    
    
RECURSIVE SUBROUTINE ElementCallTECPLOTPLOTTER(wh,lx0,ldx,limiter)
	USE, INTRINSIC :: ISO_C_BINDING
	USE TECPLOTPLOTTERmod !, only : ElementTECPLOTPLOTTER,nDOFm
	USE Parameters, only : nVar, nDim
	implicit none
	REAL, INTENT(IN) :: wh(nVar,nDOFm),lx0(nDim),ldx(nDim)
	real :: lx0_3(3),ldx_3(3)
	integer :: limiter
        lx0_3 = 0.
        ldx_3 = 0.

	lx0_3(1:nDim)=lx0(1:nDim)
	ldx_3(1:nDim)=ldx(1:nDim)
    CALL ElementTECPLOTPLOTTER(wh,lx0_3,ldx_3,limiter)
END SUBROUTINE ElementCallTECPLOTPLOTTER

RECURSIVE SUBROUTINE ElementCallTECPLOTADERDGPLOTTER(wh,lx0,ldx,limiter)
	USE, INTRINSIC :: ISO_C_BINDING
	USE TECPLOTPLOTTERmod !, only : ElementTECPLOTPLOTTER_ADERDG,nDOFm
	USE Parameters, only : nVar, nDim
	implicit none
	REAL, INTENT(IN) :: wh(nVar,nDOFm),lx0(nDim),ldx(nDim)
	real :: lx0_3(3),ldx_3(3)
	integer :: limiter
        lx0_3 = 0.
        ldx_3 = 0.

	lx0_3(1:nDim)=lx0(1:nDim)
	ldx_3(1:nDim)=ldx(1:nDim)
	CALL ElementTECPLOTPLOTTER_ADERDG(wh,lx0_3,ldx_3,limiter)
 END SUBROUTINE ElementCallTECPLOTADERDGPLOTTER

RECURSIVE SUBROUTINE ElementCallTECPLOTFVPLOTTER(wh,lx0,ldx,limiter)
	USE, INTRINSIC :: ISO_C_BINDING
	USE TECPLOTPLOTTERmod !, only : ElementTECPLOTPLOTTER_FV,nSubLim,nSubLim_GL
	USE Parameters, only : nVar, nDim
	implicit none
	REAL, INTENT(IN) :: wh(nVar,(nSubLim+2*nSubLim_GL)**nDIM),lx0(nDim),ldx(nDim)
	real :: lx0_3(3),ldx_3(3)
	integer :: limiter,i,j,k,Stencil
	lx0_3 = 0.
	ldx_3 = 0.
	lx0_3(1:nDim)=lx0(1:nDim)
	ldx_3(1:nDim)=ldx(1:nDim)
    !!PRINT *, lx0_3
    !!PRINT *, ldx_3
     !PRINT *, "zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz"
     !PRINT *, "zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz"
     !PRINT *,"nSubLim", nSublim
     !PRINT *,"nSubLim_GL", nSubLim_GL
     !PRINT *,"nDIM", nDIM
     !PRINT *,"nVar", nVar
     !PRINT *,"(nSubLim+2*nSubLim_GL)**nDIM", (nSubLim+2*nSubLim_GL)**nDIM
     !PRINT *,"(nSubLim+2*nSubLim_GL)", (nSubLim+2*nSubLim_GL)
     !PRINT *, "****************"
     !PRINT *, "lx0_3:", lx0_3
     !PRINT *, "ldx_3:", ldx_3
     !PRINT *, "****************" 
    !Stencil=(nSubLim+2*nSubLim_GL)
     !DO j=0,Stencil-1
     !    DO i=0,Stencil-1
     !        WRITE(*,'(a7,i2.2,a1,i2.2,a6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6)') ,"wh(2,:,",i+1,",",j+1,")     ", wh(2,j*Stencil**2+i*Stencil+1:j*Stencil**2+i*Stencil+Stencil)
     !    ENDDO
     !ENDDO 
     !STOP
    !WRITE(*,'(a,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6)') "wh(1)", wh(1,1:(nSubLim+2*nSubLim_GL))
    !WRITE(*,'(a,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6)') "wh(2)", wh(2,1:(nSubLim+2*nSubLim_GL))  
    !WRITE(*,'(a,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6)') "wh(1')", wh(1,(nSubLim+2*nSubLim_GL)+1:2*(nSubLim+2*nSubLim_GL))  
    !WRITE(*,'(a,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6)') "wh(2')", wh(2,(nSubLim+2*nSubLim_GL)+1:2*(nSubLim+2*nSubLim_GL)) 
    !WRITE(*,'(a,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6)') "wh(1'')", wh(1,2*(nSubLim+2*nSubLim_GL)+1:3*(nSubLim+2*nSubLim_GL))  
    !WRITE(*,'(a,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6)') "wh(2'')", wh(2,2*(nSubLim+2*nSubLim_GL)+1:3*(nSubLim+2*nSubLim_GL))  
    !WRITE(*,'(a,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6)') "wh(1*)", wh(1,(nSubLim+2*nSubLim_GL)**2:(nSubLim+2*nSubLim_GL)**2+(nSubLim+2*nSubLim_GL))    
    !WRITE(*,'(a,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6,E16.6)') "wh(2*)", wh(2,(nSubLim+2*nSubLim_GL)**2:(nSubLim+2*nSubLim_GL)**2+(nSubLim+2*nSubLim_GL))    
    !PRINT *, "zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz"
    !PRINT *, "------------------------------------------------------------"
    !PRINT *, "zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz"
	CALL ElementTECPLOTPLOTTER_FV(wh,lx0_3,ldx_3,limiter)
END SUBROUTINE ElementCallTECPLOTFVPLOTTER

RECURSIVE SUBROUTINE InitializeTECPLOTPLOTTER(time)
	USE, INTRINSIC :: ISO_C_BINDING
	USE TECPLOTPLOTTERmod !, ONLY: InitTECPLOTPLOTTER
	implicit none
	real :: time
	CALL InitTECPLOTPLOTTER(time)
END SUBROUTINE InitializeTECPLOTPLOTTER

RECURSIVE SUBROUTINE FinishTECPLOTPLOTTER(Myrank)
	USE, INTRINSIC :: ISO_C_BINDING
	USE TECPLOTPLOTTERmod !, ONLY: FinalizeTECPLOTPLOTTER
	implicit none
    !INTERFACE RECURSIVE SUBROUTINEFinalizeTECPLOTPLOTTER(Myrank)
    !integer :: Myrank
    !END INTERFACE
    !
	integer :: Myrank
	CALL FinalizeTECPLOTPLOTTER(Myrank)
   ! PRINT *, "FINALIZETECPLOTPLOTTER : done!"
END SUBROUTINE FinishTECPLOTPLOTTER


