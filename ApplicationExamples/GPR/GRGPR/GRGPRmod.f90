! DIM Parameters
  
  MODULE GRGPRmod 
    IMPLICIT NONE  

	contains

 SUBROUTINE MatrixInverse3x3(M,iM,det,whereindex) 
    !---------------
    ! compute the determinant det of the NxN-matrix M
    !---------------
    IMPLICIT NONE
    ! input variables 
    REAL, INTENT(IN)   :: M(3,3)
    ! output variables
    REAL, INTENT(OUT)    :: iM(3,3)
    REAL, INTENT(OUT)    :: det
    ! output variables
    REAL    :: ComputeDet,Id(3,3)
    INTEGER :: i,j
	INTEGER :: whereindex

	
    det = M(1,1)*M(2,2)*M(3,3)-M(1,1)*M(2,3)*M(3,2)-M(2,1)*M(1,2)*M(3,3)+M(2,1)*M(1,3)*M(3,2)+M(3,1)*M(1,2)*M(2,3)-M(3,1)*M(1,3)*M(2,2) !ComputeDet(M,3)
    IF(ABS(det).LT.1e-20) THEN
        print *, 'FATAL ERROR: det = 0'
        WRITE(*,*) M(1,1:3)
        WRITE(*,*) M(2,1:3)
        WRITE(*,*) M(3,1:3)
		print *, 'Whereindex=',whereindex
        STOP
    ENDIF 
    !
    iM(1,1) =M(2,2)*M(3,3)-M(2,3)*M(3,2)
    iM(1,2) =M(1,3)*M(3,2)-M(1,2)*M(3,3)
    iM(1,3) =M(1,2)*M(2,3)-M(1,3)*M(2,2)
    iM(2,1) =M(2,3)*M(3,1)-M(2,1)*M(3,3)
    iM(2,2) =M(1,1)*M(3,3)-M(1,3)*M(3,1)
    iM(2,3) =M(1,3)*M(2,1)-M(1,1)*M(2,3)
    iM(3,1) =M(2,1)*M(3,2)-M(2,2)*M(3,1)
    iM(3,2) =M(1,2)*M(3,1)-M(1,1)*M(3,2)
    iM(3,3) =M(1,1)*M(2,2)-M(1,2)*M(2,1)
    iM(1:3,1:3) = iM(1:3,1:3)/det
    !
    !Id = MATMUL(M(1:3,1:3),iM(1:3,1:3))
    !DO i=1,3
    !    DO j=1,3
    !        IF(i.eq.j) THEN
    !            IF((Id(i,j)-1.)**2..GT.1e-18) THEN
    !                print *, 'FATAL ERROR: iM*M !=1'
    !                STOP
    !            ENDIF
    !        ELSE
    !            IF((Id(i,j)**2).GT.1e-18) THEN
    !                print *, 'FATAL ERROR: iM*M !=1'
    !                STOP
    !            ENDIF
    !        ENDIF
    !    ENDDO
    !ENDDO
    !
    CONTINUE
    !
    END SUBROUTINE MatrixInverse3x3
    
 SUBROUTINE MatrixInverse(N,A,iA)
  IMPLICIT NONE
  INTEGER       :: N
  REAL          :: A(N,N), iA(N,N)
  !
  INTEGER       :: i,j,flag,ml(1) 
  REAL          :: piv
  REAL          :: temp(2*N)
  REAL          :: C(2*N,N)  
  !
  C(1:N,:)     = TRANSPOSE(A(:,:))
  C(N+1:2*N,:) = 0. 
  DO i = 1, N
     C(N+i,i) = 1.
  ENDDO
  !    
  ! Forward elimination and row swapping (if necessary)
  ! 
  DO i = 1, N
     ! If pivot element is zero, then swap rows 
     ml = MAXLOC(ABS(C(i,i:N))) 
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
     DO j = i+1, N 
        C(:,j) = C(:,j) - C(i,j)*C(:,i)
     ENDDO
  ENDDO
  !
  ! Back substitution
  !
  DO i = N,1,-1   
     DO j = i-1,1,-1
        C(:,j) = C(:,j) - C(i,j)*C(:,i)
     ENDDO
  ENDDO
  !
  iA = TRANSPOSE( C(N+1:2*N,:) ) 
  !
END SUBROUTINE MatrixInverse
 


SUBROUTINE RTSAFE_C2P_RMHD1(v2,X1,X2,XACC,gam,d,e,s2,b2,sb2,w,FAILED)
  IMPLICIT NONE
  INTEGER, PARAMETER    :: MAXIT=200
  INTEGER               :: J
  REAL                  :: v2 
  REAL                  :: X1,X2,XACC,gam,d,e,s2,b2,sb2,w
  REAL                  :: FL,FH,DF,XH,XL,SWAP,DXOLD,DX,F,TEMP
  LOGICAL               :: FAILED
  !
  FAILED = .FALSE.
  CALL FUNC_C2P_RMHD1(X1,FL,DF,gam,d,e,s2,b2,sb2,w)
  IF(FL.EQ.0.) THEN
     v2=X1
     RETURN
  ENDIF
  CALL FUNC_C2P_RMHD1(X2,FH,DF,gam,d,e,s2,b2,sb2,w)
  IF(FH.EQ.0.) THEN
     v2=X2
     RETURN
  ENDIF
  IF(FL*FH.GT.0.) THEN
     FAILED = .TRUE.
     v2 = 0. ! assign value even if it fails
     RETURN
  ENDIF
  IF(FL.LT.0.)THEN
     XL=X1
     XH=X2
  ELSE
     XH=X1
     XL=X2
     SWAP=FL
     FL=FH
     FH=SWAP
  ENDIF
  v2=.5*(X1+X2)
  DXOLD=ABS(X2-X1)
  DX=DXOLD
  CALL FUNC_C2P_RMHD1(v2,F,DF,gam,d,e,s2,b2,sb2,w)
  DO 11 J=1,MAXIT
     IF(((v2-XH)*DF-F)*((v2-XL)*DF-F).GE.0. &
          .OR. ABS(2.*F).GT.ABS(DXOLD*DF) ) THEN
        DXOLD=DX
        DX=0.5*(XH-XL)
        v2=XL+DX
        IF(XL.EQ.v2)RETURN
     ELSE
        DXOLD=DX
        DX=F/DF
        TEMP=v2
        v2=v2-DX
        IF(TEMP.EQ.v2)RETURN
     ENDIF
     IF(ABS(DX).LT.XACC) RETURN
     CALL FUNC_C2P_RMHD1(v2,F,DF,gam,d,e,s2,b2,sb2,w)
     IF(F.LT.0.) THEN
        XL=v2
        FL=F
     ELSE
        XH=v2
        FH=F
     ENDIF
11 CONTINUE
     FAILED = .TRUE.
     v2 = 0. ! assign value even if it fails


END SUBROUTINE RTSAFE_C2P_RMHD1

SUBROUTINE FUNC_C2P_RMHD1(x,f,df,gam,d,e,s2,b2,sb2,w)
  !
  ! This is the CONS2PRIM strategy adopted by Del Zanna et al. (2007) A&A, 473, 11-30
  ! and it corresponds to their choice 3 in Section 3.2
  !
  IMPLICIT NONE
  REAL, PARAMETER :: third=1./3.
  INTEGER    :: iter
  REAL       :: x,f,df,v2,rho,c0,c2,c3,dw,dc2,dc3,dlogw,wb,vb2
  REAL       :: gam,d,e,s2,b2,sb2,w
  INTENT(IN) :: x,gam,d,e,s2,b2,sb2
  INTENT(OUT):: f,df,w
  
  v2=x
  rho=d*sqrt(1.-v2)
  
  c3=1.-gam*(1.-v2)
  c2=gam*rho+.5*b2*(1.+v2)-e
  c0=-0.5*sb2
  
  ! For every x=v2, we solve a cubic in W of the form: 
  ! c3*W^3+c2*W^2+c0=0 (c3>0, c0<=0)
  ! W=y of the paper. If sb=0 ( meaning c0 = 0), 
  ! w = -c2/c3 > 0 and dw = 0 in the do loop below. 
  ! If -c2/c3 < 0 when sb=0, which makes w=0, 
  ! this is a signature that something was wrong before.
  
  if ( abs ( c0) < 1.0d-20) then
     w = -c2 / c3
  else
     w = max ( - c2 / c3, ( -c0 / c3)**third)
     do iter = 1,100
        dw = -((c3*w + c2)*w**2 + c0)/((3*c3*w + 2*c2)*w)
        if (abs(dw/w)<1.e-10) exit
        w = w + dw
     end do
  endif
  
  dc3   = gam
  dc2   = 0.5 * ( b2 - gam * rho / (1.0 - v2))
  dlogw = -( dc3 * w + dc2 ) / ( 3.0 * c3 * w + 2.0 * c2)
  wb    = w + b2
  vb2   = sb2 / w**2
  f     = wb**2 * v2 - ( 2.0 * w + b2) * vb2 - s2
  df    = wb * ( wb + 2.0 * dlogw * ( w * v2 + vb2))
  
END SUBROUTINE FUNC_C2P_RMHD1 
	
END MODULE GRGPRmod  
