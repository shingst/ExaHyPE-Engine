! GRMHD PDE.f90
! Trento (EQNTYPE4)

RECURSIVE SUBROUTINE PDEFlux(f,g,h,Q)
	USE Parameters, ONLY : nVar, nDim, epsilon1
	USE SpecificVarEqn99
	USE iso_c_binding
	IMPLICIT NONE
	REAL :: f(nVar), g(nVar), h(nVar), hz(nVar), Q(nVar), V(nVar)
	REAL, PARAMETER :: epsilon = 1e-14 
	INTENT(IN)  :: Q
	INTENT(OUT) :: f, g, h
	! Local variable
	REAL :: Kelast,detA, A(3,3),p,LEalpha,ialpha, u(3), sigma_vec(6), AU(3)
  
	Kelast=Q(15)+2.0/3.0*Q(16)  
	detA = Q(5) * Q(9) * Q(13) - Q(5) * Q(10) * Q(12) - Q(8) * Q(6) * Q(13) + Q(8) * Q(7) * Q(12) + Q(11) * Q(6) * Q(10) - Q(11) * Q(7) * Q(9)
	A(1,1)=Q(5)
	A(1,2)=Q(6)
	A(1,3)=Q(7)
	A(2,1)=Q(8)
	A(2,2)=Q(9)
	A(2,3)=Q(10)
	A(3,1)=Q(11)
	A(3,2)=Q(12)
	A(3,3)=Q(13)  
	p=-Kelast*(detA)**2*(1-detA)
	f = 0.0
	g = 0.0 
	h = 0.0
	
	call ComputeGPRLEstress(Kelast,sigma_vec,Q, .false.)
	
	LEalpha=Q(14)           ! GPR+LE+DI
    if(epsilon1<0) then
        ialpha=1/LEalpha
    else
        ialpha=LEalpha/(LEalpha**2+epsilon1*(1-LEalpha))
    end if
	
	u=Q(2:4)/Q(1)*ialpha
    AU = MATMUL( A, u(1:3) )
	
	f(1)=Q(1)*u(1)
	f(2)=Q(2)*u(1)    +LEalpha*p  -LEalpha*sigma_vec(1)
	f(3)=Q(3)*u(1)                -LEalpha*sigma_vec(4)
	f(4)=Q(4)*u(1)                -LEalpha*sigma_vec(6)
	f(5 )=AU(1)
	f(8 )=AU(2)
	f(11)=AU(3)
	f(19)=u(1)*Q(19)+u(1)*p-u(1)*sigma_vec(1)-u(2)*sigma_vec(4)-u(3)*sigma_vec(6)

	g(1)=Q(1)*u(2)         
	g(2)=Q(2)*u(2)          -LEalpha*sigma_vec(4)
	g(3)=Q(3)*u(2) +LEalpha*p -LEalpha*sigma_vec(2)
	g(4)=Q(4)*u(2)          -LEalpha*sigma_vec(5)
	g(6 )=AU(1)
	g(9 )=AU(2)
	g(12)=AU(3)
	g(19)=u(2)*Q(19)+u(2)*p-u(1)*sigma_vec(4)-u(2)*sigma_vec(2)-u(3)*sigma_vec(5)

	h(1)=Q(1)*u(3)
	h(2)=Q(2)*u(3)          -LEalpha*sigma_vec(6)
	h(3)=Q(3)*u(3)          -LEalpha*sigma_vec(5)
	h(4)=Q(4)*u(3) +LEalpha*p -LEalpha*sigma_vec(3)
	h(7 )=AU(1)
	h(10)=AU(2)
	h(13)=AU(3)
	h(19)=u(3)*Q(19)+u(3)*p-u(1)*sigma_vec(6)-u(2)*sigma_vec(5)-u(3)*sigma_vec(3)
	
	f(23)=Q(2)/Q(1)
	g(23)=Q(2)/Q(1)
	h(23)=Q(2)/Q(1)
	
  END SUBROUTINE PDEFlux


RECURSIVE SUBROUTINE PDENCP(BgradQ,Q,gradQ) 
	USE Parameters, ONLY :  nVar, nDim, epsilon1
	IMPLICIT NONE
	REAL, PARAMETER :: epsilon = 1e-14
	! 11. Oct 21:40: This was a matrix BGradQ(nVar, nDim) but is a vector in spaceTimePredictorNonlinear
	REAL, INTENT(OUT) :: BgradQ(nVar)
	REAL, INTENT(IN)  :: gradQ(nVar, 3)
	REAL, INTENT(IN)  :: Q(nVar)
	! Linear elasticity variables
	REAL :: lam,mu,irho, ialpha, u(3), alpha
	REAL :: AQx(nVar), BQy(nVar), CQz(nVar) , Qx(nVar), Qy(nVar), Qz(nVar) 
   
	Qx = gradQ(:,1)
	Qy = gradQ(:,2)
	IF(nDim==3) THEN
		Qz = gradQ(:,3)
	ELSE
		Qz = 0.0 
	ENDIF 
  !
    ! Linear elasticity part
	AQx = 0.   
	BQy = 0. 
	CQz = 0. 

	alpha=Q(14) ! GPR+LE+DI
    if(epsilon1<0) then
        ialpha=1/alpha
    else
        ialpha=alpha/(alpha**2+epsilon1*(1-alpha))
    end if
    u=Q(2:4)/Q(1)*ialpha
    !
	AQx(5) =  -u(2)*Qx(6) - u(3)*Qx(7)
	AQx(6) =  +u(1)*Qx(6)
	AQx(7) =  +u(1)*Qx(7)
	AQx(8) =  -u(2)*Qx(9) - u(3)*Qx(10)
	AQx(9) =  +u(1)*Qx(9)
	AQx(10) = +u(1)*Qx(10)
	AQx(11) = -u(2)*Qx(12) - u(3)*Qx(13)
	AQx(12) = +u(1)*Qx(12)
    AQx(13) = +u(1)*Qx(13)
    
    BQy(5 ) =  +u(2)*Qy(5)
    BQy(6 ) =  -u(1)*Qy(5) - u(3)*Qy(7)
    BQy(7 ) =  +u(2)*Qy(7)
    BQy(8 ) =  +u(2)*Qy(8)
    BQy(9 ) = -u(1)*Qy(8) - u(3)*Qy(10)
    BQy(10) = +u(2)*Qy(10)
    BQy(11) = +u(2)*Qy(11)
    BQy(12) = -u(1)*Qy(11) - u(3)*Qy(13)
    BQy(13) = +u(2)*Qy(13)
    
    CQz(5) =  +u(3)*Qz(5)
    CQz(6) =  +u(3)*Qz(6)
    CQz(7) =  -u(1)*Qz(5) - u(2)*Qz(6)
    CQz(8) =  +u(3)*Qz(8)
    CQz(9) = +u(3)*Qz(9) 
    CQz(10) = -u(1)*Qz(8) - u(2)*Qz(9)
    CQz(11) = +u(3)*Qz(11) 
    CQz(12) = +u(3)*Qz(12)
    CQz(13) = -u(1)*Qz(11) - u(2)*Qz(12)
	
    if( nDim .eq. 2) then
        BgradQ = AQx + BQy        
    else
        BgradQ = AQx + BQy + CQz     
    end if    
    
END SUBROUTINE PDENCP


RECURSIVE SUBROUTINE PDEEigenvalues(L,Q,n)
	USE Parameters, ONLY :  nVar, nDim
	USE iso_c_binding
	IMPLICIT NONE
	REAL :: L(nVar), n(3), Q(nVar), V(nVar)
	INTENT(IN)  :: Q,n
	INTENT(OUT) :: L 
	REAL  :: lam,mu,irho,VPR(3),uu,c0,cs

	L = 0.
  
	lam  = Q(15)   
    mu   = Q(16)
    irho = 1./Q(1)

	VPR=Q(2:4)/Q(1)  ! Skip alpha for now!
    uu = SQRT(sum(VPR(:)**2) ) 
    c0=sqrt((lam+2*mu)*irho)
    cs=sqrt(mu*irho)
    L(1)=uu+c0
    L(2)=uu-c0
    L(3)=uu+cs
    L(4)=uu-cs
    L(5)=uu+cs
    L(6)=uu-cs
    
    L(11:17)=uu
END SUBROUTINE PDEEigenvalues

RECURSIVE SUBROUTINE PDESource(S,Q) 
  USE Parameters, ONLY:  nVar, nDim
  USE SpecificVarEqn99
  USE iso_c_binding
  IMPLICIT NONE
  ! --------------------------------------------
  ! Argument list declaration
  REAL :: S(nvar), Q(nvar)
  INTENT(IN)  :: Q 
  INTENT(OUT) :: S
! Local variable declaration 
  REAL :: AM(3,3), Id(3,3), G(3,3) ,devG(3,3), detA, detA2,psiM(3,3), temp2, T
  REAL :: V(nvar), tau,ee
  ! --------------------------------------------
  ! Local variable declaration 

	S = 0.
	CALL PDECons2Prim(V,Q)
	
	AM(1,1)=V(5)
	AM(1,2)=V(6)
	AM(1,3)=V(7)
	AM(2,1)=V(8)
	AM(2,2)=V(9)
	AM(2,3)=V(10)
	AM(3,1)=V(11)
	AM(3,2)=V(12)
	AM(3,3)=V(13)
	
	G      = MATMUL( TRANSPOSE(AM), AM )
	
	Id     = 0. 
	Id(1,1) = 1.0; Id(2,2) = 1.0; Id(3,3) = 1.0 
	devG   = G - (G(1,1)+G(2,2)+G(3,3))/3.*Id  
	detA   = AM(1,1)*AM(2,2)*AM(3,3)-AM(1,1)*AM(2,3)*AM(3,2)-AM(2,1)*AM(1,2)*AM(3,3)+AM(2,1)*AM(1,3)*AM(3,2)+AM(3,1)*AM(1,2)*AM(2,3)-AM(3,1)*AM(1,3)*AM(2,2)     ! this is the determinant we have     
	psiM   = 3./(detA)*MATMUL(AM,devG)                    ! in the coefficient of the source term, we use the detA2 from the compatibility relation, to reduce nonlinearities in A 
	! 
	ee=mu2tauFL(Q(24), Q(15),Q(16),Q(1),detA)    ! Conversion from mu->tau
    tau = 1.e+16*max(0.0,Q(18))+ee*min(1.0,(1.0-Q(18)))
	if(Q(18)<1.e-2) then
        tau = ee*min(1.0,(1.0-Q(18)))   
    end if
	
	S(5)= -psiM(1,1)/tau*detA**(8./3.)
	S(6)= -psiM(1,2)/tau*detA**(8./3.)
	S(7)= -psiM(1,3)/tau*detA**(8./3.)
	S(8)= -psiM(2,1)/tau*detA**(8./3.)
	S(9)=-psiM(2,2)/tau*detA**(8./3.)
	S(10)=-psiM(2,3)/tau*detA**(8./3.)
	S(11)=-psiM(3,1)/tau*detA**(8./3.)
	S(12)=-psiM(3,2)/tau*detA**(8./3.)
	S(13)=-psiM(3,3)/tau*detA**(8./3.)
	  
	S(20)=Q(2)/Q(1)
	S(21)=Q(3)/Q(1)
	S(22)=Q(4)/Q(1)
END SUBROUTINE PDESource

RECURSIVE SUBROUTINE PDEVarName(MyNameOUT,ind) 
  USE Parameters, ONLY: nVar  
  IMPLICIT NONE     
  CHARACTER(LEN=10):: MyName(nVar),MyNameOUT
  INTEGER			:: ind


  ! EQNTYPE99
    MyName(1)  = 'rho' 
    MyName(2)  = 'u' 
    MyName(3)  = 'v'
    MyName(4)  = 'w' 
    MyName(5)  = 'A11' 
    MyName(6)  = 'A12' 
    MyName(7)  = 'A13' 
    MyName(8)  = 'A21' 
    MyName(9)  = 'A22' 
    MyName(10) = 'A23'
    MyName(11) = 'A31'
    MyName(12) = 'A32' 
    MyName(13) = 'A33'
    MyName(14) = 'alpha'
    MyName(15) = 'lambda'
    MyName(16) = 'mu'
    MyName(17) = 'Y0'
    MyName(18) = 'xi'
    MyName(19) = 'rhoE'
	
    MyName(20) = 'xDisp'
    MyName(21) = 'yDisp'
    MyName(22) = 'zDisp'
    MyName(23) = 'Disp'
    MyName(24) = 'mu_f'
	
	MyNameOUT=MyName(ind+1)
    END SUBROUTINE PDEVarName

RECURSIVE SUBROUTINE PDEAuxName(MyNameOUT,ind) 
	USE Parameters, ONLY: nVar  
	IMPLICIT NONE     
	CHARACTER(LEN=10):: AuxName(nVar),MyNameOUT
	INTEGER			:: ind


	! EQNTYPE99
	AuxName(1) = 'sxx'
	AuxName(2) = 'syy'
	AuxName(3) = 'szz'
	AuxName(4) = 'sxy'
	AuxName(5) = 'syz'
	AuxName(6) = 'sxz'
	AuxName(7) = 'YY'
	AuxName(8) = 'TT'

	MyNameOUT=AuxName(ind+1)
END SUBROUTINE PDEAuxName
	
RECURSIVE SUBROUTINE PDEAuxVar(aux,Q,x,time)
    USE Parameters, ONLY : nVar,nAux
	USE SpecificVarEqn99
	implicit none
	real :: aux(nAux),Q(nVar),x(3),time
	real :: detA
	! print *, "ENTERED HERE ------------------------------------------------------------------------------"
	aux=0.
	!return
	call computeGPRLEstress(aux(7),aux(1:6),Q,.true.)
	call ComputeAcontribution(detA,Q)
	aux(8)=1.0/Q(1)*(Q(19)-0.5*1.0/Q(1)*(sum(Q(2:4)**2))-detA)
END SUBROUTINE PDEAuxVar
	
RECURSIVE SUBROUTINE PDEJacobian(An,Q,gradQ,nv)
  USE Parameters, ONLY : nVar, nDim
  USE iso_c_binding
  IMPLICIT NONE
  ! Argument list 
  REAL :: An(nVar,nVar)
  REAL :: Q(nVar), gradQ(nVar,nDim), nv(3) 
  INTENT(IN)  :: Q, gradQ, nv
  INTENT(OUT) :: An 
  
  CALL PDEMatrixB(An,Q,nv) 
  
END SUBROUTINE PDEJacobian

RECURSIVE SUBROUTINE PDEMatrixB(An,Q,nv) 
	USE Parameters, ONLY : nVar, nDim, epsilon1
	USE iso_c_binding
	IMPLICIT NONE
	! Argument list 
	REAL :: An(nVar,nVar)
	REAL :: Q(nVar), nv(3) 
	INTENT(IN)  :: Q,nv
	INTENT(OUT) :: An  
	! Local variables
    ! Linear elasticity variables
   REAL :: lam,mu,irho, ialpha, uv(3),alpha
   REAL :: A(nVar,nVar), B(nVar,nVar), C(nVar,nVar)

    A = 0.0
    B = 0.0
    C = 0.0 

    alpha=Q(14)
    if(epsilon1<0) then
        ialpha=1/alpha
    else
        ialpha=alpha/(alpha**2+epsilon1*(1-alpha))
    end if
	uv=Q(2:4)/Q(1)*ialpha
	

    A(5,6) =  -uv(2)
    A(5,7) =  -uv(3)
	A(6,6) =  +uv(1)
	A(7,7) =  +uv(1)
	A(8,9) =  -uv(2)
    A(8,10)=  -uv(3)
	A(9,9) =  +uv(1)
	A(10,10)= +uv(1)
	A(11,12)= -uv(2)
    A(11,13)= -uv(3)
	A(12,12)= +uv(1)
    A(13,13)= +uv(1)

    B(5 ,5) =  +uv(2)
    B(6 ,5) =  -uv(1)
    B(6 ,7) =  -uv(3)
    B(7 ,7) =  +uv(2)
    B(8 ,8) =  +uv(2)
    B(9 ,8) =  -uv(1)
    B(9,10) =  -uv(3)
    B(10,10)=  +uv(2)
    B(11,11)=  +uv(2)
    B(12,11)=  -uv(1)
    B(12,13)=  -uv(3)
    B(13,13)=  +uv(2)    
    
    C(5 ,5)  =  +uv(3)
    C(6 ,6)  =  +uv(3)
    C(7 ,5)  =  -uv(1)
    C(7 ,6)  =  -uv(2)
    C(8 ,8)  =  +uv(3)
    C(9 ,9)  =  +uv(3)
    C(10,8)  =  -uv(1)
    C(10,9)  =  -uv(2)
    C(11,11) =  +uv(3) 
    C(12,12) =  +uv(3)
    C(13,11) =  -uv(1) 
    C(13,12) =  -uv(2)	
	
    if( nDim .eq. 2) then
        An = A*nv(1) + B*nv(2)         
    else
        An = A*nv(1) + B*nv(2) + C*nv(3)     
    end if
    
    

  
END SUBROUTINE PDEMatrixB
    
RECURSIVE SUBROUTINE PDEIntermediateFields(RL,LL,iRL,Q,nv) 
  USE Parameters, ONLY : nVar, nDim, nLin
  USE iso_c_binding 
  IMPLICIT NONE
  ! Argument list 
  REAL :: RL(nVar,nLin), LL(nLin,nLin), iRL(nLin,nVar)
  REAL :: Q(nVar), nv(3)
  REAL :: x(nDim), time
  INTENT(IN)  :: Q,nv
  INTENT(OUT) :: RL,LL,iRL 
  ! Local variables
  INTEGER :: i,j,k, zero,iErr, minl(1), maxl(1) 
  REAL :: R(nVar,nVar), L(nVar,nVar), Lv(nVar), iR(nVar,nVar)
  !  
  
  CALL PDEEigenvectors(R,L,iR,Q,nv)
  RL(:,1:5)  = R(:,14:18)
  iRL(1:5,:) = iR(14:18,:)
  
  RL(:,6:10)  = R(:,20:24)
  iRL(6:10,:) = iR(20:24,:)
  LL = 0. 
  LL(1,1) = L(14,14)      ! lamb=0
  LL(2,2) = L(15,15)      ! lamb=0
  LL(3,3) = L(16,16)      ! lamb=0
  LL(4,4) = L(17,17)      ! lamb=0 
  LL(5,5) = L(18,18)      ! lamb=0 
  
  LL(6,6) = L(20,20)      ! lamb=0 
  LL(7,7) = L(21,21)      ! lamb=0 
  LL(8,8) = L(22,22)      ! lamb=0 
  LL(9,9) = L(23,23)      ! lamb=0
  LL(10,10) = L(24,24)      ! lamb=0   
  !
END SUBROUTINE PDEIntermediateFields
    
    
RECURSIVE SUBROUTINE PDEEigenvectors(R,L,iR,Q,nv) 
	USE Parameters, ONLY : nVar, nDim, epsilon1
	USE iso_c_binding  
	IMPLICIT NONE
	! Argument list 
	REAL :: R(nVar,nVar), L(nVar,nVar), iR(nVar,nVar)
	REAL :: Q(nVar), nv(3)
	REAL :: x(nDim), time
	INTENT(IN)  :: Q,nv
	INTENT(OUT) :: R,L,iR 
	! Local variables
	REAL    :: Lambda(nVar)
INTEGER :: j    

	Lambda=0.
	R=0.
	iR=0.
	do j=1,nVar
		R(j,j)=1.
		iR(j,j)=1.
	end do
	L = 0.
    END SUBROUTINE PDEEigenvectors 
    
	
RECURSIVE SUBROUTINE DynamicRupture(x_in, t, Q)
	USE SpecificVarEqn99
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar, nDim, ICType
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: x_in(nDim), t     ! 
	REAL, INTENT(OUT)              :: Q(nVar)        ! 
	! Local variables
	REAL :: stressnorm, sigma_vec(6),ee,u(3),x(3)
	! Compute the normal stress using the Formula by Barton (IJNMF 2009)
	x=0.
	x(1:nDim)=x_in(1:nDim)
	call computeGPRLEstress(stressnorm,sigma_vec,Q,.true.)
	IF( stressnorm > Q(17) ) THEN
		! If the normal stress is greater than the illness stress Y0 (stored in Q(17), then broke the material
		Q(18)=0.!0.
		if(USEFrictionLaw) then
			! Compute the proper friction coefficient according to the prescribed friction law
			u=Q(2:4)/Q(1)
			ee=LEfriction_mu(x,t,u(2:4),Q(23),SSCRACKFL)
			Q(24)=ee*2.0*SSCRACKFL%dx
		end if
		! Now in the this DoF it follows the NS friction law with mu proportional to mu_f ( stored in Q(24) )
		!print *, 'Crack!'
	ELSE
		!Q(18)=1. ! Reattach 
	END IF
	! Add here some static rupture (like for the SSCRACK test case
	!
	!
	if(.not. SSCRACKFL%DynamicFL .and. USEFrictionLaw) then ! Static marker
	    u=Q(2:4)/Q(1)
		ee=LEfriction_mu(x,t,u(2:4),Q(23),SSCRACKFL)
		if(ee<SSCRACKFL%mu_s*0.999) then
			if(abs(x(1)).le.10000 .and. abs(x(2)).le.SSCRACKFL%dx) then
				Q(18)=0.
				Q(24)=ee*2.0*SSCRACKFL%dx
			end if
		end if
	end if				
				
	!
	!
	!
END SUBROUTINE DynamicRupture	
	
	
RECURSIVE SUBROUTINE getNumericalSolution(V,Q) 
  USE Parameters, ONLY: nVar  
  IMPLICIT NONE     
  REAL				:: V(nVar), Q(nVar)
  CALL PDECons2Prim(V,Q)
  !V=0.
END SUBROUTINE getNumericalSolution

RECURSIVE SUBROUTINE getExactSolution(V,pos, timeStamp) 
  USE Parameters, ONLY: nVar , nDim  
  IMPLICIT NONE     
  REAL				:: V(nVar), Q(nVar), pos(nDim), timeStamp
  call InitialData(pos, timeStamp, Q)
  CALL PDECons2Prim(V,Q)
  !V(1)=pos(1)**2!*pos(2)**2
END SUBROUTINE getExactSolution
	
RECURSIVE SUBROUTINE HLLEMRiemannSolver(basisSize,NormalNonZero,lFbndL,lFbndR,lQbndL,lQbndR,QavL,QavR) 
  USE Parameters, ONLY : nVar, nDim, nLin
  USE iso_c_binding 
  ! Local variables
  INTEGER, INTENT(IN)   :: NormalNonZero, basisSize
  REAL, INTENT(IN)     :: lQbndL(nVar,basisSize,basisSize)
  REAL, INTENT(IN)     :: lQbndR(nVar,basisSize,basisSize)
  REAL, INTENT(INOUT)  :: lFbndL(nVar,basisSize,basisSize)
  REAL, INTENT(INOUT)  :: lFbndR(nVar,basisSize,basisSize)
    ! Local variables 
	REAL :: f(nVar), g(nVar), h(nVar)
INTEGER           :: i,j,k,l, ml(1)  
REAL              :: aux(nDim), Id(nVar,nVar), smax, Qav(nVar),QavL(nVar), QavR(nVar) 
REAL              ::  xGP, yGP, xdeb, ydeb  
REAL              :: Bn(nVar,nVar), DM(nVar,nVar), ncp(nVar), nv(3)
REAL              :: gradQ(nVar,3), src(nVar),flattener(nLin)
  REAL    :: absA(nVar,nVar), amax  
  REAL    :: QM(nVar),LL(nVar),LR(nVar),LM(nVar)
  REAL    :: deltaL(nLin,nLin),Lam(nLin,nLin),Lap(nLin,nLin) 
  REAL    :: RL(nVar,nLin),iRL(nLin,nVar),LLin(nLin,nLin) 
  
  nv(:)=0.
  nv(NormalNonZero+1)=1.
  !print *, "Normal non zero in fortran=", NormalNonZero
  !print *, "basisSize=", basisSize
  !print *, "NormalNonZero=", NormalNonZero
  !print *, "QavR=",QavR(1)
  !return
  !nv(NormalNonZero)=1.;
  !FL=0.
  !FR=0.
	flattener=1.
	lFbndL=0.
	lFbndR=0.
	CALL PDEflux(f,g,h,QavL)
	lFbndL(:,1,1)=f*nv(1)+g*nv(2)*h*nv(3)

	CALL PDEflux(f,g,h,QavR)
	lFbndR(:,1,1)=f*nv(1)+g*nv(2)*h*nv(3)	

    CALL PDEEigenvalues(LL,QavL,nv) 
    CALL PDEEigenvalues(LR,QavR,nv) 
    smax = MAX( MAXVAL(ABS(LL)), MAXVAL(ABS(LR)) )
    ! Identity matrix 
    Id = 0.0 
    DO i=1,nVar
        Id(i,i)=1.0
    ENDDO
    gradQ = 0.0      
    ! HLLEM
    Qav = 0.5*(QavL+QavR) 
    CALL PDEIntermediateFields(RL,LLin,iRL,Qav,nv) 
    Lam = 0.5*(LLin-ABS(LLin))
    Lap = 0.5*(LLin+ABS(LLin)) 
    deltaL = 0.0
    DO i = 1, nLin
        deltaL(i,i) = ( 1. - Lam(i,i)/(-smax-1e-14) - Lap(i,i)/(smax+1e-14) )*flattener(i)  
    ENDDO    
    absA = 0. 
    DO i = 1, nVar
        absA(i,i) = -0.5*smax  ! regular Rusanov diffusion, only on the diagonal 
    ENDDO  
    absA = absA + 0.5*smax*MATMUL( RL, MATMUL(deltaL, iRL) )  ! HLLEM anti-diffusion

        DO k = 1, basisSize
          DO j = 1, basisSize
                Qav = 0.5*(lQbndR(:,j,k)+lQbndL(:,j,k)) 
                gradQ(:,NormalNonZero+1) = lQbndR(:,j,k) - lQbndL(:,j,k) 
                CALL PDENCP(ncp,Qav,gradQ)
				lFbndL(:,j,k) = 0.5*( lFbndR(:,j,k) + lFbndL(:,j,k) ) + MATMUL(absA, lQbndR(:,j,k) - lQbndL(:,j,k) )    ! purely conservative flux 
				!lFbndL(:,j,k) = 0.5*( lFbndR(:,j,k) + lFbndL(:,j,k) ) - 0.5*smax*( lQbndR(:,j,k) - lQbndL(:,j,k) )      ! purely conservative flux 
				lFbndR(:,j,k) = lFbndL(:,j,k) - 0.5*ncp(:)                                                              ! subtract the non-conservative product 
                lFbndL(:,j,k) = lFbndL(:,j,k) + 0.5*ncp(:)
            ENDDO
        ENDDO				
END SUBROUTINE HLLEMRiemannSolver


RECURSIVE SUBROUTINE TESTRiemannSolver(basisSize,NormalNonZero,lFbndL,lFbndR,lQbndL,lQbndR,QavL,QavR) 
  USE Parameters, ONLY : nVar, nDim, nLin
  USE iso_c_binding 
  ! Local variables
  INTEGER, INTENT(IN)   :: NormalNonZero, basisSize
  REAL, INTENT(IN)     :: lQbndL(nVar,basisSize,basisSize)
  REAL, INTENT(IN)     :: lQbndR(nVar,basisSize,basisSize)
  REAL, INTENT(INOUT)  :: lFbndL(nVar,basisSize,basisSize)
  REAL, INTENT(INOUT)  :: lFbndR(nVar,basisSize,basisSize)
    ! Local variables 
INTEGER           :: i,j,k,l, ml(1)  
REAL              :: aux(nDim), Id(nVar,nVar), smax, Qav(nVar),QavL(nVar), QavR(nVar) 
REAL              ::  xGP, yGP, xdeb, ydeb  
REAL              :: Bn(nVar,nVar), DM(nVar,nVar), ncp(nVar), nv(3)
REAL              :: gradQ(nVar,3), src(nVar),flattener(nLin)
  REAL    :: absA(nVar,nVar), amax  
  REAL    :: QM(nVar),LL(nVar),LR(nVar),LM(nVar)
  REAL    :: deltaL(nLin,nLin),Lam(nLin,nLin),Lap(nLin,nLin) 
  REAL    :: RL(nVar,nLin),iRL(nLin,nVar),LLin(nLin,nLin) 
  nv(:)=0.
  !if(NormalNonZero .gt. 0) then
	nv(NormalNonZero+1)=1.
  !end if
	Qav = 0.5*(QavL+QavR) 
	!CALL PDEMatrixB(Bn,0.5*(QavL+QavR),nv)   ! evaluate the system matrix just once in the averaged state 
    !CALL DissipationMatrix(DM,QavL,QavR,nv)        ! according to the RIemann solver, compute the dissipation matrix 
    !DO k = 1, basisSize
    !  DO j = 1, basisSize
    !        !CALL PDEMatrixB(Bn,0.5*(lQbndL(:,j,k)+lQbndR(:,j,k)),nv,0.5*(lparbndL(:,j,k)+lparbndR(:,j,k))) ! evaluate the system matrix in each Gaussian point again (slow) 
    !        lFbndL(:,j,k) = 0.5*MATMUL( Bn + DM, lQbndR(:,j,k) - lQbndL(:,j,k) ) ! WAS -
    !        lFbndR(:,j,k) = 0.5*MATMUL( Bn - DM, lQbndR(:,j,k) - lQbndL(:,j,k) )! WAS + 	
    !    ENDDO
    !ENDDO
	!print *, NormalNonZero
	!print *, maxval(abs(lFbndL(:,1,1)-lFbndR(:,1,1)))
	!if(maxval(abs(Bn)) .gt. 0) then
	!	print *, 'OK!', maxval(abs(Bn))
	!end if
	
	
	
	CALL PDEEigenvalues(LL,QavL,nv) 
	CALL PDEEigenvalues(LR,QavR,nv) 
	CALL PDEEigenvalues(LM,Qav,nv) 
	sL = MIN(0.0, MINVAL(LL), MINVAL(LM) ) 
	sR = MAX(0.0, MAXVAL(LR), MAXVAL(LM) )
    ! Identity matrix 
    Id = 0.0 
    DO i=1,nVar
        Id(i,i)=1.0
    ENDDO
	lFbndL=0;
	lFbndR=0;
	gradQ = 0.0  
        DO k = 1, basisSize
          DO j = 1, basisSize
                Qav = 0.5*(lQbndR(:,j,k)+lQbndL(:,j,k)) 
                gradQ(:,NormalNonZero+1) = lQbndR(:,j,k) - lQbndL(:,j,k) 
                CALL PDENCP(ncp,Qav,gradQ)
				lFbndL(:,j,k) = sL*sR/(sR-sL)*( lQbndR(:,j,k) - lQbndL(:,j,k) )
                !lFbndL(:,j,k) = ( sR*lFbndL(:,j,k) - sL*lFbndR(:,j,k) )/(sR-sL) + sR*sL/(sR-sL)*( lQbndR(:,j,k) - lQbndL(:,j,k) )   ! purely conservative flux 
                lFbndR(:,j,k) = lFbndL(:,j,k) - sR/(sR-sL)*ncp(:)                                                                   ! subtract the non-conservative product 
                lFbndL(:,j,k) = lFbndL(:,j,k) - sL/(sR-sL)*ncp(:)                                                                   ! add the non-conservative product 
            ENDDO
        ENDDO	
END SUBROUTINE TESTRiemannSolver

RECURSIVE SUBROUTINE InitTECPLOT(N_in,M_in)
	USE TECPLOTPLOTTERmod
	implicit none
	INTEGER :: N_in,M_in
	CALL SetMainParameters(N_in,M_in)
END SUBROUTINE InitTECPLOT





