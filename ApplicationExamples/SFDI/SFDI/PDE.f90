! GRMHD PDE.f90
! Trento (EQNTYPE4)

RECURSIVE SUBROUTINE PDEFlux(f,g,hz,Q)
  USE MainVariables, ONLY : nVar, nDim, EQN
  USE iso_c_binding
  IMPLICIT NONE
  REAL :: f(nVar), g(nVar), h(nVar), hz(nVar), Q(nVar), V(nVar)
  INTENT(IN)  :: Q
  INTENT(OUT) :: f, g, hz
  ! LOCAL VARIABLES
  REAL :: p


	f=0
	g=0
	h=0
	hz=0.
	CALL PDECons2Prim(V,Q)

  p = V(5) 
  ! 
  f(1) = V(6)*V(1)*V(2) 
  f(2) = V(6)*( V(1)*V(2)*V(2) + p ) 
  f(3) = V(6)*V(1)*V(2)*V(3) 
  f(4) = V(6)*V(1)*V(2)*V(4) 
  f(5) = V(2)*(Q(5) + V(6)*p) 
  f(6) = 0. 
  f(7) = 0. 
  f(8) = 0. 
  f(9) = 0. 
  !
  g(1) = V(6)*V(1)*V(3) 
  g(2) = V(6)*V(1)*V(3)*V(2) 
  g(3) = V(6)*( V(1)*V(3)*V(3) + p ) 
  g(4) = V(6)*V(1)*V(3)*V(4) 
  g(5) = V(3)*(Q(5) + V(6)*p) 
  g(6) = 0. 
  g(7) = 0. 
  g(8) = 0. 
  g(9) = 0. 
  ! 
  h(1) = V(6)*V(1)*V(4) 
  h(2) = V(6)*V(1)*V(4)*V(2) 
  h(3) = V(6)*V(1)*V(4)*V(3)    
  h(4) = V(6)*( V(1)*V(4)*V(4) + p ) 
  h(5) = V(4)*( Q(5) + V(6)*p) 
  h(6) = 0. 
  h(7) = 0. 
  h(8) = 0. 
  h(9) = 0. 
  !
  IF(nDim==3) THEN
    hz=h
  END IF  
  END SUBROUTINE PDEFlux


RECURSIVE SUBROUTINE PDENCP(BgradQ,Q,gradQ) 
   USE MainVariables, ONLY :  nVar, nDim, EQN
   IMPLICIT NONE
   ! 11. Oct 21:40: This was a matrix BGradQ(nVar, nDim) but is a vector in spaceTimePredictorNonlinear
   REAL, INTENT(OUT) :: BgradQ(nVar)
   REAL, INTENT(IN)  :: gradQ(nVar, 3)
   REAL, INTENT(IN)  :: Q(nVar)
   ! LOCAL VARIABLES
   REAL				 :: VP(nVar),Qx(nVar),Qy(nVar),Qz(nVar)
   REAL 			 :: AQx(nVar),BQy(nVar),CQz(nVar)
   CALL PDECons2Prim(VP,Q)
	Qx = gradQ(:,1)
	Qy = gradQ(:,2)
	IF(nDim==3) THEN
		Qz = gradQ(:,3)
	ELSE
		Qz = 0.0 
	ENDIF 
     ! Initialize NCP ---!
    AQx=0.              !
    BQy=0.              !
    CQz=0.              !
    ! ------------------!
 
    ! 
   AQx(2) = -Vp(5)*Qx(6) 
   AQx(5) = -Vp(5)*Vp(7)*Qx(6) 
   AQx(6) =  Vp(7)*Qx(6) 
   BQy(3) = -Vp(5)*Qy(6) 
   BQy(5) = -Vp(5)*Vp(8)*Qy(6) 
   BQy(6) =  Vp(8)*Qy(6) 
 
    if( nDim .eq. 2) then
        BgradQ = AQx + BQy         
    else
        BgradQ = AQx + BQy + CQz     
    end if
END SUBROUTINE PDENCP


RECURSIVE SUBROUTINE PDEEigenvalues(L,Q,n)
  USE MainVariables, ONLY :  nVar, nDim, EQN
  USE iso_c_binding
  IMPLICIT NONE
  REAL :: L(nVar), n(3), Q(nVar), V(nVar)
  INTENT(IN)  :: Q,n
  INTENT(OUT) :: L 
  ! Local Variables 
  REAL :: u,vs,c

	IF(MAXVAL(ABS(Q))<1.E-10) THEN
		L=1
		
		RETURN
	END IF
    IF(MAXVAL(ABS(n))<1.E-10) THEN
		print *, 'Impossible error', n
		!stop
	end if
   CALL PDECons2Prim(V,Q)   
   !    
   
   IF(SUM(n**2).EQ.0.) THEN  
     u  = SQRT( SUM(V(2:4)**2) )
     vs = SQRT( SUM(V(7:9)**2) )
   ELSE
     u  = ( V(2)*n(1) + V(3)*n(2) + V(4)*n(3) )  
     vs = ( V(7)*n(1) + V(8)*n(2) + V(9)*n(3) )  
   ENDIF
   c = SQRT( EQN%gamma*V(5)/V(1) ) 
   !IF(V(6)<1e-2) THEN
   !    c = 1.5*c
   !ENDIF
   ! 
   L(1)   = u-c 
   L(2:4) = u
   L(5)   = u+c
   L(6)   = vs 
   L(7:9) = 0. 
   
   
   if(any(isnan(L))) then
	!print *,'Impossible error 3'
	!print *, Q
	!print *, V
	!print *, c,V(5),V(1)
	L=1
	!stop
   end if
   if(maxval(abs(L))<1.e-12) then
	print *, 'Another impossible error', L,u,vs,n
   end if
END SUBROUTINE PDEEigenvalues

RECURSIVE SUBROUTINE PDESource(S,Q) 
  USE MainVariables, ONLY:  nVar, nDim,EQN
  USE iso_c_binding
  IMPLICIT NONE
  ! --------------------------------------------
  ! Argument list declaration
  REAL :: S(nvar), Q(nvar), V(nvar)
  INTENT(IN)  :: Q 
  INTENT(OUT) :: S
  ! --------------------------------------------
  S = 0.
  CALL PDECons2Prim(V,Q) 
  S(2:4) = - EQN%nu* SIGN(1.0, 0.99 - V(6) ) * ( V(2:4) - V(7:9) ) 
  S(5)   = DOT_PRODUCT( V(7:9),S(2:4) ) 
END SUBROUTINE PDESource

RECURSIVE SUBROUTINE PDEVarName(MyNameOUT,ind) 
  USE MainVariables, ONLY: nVar  
  IMPLICIT NONE     
  CHARACTER(LEN=10):: MyName(nVar),MyNameOUT
  INTEGER			:: ind

    MyName(1) = 'rho' 
    MyName(2) = 'u' 
    MyName(3) = 'v'
    MyName(4) = 'w'
    MyName(5) = 'p'   
    MyName(6) = 'alpha' 
    MyName(7) = 'us' 
    MyName(8) = 'vs' 
    MyName(9) = 'ws'
	
	MyNameOUT=MyName(ind+1)
    END SUBROUTINE PDEVarName

RECURSIVE SUBROUTINE PDEAuxName(MyNameOUT,ind) 
	USE MainVariables, ONLY: nVar  
	IMPLICIT NONE     
	CHARACTER(LEN=10):: AuxName(nVar),MyNameOUT
	INTEGER			:: ind


	  ! EQNTYPE99
	  !AuxName(1) = 'sxx'
	  !AuxName(2) = 'syy'
	  !AuxName(3) = 'szz'
	  !AuxName(4) = 'sxy'
	  !AuxName(5) = 'syz'
	  !AuxName(6) = 'sxz'
	  !AuxName(7) = 'YY'
	  !AuxName(8) = 'YYsh'
	  !AuxName(9) = 'YYp'
	  !AuxName(10) = 'tau1'
	  !AuxName(11) = 'tau2'
	  !AuxName(12) = 'LL'
	  !AuxName(13) = 'MM'
	  !AuxName(14) = 'T'
	  !AuxName(15) = 'diff|ADG|'

	MyNameOUT=AuxName(ind+1)
END SUBROUTINE PDEAuxName
	
RECURSIVE SUBROUTINE PDEAuxVar(aux,Q,x,time)
    USE MainVariables, ONLY : nVar,nAux,EQN
	implicit none
	real :: aux(nAux),Q(nVar),x(3),time
	
    
END SUBROUTINE PDEAuxVar
	
	
RECURSIVE SUBROUTINE PDEIntermediateFields(RL,LL,iRL,Q,nv) 
  USE MainVariables, ONLY : nVar, nDim, nLin
  USE iso_c_binding 
  IMPLICIT NONE
  ! Argument list 
  REAL :: RL(nVar,nLin), LL(nLin,nLin), iRL(nLin,nVar)
  REAL :: Q(nVar), nv(3), R(nVar,nVar), iR(nVar,nVar), L(nVar,nVar)
  REAL :: x(nDim), time
  INTENT(IN)  :: Q,nv
  INTENT(OUT) :: RL,LL,iRL 
  ! Local variables

  !  
	  CALL PDEEigenvectors(R,L,iR,Q,nv)
	  RL  = R(:,  (/1,2,3,4,5,6,7/) ) 
	  iRL = iR((/1,2,3,4,5,6,7/), : ) 
	  LL = 0. 
	  LL(1,1) = L(1,1) 
	  LL(2,2) = L(2,2) 
	  LL(3,3) = L(3,3) 
	  LL(4,4) = L(4,4) 
	  LL(5,5) = L(5,5) 
	  LL(6,6) = L(6,6) 
	  LL(7,7) = L(7,7) 
END SUBROUTINE PDEIntermediateFields
    
    
RECURSIVE SUBROUTINE PDEEigenvectors(R,L,iR,Q,nv) 
	USE MainVariables, ONLY : nVar, nDim, EQN, EPSalpha
	USE iso_c_binding  
	IMPLICIT NONE
	! Argument list 
	REAL :: R(nVar,nVar), L(nVar,nVar), iR(nVar,nVar)
	REAL :: Q(nVar), nv(3),VP(nVar)
	REAL :: x(nDim), time
	INTENT(IN)  :: Q,nv
	INTENT(OUT) :: R,L,iR 
	! Local variables
	REAL :: sv(3),tv(3),QPR(nVar),ialpha,dQdV(nVar,nVar),dVdQ(nVar,nVar),RM(nVar,nVar),iRM(nVar,nVar),TM(nVar,nVar)
	integer :: j,k

  IF(ABS(ABS(nv(1))-1.0).LE.1e-14) THEN 
        sv = (/ 0., 1., 0. /) 
        tv = (/ 0., 0., 1. /) 
  ENDIF 
  IF(ABS(ABS(nv(2))-1.0).LE.1e-14) THEN 
        sv = (/ 0., 0., 1. /) 
        tv = (/ 1., 0., 0. /) 
  ENDIF 
  IF(ABS(ABS(nv(3))-1.0).LE.1e-14) THEN 
        sv = (/ 1., 0., 0. /) 
        tv = (/ 0., 1., 0. /) 
  ENDIF

   QPR =  Q 
   QPR(2) =  DOT_PRODUCT(nv, Q(2:4) ) 
   QPR(3) =  DOT_PRODUCT(sv, Q(2:4) ) 
   QPR(4) =  DOT_PRODUCT(tv, Q(2:4) )
   QPR(7) =  DOT_PRODUCT(nv, Q(7:9) ) 
   QPR(8) =  DOT_PRODUCT(sv, Q(7:9) ) 
   QPR(9) =  DOT_PRODUCT(tv, Q(7:9) )


  CALL PDECons2Prim(VP,QPR) 

	  ialpha = VP(6)/(VP(6)**2 + EPSalpha) 

      dQdV(1,1) = VP(6)
      dQdV(1,2) = 0
      dQdV(1,3) = 0
      dQdV(1,4) = 0
      dQdV(1,5) = 0
      dQdV(1,6) = VP(1)
      dQdV(1,7) = 0
      dQdV(1,8) = 0
      dQdV(1,9) = 0
      dQdV(2,1) = VP(6)*VP(2)
      dQdV(2,2) = VP(1)*VP(6)
      dQdV(2,3) = 0
      dQdV(2,4) = 0
      dQdV(2,5) = 0
      dQdV(2,6) = VP(1)*VP(2)
      dQdV(2,7) = 0
      dQdV(2,8) = 0
      dQdV(2,9) = 0
      dQdV(3,1) = VP(6)*VP(3)
      dQdV(3,2) = 0
      dQdV(3,3) = VP(1)*VP(6)
      dQdV(3,4) = 0
      dQdV(3,5) = 0
      dQdV(3,6) = VP(1)*VP(3)
      dQdV(3,7) = 0
      dQdV(3,8) = 0
      dQdV(3,9) = 0
      dQdV(4,1) = VP(6)*VP(4)
      dQdV(4,2) = 0
      dQdV(4,3) = 0
      dQdV(4,4) = VP(1)*VP(6)
      dQdV(4,5) = 0
      dQdV(4,6) = VP(1)*VP(4)
      dQdV(4,7) = 0
      dQdV(4,8) = 0
      dQdV(4,9) = 0
      dQdV(5,1) = VP(6)*(VP(2)**2+VP(3)**2+VP(4)**2)/2
      dQdV(5,2) = VP(1)*VP(6)*VP(2)
      dQdV(5,3) = VP(1)*VP(6)*VP(3)
      dQdV(5,4) = VP(1)*VP(6)*VP(4)
      dQdV(5,5) = VP(6)/(EQN%gamma-1)
      dQdV(5,6) = (EQN%gamma*VP(1)*VP(2)**2+EQN%gamma*VP(1)*VP(3)**2+EQN%gamma*VP(1)*VP(4)**2-VP(1)*VP(2)**2-VP(1)*VP(3)**2-VP(1)*VP(4)**2+2*VP(5))/(EQN%gamma-1)/2
      dQdV(5,7) = 0
      dQdV(5,8) = 0
      dQdV(5,9) = 0
      dQdV(6,1) = 0
      dQdV(6,2) = 0
      dQdV(6,3) = 0
      dQdV(6,4) = 0
      dQdV(6,5) = 0
      dQdV(6,6) = 1
      dQdV(6,7) = 0
      dQdV(6,8) = 0
      dQdV(6,9) = 0
      dQdV(7,1) = 0
      dQdV(7,2) = 0
      dQdV(7,3) = 0
      dQdV(7,4) = 0
      dQdV(7,5) = 0
      dQdV(7,6) = 0
      dQdV(7,7) = 1
      dQdV(7,8) = 0
      dQdV(7,9) = 0
      dQdV(8,1) = 0
      dQdV(8,2) = 0
      dQdV(8,3) = 0
      dQdV(8,4) = 0
      dQdV(8,5) = 0
      dQdV(8,6) = 0
      dQdV(8,7) = 0
      dQdV(8,8) = 1
      dQdV(8,9) = 0
      dQdV(9,1) = 0
      dQdV(9,2) = 0
      dQdV(9,3) = 0
      dQdV(9,4) = 0
      dQdV(9,5) = 0
      dQdV(9,6) = 0
      dQdV(9,7) = 0
      dQdV(9,8) = 0
      dQdV(9,9) = 1

      dVdQ(1,1) = 1*ialpha
      dVdQ(1,2) = 0
      dVdQ(1,3) = 0
      dVdQ(1,4) = 0
      dVdQ(1,5) = 0
      dVdQ(1,6) = -VP(1)*ialpha
      dVdQ(1,7) = 0
      dVdQ(1,8) = 0
      dVdQ(1,9) = 0
      dVdQ(2,1) = -1*ialpha*VP(2)/VP(1)
      dVdQ(2,2) = 1/VP(1)*ialpha
      dVdQ(2,3) = 0
      dVdQ(2,4) = 0
      dVdQ(2,5) = 0
      dVdQ(2,6) = 0
      dVdQ(2,7) = 0
      dVdQ(2,8) = 0
      dVdQ(2,9) = 0
      dVdQ(3,1) = -1*ialpha*VP(3)/VP(1)
      dVdQ(3,2) = 0
      dVdQ(3,3) = 1/VP(1)*ialpha
      dVdQ(3,4) = 0
      dVdQ(3,5) = 0
      dVdQ(3,6) = 0
      dVdQ(3,7) = 0
      dVdQ(3,8) = 0
      dVdQ(3,9) = 0
      dVdQ(4,1) = -1*ialpha*VP(4)/VP(1)
      dVdQ(4,2) = 0
      dVdQ(4,3) = 0
      dVdQ(4,4) = 1/VP(1)*ialpha
      dVdQ(4,5) = 0
      dVdQ(4,6) = 0
      dVdQ(4,7) = 0
      dVdQ(4,8) = 0
      dVdQ(4,9) = 0
      dVdQ(5,1) = 1*ialpha*(VP(2)**2+VP(3)**2+VP(4)**2)*(EQN%gamma-1)/2
      dVdQ(5,2) = -1*ialpha*VP(2)*(EQN%gamma-1)
      dVdQ(5,3) = -1*ialpha*VP(3)*(EQN%gamma-1)
      dVdQ(5,4) = -1*ialpha*VP(4)*(EQN%gamma-1)
      dVdQ(5,5) = 1*ialpha*(EQN%gamma-1)
      dVdQ(5,6) = -1*ialpha*VP(5)
      dVdQ(5,7) = 0
      dVdQ(5,8) = 0
      dVdQ(5,9) = 0
      dVdQ(6,1) = 0
      dVdQ(6,2) = 0
      dVdQ(6,3) = 0
      dVdQ(6,4) = 0
      dVdQ(6,5) = 0
      dVdQ(6,6) = 1
      dVdQ(6,7) = 0
      dVdQ(6,8) = 0
      dVdQ(6,9) = 0
      dVdQ(7,1) = 0
      dVdQ(7,2) = 0
      dVdQ(7,3) = 0
      dVdQ(7,4) = 0
      dVdQ(7,5) = 0
      dVdQ(7,6) = 0
      dVdQ(7,7) = 1
      dVdQ(7,8) = 0
      dVdQ(7,9) = 0
      dVdQ(8,1) = 0
      dVdQ(8,2) = 0
      dVdQ(8,3) = 0
      dVdQ(8,4) = 0
      dVdQ(8,5) = 0
      dVdQ(8,6) = 0
      dVdQ(8,7) = 0
      dVdQ(8,8) = 1
      dVdQ(8,9) = 0
      dVdQ(9,1) = 0
      dVdQ(9,2) = 0
      dVdQ(9,3) = 0
      dVdQ(9,4) = 0
      dVdQ(9,5) = 0
      dVdQ(9,6) = 0
      dVdQ(9,7) = 0
      dVdQ(9,8) = 0
      dVdQ(9,9) = 1

      R(1,1) = 0
      R(1,2) = 0
      R(1,3) = 0
      R(1,4) = (VP(2)-VP(7))**2
      R(1,5) = 1
      R(1,6) = 0
      R(1,7) = 0
      R(1,8) = 1
      R(1,9) = 1
      R(2,1) = 0
      R(2,2) = 0
      R(2,3) = 0
      R(2,4) = -(VP(2)-VP(7))*EQN%gamma*VP(5)/VP(1)**2
      R(2,5) = 0
      R(2,6) = 0
      R(2,7) = 0
      R(2,8) = sqrt(EQN%gamma*VP(1)*VP(5))/VP(1)**2
      R(2,9) = -sqrt(EQN%gamma*VP(1)*VP(5))/VP(1)**2
      R(3,1) = 0
      R(3,2) = 0
      R(3,3) = 0
      R(3,4) = 0
      R(3,5) = 0
      R(3,6) = 1
      R(3,7) = 0
      R(3,8) = 0
      R(3,9) = 0
      R(4,1) = 0
      R(4,2) = 0
      R(4,3) = 0
      R(4,4) = 0
      R(4,5) = 0
      R(4,6) = 0
      R(4,7) = 1
      R(4,8) = 0
      R(4,9) = 0
      R(5,1) = 0
      R(5,2) = 0
      R(5,3) = 0
      R(5,4) = (VP(2)-VP(7))**2*VP(5)*EQN%gamma/VP(1)
      R(5,5) = 0
      R(5,6) = 0
      R(5,7) = 0
      R(5,8) = VP(5)*EQN%gamma/VP(1)
      R(5,9) = VP(5)*EQN%gamma/VP(1)
      R(6,1) = 0
      R(6,2) = 0
      R(6,3) = 0
      R(6,4) = (-VP(1)*VP(2)**2+2*VP(1)*VP(2)*VP(7)-VP(1)*VP(7)**2+VP(5)*EQN%gamma)*VP(6)/VP(1)**2
      R(6,5) = 0
      R(6,6) = 0
      R(6,7) = 0
      R(6,8) = 0
      R(6,9) = 0
      R(7,1) = 1
      R(7,2) = 0
      R(7,3) = 0
      R(7,4) = 0
      R(7,5) = 0
      R(7,6) = 0
      R(7,7) = 0
      R(7,8) = 0
      R(7,9) = 0
      R(8,1) = 0
      R(8,2) = 1
      R(8,3) = 0
      R(8,4) = 0
      R(8,5) = 0
      R(8,6) = 0
      R(8,7) = 0
      R(8,8) = 0
      R(8,9) = 0
      R(9,1) = 0
      R(9,2) = 0
      R(9,3) = 1
      R(9,4) = 0
      R(9,5) = 0
      R(9,6) = 0
      R(9,7) = 0
      R(9,8) = 0
      R(9,9) = 0

      iR(1,1) = 0
      iR(1,2) = 0
      iR(1,3) = 0
      iR(1,4) = 0
      iR(1,5) = 0
      iR(1,6) = 0
      iR(1,7) = 1
      iR(1,8) = 0
      iR(1,9) = 0
      iR(2,1) = 0
      iR(2,2) = 0
      iR(2,3) = 0
      iR(2,4) = 0
      iR(2,5) = 0
      iR(2,6) = 0
      iR(2,7) = 0
      iR(2,8) = 1
      iR(2,9) = 0
      iR(3,1) = 0
      iR(3,2) = 0
      iR(3,3) = 0
      iR(3,4) = 0
      iR(3,5) = 0
      iR(3,6) = 0
      iR(3,7) = 0
      iR(3,8) = 0
      iR(3,9) = 1
      iR(4,1) = 0
      iR(4,2) = 0
      iR(4,3) = 0
      iR(4,4) = 0
      iR(4,5) = 0
      iR(4,6) = VP(1)**2/(-VP(1)*VP(2)**2+2*VP(1)*VP(2)*VP(7)-VP(1)*VP(7)**2+VP(5)*EQN%gamma)*ialpha
      iR(4,7) = 0
      iR(4,8) = 0
      iR(4,9) = 0
      iR(5,1) = 1
      iR(5,2) = 0
      iR(5,3) = 0
      iR(5,4) = 0
      iR(5,5) = -VP(1)/VP(5)/EQN%gamma
      iR(5,6) = 0
      iR(5,7) = 0
      iR(5,8) = 0
      iR(5,9) = 0
      iR(6,1) = 0
      iR(6,2) = 0
      iR(6,3) = 1
      iR(6,4) = 0
      iR(6,5) = 0
      iR(6,6) = 0
      iR(6,7) = 0
      iR(6,8) = 0
      iR(6,9) = 0
      iR(7,1) = 0
      iR(7,2) = 0
      iR(7,3) = 0
      iR(7,4) = 1
      iR(7,5) = 0
      iR(7,6) = 0
      iR(7,7) = 0
      iR(7,8) = 0
      iR(7,9) = 0
      iR(8,1) = 0
      iR(8,2) = VP(1)**2/sqrt(EQN%gamma*VP(1)*VP(5))/2
      iR(8,3) = 0
      iR(8,4) = 0
      iR(8,5) = VP(1)/VP(5)/EQN%gamma/2
      iR(8,6) = -(VP(2)-VP(7))*(sqrt(EQN%gamma*VP(1)*VP(5))*VP(2)-sqrt(EQN%gamma*VP(1)*VP(5))*VP(7)-VP(5)*EQN%gamma)*VP(1)**2/(-VP(1)*VP(2)**2+2*VP(1)*VP(2)*VP(7)-VP(1)*VP(7)**2+VP(5)*EQN%gamma)*ialpha/sqrt(EQN%gamma*VP(1)*VP(5))/2
      iR(8,7) = 0
      iR(8,8) = 0
      iR(8,9) = 0
      iR(9,1) = 0
      iR(9,2) = -VP(1)**2/sqrt(EQN%gamma*VP(1)*VP(5))/2
      iR(9,3) = 0
      iR(9,4) = 0
      iR(9,5) = VP(1)/VP(5)/EQN%gamma/2
      iR(9,6) = -(VP(2)-VP(7))*(sqrt(EQN%gamma*VP(1)*VP(5))*VP(2)-sqrt(EQN%gamma*VP(1)*VP(5))*VP(7)+VP(5)*EQN%gamma)*VP(1)**2/(-VP(1)*VP(2)**2+2*VP(1)*VP(2)*VP(7)-VP(1)*VP(7)**2+VP(5)*EQN%gamma)*ialpha/sqrt(EQN%gamma*VP(1)*VP(5))/2
      iR(9,7) = 0
      iR(9,8) = 0
      iR(9,9) = 0

      L(1,1) = 0
      L(1,2) = 0
      L(1,3) = 0
      L(1,4) = 0
      L(1,5) = 0
      L(1,6) = 0
      L(1,7) = 0
      L(1,8) = 0
      L(1,9) = 0
      L(2,1) = 0
      L(2,2) = 0
      L(2,3) = 0
      L(2,4) = 0
      L(2,5) = 0
      L(2,6) = 0
      L(2,7) = 0
      L(2,8) = 0
      L(2,9) = 0
      L(3,1) = 0
      L(3,2) = 0
      L(3,3) = 0
      L(3,4) = 0
      L(3,5) = 0
      L(3,6) = 0
      L(3,7) = 0
      L(3,8) = 0
      L(3,9) = 0
      L(4,1) = 0
      L(4,2) = 0
      L(4,3) = 0
      L(4,4) = VP(7)
      L(4,5) = 0
      L(4,6) = 0
      L(4,7) = 0
      L(4,8) = 0
      L(4,9) = 0
      L(5,1) = 0
      L(5,2) = 0
      L(5,3) = 0
      L(5,4) = 0
      L(5,5) = VP(2)
      L(5,6) = 0
      L(5,7) = 0
      L(5,8) = 0
      L(5,9) = 0
      L(6,1) = 0
      L(6,2) = 0
      L(6,3) = 0
      L(6,4) = 0
      L(6,5) = 0
      L(6,6) = VP(2)
      L(6,7) = 0
      L(6,8) = 0
      L(6,9) = 0
      L(7,1) = 0
      L(7,2) = 0
      L(7,3) = 0
      L(7,4) = 0
      L(7,5) = 0
      L(7,6) = 0
      L(7,7) = VP(2)
      L(7,8) = 0
      L(7,9) = 0
      L(8,1) = 0
      L(8,2) = 0
      L(8,3) = 0
      L(8,4) = 0
      L(8,5) = 0
      L(8,6) = 0
      L(8,7) = 0
      L(8,8) = (VP(1)*VP(2)+sqrt(EQN%gamma*VP(1)*VP(5)))/VP(1)
      L(8,9) = 0
      L(9,1) = 0
      L(9,2) = 0
      L(9,3) = 0
      L(9,4) = 0
      L(9,5) = 0
      L(9,6) = 0
      L(9,7) = 0
      L(9,8) = 0
      L(9,9) = -(-VP(1)*VP(2)+sqrt(EQN%gamma*VP(1)*VP(5)))/VP(1)

      
    RM  = MATMUL(dQdV,R) 
    iRM = MATMUL(iR,dVdQ) 
      
     
    TM(:,:)    = 0.
    TM(1,1)    = 1.
    TM(5,5)    = 1.
    TM(6,6)    = 1.
    TM(2:4,2)  = nv(:)
    TM(2:4,3)  = sv(:)
    TM(2:4,4)  = tv(:)
    TM(7:9,7)  = nv(:)
    TM(7:9,8)  = sv(:)
    TM(7:9,9)  = tv(:)
      
      
    R  = 0.
    iR = 0. 
    DO j = 1, nVar
     DO k = 1, nVar
        R(:,j) = R(:,j) + TM(:,k)*RM(k,j)
     ENDDO
    ENDDO
    DO k = 1, nVar
     DO j = 1, nVar
       iR(:,j) = iR(:,j) + iRM(:,k)*TM(j,k)
     ENDDO
    ENDDO

 
    END SUBROUTINE PDEEigenvectors 
	
RECURSIVE SUBROUTINE PDEJacobian(An,Q,gradQ,nv)
  USE MainVariables, ONLY : nVar, nDim
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
  USE MainVariables, ONLY : nVar, nDim, EQN
  USE iso_c_binding
  IMPLICIT NONE
  ! Argument list 
  REAL :: An(nVar,nVar)
  REAL :: Q(nVar), nv(3) , Vp(nVar)
  REAL :: A(nVar,nVar),B(nVar,nVar),C(nVar,nVar)
  INTENT(IN)  :: Q,nv
  INTENT(OUT) :: An  
  ! Local variables
  ! Linear elasticity variables

   
    A=0.
    B=0.
    C=0.
	
	CALL PDECons2Prim(Vp,Q)  
    
	A(2,6)    =  -Vp(5)
  A(5,6)    =  -Vp(5)*Vp(7) 
  A(6,6)    =  +Vp(7)  
  B(3,6)    =  -Vp(5)
  B(5,6)    =  -Vp(5)*Vp(8) 
  B(6,6)    =  +Vp(8)  
  C(4,6)    =  -Vp(5) 
  C(5,6)    =  -Vp(5)*Vp(9) 
  C(6,6)    =  +Vp(9)  
  !
  An = A*nv(1) + B*nv(2) + C*nv(3) 	
	
    if( nDim .eq. 2) then
        An = A*nv(1) + B*nv(2)         
    else
        An = A*nv(1) + B*nv(2) + C*nv(3)     
    end if

END SUBROUTINE PDEMatrixB

RECURSIVE SUBROUTINE RoeMatrix(ARoe,QL,QR,nv) 
  USE MainVariables, ONLY : nVar, nDim,sGP3,wGP3,nGP3
  IMPLICIT NONE
  ! Argument list 
  REAL        :: ARoe(nVar,nVar), QL(nVar), QR(nVar), nv(3) 
  INTENT(IN)  :: QL, QR, nv  
  INTENT(OUT) :: ARoe 
  ! Local variables 
  INTEGER             :: iGP  
  REAL                :: psi(nVar) 
  REAL                :: A(nVar,nVar)
  REAL                :: gradQ(nVar,3) 
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

RECURSIVE SUBROUTINE HLLEMFluxFV(smax,FL,FR,QL,QR,NormalNonZero) 
  USE MainVariables, ONLY : nVar, nDim, nLin
  USE iso_c_binding 
  ! Local variables
  INTEGER, INTENT(IN)   :: NormalNonZero
  REAL, INTENT(IN)     :: QL(nVar)
  REAL, INTENT(IN)     :: QR(nVar)
  REAL, INTENT(INOUT)  :: FL(nVar)
  REAL, INTENT(INOUT)  :: FR(nVar)
  REAL, INTENT(OUT)    :: smax
  REAL    :: QavL(nVar), QavR(nVar)  
    ! Local variables 
  INTEGER           :: i,j,k,l, ml(1)  ,iErr
  REAL              :: Qav(nVar)!,LambdaOut
  REAL              ::  nv(3), flattener(nLin)
  REAL    :: absA(nVar,nVar), amax  ,gradQ(nVar,3), ncp(nVar)
  REAL    :: QM(nVar),LL(nVar),LR(nVar),LM(nVar)
  REAL    :: deltaL(nLin,nLin),Lam(nLin,nLin),Lap(nLin,nLin) 
  REAL    :: RL(nVar,nLin),iRL(nLin,nVar),LLin(nLin,nLin) , TMPM(nLin, nVar),TMPM2(nVar,nVar)
  REAL    :: Aroe(nVar,nVar),Aroep(nVar,nVar), Aroem(nVar,nVar), Dm(nVar), Dp(nVar), dQ(nVar)
  REAL :: f1R(nVar), g1R(nVar), h1R(nVar) ,flux(nVar),flux2(nVar)
  REAL :: f1L(nVar), g1L(nVar), h1L(nVar) , VM(nVar) 
  REAL :: minLL,minLM,maxLR,maxLM
  REAL  :: QRL(nVar),fR2(nVar),fL2(nVar)
  
  REAL :: XX0(3),TIME0
  XX0=0.
  TIME0=0.
  !  
  nv(:)=0.
  nv(NormalNonZero+1)=1.
  !
  flattener=1.!0.8
  if(min(QL(6),QR(6))<0.02) then
	flattener(6)=0.5
  end if
  !
CALL PDEFlux(f1L,g1L,h1L,QL)
CALL PDEFlux(f1R,g1R,h1R,QR)
!
fR = f1R*nv(1)+g1R*nv(2)+h1R*nv(3)
fL = f1L*nv(1)+g1L*nv(2)+h1L*nv(3)
  
  QM=0.5*(QL+QR)
  
   CALL PDEEigenvalues(LL,QL,nv)  
  CALL PDEEigenvalues(LR,QR,nv)  
  CALL PDEEigenvalues(LM,QM,nv) 
  minLL=MINVAL(LL(:))
  minLM=MINVAL(LM(:))
  maxLR=MAXVAL(LR(:))
  maxLM=MAXVAL(LM(:))  
  smax = MAX( MAXVAL(ABS(LL)), MAXVAL(ABS(LR)) ) 
  sL = -smax
  sR = +smax 
  !sL = MIN( 0., minLL , minLM ) 
  !sR = MAX( 0., maxLR , maxLM )  
  !LAMBDA = smax ! max(abs(sL),abs(sR))
  CALL PDEIntermediateFields(RL,LLin,iRL,QM,nv) 
  Lam = 0.5*(LLin-ABS(LLin))
  Lap = 0.5*(LLin+ABS(LLin)) 

  deltaL = 0.0
  DO i = 1, nLin
      deltaL(i,i) = ( 1. - Lam(i,i)/(sL-1e-14) - Lap(i,i)/(sR+1e-14) )*flattener(i)  
  ENDDO  
  amax = 0. 

  absA = 0. 
  DO i = 1, nVar
      absA(i,i) = sR*sL/(sR-sL)  - 0.5*amax ! regular HLL diffusion, only on the diagonal 
  ENDDO  
  TMPM=MATMUL(deltaL, iRL)
  TMPM2=MATMUL( RL,TMPM)
  !absA = absA - sR*sL/(sR-sL)*TMPM2 ! HLLEM anti-diffusion  
  !    
  !PRINT *, "RoeMatrix"
  CALL RoeMatrix(ARoe,QL,QR,nv)
  !PRINT *, "RoeMatrix done!"
  !
  ARoem = -sL*ARoe/(sR-sL)
  ARoep = +sR*ARoe/(sR-sL)
  ! 

  dQ = QR - QL
  flux = (sR*fL - sL*fR)/(sR-sL) + MATMUL( absA, QR - QL ) 
  !
  !Dp = -MATMUL(Aroep,dQ)
  !Dm = -MATMUL(Aroem,dQ)        ! these are the path integral of the MatrixB from QL to QR. (the NCP as a first approximation)
  !
  gradQ=0.
  gradQ(:,NormalNonZero+1) = QR(:) - QL(:) 
  CALL PDENCP(ncp,QM,gradQ)
  Dp = MATMUL(Aroep,dQ)
  Dm = MATMUL(Aroem,dQ)        ! these are the path integral of the MatrixB from QL to QR. (the NCP as a first approximation)
  !
  !if(abs(flux(21))>1.e-0) then
  !   print *, '*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
	!	print *,flux
	!	print *,'---------------------'
	!	print *,fl
	!	print *,'---------------------'
	!	print *,fr
	!	print *,'---------------------'
	!	print *,Ql
	!	print *,'---------------------'
	!	print *,qr
  !     print *, '====================='
  !     print *, sR
  !     print *, '====================='
  !     print *, sL
	! print *, '====================='
	!	print *,'---------------------'
	!	print *,deltaL
	!	print *,'---------------------'
	!	print *, MATMUL( absA, QR - QL ) 
	!	print *,'---------------------'
	!	print *, sR*sL/(sR-sL)
	!	print *,'---------------------'
	!	print *, absA(18:nVar,:)
	!	print *,'---------------------'
	!	print *, absA(:,:)
	!	print *, '*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
	!	pause
  !end if
  !flux = 0.5*( fR + fL ) + MATMUL(absA, QR - QL)
  
  !fR2=fR
  !fL2=fL
  
  fR = flux - Dp
  fL = flux + Dm
  

  
 !fR(18:nVar) = 0.- 0.5*ncp(18:nVar)
 !fL(18:nVar) = 0.+ 0.5*ncp(18:nVar) 
  if(any(isnan(fR))) then
	print *, 'Issue here'
    print *,'---------------------'
    print *,QL
    print *, '====================='
    print *,QR
    print *, '====================='
    print *, fR
    print *, '====================='
    print *, h1R
    print *, '====================='
    print *, sR
    print *, '====================='
    print *, sL
	 print *, '====================='
	print *, deltaL
    print *, '***********************'
    print *, flux2
    print *, '***********************'
    print *, flux
	print *, '***********************'
    !print *, Lam
	!    print *, '***********************'
    !print *, Lam
    print *,'---------------------'	
	pause
  end if
    END SUBROUTINE HLLEMFluxFV	



