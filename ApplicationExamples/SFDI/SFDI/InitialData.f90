! DIM Initial Data


RECURSIVE SUBROUTINE InitParameters(STRLEN,PARSETUP) 
	USE MainVariables, ONLY : nVar, nDim, EQN, ICType, MHDRotomega, EPRadius,ICdelta, ICRot,ICuL,ICA,EPCenter, NACA,nBlades,pitcha
	IMPLICIT NONE  
	integer          :: STRLEN
	character(len=STRLEN) :: PARSETUP
	
	EQN%Pi=acos(-1.0)
	ICType=trim(parsetup)
	print *, "****************************************************************"
	print *, 'Chosen setup: 	',ICType
	print *, "****************************************************************"
	
	select case(ICType)
		case('NBodies')
			EQN%gamma        =  1.4 
			EQN%nu           =  0.0 
			EQN%mu           =  0.0 
			EQN%cv           =  2.5 
			EQN%Pr           =  0.75
			EQN%R     = EQN%cv*(EQN%gamma-1.)
			EQN%kappa = EQN%gamma*EQN%cv*EQN%mu/EQN%Pr
			
			ICuL(1:5)        =   (/ 1.4 , 0.0 , 0.0 , 0.0 , 1.0 /)
			ICrot            =   3
			EPCenter(1:2) =   (/ 0.0 , 0.0 /)
			EPRadius         =   1.0   
			ICA              =   0.2   
			ICdelta          =   0.01 
			MHDRotomega      =   3.0   
		case('DoubleRotor3D')
			EQN%gamma        =  1.4 
			EQN%nu           =  0.0 
			EQN%mu           =  0.0 
			EQN%cv           =  2.5 
			EQN%Pr           =  0.75
			EQN%R     = EQN%cv*(EQN%gamma-1.)
			EQN%kappa = EQN%gamma*EQN%cv*EQN%mu/EQN%Pr
			
			ICuL(1:5)        =   (/ 1.4 , 1.0 , 0.0 , 0.0 , 1.0 /)
			NACA(1:4)        =   (/ 2, 4, 1, 2 /)
			nBlades  		 =   5
			pitcha           =30.00
			EPCenter         =   0.
			EPRadius         =   1.0   
			ICA              =   0.2   
			ICdelta          =   0.001 
			MHDRotomega     =   1.0 
			!RotOmegaV(1:3)   = (/ 1.0, 0.0, 0.0 /)			
	END SELECT	
	
END SUBROUTINE InitParameters

function SmoothInterface(r,ICsig,epsilon,smooth_order_in)
        implicit none
        real    :: r
        real    :: SmoothInterface,smooth_order
        real    :: eta,ICsig,xi 
        real    :: epsilon
        real    :: smooth_order_in
        
        !if(.not. present(epsilon)) then
        !    epsilon=1.e-9    
        !end if
        !if(.not. present(smooth_order_in)) then
        !    smooth_order=4   
        !end if
        smooth_order=smooth_order_in

        eta=0.0
        smooth_order=2.0
        ! =============== WAY 1 ================================
        if(r>(1+eta)*ICsig) then
            xi=1    
        elseif(r<-(1-eta)*ICsig) then
            xi=0
        else
            xi = (r+(1-eta)*ICsig)/(2.0*ICsig) 
        end if
        ! =============== WAY 2 ================================
        SmoothInterface  = (1.0)*(1-xi)**smooth_order + (0.0+epsilon)*xi**smooth_order 
        if(abs(smooth_order) .le. 1.e-3) then
            xi = 0.5+0.5*ERF(r/(2*ICsig))  
            SmoothInterface  = (1.0-epsilon)*(1-xi) + (0.0+epsilon)*xi
        end if
    end function SmoothInterface

RECURSIVE SUBROUTINE InitialData(xGP, tGp, Q)
	USE, INTRINSIC :: ISO_C_BINDING
	USE MainVariables, ONLY : nVar, nDim, EQN, ICType, MHDRotomega, EPRadius,ICdelta, ICRot, ICuL,ICA,EPCenter, NACA,nBlades,pitcha
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: xGP(3), tGp        ! 
	REAL, INTENT(OUT)              :: Q(nVar)        ! 
	REAL :: up(nVar), r,tmp
	REAL :: lm,lp,lt,EPS,theta,XX,YY,ZZ,chord,yt,yc,dydx,xxU,xxL,yyU,yyL,locDISIZE
	REAL :: diX1,diX2,diY1,diY2,diZ1,diZ2
	real, external :: SmoothInterface
	INTEGER :: i
	
	select case(ICType)
		case('NBodies')
			up(1:5) = ICuL(1:5)
			up(6)   = 1.0 - ICdelta
			DO i = 1, ICrot     
				r = SQRT( (xGP(1)-EPRadius*COS(i*2*EQN%Pi/ICrot))**2+(xGP(2)-EPRadius*SIN(i*2*EQN%Pi/ICRot))**2) 
				tmp=SmoothInterface(-r+ICA,0.02,ICdelta,1.0)
				up(6)=min(up(6),tmp)	
				!IF(r<=2*ICA) THEN
                 !   up(6)=SmoothInterface(-r+ICA,0.02,ICdelta,1.0)				
				!	!up(6) = ICdelta
				!ENDIF 
			ENDDO                           
			up(7)   =  MHDRotomega/EPRadius*xGP(2) 
			up(8)   = -MHDRotomega/EPRadius*xGP(1) 
			up(9)   = 0.0 
			!up(6)   = 1.0
		case('DoubleRotor3D')
		    locDISIZE=0.05
			up(1:5) = ICuL(1:5)
			up(6)   = 1.0 - ICdelta
			r = SQRT( (xGP(2)-EPCenter(2))**2 + (xGP(3)-EPCenter(3))**2 )
			eps = 1e-20
			lm = NACA(1)/100.
			lp = NACA(2)/10.
			lt = (10*NACA(3)+NACA(4))/100. 
			IF(xGP(1)>=-1.25 .and. xGP(1)<=-1.0) THEN
				IF(r < ICA*SQRT(xGP(1)+1.25)/SQRT(0.25)) THEN
					up(6) = ICdelta
				ENDIF            
			ENDIF
			IF( xGP(1)>=-1.0 .and. xGP(1)<=1.0 ) THEN
				!tmp=SmoothInterface(-r+ICA,locDISIZE,ICdelta,1.0)
				!up(6)=min(up(6),tmp)
				IF(r < ICA) THEN
					up(6) = ICdelta
				ENDIF            
			ENDIF
			IF(xGP(1)>=1.0 .and. xGP(1)<=1.25) THEN
			    !tmp=SmoothInterface(-r+ICA*SQRT(1.25-xGP(1))/SQRT(0.25),locDISIZE,ICdelta,1.0)
				!up(6)=min(up(6),tmp)
				IF(r < ICA*SQRT(1.25-xGP(1))/SQRT(0.25)) THEN
					up(6) = ICdelta
				ENDIF            
			ENDIF
			!
			! coordinate transformation to the blade-relative coordinate system  
			!
			IF(xGP(1)>=0) THEN
			    
					DO i = 1, nBlades
						!
						theta = 2*EQN%Pi*REAL(i-1)/REAL(nBlades)  
						!
						xx =  xGP(1) - 0.3  
						yy =  xGP(2)*COS(theta) + xGP(3)*SIN(theta)  
						zz = -xGP(2)*SIN(theta) + xGP(3)*COS(theta)    
						!
						chord = 0.5 - (0.5-0.2)/EPRadius*zz            
						!
						if(xx>=0.0 .and. xx<chord ) then
							yt = 5*lt*chord*(0.2969*sqrt(xx/chord)-0.126*xx/chord-0.3516*(xx/chord)**2+0.2843*(xx/chord)**3-0.1015*(xx/chord)**4)
							IF(xx <= lp*chord) THEN
								yc = lm*(2*lp*xx/chord-(xx/chord)**2)/(lp**2+eps)
								dydx = 2*lm*(lp-xx/chord)/(lp**2+eps)  
							ELSE
								yc = lm*(1-2*lp+2*lp*xx/chord-(xx/chord)**2)/(1-lp)**2
								dydx = 2*lm*(lp-xx/chord)/(1-lp)**2  
							ENDIF
							theta = ATAN(dydx)
							xxU = xx-yt*SIN(theta)
							xxL = xx+yt*SIN(theta)
							yyU = yc+yt*COS(theta)
							yyL = yc-yt*COS(theta)
							
							diX1=SmoothInterface(-zz+EPRadius,locDISIZE,ICdelta,1.0)
							diX2=SmoothInterface(+zz-0,locDISIZE,ICdelta,1.0)
							diY1=SmoothInterface(-xx+chord,locDISIZE,ICdelta,1.0)
							diY2=SmoothInterface(+xx-0,locDISIZE,ICdelta,1.0)
							diZ1=SmoothInterface(-yy+yyU,locDISIZE,ICdelta,1.0)
							diZ2=SmoothInterface(+yy-yyL,locDISIZE,ICdelta,1.0)
							tmp=max(diX1,diX2)
							tmp=max(diY1,tmp)
							tmp=max(diY2,tmp)
							tmp=max(diZ1,tmp)
							tmp=max(diZ2,tmp)
							!up(6) = min(up(6), tmp)
							if(isnan(up(6))) then
								print *, '(NaN) up6',yt,lt,chord,xx
								!pause
							end if
							if(.not. (zz > 0 .AND. zz < EPRadius .AND. xx>=0.0 .AND. xx<=chord .AND. yy<=yyU .AND. yy>=yyL) .and. tmp<0.02) then
								print *, '(tmp issue) up6',tmp,xx
								
								!pause
							end if
							IF(zz > 0 .AND. zz < EPRadius) THEN
								IF( xx>=0.0 .AND. xx<=chord .AND. yy<=yyU .AND. yy>=yyL ) THEN
										up(6) = ICdelta
								ENDIF      
							ENDIF
						end if
					ENDDO 
				
				up(7)   = 0.0
				up(8)   =  MHDRotomega/EPRadius*xGP(3)
				up(9)   = -MHDRotomega/EPRadius*xGP(2)
           
      ELSE
           DO i = 1, nBlades
               !
               theta = 2*EQN%Pi*REAL(i-1)/REAL(nBlades)  
               !
               xx = -xGP(1) - 0.3  
               yy =  xGP(2)*COS(theta) + xGP(3)*SIN(theta)  
               zz = -xGP(2)*SIN(theta) + xGP(3)*COS(theta)    
               !
               chord = 0.5 - (0.5-0.2)/EPRadius*zz            
               !
			   if(xx>=0.0 .and. xx<chord ) then
               yt = 5*lt*chord*(0.2969*sqrt(xx/chord)-0.126*xx/chord-0.3516*(xx/chord)**2+0.2843*(xx/chord)**3-0.1015*(xx/chord)**4)
               IF(xx <= lp*chord) THEN
                   yc = lm*(2*lp*xx/chord-(xx/chord)**2)/(lp**2+eps)
                   dydx = 2*lm*(lp-xx/chord)/(lp**2+eps)  
               ELSE
                   yc = lm*(1-2*lp+2*lp*xx/chord-(xx/chord)**2)/(1-lp)**2
                   dydx = 2*lm*(lp-xx/chord)/(1-lp)**2  
               ENDIF
               theta = ATAN(dydx)
               xxU = xx-yt*SIN(theta)
               xxL = xx+yt*SIN(theta)
               yyU = yc+yt*COS(theta)
               yyL = yc-yt*COS(theta)
               IF(zz > 0 .AND. zz < EPRadius) THEN
                   IF( xx>=0.0 .AND. xx<=chord .AND. yy<=yyU .AND. yy>=yyL ) THEN
                        up(6) = ICdelta
                   ENDIF      
               ENDIF
			   end if
           ENDDO    
           up(7)   = 0.0
           up(8)   = -MHDRotomega/EPRadius*xGP(3)
           up(9)   = +MHDRotomega/EPRadius*xGP(2)
           
      ENDIF 
	END SELECT
	CALL PDEPrim2Cons(Q,up)
END SUBROUTINE InitialData

RECURSIVE SUBROUTINE PDElimitervalue(limiter_value,xx,numberOfObservables, observablesMin, observablesMax)
	USE, INTRINSIC :: ISO_C_BINDING
	USE MainVariables, ONLY : nVar, nDim
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: xx(3)        ! 
	INTEGER, INTENT(IN)					:: numberOfObservables
	INTEGER, INTENT(OUT)              :: limiter_value        !
	REAL, INTENT(IN)					:: observablesMin(numberOfObservables), observablesMax(numberOfObservables)
	real	:: rr	

	limiter_value=0
   IF(observablesMin(1)<1.e-2) then
      limiter_value = 1
   end if

   !print *, observablesMax(6)
	
   IF((observablesMax(6)<0.985 .and. observablesMin(6)>0.015) .or. observablesMax(6)>1.005) THEN  ! was -0.001   .or. levelmin<0.001
		limiter_value = 1 
   ENDIF 
   !rr=sqrt(sum(xx**2))
   !if(rr<1.8) then
   !limiter_value = 1
   !end if
   !
	!limiter_value = 0
END SUBROUTINE PDElimitervalue

  


