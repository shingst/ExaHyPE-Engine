! DIM Initial Data


RECURSIVE SUBROUTINE InitialData(xGP, t, Q,maxAMR,CoarseLen, order)
	USE SpecificVarEqn99
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar, nDim, ICType
	IMPLICIT NONE 
	! Argument list 
	INTEGER, INTENT(IN)			   :: maxAMR, order
	REAL, INTENT(IN)               :: xGP(nDim), t,CoarseLen        ! 
	REAL, INTENT(OUT)              :: Q(nVar)        ! 
	! Local variables
	REAL :: up(nVar),LEsigma(6)
	REAL :: ICA, ICsig, ICxd, r, DIsize,c0
	REAL :: rho0, lambda0, mu0, ICx0(3)
	select case(ICType)
		case('SSCRACK')
			! ============================ SELF-SIMILAR RUPTURE TEST =====================================
			
			! ---------- DEFINE THE PARAMETERS -----------
			lambda0=1.3333333333333333333e+10
			mu0=1.333333333333333333333e+10
			rho0=2500.000000
			ICsig=0.5 ! Velocity
			! --------------------------------------------
			
			
			up=0.
			! Set the lamè constants ------ !
			up(15)=lambda0   ! Lambda    !
			up(16)=mu0      ! mu        !
			up(1)= rho0     ! rho       !
			! ----------------------------- !
			! Initial velocity vector ----- !
			up(2) = 0.0                     !
			up(3) = 0.0                     ! 
			up(4) = 0.0                     !
			! ----------------------------- ! 
			LEsigma=0.
			LEsigma(1)=0.!40.0*1.e+6
			LEsigma(2)=40.0*1.e+6!/1.e+6
			LEsigma(3)=0.!40.0*1.e+6
			LEsigma(4)=20.0*1.e+6!/1.e+6
			up(5:13)= GPRsigma2A(LEsigma,0.0,up(15),up(16),up(1),1.e-3)
			! ----------------------------- !
			up(14)=1     ! alpha
			up(17)=1.e+16! Y0
			up(18)=1     ! xi
			! ----------------------------- !
			! Friction parameters -----------
			if(t .eq. 0) then
			SSCRACKFL%mu_s=0.5
			SSCRACKFL%mu_d=0.25
			SSCRACKFL%DynamicFL= .false.
			SSCRACKFL%V=2000.0  ! m/s
			SSCRACKFL%t0=0.0    !s
			SSCRACKFL%L=250     ! m
			
			SSCRACKFL%dx=CoarseLen/(3.0**maxAMR)/(order+1.0)*2.0	! m
			end if
			! ------------------------------
			if(abs(xGP(1)).lt.10000.0) then
				!up(18)=1.0-EXP(-0.5*xGP(2)**2/100.0**2)
				if(abs(xGP(2))<200) then
				!   up(18)=0.0 
				end if
			end if
			!up(20:22)=xGP(1:3)
			USEFrictionLaw=.true.  ! Do not use friction law
			up(24)=LEfriction_mu(xGP(1:3),0.0,up(2:4),0.0,SSCRACKFL)
		case('DRupture')
			! ============================ DYNAMIC RUPTURE TEST CASE =====================================
			! Debug test for dynamic rupture. We introduce two zones with different velocities that generate
			! a stress sufficient to break the material
			

			
			up=0.
			up(2:4) = 0.0	
			if(xGP(1)<0) then! .and. xGP(1)<0.6) then
				up(3) = ICsig  
			else
				up(3) = -ICsig	
			end if
			
			up(5)=1.            ! A_11
			up(9)=1.            ! A_22
			up(13)=1.           ! A_33
			
			up(14)=1.           ! alpha
			
			up(15)=lambda0   ! Lambda
			up(16)=mu0       ! mu
			up(1)=rho0       ! mu
			
			up(18)=1        ! xi, 1 -> not broken, 0 -> broken with friction law
			
			up(17)=1!.e+6     ! Y0
			if(xGP(2)>0) then
				up(17)=1-1.0*EXP(-0.5*(xGP(2)-2.0*xGP(1))**2/0.1**2)
				up(17)=min(up(17),1-1.0*EXP(-0.5*(xGP(2)+2.0*xGP(1))**2/0.1**2))
			else
				up(17)=min(up(17),1-1.0*EXP(-0.5*(xGP(1))**2/0.05**2))   
			end if
			
			! Friction parametersUSEFrictionLaw=.FALSE.  ! Do not use friction law
			up(24)=0.001            ! Set the friction value mu to be constant
			case('StiffInclusion')
			! ============================ STIFF INCLUSION TEST CASE =====================================
			! p-wave travelling through a still material with p-wave speed 10 times larger. This test was
			! taken from gij1 and was modified so the wall is identify by a DI instead of physical boundaries
			!
			! ---------- DEFINE THE PARAMETERS -----------
			ICsig=-0.8
			ICxd=0.01
			ICA=0.0001
			! --------------------------------------------
			
			
			
			
			up=0.
			up(15)=2.0   ! Lambda
			up(16)=1.0   ! Mu
			up(1) =1.0   ! Rho
			
			c0=2.0
			
			LEsigma=(/ 4.0,2.0,2.0,0.0,0.0,0.0 /) ! p-wave sigma
			LEsigma=ICA*LEsigma*exp(-0.5*(xGP(1)-ICsig-c0*t)**2/ICxd**2)
			
			!if(abs(xGP(1))<0.5 .and. abs(xGP(2))<0.1) then
			!	up(15)=200.0
			!	up(16)=100.0
			!end if
			if(abs(xGP(1))<0.49794239582487282 .and. abs(xGP(2))<0.09958848025146097) then
				up(15)=200.0
				up(16)=100.0
			end if
			
			up(5:13)= GPRsigma2A(LEsigma,0.0,up(15),up(16),up(1),1.e-8) ! Assign Aij from sigma and theta
			up(2)	= -2.0*ICA*exp(-0.5*(xGP(1)-ICsig-c0*t)**2/ICxd**2)     	! Assign the velocity
			
			up(14)=1    	! alpha
			up(17)=1.e+16   ! Y0
			up(18)=1     	! xi
			
			r = abs(xGP(2))
			DIsize=0.01 !(ymax-ymin)/JMAX*1.0/(MaxRefFactor**MAXREFLEVEL)/2.0
			up(14)=SmoothInterface(-0.5+r,DIsize,0.0,1.0)
			
			USEFrictionLaw=.FALSE.  ! Do not use friction law
			up(24)=1.0            ! Set the friction value mu to be constant
			! =============================================================================================
		case('TPV3')
			! ============================ SELF-SIMILAR RUPTURE TEST =====================================
			
			! ---------- DEFINE THE PARAMETERS -----------
			lambda0=3.204376e+10
			mu0=3.203812e+10
			rho0=2670.0
			ICsig=0.5 ! Velocity
			! --------------------------------------------
			
			
			up=0.
			! Set the lamè constants ------ !
			up(15)=lambda0   ! Lambda    !
			up(16)=mu0      ! mu        !
			up(1)= rho0     ! rho       !
			! ----------------------------- !
			! Initial velocity vector ----- !
			up(2) = 0.0                     !
			up(3) = 0.0                     ! 
			up(4) = 0.0                     !
			! ----------------------------- ! 
			LEsigma=0.
			LEsigma(1)=120.0*1.e+6
			LEsigma(2)=0.
			LEsigma(3)=0.
			!LEsigma(6)=70.0
			!if(abs(xGP(1)) .le. 1500.0 .and. abs(xGP(2)) .le. 200.0 .and. abs(xGP(3)) .le. 100.0) then
			!	LEsigma(4)=81.6*1.e+6 
			!else
			!	LEsigma(4)=70.0*1e+6
			!end if
			LEsigma(6)=70.0*1.e+6
			r=sqrt((xGP(2)-0.0)**2+(xGP(3)-0.0)**2)
			if(abs(xGP(1))<500.0) then
				if(abs(xGP(2))<3000 .and. abs(xGP(3))<3000) then
					LEsigma(6)=81.6*1.e+6;
				else
						LEsigma(6)=70.0*1.e+6
				end if
				!LEsigma(6)=(70.0+11.6*exp(-r**2/(2.0*6000.0)))*1.e+6
			end if
			up(5:13)= GPRsigma2A(LEsigma,0.0,up(15),up(16),up(1),1.e-3)
			! ----------------------------- !
			up(14)=1     ! alpha
			up(18)=1     ! xi
			! ----------------------------- !
			r=xGP(1)
			up(17)=1.0e+8 - (1.0e+8-75.24*1.e+6)*exp(-r**2/(2.0*1000.0**2))
			
			USEFrictionLaw=.false.  ! Do not use friction law
			up(24)=100.0
			up(25)=0.0 ! Alpha friction

			case('TPV3_2D')
			
				! Set the lamè constants ------ !
				lambda0=3.204376e+10
				mu0=3.203812e+10
				rho0=2670.0
				! ----------------------------- !
				! Set the lamè constants ------ !
				up(15)=lambda0   ! Lambda    !
				up(16)=mu0      ! mu        !
				up(1)= rho0     ! rho       !
				! ----------------------------- !				
				! Initial velocity vector ----- !
				up(2) = 0.0                     !
				up(3) = 0.0                     ! 
				up(4) = 0.0                     !
				! ----------------------------- !
				LEsigma=0.
				LEsigma(1)=0.
				LEsigma(2)=120.0*1.e+6!/1.e+6!0.
				LEsigma(3)=0.!120.0*1.e+6
				if(abs(xGP(1)) .le. 1500.0 .and. abs(xGP(2)) .le. 200.0 .and. abs(xGP(3)) .le. 100.0) then
					LEsigma(4)=81.6*1.e+6 
				else
					LEsigma(4)=70.0*1e+6
				end if
				LEsigma(4)=70.0*1e+6!/1.e+6
				if(abs(xGP(1)) .le. 1500.0) then
				r=xGP(2)
				LEsigma(4)=(70.0*1e+6 + (11.6*1.e+6)*exp(-r**2/(2.0*20.0**2)))!/1.e+6
				end if
				up(5:13)= GPRsigma2A(LEsigma,0.0,up(15),up(16),up(1),1.e-3)
				! ----------------------------- !
				up(14)=1     ! alpha
				if(abs(xGP(2)) .le. 300.0) then
					!up(17)=1.83e+8!1.e+16! Y0
					up(17)=81.24*1.e+6
				else
					up(17)=2.2e+8  
					up(17)=2.2e+9
				end if
				r=xGP(2)
				up(17)=(1.0e+8 - (1.0e+8-75.24*1.e+6)*exp(-r**2/(2.0*200.0**2)))!/1.e+6
				
				up(18)=1     ! xi
				! ----------------------------- !
				USEFrictionLaw=.false.  ! Do not use friction law
				up(24)=10000.0
				up(25)=0.0 ! Alpha friction
				
				up(19:23)=0.0
		case default
			print *, ICType, ' Not implemented'
			stop
	end select
		CALL PDEPrim2Cons(Q,up)
END SUBROUTINE InitialData

RECURSIVE SUBROUTINE PDElimitervalue(limiter_value,xx_in,numberOfObservables, observablesMin, observablesMax)
	USE SpecificVarEqn99
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar, nDim,ICType
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: xx_in(nDim)        ! 
	INTEGER, INTENT(IN)					:: numberOfObservables
	INTEGER, INTENT(OUT)              :: limiter_value        !
	REAL, INTENT(IN)					:: observablesMin(numberOfObservables), observablesMax(numberOfObservables)
	LOGICAL :: dmpresult
	real	:: rr,ldx(3),xx(3)
	xx=0.
	xx(1:nDim)=xx_in(1:nDim)
  limiter_value=0;
return 
  if((observablesMin(1)<0.999 .and. observablesMax(1)>0.001) .or. observablesMax(1)>1.001 .or. observablesMin(1)<-0.001) THEN 
		dmpresult=.FALSE.
   else
		dmpresult=.TRUE.
   ENDIF 
   
   !if(abs(observablesMax(2)-observablesMin(2)) .gt. 1.e-2) then
	!dmpresult=.FALSE.
   !end if
	if(abs(observablesMin(2)) .lt. 1.e-2) then
		dmpresult=.FALSE.
	end if   
   ldx=0.01
   call StaticLimiterEQ99(dmpresult,xx,ldx)

   
   if(dmpresult) then
	limiter_value=0
   else
	limiter_value=1
   end if
END SUBROUTINE PDElimitervalue



RECURSIVE SUBROUTINE PDEGeometriclimitervalue(limiter_value,xx)
	USE SpecificVarEqn99
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar, nDim,ICType
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: xx(nDim)        ! 
	INTEGER, INTENT(OUT)              :: limiter_value        !
	real	:: rr	

	!if(ICType .eq. 'CGeom') then
	!	rr = RadiusFromCG(xx(1),xx(2),xx(3))
	!	!rr = xx(3)-1000.0
	!	if(abs(rr)<200) then
	!		limiter_value=1
	!	else
	!		limiter_value=0
	!	end if
	!	return
	!else
	!	print *, 'Error, geometric limiter not implemented for this test case'
	!	stop
	!end if

END SUBROUTINE PDEGeometriclimitervalue

RECURSIVE SUBROUTINE InitialCG3D(x, t, Q)
    USE, INTRINSIC :: ISO_C_BINDING
	USE SpecificVarEqn99
    USE Parameters, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(3)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 

	INTEGER	:: iErr
    REAL    :: up(nVar), Pi = ACOS(-1.0)
    REAL    :: xi, ICQR(2), ICxd, ICsig
    REAL    :: ICuL(nVar), ICuR(nVar),ICA,r, ICx0(3),nv(3)
    ! Initialize parameters

        up = 0.0
        ! Define the material properties
!        if(x(3) .gt. -1000.0) then
!            ! First medium
!            up(10) = 2.080000e+10  
!            up(11) = 1.040000e+10  
!            up(12) = 2600.0 
!        else
            ! Second medium
            up(10) = 3.240380160000000e+10  
            up(11) = 3.239809920000000e+10
            up(12) = 2700.0
!        end if
        ICsig=200.0
        r = (x(3)+1000)
        !ICsig=10
!        up(10) =0.5*(3.24038016e+10 +2.080000e+10 )+0.5*(3.24038016e+10 -2.080000e+10 )*erf(-r/ICsig)
!        up(11) =0.5*(3.239809920000000e+10 +1.040000e+10)+0.5*(3.239809920000000e+10 -1.040000e+10)*erf(-r/ICsig)
!        up(12) =0.5*(2700.0 +2600.0 )+0.5*(2700.0 -2600.0 )*erf(-r/ICsig)
        up(14)  = 1.0   ! No Fracture everywhere
        r = x(3)
        !r = RadiusFromCG(x(1),x(2),x(3))
		ICsig=50.0
		r = DistanceFromSurfaceCG(x(1),x(2),x(3),ICsig)
		!r = RadiusFromCG(x(1),x(2),x(3)) 
        !r = x(3)-(1000.0+500.0*sin(x(1)))

        up(13)=SmoothInterface(r,ICsig,0.0,0.5)                 ! Get the smooth value of alpha for that value of r
		
		!up(13)=1.0
        ICx0=[0.0, 0.0, -1000.0]
        r = SQRT((x(1)-ICx0(1))**2 + (x(2)-ICx0(2))**2+(x(3)-ICx0(3))**2 ) 
        nv(1) = 0.0
        nv(2) = 0.0
        nv(3) = 1.0
        up(7) = -0.01*EXP(-0.5*r**2/300.0**2) * nv(1)
        up(8) = -0.01*EXP(-0.5*r**2/300.0**2) * nv(2)   
        up(9) = -0.01*EXP(-0.5*r**2/300.0**2) * nv(3) 
		
    CALL PDEPrim2Cons(Q,up)
    !Q=up
    END SUBROUTINE InitialCG3D


	! ================ Specific complex geometries routines =======================
RECURSIVE subroutine ReadCGFile(MyOffset,MyDomain)
		USE	:: SpecificVarEqn99
		USE :: Parameters, only : ICType, ndim
		USE, INTRINSIC :: ISO_C_BINDING
        implicit none
        CHARACTER(LEN=100)         :: CGEOMFile
        integer     :: jcg,nx_cg_new,ny_cg_new,i,j, ix(2)
        real        :: leng(2,2),center(2)
        integer     :: n_new_in(2)
        real, allocatable   :: x_cg_new(:),y_cg_new(:),z_cg_new(:,:)
        real            :: h, phi(4), xi,gamma
        real            :: minx,maxx,miny,maxy
		logical			:: invert_coordinates, binary_input
		real			:: MyOffset(3),MyDomain(3), scalefactor(3)
		! Input parameters
		leng(1:2,1)=MyDomain(1:2)+MyOffset(1:2)
		leng(1:2,2)=-MyOffset(1:2)
		print *, "**********************************************************"
		print *, MyDomain
		print *, MyOffset
		print *, "**********************************************************"
	if(ICType .eq. 'CGeom' .and. ndim .eq. 3) then
		!leng=15000.0
		n_new_in=(/200, 200/)			! Number of elements for the min sub tri function
		
		CGEOMFile="CG.dat"			! DTM file
		center=(/0.0, 0.0/)			! UTM coordinates of the center (with respect to the DTM data file)
		binary_input=.false.
		
		!CGEOMFile="trient_003_44_48_9_13.bin"			! DTM file
		!center=(/4405.905971174,2551.552691730/)			! UTM coordinates of the center (with respect to the DTM data file)
		!binary_input=.true.
		
		!leng=15000.0
		!center=(/600.0, 5110.0/)			! UTM coordinates of the center (with respect to the DTM data file)
		!n_new_in=(/200, 200/)			! Number of elements for the min sub tri function
		!CGEOMFile="alps_01.txt"			! DTM file
		
		if(binary_input) then
			open(8, file=trim(CGEOMFile) ,form='unformatted')
				read(8) nx_cg
				read(8) ny_cg
				allocate(x_cg(nx_cg),y_cg(ny_cg),z_cg(nx_cg,ny_cg))
				read(8) scalefactor(1:3)
				read(8) x_cg(1:nx_cg)
				read(8) y_cg(1:ny_cg)
				read(8) z_cg      
			close(8)			
		else
			open(8, file=trim(CGEOMFile), action='read')
				read(8,*) nx_cg
				read(8,*) ny_cg
				allocate(x_cg(nx_cg),y_cg(ny_cg),z_cg(nx_cg,ny_cg))
				read(8,*) scalefactor(1:3)
				read(8,*) x_cg(1:nx_cg)
				read(8,*) y_cg(1:ny_cg)
				do jcg=1,ny_cg
					read(8,*) z_cg(1:nx_cg,jcg)       
				end do
			close(8)
		end if
		print *, 'Min-Max of z (DTM)=',minval(z_cg), maxval(z_cg)
		center(1)=center(1)*scalefactor(1);
		center(2)=center(2)*scalefactor(2);
		x_cg=x_cg*scalefactor(1)
		y_cg=y_cg*scalefactor(2)
		z_cg=z_cg*scalefactor(3)
        if(y_cg(2)<y_cg(1)) then
			invert_coordinates = .true.
			y_cg(1:nx_cg)=y_cg(nx_cg:1:-1)
			z_cg(nx_cg,1:ny_cg)=z_cg(nx_cg,ny_cg:1:-1)
		else
			invert_coordinates = .false.
		end if        
        ! Connect the input parameters 
        maxx=center(1)+leng(1,1)
        minx=center(1)-leng(1,2)
        maxy=center(2)+leng(2,1)
        miny=center(2)-leng(2,2)
        
        nx_cg_new=n_new_in(1);
        ny_cg_new=n_new_in(2);
        allocate(x_cg_new(nx_cg_new),y_cg_new(ny_cg_new),z_cg_new(nx_cg_new,ny_cg_new))
        h=(maxx-minx)/(nx_cg_new-1)
        x_cg_new(1)=minx
        do i=2,nx_cg_new
            x_cg_new(i)=x_cg_new(i-1)+h   
        end do
        h=(maxy-miny)/(ny_cg_new-1)
        y_cg_new(1)=miny
        do i=2,ny_cg_new
            y_cg_new(i)=y_cg_new(i-1)+h   
        end do

        do i=1,nx_cg_new
            do j=1,ny_cg_new
                ix=lookatindex_cg(x_cg_new(i),y_cg_new(j))
                phi(1)=z_cg(ix(1),ix(2))
                phi(2)=z_cg(ix(1)+1,ix(2))
                phi(3)=z_cg(ix(1),ix(2)+1)
                phi(4)=z_cg(ix(1)+1,ix(2)+1)
                xi=(x_cg_new(i)-x_cg(ix(1)))/(x_cg(ix(1)+1)-x_cg(ix(1)))
                gamma=(y_cg_new(j)-y_cg(ix(2)))/(y_cg(ix(2)+1)-y_cg(ix(2)))
                z_cg_new(i,j)=(1-xi)*(1-gamma)*phi(1)+xi*(1-gamma)*phi(2)+gamma*(1-xi)*phi(3)+xi*gamma*phi(4)
            end do
        end do
        x_cg_new=x_cg_new-center(1) ! The center for the DTM is (0,0) in the fortran code
        y_cg_new=y_cg_new-center(2) ! The center for the DTM is (0,0) in the fortran code
        
        deallocate(x_cg,y_cg,z_cg)
        nx_cg=nx_cg_new;
        ny_cg=ny_cg_new;        
        allocate(x_cg(nx_cg),y_cg(ny_cg),z_cg(nx_cg,ny_cg))
        x_cg=x_cg_new
        y_cg=y_cg_new-243.9     ! Move by 1 element since CG is shifted with respect to the real DTM
		if(invert_coordinates) then
			do i=1,nx_cg
				z_cg(i,1:ny_cg)=z_cg_new(i,ny_cg:1:-1)
			end do
		else
			do i=1,nx_cg
				z_cg(i,1:ny_cg)=z_cg_new(i,1:ny_cg)
			end do
		end if
		
        xL_cg=x_cg(1)
        xR_cg=x_cg(nx_cg)
        yL_cg=y_cg(1)
        yR_cg=y_cg(ny_cg)
        dx_cg=(xR_cg-xL_cg)/nx_cg
        dy_cg=(yR_cg-yL_cg)/ny_cg
        zMin_cg=minval(z_cg)
        zMax_cg=maxval(z_cg)
		deallocate(x_cg_new,y_cg_new,z_cg_new)
		print *, 'Min-Max of x=',minval(x_cg), maxval(x_cg)
		print *, 'Min-Max of y=',minval(y_cg), maxval(y_cg)
		print *, 'Min-Max of z=',minval(z_cg), maxval(z_cg)
		!print *, maxy,miny
		!stop
	end if
end subroutine ReadCGFile



