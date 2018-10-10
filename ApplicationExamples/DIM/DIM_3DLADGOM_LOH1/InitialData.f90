! DIM Initial Data


RECURSIVE SUBROUTINE InitialData(x, t, Q)
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar, nDim, ICFlag
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: x(nDim), t        ! 
	REAL, INTENT(OUT)              :: Q(nVar)        ! 

	!Call InitialPlaneWave(x, t, Q)
	!call GaussianBubble(x,t,Q)
	select case(ICFlag)
		case('CGeom')
			if(nDim .eq. 3) then
				call InitialCG3D(x, t, Q);
			else
				print *, ICFlag, ' Not implemented'
				stop
			end if
		case('LOH1')
			if(nDim .eq. 3) then
				call InitialLOH13D(x, t, Q);
			else
				print *, ICFlag, ' Not implemented'
				stop
			end if
		case('PlaneWave')
			call InitialPlaneWave(x, t, Q)
		case('SmoothTC')
			call GaussianBubble(x,t,Q)
		case default
			print *, ICFlag, ' Not implemented'
	end select
		
END SUBROUTINE InitialData

RECURSIVE SUBROUTINE PDElimitervalue(limiter_value,xx,numberOfObservables, observablesMin, observablesMax)
	USE SpecificVarEqn99
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar, nDim,ICFlag
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: xx(nDim)        ! 
	INTEGER, INTENT(IN)					:: numberOfObservables
	INTEGER, INTENT(OUT)              :: limiter_value        !
	REAL, INTENT(IN)					:: observablesMin(numberOfObservables), observablesMax(numberOfObservables)
	real	:: rr	

	! Plane Wave limiter_value
	!if(ICFlag .eq. 'PlaneWave') then
	!rr =sqrt(sum(xx(:)**2))
	!rr=0.1-rr
	!if(abs(rr).le. 0.25) then
	!	limiter_value=1
	!else
	!	limiter_value=0
	!end if
	!limiter_value=0
	!return
	!end if
	!
	
	! Exclude the boundaries
	!if(abs(xx(1))>2000 .or. abs(xx(2))>2000) then
	!	limiter_value=0
	!	return
	!end if
	!rr = RadiusFromCG(xx(1),xx(2),xx(3))
	!rr = DistanceFromSurfaceCG(xx(1),xx(2),xx(3),50.0)
	!rr=1.00
	!rr=xx(3)-2000
	!limiter_value=0
	!if(abs(rr)<500) then
	!    limiter_value=1
	!	return
	!end if
	if(ICFlag .eq. 'CGeom') then
		rr = RadiusFromCG(xx(1),xx(2),xx(3))
		!rr = xx(3)-1000.0
		if(abs(rr)<1000) then
			limiter_value=1
		else
			limiter_value=0
		end if
		return
	end if
	
   if((observablesMin(1)<0.999 .and. observablesMax(1)>0.001) .or. observablesMax(1)>1.001 .or. observablesMin(1)<-0.001) THEN 
       limiter_value=1
   else
		limiter_value=0
   ENDIF 
	!limiter_value=0
END SUBROUTINE PDElimitervalue



RECURSIVE SUBROUTINE PDEGeometriclimitervalue(limiter_value,xx)
	USE SpecificVarEqn99
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar, nDim,ICFlag
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: xx(nDim)        ! 
	INTEGER, INTENT(OUT)              :: limiter_value        !
	real	:: rr	

	if(ICFlag .eq. 'CGeom') then
		rr = RadiusFromCG(xx(1),xx(2),xx(3))
		!rr = xx(3)-1000.0
		if(abs(rr)<200) then
			limiter_value=1
		else
			limiter_value=0
		end if
		return
	elseif(ICFlag .eq. 'LOH1') then
		if(xx(2).le. 0 .and. xx(2)>-0.5) then
			limiter_value=1
		else
			limiter_value=0
		end if
		!limiter_value=0
		return
	elseif(ICFlag .eq. 'PlaneWave') then
		rr =sqrt(sum(xx(1:nDim)**2))
		if(abs(rr)<0.25) then
			limiter_value=1
		else
			limiter_value=0
		end if
	else
		print *, 'Error, geometric limiter not implemented for this test case'
		stop
	end if

END SUBROUTINE PDEGeometriclimitervalue



RECURSIVE SUBROUTINE InitialLOH13D(x, t, Q)
    USE, INTRINSIC :: ISO_C_BINDING
	USE SpecificVarEqn99
    USE Parameters, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 

	INTEGER	:: iErr
    REAL    :: up(nVar), Pi = ACOS(-1.0)
    REAL    :: xi, ICQR(2), ICxd, ICsig
    REAL    :: ICuL(nVar), ICuR(nVar),ICA,r, ICx0(3),nv(3),parL(3),parR(3)
    ! Initialize parameters

        up = 0.0
        ! Define the material properties
        if(x(1) < 1.0) then
            ! First medium
            up(10) = 2.080000e+01  
            up(11) = 1.040000e+01  
            up(12) = 2.600000
			
			up(10) = 3.685150e+01  
            up(11) = 3.017425e+01
            up(12) = 2.700000
        else
            ! Second medium
            up(10) = 3.685150e+01  
            up(11) = 3.017425e+01
            up(12) = 2.700000
        end if
		
            parR(1) = 2.080000e+01  
            parR(2) = 1.040000e+01  
            parR(3) = 2.600000
            parL(1) = 3.685150e+01  
            parL(2) = 3.017425e+01
            parL(3) = 2.700000
			!parL=parR
        r = (x(2)-1)
        ICsig=0.1
        !up(10) =0.5*(3.24038016e+10 +2.080000e+10 )+0.5*(3.24038016e+10 -2.080000e+10 )*erf(-r/ICsig)
!        up(11) =0.5*(3.239809920000000e+10 +1.040000e+10)+0.5*(3.239809920000000e+10 -1.040000e+10)*erf(-r/ICsig)
        up(10:12) =0.5*(parR(:) +parL(:) )+0.5*(parR(:) -parL(:) )*erf(-r/ICsig)
        !up(14)  = 1.0   ! No Fracture everywhere
        
		r = x(2)-0.0
        !r = RadiusFromCG(x(1),x(2),x(3))
		ICsig=0.2
		!r = RadiusFromCG(x(1),x(2),x(3)) 
        r = x(3)-12.5
		!up(13)=1.0
        !up(13)=SmoothInterface(r,ICsig,0.0,1)                 ! Get the smooth value of alpha for that value of r
        if(x(2)<0) then
			up(13)=0.0
			up(1:11)=up(1:11)*(up(13)+1.e-14)
		else
			up(13)=1.0
		end if
		!ICx0=[7.5, 7.5, 7.5]
        !r = SQRT((x(1)-ICx0(1))**2 + (x(2)-ICx0(2))**2+(x(3)-ICx0(3))**2 ) 
        !nv(1) = 0.0
        !nv(2) = 0.0
        !nv(3) = 1.0
        !up(7) = -0.01*EXP(-0.5*r**2/1.0**2) * nv(1)
        !up(8) = -0.01*EXP(-0.5*r**2/1.0**2) * nv(2)   
        !up(9) = -0.01*EXP(-0.5*r**2/1.0**2) * nv(3) 
		
		
    CALL PDEPrim2Cons(Q,up)
    !Q=up
    END SUBROUTINE InitialLOH13D


RECURSIVE SUBROUTINE InitialCG3D(x, t, Q)
    USE, INTRINSIC :: ISO_C_BINDING
	USE SpecificVarEqn99
    USE Parameters, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
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

        up(13)=SmoothInterface(r,ICsig,0.0,1)                 ! Get the smooth value of alpha for that value of r
		
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

RECURSIVE SUBROUTINE InitialPlaneWave(x, t, Q)
    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim
	USE	:: SpecificVarEqn99
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 

	INTEGER	:: iErr
    REAL    :: up(nVar), Pi = ACOS(-1.0)
    REAL    :: xi, ICQR(2), ICxd, ICsig
    REAL    :: ICuL(nVar), ICuR(nVar),ICA,r
    ! Initialize parameters
    ICuL(:)=(/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0, 1.0, 1.0, 1.0 /)
    ICuR(:)=(/ 0.4, 0.2, 0.2, 0.0, 0.0, 0.0, -0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
    ICA=1.0     ! wave length
    ICxd=0.25   ! radius of the circle
    ICsig=1e-2  ! smoothing parameter
    ICQR(:)= (/0.0, 1.0 /)
    ! Initialize the variables vector V
    up = ICuL + ICuR*SIN(2*Pi*(x(1)-2.0*t)/ICA)
    r = SQRT(sum(x(1:nDim)**2)) 
	! Way 0 to compute alpha
    xi = 0.5+0.5*ERF((r-ICxd)/ICsig)
    up(13)  = ICQR(1)*(1-xi) + ICQR(2)*xi 
	up(13)  = 1.0
	! Way one to compute alpha
	ICsig=1.e-2  ! smoothing parameter
	r=ICxd-r
	up(13)=SmoothInterface(r,ICsig,0.0,0)
    up(13)=1.0
	
    !up(1:9) = up(1:9)*up(13)
    up(14)=1.0
	up(1:6)=up(1:6)*up(13)

    CALL PDEPrim2Cons(Q,up)
    !Q=up
    END SUBROUTINE InitialPlaneWave
    
RECURSIVE SUBROUTINE GaussianBubble(x, t, Q)
    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 

    REAL    :: up(nVar), Pi = ACOS(-1.0)
    REAL    :: xi, ICQR(2), ICxd, ICsig
    REAL    :: ICuL(nVar), ICuR(nVar),ICA,r

        up=0;
        r = SQRT((x(1))**2 + (x(2))**2)
        up(1:2)=exp(-40.0*r**2)
        up(13)=1.0
        up(14)=1.0
        up(10) = 2.0  
        up(11) = 1.0  
        up(12) = 1.0 
    
    !CALL PDEPrim2Cons(Q,up)
        Q=up
END SUBROUTINE GaussianBubble
	! ================ Specific complex geometries routines =======================
RECURSIVE subroutine ReadCGFile(MyOffset,MyDomain)
		USE	:: SpecificVarEqn99
		USE :: Parameters, only : ICFlag, ndim
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
	if(ICFlag .eq. 'CGeom' .and. ndim .eq. 3) then
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

!RECURSIVE SUBROUTINE SmoothInterface(alpha,r,ICsig,epsilon,smooth_order)
!    USE, INTRINSIC :: ISO_C_BINDING
!    USE Parameters, ONLY : nVar, nDim
!    IMPLICIT NONE 
!    ! Argument list 
!    REAL, INTENT(IN)               :: r,ICsig
!    REAL, INTENT(OUT)              :: alpha        ! 
!
!    REAL    :: eta,xi
!        real :: epsilon
!        integer :: smooth_order
!        
!        if(epsilon<0) then
!            epsilon=1.e-9    
!        end if
 !       if(smooth_order<0) then
 !           smooth_order=4   
 !       end if
 !       eta=0.8 ! Optimal for v1
 !       ! =============== WAY 1 ================================
 !       if(r>(1+eta)*ICsig) then
 !           xi=1    
 !       elseif(r<-(1-eta)*ICsig) then
 !           xi=0
!        else
 !           xi = (r+(1-eta)*ICsig)/(2.0*ICsig) 
 !       end if
 !       ! =============== WAY 2 ================================
 !       alpha  = (1.0-epsilon)*(1-xi)**smooth_order + (0.0+epsilon)*xi**smooth_order 
 !       if(smooth_order .eq. 0) then
 !           xi = 0.5+0.5*ERF(r/(2*ICsig))  
 !           alpha  = (1.0-epsilon)*(1-xi) + (0.0+epsilon)*xi
 !       end if
!END SUBROUTINE SmoothInterface

