! DIM Initial Data


RECURSIVE SUBROUTINE InitialData(x, t, Q)
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar, nDim
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: x(nDim), t        ! 
	REAL, INTENT(OUT)              :: Q(nVar)        ! 

	!Call InitialPlaneWave(x, t, Q)
	!call GaussianBubble(x,t,Q)
		
	call InitialCG3D(x, t, Q);
END SUBROUTINE InitialData

RECURSIVE SUBROUTINE PDElimitervalue(limiter_value,xx)
	USE, INTRINSIC :: ISO_C_BINDING
	USE SpecificVarEqn99
	USE Parameters, ONLY : nVar, nDim
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: xx(nDim)        ! 
	INTEGER, INTENT(OUT)              :: limiter_value        !
	real	:: rr	

	! Plane Wave limiter_value
	!rr =sqrt(sum(xx(:)**2))
	!rr=0.25-rr
	!if(abs(rr).le. 0.25) then
	!	limiter_value=1
	!else
	!	limiter_value=0
	!end if
	!
	
	! Exclude the boundaries
	if(abs(xx(1))>2000 .or. abs(xx(2))>2000) then
		limiter_value=0
		return
	end if
	!rr = RadiusFromCG(xx(1),xx(2),xx(3))
	rr = DistanceFromSurfaceCG(xx(1),xx(2),xx(3))
	!rr=xx(3)-2000
	if(abs(rr)<500) then
		limiter_value=1
	else
		limiter_value=0
	end if	
	limiter_value=0
END SUBROUTINE PDElimitervalue

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
        if(x(3) .gt. -1000.0) then
            ! First medium
            up(10) = 2.080000e+10  
            up(11) = 1.040000e+10  
            up(12) = 2600.0 
        else
            ! Second medium
            up(10) = 3.240380160000000e+10  
            up(11) = 3.239809920000000e+10
            up(12) = 2700.0
        end if
        ICsig=200.0
        r = (x(3)+1000)
        !ICsig=10
        up(10) =0.5*(3.24038016e+10 +2.080000e+10 )+0.5*(3.24038016e+10 -2.080000e+10 )*erf(-r/ICsig)
        up(11) =0.5*(3.239809920000000e+10 +1.040000e+10)+0.5*(3.239809920000000e+10 -1.040000e+10)*erf(-r/ICsig)
        up(12) =0.5*(2700.0 +2600.0 )+0.5*(2700.0 -2600.0 )*erf(-r/ICsig)
        up(14)  = 1.0-1e-7   ! No Fracture everywhere
        r = x(3)
        !r = RadiusFromCG(x(1),x(2),x(3))
		r = DistanceFromSurfaceCG(x(1),x(2),x(3))
        ICsig=400.0

        up(13)=SmoothInterface(r,ICsig,1.e-4,1)                 ! Get the smooth value of alpha for that value of r
		!up(13)=1.0
        ICx0=[0.0, 0.0, 0.0]
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
    ICsig=3e-3  ! smoothing parameter
    ICQR(:)= (/1e-7, 0.9999999 /)
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
	CALL SmoothInterface(up(13),r,ICsig,1.e-4,0)
    !up(13)=1.0
	
    !up(1:9) = up(1:9)*up(13)
    up(14)=1.0-1.e-10

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
    subroutine ReadCGFile()
		USE	:: SpecificVarEqn99
		USE, INTRINSIC :: ISO_C_BINDING
        implicit none
        CHARACTER(LEN=100)         :: CGEOMFile
        integer     :: jcg,nx_cg_new,ny_cg_new,i,j, ix(2)
        real        :: leng(2),center(2)
        integer     :: n_new_in(2)
        real, allocatable   :: x_cg_new(:),y_cg_new(:),z_cg_new(:,:)
        real            :: h, phi(4), xi,gamma
        real            :: minx,maxx,miny,maxy
		! Input parameters
		leng=(/15000.0, 15000.0/)	! Dimension of the area (in the x- and y- direction)
		center=(/0.0, 0.0/)			! UTM coordinates of the center (with respect to the DTM data file)
		n_new_in=(/50, 50/)			! Numer of elements for the min sub tri function
		CGEOMFile="CG.dat"			! DTM file
		
		
        open(8, file=trim(CGEOMFile), action='read')
            read(8,*) nx_cg
            read(8,*) ny_cg
            allocate(x_cg(nx_cg),y_cg(ny_cg),z_cg(nx_cg,ny_cg))
            read(8,*) x_cg(1:nx_cg)
            read(8,*) y_cg(1:ny_cg)
            do jcg=1,ny_cg
                read(8,*) z_cg(1:nx_cg,jcg)       
            end do
        close(8)
        
        y_cg(1:nx_cg)=y_cg(nx_cg:1:-1)
        z_cg(nx_cg,1:ny_cg)=z_cg(nx_cg,ny_cg:1:-1)
        
        ! Connect the input parameters 
        maxx=center(1)+leng(1)
        minx=center(1)-leng(1)
        maxy=center(2)+leng(2)
        miny=center(2)-leng(2)
        
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
                !z_cg_new(i,j)=0.2*y_cg_new(j)
            end do
        end do
        x_cg_new=x_cg_new-center(1) ! The center for the DTM is (0,0) in the fortran code
        y_cg_new=y_cg_new-center(2) ! The center for the DTM is (0,0) in the fortran code
        
        deallocate(x_cg,y_cg,z_cg)
        nx_cg=nx_cg_new;
        ny_cg=ny_cg_new;        
        allocate(x_cg(nx_cg),y_cg(ny_cg),z_cg(nx_cg,ny_cg))
        x_cg=x_cg_new
        y_cg=y_cg_new
        z_cg=z_cg_new
		deallocate(x_cg_new,y_cg_new,z_cg_new)
		print *, 'Min-Max of x=',minval(x_cg), maxval(x_cg)
		print *, 'Min-Max of y=',minval(y_cg), maxval(y_cg)
		print *, 'Min-Max of z=',minval(z_cg), maxval(z_cg)
        
    end subroutine ReadCGFile

RECURSIVE SUBROUTINE SmoothInterface(alpha,r,ICsig,epsilon,smooth_order)
    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: r,ICsig
    REAL, INTENT(OUT)              :: alpha        ! 

    REAL    :: eta,xi
        real :: epsilon
        integer :: smooth_order
        
        if(epsilon<0) then
            epsilon=1.e-9    
        end if
        if(smooth_order<0) then
            smooth_order=4   
        end if
        eta=0.8 ! Optimal for v1
        ! =============== WAY 1 ================================
        if(r>(1+eta)*ICsig) then
            xi=1    
        elseif(r<-(1-eta)*ICsig) then
            xi=0
        else
            xi = (r+(1-eta)*ICsig)/(2.0*ICsig) 
        end if
        ! =============== WAY 2 ================================
        alpha  = (1.0-epsilon)*(1-xi)**smooth_order + (0.0+epsilon)*xi**smooth_order 
        if(smooth_order .eq. 0) then
            xi = 0.5+0.5*ERF(r/(2*ICsig))  
            alpha  = (1.0-epsilon)*(1-xi) + (0.0+epsilon)*xi
        end if
END SUBROUTINE SmoothInterface

