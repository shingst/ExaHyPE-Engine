module SpecificVarEqn99
    IMPLICIT NONE  
    PUBLIC  
    integer                             :: nLayers
    real, allocatable, dimension(:,:)    :: LayerProperties
    ! Data for the complex 3D profile
    real, allocatable                   :: x_cg(:), y_cg(:),z_cg(:,:)
	real    :: xL_cg, xR_cg,yL_cg,yR_cg,dx_cg,dy_cg,zMax_cg,zMin_cg
    integer                             :: nx_cg,ny_cg
    contains
    !*********************************************************************************
    ! Three dimensional complex geometry
    !*********************************************************************************
     real function DistanceFromSurfaceCG(x_in,y_in,z_in,di_size) ! Look at the minimum on the loaded configuration
        implicit none
        real    :: x_in,y_in,z_in
        ! Local variables
        real    :: minvals(3), rr, np(3), np_prj(3), point(3), xx(3,3), nv(3),u(3),v(3),cv(3),tv(3), sign, np_norm,np_nv_norm,di_size
        integer :: minpos(3,2),i,j,k,l,shift(8,2,2)
        integer :: ix_mid,iy_mid,minpos_deb(2)
        point=(/x_in,y_in,z_in /)
        

		
        if(z_in<zMin_cg-10.0*di_size) then
            DistanceFromSurfaceCG=-1.e+14 
            return
        end if
        if(z_in>zMax_cg+10.0*di_size) then
            DistanceFromSurfaceCG=1.e+14 
            return
        end if
		
		ix_mid=floor((x_in-xL_cg)/dx_cg)
        iy_mid=floor((y_in-yL_cg)/dy_cg)
		! Second reduction
        !if(z_in<minval(z_cg(max(ix_mid-40,1):min(nx_cg,ix_mid+40),max(iy_mid-40,1):min(ny_cg,iy_mid+40)))-10.0*di_size) then
        !    DistanceFromSurfaceCG=-1.e+14 
        !    return
        !end if
        !if(z_in>maxval(z_cg(max(ix_mid-40,1):min(nx_cg,ix_mid+40),max(iy_mid-40,1):min(ny_cg,iy_mid+40)))+10.0*di_size) then
        !    DistanceFromSurfaceCG=1.e+14 
        !    return
        !end if
        ! Fast way
        
        minvals=1.e+14
        minpos=0
        do i=max(ix_mid-40,1),min(nx_cg,ix_mid+40)
            do j=max(iy_mid-40,1),min(ny_cg,iy_mid+40)
               rr=sqrt((point(1)-x_cg(i))**2+(point(2)-y_cg(j))**2+(point(3)-z_cg(i,j))**2)
               if(rr<minvals(1)) then
                   minvals(1)=rr
                   minpos(1,1)= i!(/i,j/);
				   minpos(1,2)= j
                end if
            end do
        end do  
        !if(maxval(abs(minpos_deb-minpos(1,:)))>0) then
        !    print *, 'Warning Distance search!'
        !    print *, minpos(1,:),minpos_deb
        !    print *, 'ix=',ix_mid,'i_est=',(x_in-xL_cg)/dx_cg
        !    print *, x_in,xL_cg,dx_cg
        !    print *, 'ix=',iy_mid,'i_est=',(y_in-yL_cg)/dy_cg
        !    print *, y_in,yL_cg,dy_cg
        !    pause
        !end if
        DistanceFromSurfaceCG=1.e+14
        shift(1,1,:)=(/ 0,1  /)
        shift(1,2,:)=(/ -1,1  /)
        shift(2,1,:)=(/ 1,0  /)
        shift(2,2,:)=(/ 0,1  /) 
        shift(3,1,:)=(/ 1,-1  /)
        shift(3,2,:)=(/ 1,0  /)
        shift(4,1,:)=(/ 0,-1  /)
        shift(4,2,:)=(/ 1,-1  /)
        shift(5,1,:)=(/ -1,0  /)
        shift(5,2,:)=(/ 0,-1  /)
        shift(6,1,:)=(/ -1,1  /)
        shift(6,2,:)=(/ -1,0  /)
		np_nv_norm=100;
        do l=1,6
            minpos(2,:)=minpos(1,:)+shift(l,1,:)
            minpos(3,:)=minpos(1,:)+shift(l,2,:)
            if(minval(minpos) .eq. 0 .or. maxval(minpos(:,1))>nx_cg  .or. maxval(minpos(:,2))>ny_cg  ) then
                cycle
            end if    

            ! COnstruct the triangle composed by the three minimum
            do k=1,3
                xx(k,:)=[x_cg(minpos(k,1)),y_cg(minpos(k,2)),z_cg(minpos(k,1),minpos(k,2))];
            end do
            ! COmpute the normal and the tangential vectors
            u=xx(3,:)-xx(1,:)
            v=xx(2,:)-xx(1,:)
            nv(1)=u(2)*v(3)-u(3)*v(2)
            nv(2)=u(3)*v(1)-u(1)*v(3)
            nv(3)=u(1)*v(2)-u(2)*v(1)
            nv=-nv/sqrt(sum(nv**2))  ! Normalize nv
        
            np_norm=dot_product((point-xx(1,:)),nv)
            np_prj=point-np_norm*nv
            
            if(np_norm>0) then
                sign=+1.0    
            else
                sign=-1.0     
            end if
            np=np_prj
            call GetTirangleMinPoint(xx(:,1),xx(:,2),xx(:,3),np)
            if(abs(np_norm)<1.e-1 .and. maxval(np-np_prj)>1.e-6) then ! Normal vector perturbation problem
                cycle    
            end if
            
            rr=sqrt(sum(  (np-point)**2  ))
            if(rr.le. abs(DistanceFromSurfaceCG)) then
                if(abs(rr-abs(DistanceFromSurfaceCG))<1.e-6 .and. abs(np_norm)<abs(np_nv_norm)) then
                    cycle    
                else
                    DistanceFromSurfaceCG=sign*rr
                    np_nv_norm=np_norm
                end if
            end if
        end do
    end function DistanceFromSurfaceCG
    
    
    subroutine GetTirangleMinPoint(XX,YY,ZZ,PP)
        implicit none
        real    :: XX(3),YY(3),ZZ(3),PP(3)
        real    :: xi, gamma, xi_e, gamma_e, x_out(3), alpha
        
    gamma=((PP(1)-XX(1))*(YY(2)-YY(1))-(XX(2)-XX(1))*(PP(2)-YY(1)))/((XX(3)-XX(1))*(YY(2)-YY(1))+(YY(1)-YY(3))*(XX(2)-XX(1)))
    xi=((PP(1)-XX(1))*(YY(3)-YY(1))-(XX(3)-XX(1))*(PP(2)-YY(1)))/((XX(2)-XX(1))*(YY(3)-YY(1))+(YY(1)-YY(2))*(XX(3)-XX(1)))
    
    if(xi>0 .and. gamma > 0 .and. xi+gamma<1) then
        xi_e=xi
        gamma_e=gamma
    elseif(xi<=0) then
      if(gamma<=0) then
          xi_e=0
          gamma_e=0
      elseif(gamma <=1) then
          xi_e=0
          gamma_e=gamma      
      else
          xi_e=0
          gamma_e=1          
      end if
    elseif(gamma<=0) then
      if(xi <=1) then
          xi_e=xi
          gamma_e=0         
      else
          xi_e=1
          gamma_e=0;        
      end if       
    else
        alpha=0.5*(gamma-xi+1)
        if(alpha<=0) then
          xi_e=1
          gamma_e=0              
        elseif(alpha>=1) then
          xi_e=0
          gamma_e=1             
        else
          xi_e=0.5*(1-alpha)
          gamma_e=0.5*(1+alpha)              
        end if
    end if
    
    x_out(1)=XX(1)+(XX(2)-XX(1))*xi_e+(XX(3)-XX(1))*gamma_e
    x_out(2)=YY(1)+(YY(2)-YY(1))*xi_e+(YY(3)-YY(1))*gamma_e
    x_out(3)=ZZ(1)+(ZZ(2)-ZZ(1))*xi_e+(ZZ(3)-ZZ(1))*gamma_e
    
    PP=x_out
    end subroutine GetTirangleMinPoint
    
	
	
    real function RadiusFromCG(x_in,y_in,z_in)
		USE, INTRINSIC :: ISO_C_BINDING
        implicit none
        real    :: x_in,y_in,z_in,z_out
        real    :: i,j,phi(4),xi,gamma
        integer :: ix(2)
        
        ix=lookatindex_cg(x_in,y_in)
        phi(1)=z_cg(ix(1),ix(2))
        phi(2)=z_cg(ix(1)+1,ix(2))
        phi(3)=z_cg(ix(1),ix(2)+1)
        phi(4)=z_cg(ix(1)+1,ix(2)+1)
        xi=(x_in-x_cg(ix(1)))/(x_cg(ix(1)+1)-x_cg(ix(1)))
        gamma=(y_in-y_cg(ix(2)))/(y_cg(ix(2)+1)-y_cg(ix(2)))
        z_out=(1-xi)*(1-gamma)*phi(1)+xi*(1-gamma)*phi(2)+gamma*(1-xi)*phi(3)+xi*gamma*phi(4)
        !z_out=z_cg(ix(1),ix(2))+z_cg(ix(1)+1,ix(2))+z_cg(ix(1),ix(2)+1)+z_cg(ix(1)+1,ix(2)+1)
        ! Reverse for negative representation (z down)
        !z_out=-z_out
        !RadiusFromCG=z_out-z_in
        RadiusFromCG=-z_out+z_in
    end function RadiusFromCG

    
    function lookatindex_cg(x_in,y_in)
        implicit none
        real    :: x_in,y_in
        integer    :: lookatindex_cg(2)
        integer :: i,j
        
        lookatindex_cg=-1
        do i=1,nx_cg-1
            if(x_cg(i).le.x_in .and. x_in .le. x_cg(i+1)) then
               lookatindex_cg(1)=i
               exit
            end if
        end do
        do j=1,ny_cg-1
            if(y_cg(j).le.y_in .and. y_in .le. y_cg(j+1)) then
               lookatindex_cg(2)=j
               exit
            end if
        end do   
        if(minval(lookatindex_cg(:))<0) then
            print *, 'LookIndex error for x_in=',x_in, ' and y_in=',y_in, '! Please choose a larger CG domain'
            stop
        end if
    end function lookatindex_cg
    
    function lookatindex_cg_fast(x_in,y_in)
        implicit none
        real    :: x_in,y_in
        integer    :: lookatindex_cg_fast(2)
        integer :: i,j
        
		
		
		lookatindex_cg_fast(1)=floor((x_in-xL_cg)/dx_cg)+1
        lookatindex_cg_fast(2)=floor((y_in-yL_cg)/dy_cg)+1
		!print *, 'DyL_cg=',y_in-yL_cg,x_in-xL_cg,lookatindex_cg_fast
        if(minval(lookatindex_cg_fast(1:2))<0) then
            print *, 'LookIndex error for x_in=',x_in, ' and y_in=',y_in, '! Please choose a larger CG domain'
            stop
        end if
    end function lookatindex_cg_fast
    !*********************************************************************************
    !           Use constant rings in the parameter file
    !*********************************************************************************
    subroutine ComputeMaterialProperty(Prop,r,ICsig)
		USE, INTRINSIC :: ISO_C_BINDING
        implicit none
        real    :: Prop(1:4),r, ICsig
        intent(OUT) :: Prop
        ! Local varialbes
        real    :: xi, rcrit,MPL(1:4), MPR(1:4) ! Material property left and right
        real    :: Drin, Drout,r_layer
        integer :: i,iL,iR
        
        i=getmaterialring(r)
        if(r<4) then
            continue    
        end if
        if(abs(r-LayerProperties(1,i))<abs(r-LayerProperties(2,i))) then
           iR=i
           iL=getmaterialring(max(0.0,LayerProperties(1,i)-1.e-12))
           r_layer=LayerProperties(1,i)
        else
           iR=getmaterialring(LayerProperties(2,i))
           iL=i
           r_layer=LayerProperties(2,i)
        end if
        
        xi = 0.5+0.5*ERF((r-r_layer)/ICsig) 
        MPL=LayerProperties(3:6,iL)
        MPR=LayerProperties(3:6,iR)
        Prop(1:4)  = MPL(1:4)*(1-xi) + MPR(1:4)*xi  ! Smooth transformation of material parameter
        if(r-r_layer.gt.0 ) then
            Prop(1:3)  = MPR(1:3)                   ! Discontinuous transformation of material parameter  
        else
            Prop(1:3)  = MPL(1:3)                   ! Discontinuous transformation of material parameter  
        end if
    end subroutine ComputeMaterialProperty
    ! ********************************************************************************
    function getmaterialring(r)
		USE, INTRINSIC :: ISO_C_BINDING
        implicit none
        real        :: r
        integer     :: i, getmaterialring
        
        do i=1,nLayers
            if(r-LayerProperties(1,i).ge.0 .and. r-LayerProperties(2,i).lt.0) then
                getmaterialring=i
                return
            end if
        end do  
        getmaterialring=0
    end function getmaterialring    
    !*********************************************************************************
    subroutine PREM(Prop,rin,ICsig)
		USE, INTRINSIC :: ISO_C_BINDING
        implicit none
        real    :: Prop(1:4),r,rin, ICsig
        intent(OUT) :: Prop
        ! Local varialbes
        real    :: rnorm,Rmax,Rring(2)
        
        real    :: xi, rcrit,MPL(1:4), MPR(1:4),QV(4) ! Material property left and right
        real    :: Drin, Drout,r_layer,factor;
        integer :: i,iL,iR
        factor=1000.0;
        r=rin/factor; ! m-> Km since the datas are in kilometers
        Rmax=6371       ! Earth radius
        rnorm=r/RMAX    ! Normalized radius
        QV(4)=1.0       ! Material inside
        if(r .ge. 0 .and. r .lt. 1221.5) then           !++++++ Inner core +++++++
            Rring(:)=(/ 0.0,    1221.5 /)
            QV(1)= 13.0885-8.8381*rnorm**2  ! Density
            QV(2)= 11.2622-6.3640*rnorm**2  ! P-wave speed
            QV(3)= 3.6678 -4.4475*rnorm**2  ! S-wave speed
        elseif(r .ge. 1221.5 .and. r .lt. 3480.0) then  !++++++ Outer core +++++++
            Rring(:)=(/ 1221.5,    3480.0 /)
            QV(1)= 12.5815-1.2638*rnorm-3.6426*rnorm**2-5.5281*rnorm**3  ! Density
            QV(2)= 11.0487-4.0362*rnorm+4.8023*rnorm**2-13.5732*rnorm**3 ! P-wave speed
            QV(3)= 0.                                                    ! S-wave speed            
        elseif(r .ge. 3480.0 .and. r .lt. 3630.0) then  !++++++ Lower mantle (1) +++++++  
            Rring(:)=(/ 3480.0,    3630.0 /)
            QV(1)= 7.9565-6.4761*rnorm+5.5283*rnorm**2-3.0807*rnorm**3  ! Density
            QV(2)= 15.3891-5.3181*rnorm+5.5242*rnorm**2-2.5514*rnorm**3 ! P-wave speed
            QV(3)= 6.9254+1.4672*rnorm-2.0834*rnorm**2+0.9783*rnorm**3  ! S-wave speed      
        elseif(r .ge. 3630.0 .and. r .lt. 5600.0) then  !++++++ Lower mantle (2) +++++++ 
            Rring(:)=(/ 3630.0,    5600.0 /)
            QV(1)= 7.9565-6.4761*rnorm+5.5283*rnorm**2-3.0807*rnorm**3  ! Density
            QV(2)= 24.9520-40.4673*rnorm+51.4832*rnorm**2-26.6419*rnorm**3 ! P-wave speed
            QV(3)= 11.1671-13.7818*rnorm+17.4575*rnorm**2-9.2777*rnorm**3  ! S-wave speed   
        elseif(r .ge. 5600.0 .and. r .lt. 5701.0) then  !++++++ Lower mantle (3) +++++++  
            Rring(:)=(/ 5600.0 ,    5701.0/)
            QV(1)= 7.9565-6.4761*rnorm+5.5283*rnorm**2-3.0807*rnorm**3  ! Density
            QV(2)= 29.2766-23.6027*rnorm+5.5242*rnorm**2-2.5514*rnorm**3 ! P-wave speed
            QV(3)= 22.3459-17.2473*rnorm-2.0834*rnorm**2+0.9783*rnorm**3  ! S-wave speed   
        elseif(r .ge. 5701.0 .and. r .lt. 5771.0) then  !++++++ Transition zone (1) +++++++ 
            Rring(:)=(/ 5701.0 ,    5771.0 /)
            QV(1)= 5.3197-1.4836*rnorm  ! Density
            QV(2)= 19.0957-9.8672*rnorm ! P-wave speed
            QV(3)= 9.9839-4.9324*rnorm ! S-wave speed   
        elseif(r .ge. 5771.0 .and. r .lt. 5971.0) then  !++++++ Transition zone (2) +++++++ 
            Rring(:)=(/ 5771.0,    5971.0 /)
            QV(1)= 11.2494-8.0298*rnorm  ! Density
            QV(2)= 39.7027-32.6166*rnorm ! P-wave speed
            QV(3)= 22.3512-18.5856*rnorm ! S-wave speed  
        elseif(r .ge. 5971.0 .and. r .lt. 6151.0) then  !++++++ Transition zone (3) +++++++  
            Rring(:)=(/ 5971.0,    6151.0 /)
            QV(1)= 7.1089-3.8045*rnorm  ! Density
            QV(2)= 20.3926-12.2569*rnorm ! P-wave speed
            QV(3)= 8.9496-4.4597*rnorm ! S-wave speed  
        elseif(r .ge. 6151.0 .and. r .lt. 6346.6) then  !++++++ LVZ and LID isotropic approximation +++++++  
            Rring(:)=(/ 6151.0,    6346.6 /)
            QV(1)= 2.6910+0.6924*rnorm  ! Density
            QV(2)= 4.1875+3.9382*rnorm ! P-wave speed
            QV(3)= 2.1519+2.3481*rnorm ! S-wave speed  
        elseif(r .ge. 6346.6 .and. r .lt.6356.0) then  !++++++ Crust (1) +++++++ 
            Rring(:)=(/ 6346.6,    6356.0 /)
            QV(1)= 2.900  ! Density
            QV(2)= 6.800 ! P-wave speed
            QV(3)= 3.900 ! S-wave speed  
        elseif(r .ge. 6356.0 .and. r .lt.6368.0) then  !++++++ Crust (2) +++++++ 
            Rring(:)=(/ 6356.0,    6368.0 /)
            QV(1)= 2.600  ! Density
            QV(2)= 5.800 ! P-wave speed
            QV(3)= 3.200 ! S-wave speed 
        elseif(r .ge. 6368.0 .and. r .lt.6371.0) then  !++++++ Ocean +++++++  
            Rring(:)=(/ 6368.0,    6371.0 /)
            QV(1)= 1.020  ! Density
            QV(2)= 1.450 ! P-wave speed
            QV(3)= 0.0 ! S-wave speed 
        else
            Rring(:)=(/ 6371.0,    1.e16 /)
            QV(1)= 1.020  ! Last one
            QV(2)= 1.450 ! Last one
            QV(3)= 0.0 ! Last one
            QV(4)=0.0
        end if
        ! Sismic Unit conversion
        !QV(1)=QV(1)*1.e+12
        ! SI conversion
        QV(1)=QV(1)/1000.0 ! g/cm3-> Kg/m3
        QV(2)=QV(2)*1000.0 ! km/s-> m->s
        QV(3)=QV(3)*1000.0 ! km/s-> m->s
        
        if(abs(r-Rmax)<1.0) then
            xi = 0.5+0.5*ERF((r-Rmax)/ICsig) 
            QV(4)=(1.0-1.e-7)*(1-xi) + (0.0+1.e-7)*xi           ! Smooth transation in the external ring for alpha
        end if
        Prop(1)=  QV(1)*(QV(2)**2-2.0*QV(3)**2)         ! mu
        Prop(2)=  QV(1)*QV(3)**2         ! mu
        Prop(3)=  QV(1)
        Prop(4)=  QV(4)
    end subroutine PREM
    !*********************************************************************************
    !*********************************************************************************
    subroutine GetDistSurfaceInterface(dist,nv,x0,ind)
		USE, INTRINSIC :: ISO_C_BINDING
        implicit none
        real        :: f,x0(2),nv(3),dist
        integer     :: ind,i
        real        :: dd,ddd,xx(3),dx(3),ddx(3),t0(2)
        real        :: xx0(3)
        select case(ind)
        case(0) ! 2D complex geometry
                t0=xy2t(x0,ind)
                call FreeSurfaceInterface(xx0,dx,ddx,t0,ind)
                do i=1,100   ! Newton method to look at the local minimum
                    call FreeSurfaceInterface(xx,dx,ddx,t0,ind)  
                    dd=(xx(1)-x0(1))*dx(1)+(xx(2)-x0(2))*dx(2)
                    ddd=dx(1)**2+(xx(1)-x0(1))*ddx(1)+dx(2)**2+(xx(2)-x0(2))*ddx(2)
                    t0(1)=t0(1)-dd/ddd
                    if(abs(dd).lt.1.e-5 .or. t0(1).gt.1 .or. t0(1).lt.0) then
                        exit    
                    end if
                end do
                if(i.ge. 100) then
                    t0=xy2t(x0,ind)  
                end if
                call FreeSurfaceInterface(xx,dx,ddx,t0,ind) 
                nv=(/  -dx(2),dx(1), 0.0 /)
                nv=nv/sqrt(nv(1)**2+nv(2)**2)
                dist  = (x0(1)-xx(1))*nv(1) + (x0(2)-xx(2))*nv(2)
                if(x0(2)>xx0(2)) then
                    dist=sqrt(sum((x0(1:2)-xx(1:2))**2))
                else
                    dist=-sqrt(sum((x0(1:2)-xx(1:2))**2))   
                end if
            case default
                print *, 'Free Surface function Not implemented!'
                stop
        end select
            
    end subroutine GetDistSurfaceInterface
    !*********************************************************************************
    function SmoothInterface(r,ICsig,epsilon,smooth_order_in)
		USE, INTRINSIC :: ISO_C_BINDING
        implicit none
        real    :: r
        real    :: SmoothInterface
        real    :: eta,ICsig,xi 
        real, optional :: epsilon
        integer, optional :: smooth_order_in
		real	:: smooth_order
        
        if(.not. present(epsilon)) then
            epsilon=1.e-9    
        end if

		smooth_order=0.5
        eta=-0.6
        ! =============== WAY 1 ================================
        if(r>(1+eta)*ICsig) then
            xi=1    
        elseif(r<-(1-eta)*ICsig) then
            xi=0
        else
            xi = (r+(1-eta)*ICsig)/(2.0*ICsig) 
        end if
        ! =============== WAY 2 ================================
        SmoothInterface  = (1.0-epsilon)*(1-xi)**smooth_order + (0.0+epsilon)*xi**smooth_order 
        !if(smooth_order .eq. 0) then
        !    xi = 0.5+0.5*ERF(r/(2*ICsig))  
        !    SmoothInterface  = (1.0-epsilon)*(1-xi) + (0.0+epsilon)*xi
        !end if
    end function SmoothInterface
    
    function xy2t(x,ind)
		USE, INTRINSIC :: ISO_C_BINDING
        implicit none
        real :: x(2),xy2t(2)
        integer :: ind
        select case(ind)
            case(0) ! 2D complex geometry
              xy2t=(/ x(1)/2000.0,0.0 /)  
            case default
                print *, 'Free Surface function Not implemented!'
                stop
        end select            
    end function xy2t
    !*********************************************************************************
    subroutine FreeSurfaceInterface(xx,dx,ddx,t,ind)
		USE, INTRINSIC :: ISO_C_BINDING
        implicit none
        real    :: xx(3),dx(3),ddx(3)   ! Value of the function, I and II derivative (only 2d for the moment)
        real    :: t(2),x_t,y_t
        integer :: ind
        ! t in [0,1]^2 is the parametrization of the surface
        select case(ind)
        case(0) ! 2D complex geometry
                x_t=2000*t(1)
                y_t=2000.0+100.0*sin(3.0*x_t/200.0)+100.0*sin(2.0*x_t/200.0)
                xx=(/  x_t,y_t, 0.0 /)
                dx=2000.0*(/  1.0,3.0/2.0*cos(3.0*x_t/200.0)+1.0*cos(2*x_t/200.0), 0.0 /)
                ddx=2000.0**2*(/  0.0,-9.0/400.0*sin(3.0*x_t/200.0)-1.0/100.0*sin(2*x_t/200.0), 0.0 /)
            case default
                print *, 'Free Surface function Not implemented!'
                stop
        end select
        
    end subroutine FreeSurfaceInterface
    !*********************************************************************************
end module SpecificVarEqn99