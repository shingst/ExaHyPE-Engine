    module SpecificVarEqn99
    integer                             :: nLayers
    real, allocatable, dimension(:,:)    :: LayerProperties
    ! Datas for the complex 3D profile
    real, allocatable                   :: x_cg(:), y_cg(:),z_cg(:,:)
    integer                             :: nx_cg,ny_cg
    real    :: xL_cg, xR_cg,yL_cg,yR_cg,dx_cg,dy_cg,zMax_cg,zMin_cg
    ! Triangulation 
    integer, allocatable                :: CGfaces(:,:)
    real, allocatable                   :: CGvertices(:,:)
    integer                             :: NCGTri, NCGvertices
    
    integer, parameter                  :: GPRLEversion=2       ! 0 -> sigma=-rho0*|A|*(1-eps)*E_eps    1 -> sigma=-rho0*|A|*G*E_A      2-> sigma=-rho A^top E_A
    
    TYPE tFriction
        REAL :: mu_s    ! Static friction coefficient
        REAL :: mu_d    ! Dynamic friction coefficient
        REAL :: Dc      ! Critical displacement
        ! ------------------------------------------------------------------
        LOGICAL :: DynamicFL ! Set if it is a dynamic or static friction law
        REAL    :: V         ! Rupture velocity
        REAL    :: t0        ! Rupture initial time
        REAL    :: L        ! Nucleation size
        ! ------------------------------------------------------------------
    END TYPE tFriction
    TYPE(tFriction) :: SSCRACKFL
    LOGICAL         :: USEFrictionLaw
    contains
    !*********************************************************************************
    ! Three dimensional complex geometry
    !*********************************************************************************
    !subroutine ReadCGFile(minx,maxx,miny,maxy)
    
    subroutine GetCGTriangulation()
        implicit none
        real, allocatable   :: Vertices(:,:)
        integer                 :: nVert,i,j, nFaces, ind1,ind2,ind3,ind4,ind5,ind6, stopind
        
        allocate(CGfaces((nx_cg-1)*(ny_cg-1)*2, 3   ))
        allocate(Vertices((nx_cg)*(ny_cg), 3   ))
        nFaces=0
        nVert=0
        stopind=nx_cg*4;
        do i=1,nx_cg-1
            do j=1,ny_cg-1
                nFaces=nFaces+1
                call GetVertexIndex(ind1,Vertices,nVert,(/ x_cg(i),y_cg(j),z_cg(i,j) /),stopind,0    )
                call GetVertexIndex(ind2,Vertices,nVert,(/ x_cg(i+1),y_cg(j),z_cg(i+1,j) /),stopind,0    )
                call GetVertexIndex(ind3,Vertices,nVert,(/ x_cg(i),y_cg(j+1),z_cg(i,j+1) /),stopind,0    )
                CGfaces(nfaces,:)=(/ ind1,ind2,ind3 /)
                nFaces=nFaces+1
                call GetVertexIndex(ind4,Vertices,nVert,(/ x_cg(i+1),y_cg(j+1),z_cg(i+1,j+1) /),stopind,0    )
                ind5=ind3
                ind6=ind2
                CGfaces(nfaces,:)=(/ ind4,ind5,ind6 /)
            end do
        end do
        allocate(CGvertices(nVert,3))
        CGvertices=Vertices(1:nVert,:)
        deallocate(Vertices)
    end subroutine GetCGTriangulation
    
    subroutine GetVertexIndex(ind,vertices,nvertices,xv,stopind,skip   )
        implicit none
        integer :: skip,stopind,nvertices,ind
        real :: vertices((nx_cg)*(ny_cg), 3)
        integer :: k
        real    :: xv(3)
        
        if(skip .eq. 0) then
            do k=nvertices,max(1,nvertices-stopind),-1   
                if(sum((xv-vertices(k,:))**2 )<1.e-10) then
                    ind=k
                    return
                end if
            end do
        end if
        nvertices=nvertices+1
        ind=nvertices
        vertices(ind,:)=xv
    end subroutine GetVertexIndex
    
    !real function DistanceFromSurfaceCG(x_in,y_in,z_in) ! Look at the minimum on the loaded configuration
    !    implicit none
    !    real    :: x_in,y_in,z_in
    !    ! Local variables
    !    real    :: minvals(3), rr, np(3), point(3), xx(3,3), nv(3),u(3),v(3),cv(3),tv(3)
    !    integer :: minpos(3,2),i,j,k
    !    point=(/x_in,y_in,z_in /)
    !    minvals=1.e+14
    !    minpos=0
    !    do i=1,nx_cg
    !        do j=1,ny_cg
    !           rr=sqrt((point(1)-x_cg(i))**2+(point(2)-y_cg(j))**2+(point(3)-z_cg(i,j))**2)
    !           if(rr<minvals(1)) then
    !               minvals(3)=minvals(2)
    !               minpos(3,:)=minpos(2,:) 
    !               minvals(2)=minvals(1)
    !               minpos(2,:)=minpos(1,:)
    !               minvals(1)=rr
    !               minpos(1,:)=(/i,j/);              
    !           elseif(rr<minvals(2)) then
    !               minvals(3)=minvals(2)
    !               minpos(3,:)=minpos(2,:)  
    !               minvals(2)=rr
    !               minpos(2,:)=(/i,j/)               
    !           elseif(rr<minvals(3)) then
    !               minvals(3)=rr
    !               minpos(3,:)=(/i,j/)
    !           end if
    !        end do
    !    end do
    !    ! COnstruct the triangle composed by the three minimum
    !    do k=1,3
    !        xx(k,:)=[x_cg(minpos(k,1)),y_cg(minpos(k,2)),z_cg(minpos(k,1),minpos(k,2))];
    !    end do
    !    ! COmpute the normal and the tangential vectors
    !    u=xx(3,:)-xx(1,:)
    !    v=xx(2,:)-xx(1,:)
    !    nv(1)=u(2)*v(3)-u(3)*v(2)
    !    nv(2)=u(3)*v(1)-u(1)*v(3)
    !    nv(3)=u(1)*v(2)-u(2)*v(1)
    !    nv=nv/sqrt(sum(nv**2))  ! Normalize nv
    !    
    !    np=point-dot_product((point-xx(1,:)),nv)*nv
    !    call GetTirangleMinPoint(xx(:,1),xx(:,2),xx(:,3),np)
    !    DistanceFromSurfaceCG=sqrt(sum(  (np-point)**2  ))
    !end function DistanceFromSurfaceCG
    !
    real function DistanceFromSurfaceCG(x_in,y_in,z_in,di_size) ! Look at the minimum on the loaded configuration
        implicit none
        real    :: x_in,y_in,z_in
        ! Local variables
        real    :: minvals(3), rr, np(3), np_prj(3), point(3), xx(3,3), nv(3),u(3),v(3),cv(3),tv(3), sign, np_norm,np_nv_norm,di_size
        integer :: minpos(3,2),i,j,k,l,shift(8,2,2)
        integer :: ix_mid,iy_mid,minpos_deb(2),inside
        point=(/x_in,y_in,z_in /)
        
        if(z_in<zMin_cg-10.0*di_size) then
            DistanceFromSurfaceCG=-1.e+14 
            return
        end if
        if(z_in>zMax_cg+10.0*di_size) then
            DistanceFromSurfaceCG=1.e+14 
            return
        end if
        ! Fast way
        ix_mid=floor((x_in-xL_cg)/dx_cg)
        iy_mid=floor((y_in-yL_cg)/dy_cg)
        minvals=1.e+14
        minpos=0
        do i=max(ix_mid-40,1),min(nx_cg,ix_mid+40)
            do j=max(iy_mid-40,1),min(ny_cg,iy_mid+40)
               rr=sqrt((point(1)-x_cg(i))**2+(point(2)-y_cg(j))**2+(point(3)-z_cg(i,j))**2)
               if(rr<minvals(1)) then
                   minvals(1)=rr
                   minpos(1,:)=(/i,j/);
                end if
            end do
        end do 
		if(point(3)-z_cg(minpos(1,1),minpos(1,2))>0) then
			DistanceFromSurfaceCG=minvals(1)
		else
			DistanceFromSurfaceCG=-minvals(1)
		end if
		if(minvals(1)>2.0*di_size) then
			!rr=RadiusFromCG(point(1),point(2),point(3))
			return
		end if
        !if(maxval(abs(minpos_deb-minpos(1,:)))>0) then
        !    print *, 'Warning Distance search!'
        !    print *, minpos(1,:),minpos_deb
        !    print *, 'ix=',ix_mid,'i_est=',(x_in-xL_cg)/dx_cg
        !    print *, x_in,xL_cg,dx_cg
        !    print *, 'ix=',iy_mid,'i_est=',(y_in-yL_cg)/dy_cg
        !    print *, y_in,yL_cg,dy_cg
        !    pause
        !end if
		
		!rr=RadiusFromCG(point(1),point(2),point(3))
		!if(rr>0) then
		!	sign=+1.0    
		!else
		!	sign=-1.0     
		!end if
		
        !DistanceFromSurfaceCG=1.e+14
        !shift(1,1,:)=(/ 0,1  /)
        !shift(1,2,:)=(/ -1,0  /)
        !shift(2,1,:)=(/ 1,0  /)
        !shift(2,2,:)=(/ 0,1  /) 
        !shift(3,1,:)=(/ 0,-1  /)
        !shift(3,2,:)=(/ 1,0  /)
        !shift(4,1,:)=(/ -1,0  /)
        !shift(4,2,:)=(/ 0,-1  /)
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
		
		np_nv_norm=100.0
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
            if(abs(np_norm)<5.e-1*di_size .and. maxval(np-np_prj)>1.e-6) then ! Normal vector perturbation problem
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
        implicit none
        real    :: x_in,y_in,z_in,z_out
        real    :: i,j,phi(4),xi,gamma
        integer :: ix(2)
        y_in=-y_in ! Valid only for comparison with PDESOL
        
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
            !pause
            stop
        end if
    end function lookatindex_cg
    !*********************************************************************************
    !           Use constant rings in the parameter file
    !*********************************************************************************
    subroutine ComputeMaterialProperty(Prop,r,ICsig)
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

    subroutine computeGPRLEstress(stressnorm,sigma_vec,Q, addp)
        use Parameters, only : nVar

        implicit none
        real, intent(in) :: Q(nVar)
        real, intent(out) :: sigma_vec(6)
        real :: detA,stressnorm
        logical :: addp
        select case(GPRLEversion)
            case(0)
                sigma_vec(1)=(Q(5) * Q(9) * Q(13) - Q(5) * Q(10) * Q(12) - Q(8) * Q(6) * Q(13) + Q(8) * Q(7) * Q(12) + Q(11) * Q(6) * Q(10) - Q(11) * Q(7) * Q(9)) * (Q(5) ** 2 + Q(8) ** 2 + Q(11) ** 2) * Q(16) * (-0.2D1 / 0.3D1 * Q(5) ** 2 - 0.2D1 / 0.3D1 * Q(8) ** 2 - 0.2D1 / 0.3D1 * Q(11) ** 2 + Q(6) ** 2 / 0.3D1 + Q(9) ** 2 / 0.3D1 + Q(12) ** 2 / 0.3D1 + Q(7) ** 2 / 0.3D1 + Q(10) ** 2 / 0.3D1 + Q(13) ** 2 / 0.3D1)
                sigma_vec(2)=(Q(5) * Q(9) * Q(13) - Q(5) * Q(10) * Q(12) - Q(8) * Q(6) * Q(13) + Q(8) * Q(7) * Q(12) + Q(11) * Q(6) * Q(10) - Q(11) * Q(7) * Q(9)) * (Q(6) ** 2 + Q(9) ** 2 + Q(12) ** 2) * Q(16) * (Q(5) ** 2 / 0.3D1 + Q(8) ** 2 / 0.3D1 + Q(11) ** 2 / 0.3D1 - 0.2D1 / 0.3D1 * Q(6) ** 2 - 0.2D1 / 0.3D1 * Q(9) ** 2 - 0.2D1 / 0.3D1 * Q(12) ** 2 + Q(7) ** 2 / 0.3D1 + Q(10) ** 2 / 0.3D1 + Q(13) ** 2 / 0.3D1)
                sigma_vec(3)=(Q(5) * Q(9) * Q(13) - Q(5) * Q(10) * Q(12) - Q(8) * Q(6) * Q(13) + Q(8) * Q(7) * Q(12) + Q(11) * Q(6) * Q(10) - Q(11) * Q(7) * Q(9)) * (Q(7) ** 2 + Q(10) ** 2 + Q(13) ** 2) * Q(16) * (Q(5) ** 2 / 0.3D1 + Q(8) ** 2 / 0.3D1 + Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.3D1 + Q(9) ** 2 / 0.3D1 + Q(12) ** 2 / 0.3D1 - 0.2D1 / 0.3D1 * Q(7) ** 2 - 0.2D1 / 0.3D1 * Q(10) ** 2 - 0.2D1 / 0.3D1 * Q(13) ** 2)
                sigma_vec(4)=0.2D1 * (Q(5) * Q(9) * Q(13) - Q(5) * Q(10) * Q(12) - Q(8) * Q(6) * Q(13) + Q(8) * Q(7) * Q(12) + Q(11) * Q(6) * Q(10) - Q(11) * Q(7) * Q(9)) * (1 + Q(5) * Q(6) + Q(8) * Q(9) + Q(11) * Q(12)) * (Q(16)) * (-(Q(5) * Q(6)) / 0.2D1 - (Q(8) * Q(9)) / 0.2D1 - (Q(11) * Q(12)) / 0.2D1)
                sigma_vec(5)=0.2D1 * (Q(5) * Q(9) * Q(13) - Q(5) * Q(10) * Q(12) - Q(8) * Q(6) * Q(13) + Q(8) * Q(7) * Q(12) + Q(11) * Q(6) * Q(10) - Q(11) * Q(7) * Q(9)) * (1 + Q(6) * Q(7) + Q(9) * Q(10) + Q(12) * Q(13)) * (Q(16)) * (-(Q(6) * Q(7)) / 0.2D1 - (Q(9) * Q(10)) / 0.2D1 - (Q(12) * Q(13)) / 0.2D1)
                sigma_vec(6)=0.2D1 * (Q(5) * Q(9) * Q(13) - Q(5) * Q(10) * Q(12) - Q(8) * Q(6) * Q(13) + Q(8) * Q(7) * Q(12) + Q(11) * Q(6) * Q(10) - Q(11) * Q(7) * Q(9)) * (1 + Q(5) * Q(7) + Q(8) * Q(10) + Q(11) * Q(13)) * (Q(16)) * (-(Q(5) * Q(7)) / 0.2D1 - (Q(8) * Q(10)) / 0.2D1 - (Q(11) * Q(13)) / 0.2D1)
            case(1)
                sigma_vec(1) = -Q(1) * (Q(5) * Q(9) * Q(13) - Q(5) * Q(10) * Q(12) - Q(6&
                     &) * Q(8) * Q(13) + Q(6) * Q(10) * Q(11) + Q(7) * Q(8) * Q(12) - Q(&
                     &7) * Q(9) * Q(11)) * ((Q(5) ** 2 + Q(8) ** 2 + Q(11) ** 2) * Q(16)&
                     & * (-0.4D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** 2 / 0.3D1 - Q(1&
                     &1) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) **&
                     & 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 /&
                     & 0.6D1) * Q(5) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D1 - Q(9) ** 2 /&
                     & 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D&
                     &1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + &
                     &Q(13) ** 2 / 0.6D1) * Q(5) + 0.2D1 / 0.3D1 * (-Q(7) ** 2 / 0.3D1 -&
                     & Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8&
                     &) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** &
                     &2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(5) - 0.2D1 * (-Q(5) * Q(6) / 0&
                     &.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.2D1) * Q(6) - 0.2D1&
                     & * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D1 - Q(11) * Q(13) / &
                     &0.2D1) * Q(7)) / Q(1) + (Q(5) * Q(6) + Q(8) * Q(9) + Q(11) * Q(12&
                     &)) * Q(16) * (-0.4D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** 2 / 0&
                     &.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 &
                     &+ Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(&
                     &13) ** 2 / 0.6D1) * Q(8) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D1 - Q&
                     &(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) *&
                     &* 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 &
                     &/ 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(8) + 0.2D1 / 0.3D1 * (-Q(7) ** 2&
                     & / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** 2 / 0&
                     &.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0.6D1 &
                     &+ Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(8) - 0.2D1 * (-Q(5) &
                     &* Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.2D1) * Q(&
                     &9) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D1 - Q(11) &
                     &* Q(13) / 0.2D1) * Q(10)) / Q(1) + (Q(5) * Q(7) + Q(8) * Q(10) + &
                     &Q(11) * Q(13)) * Q(16) * (-0.4D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q&
                     &(8) ** 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) *&
                     &* 2 / 0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 &
                     &/ 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(11) + 0.2D1 / 0.3D1 * (-Q(6) ** &
                     &2 / 0.3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0&
                     &.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 &
                     &+ Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(11) + 0.2D1 / 0.3D1&
                     & * (-Q(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 +&
                     & Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6)&
                     & ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(11) - &
                     &0.2D1 * (-Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12&
                     &) / 0.2D1) * Q(12) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) &
                     &/ 0.2D1 - Q(11) * Q(13) / 0.2D1) * Q(13)) / Q(1))
                      sigma_vec(2) = -Q(1) * (Q(5) * Q(9) * Q(13) - Q(5) * Q(10) * Q(12) - Q(6&
                     &) * Q(8) * Q(13) + Q(6) * Q(10) * Q(11) + Q(7) * Q(8) * Q(12) - Q(&
                     &7) * Q(9) * Q(11)) * ((Q(5) * Q(6) + Q(8) * Q(9) + Q(11) * Q(12)) &
                     &* Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** 2 / 0.3D1&
                     & - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(&
                     &12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) &
                     &** 2 / 0.6D1) * Q(6) - 0.4D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D1 - Q(9) &
                     &** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 &
                     &/ 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.&
                     &6D1 + Q(13) ** 2 / 0.6D1) * Q(6) + 0.2D1 / 0.3D1 * (-Q(7) ** 2 / 0&
                     &.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1&
                     & + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0.6D1 + Q(&
                     &9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(6) - 0.2D1 * (-Q(5) * Q(&
                     &6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.2D1) * Q(5) -&
                     & 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(&
                     &13) / 0.2D1) * Q(7)) / Q(1) + (Q(6) ** 2 + Q(9) ** 2 + Q(12) ** 2&
                     &) * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** 2 / 0.3&
                     &D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + &
                     &Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13&
                     &) ** 2 / 0.6D1) * Q(9) - 0.4D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D1 - Q(9&
                     &) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** &
                     &2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / &
                     &0.6D1 + Q(13) ** 2 / 0.6D1) * Q(9) + 0.2D1 / 0.3D1 * (-Q(7) ** 2 /&
                     & 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6&
                     &D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0.6D1 + &
                     &Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(9) - 0.2D1 * (-Q(5) * &
                     &Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.2D1) * Q(8)&
                     & - 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * &
                     &Q(13) / 0.2D1) * Q(10)) / Q(1) + (Q(6) * Q(7) + Q(9) * Q(10) + Q(&
                     &12) * Q(13)) * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8)&
                     & ** 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2&
                     & / 0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0&
                     &.6D1 + Q(13) ** 2 / 0.6D1) * Q(12) - 0.4D1 / 0.3D1 * (-Q(6) ** 2 /&
                     & 0.3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D&
                     &1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q&
                     &(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(12) + 0.2D1 / 0.3D1 * &
                     &(-Q(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(&
                     &5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) **&
                     & 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(12) - 0.2&
                     &D1 * (-Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) /&
                     & 0.2D1) * Q(11) - 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0&
                     &.2D1 - Q(12) * Q(13) / 0.2D1) * Q(13)) / Q(1))
                      sigma_vec(3) = -Q(1) * (Q(5) * Q(9) * Q(13) - Q(5) * Q(10) * Q(12) - Q(6&
                     &) * Q(8) * Q(13) + Q(6) * Q(10) * Q(11) + Q(7) * Q(8) * Q(12) - Q(&
                     &7) * Q(9) * Q(11)) * ((Q(5) * Q(7) + Q(8) * Q(10) + Q(11) * Q(13))&
                     & * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** 2 / 0.3D&
                     &1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q&
                     &(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13)&
                     & ** 2 / 0.6D1) * Q(7) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D1 - Q(9)&
                     & ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2&
                     & / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0&
                     &.6D1 + Q(13) ** 2 / 0.6D1) * Q(7) - 0.4D1 / 0.3D1 * (-Q(7) ** 2 / &
                     &0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D&
                     &1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0.6D1 + Q&
                     &(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(7) - 0.2D1 * (-Q(6) * Q&
                     &(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(13) / 0.2D1) * Q(6)&
                     & - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D1 - Q(11) * &
                     &Q(13) / 0.2D1) * Q(5)) / Q(1) + (Q(6) * Q(7) + Q(9) * Q(10) + Q(1&
                     &2) * Q(13)) * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) &
                     &** 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 &
                     &/ 0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.&
                     &6D1 + Q(13) ** 2 / 0.6D1) * Q(10) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / &
                     &0.3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1&
                     & + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(&
                     &10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(10) - 0.4D1 / 0.3D1 * (&
                     &-Q(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5&
                     &) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** &
                     &2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(10) - 0.2D&
                     &1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(13) /&
                     & 0.2D1) * Q(9) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.&
                     &2D1 - Q(11) * Q(13) / 0.2D1) * Q(8)) / Q(1) + (Q(7) ** 2 + Q(10) &
                     &** 2 + Q(13) ** 2) * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 &
                     &- Q(8) ** 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9&
                     &) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) **&
                     & 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(13) + 0.2D1 / 0.3D1 * (-Q(6) &
                     &** 2 / 0.3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 &
                     &/ 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6&
                     &D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(13) - 0.4D1 / 0.&
                     &3D1 * (-Q(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D&
                     &1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q&
                     &(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(13)&
                     & - 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * &
                     &Q(13) / 0.2D1) * Q(12) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(&
                     &10) / 0.2D1 - Q(11) * Q(13) / 0.2D1) * Q(11)) / Q(1))
                      sigma_vec(4) = -Q(1) * (Q(5) * Q(9) * Q(13) - Q(5) * Q(10) * Q(12) - Q(6&
                     &) * Q(8) * Q(13) + Q(6) * Q(10) * Q(11) + Q(7) * Q(8) * Q(12) - Q(&
                     &7) * Q(9) * Q(11)) * ((Q(5) ** 2 + Q(8) ** 2 + Q(11) ** 2) * Q(16)&
                     & * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** 2 / 0.3D1 - Q(11&
                     &) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** &
                     &2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / &
                     &0.6D1) * Q(6) - 0.4D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D1 - Q(9) ** 2 / &
                     &0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1&
                     & + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q&
                     &(13) ** 2 / 0.6D1) * Q(6) + 0.2D1 / 0.3D1 * (-Q(7) ** 2 / 0.3D1 - &
                     &Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8)&
                     & ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2&
                     & / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(6) - 0.2D1 * (-Q(5) * Q(6) / 0.&
                     &2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.2D1) * Q(5) - 0.2D1 &
                     &* (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(13) / 0&
                     &.2D1) * Q(7)) / Q(1) + (Q(5) * Q(6) + Q(8) * Q(9) + Q(11) * Q(12)&
                     &) * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** 2 / 0.3&
                     &D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + &
                     &Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13&
                     &) ** 2 / 0.6D1) * Q(9) - 0.4D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D1 - Q(9&
                     &) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** &
                     &2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / &
                     &0.6D1 + Q(13) ** 2 / 0.6D1) * Q(9) + 0.2D1 / 0.3D1 * (-Q(7) ** 2 /&
                     & 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6&
                     &D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0.6D1 + &
                     &Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(9) - 0.2D1 * (-Q(5) * &
                     &Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.2D1) * Q(8)&
                     & - 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * &
                     &Q(13) / 0.2D1) * Q(10)) / Q(1) + (Q(5) * Q(7) + Q(8) * Q(10) + Q(&
                     &11) * Q(13)) * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8)&
                     & ** 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2&
                     & / 0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0&
                     &.6D1 + Q(13) ** 2 / 0.6D1) * Q(12) - 0.4D1 / 0.3D1 * (-Q(6) ** 2 /&
                     & 0.3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D&
                     &1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q&
                     &(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(12) + 0.2D1 / 0.3D1 * &
                     &(-Q(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(&
                     &5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) **&
                     & 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(12) - 0.2&
                     &D1 * (-Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) /&
                     & 0.2D1) * Q(11) - 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0&
                     &.2D1 - Q(12) * Q(13) / 0.2D1) * Q(13)) / Q(1))
                      sigma_vec(5) = -Q(1) * (Q(5) * Q(9) * Q(13) - Q(5) * Q(10) * Q(12) - Q(6&
                     &) * Q(8) * Q(13) + Q(6) * Q(10) * Q(11) + Q(7) * Q(8) * Q(12) - Q(&
                     &7) * Q(9) * Q(11)) * ((Q(5) * Q(6) + Q(8) * Q(9) + Q(11) * Q(12)) &
                     &* Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** 2 / 0.3D1&
                     & - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(&
                     &12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) &
                     &** 2 / 0.6D1) * Q(7) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D1 - Q(9) &
                     &** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 &
                     &/ 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.&
                     &6D1 + Q(13) ** 2 / 0.6D1) * Q(7) - 0.4D1 / 0.3D1 * (-Q(7) ** 2 / 0&
                     &.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1&
                     & + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0.6D1 + Q(&
                     &9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(7) - 0.2D1 * (-Q(6) * Q(&
                     &7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(13) / 0.2D1) * Q(6) &
                     &- 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D1 - Q(11) * Q&
                     &(13) / 0.2D1) * Q(5)) / Q(1) + (Q(6) ** 2 + Q(9) ** 2 + Q(12) ** &
                     &2) * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** 2 / 0.&
                     &3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 +&
                     & Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(1&
                     &3) ** 2 / 0.6D1) * Q(10) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D1 - Q&
                     &(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) *&
                     &* 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 &
                     &/ 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(10) - 0.4D1 / 0.3D1 * (-Q(7) ** &
                     &2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** 2 / &
                     &0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0.6D1&
                     & + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(10) - 0.2D1 * (-Q(6&
                     &) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(13) / 0.2D1) *&
                     & Q(9) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D1 - Q(1&
                     &1) * Q(13) / 0.2D1) * Q(8)) / Q(1) + (Q(6) * Q(7) + Q(9) * Q(10) &
                     &+ Q(12) * Q(13)) * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - &
                     &Q(8) ** 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) &
                     &** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2&
                     & / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(13) + 0.2D1 / 0.3D1 * (-Q(6) **&
                     & 2 / 0.3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / &
                     &0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1&
                     & + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(13) - 0.4D1 / 0.3D&
                     &1 * (-Q(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 &
                     &+ Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6&
                     &) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(13) -&
                     & 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(&
                     &13) / 0.2D1) * Q(12) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10&
                     &) / 0.2D1 - Q(11) * Q(13) / 0.2D1) * Q(11)) / Q(1))
                      sigma_vec(6) = -Q(1) * (Q(5) * Q(9) * Q(13) - Q(5) * Q(10) * Q(12) - Q(6&
                     &) * Q(8) * Q(13) + Q(6) * Q(10) * Q(11) + Q(7) * Q(8) * Q(12) - Q(&
                     &7) * Q(9) * Q(11)) * ((Q(5) ** 2 + Q(8) ** 2 + Q(11) ** 2) * Q(16)&
                     & * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** 2 / 0.3D1 - Q(11&
                     &) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** &
                     &2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / &
                     &0.6D1) * Q(7) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D1 - Q(9) ** 2 / &
                     &0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1&
                     & + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q&
                     &(13) ** 2 / 0.6D1) * Q(7) - 0.4D1 / 0.3D1 * (-Q(7) ** 2 / 0.3D1 - &
                     &Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8)&
                     & ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2&
                     & / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(7) - 0.2D1 * (-Q(6) * Q(7) / 0.&
                     &2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(13) / 0.2D1) * Q(6) - 0.2D1&
                     & * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D1 - Q(11) * Q(13) / &
                     &0.2D1) * Q(5)) / Q(1) + (Q(5) * Q(6) + Q(8) * Q(9) + Q(11) * Q(12&
                     &)) * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** 2 / 0.&
                     &3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 +&
                     & Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(1&
                     &3) ** 2 / 0.6D1) * Q(10) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D1 - Q&
                     &(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) *&
                     &* 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 &
                     &/ 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(10) - 0.4D1 / 0.3D1 * (-Q(7) ** &
                     &2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** 2 / &
                     &0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0.6D1&
                     & + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(10) - 0.2D1 * (-Q(6&
                     &) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(13) / 0.2D1) *&
                     & Q(9) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D1 - Q(1&
                     &1) * Q(13) / 0.2D1) * Q(8)) / Q(1) + (Q(5) * Q(7) + Q(8) * Q(10) &
                     &+ Q(11) * Q(13)) * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - &
                     &Q(8) ** 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) &
                     &** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2&
                     & / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(13) + 0.2D1 / 0.3D1 * (-Q(6) **&
                     & 2 / 0.3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / &
                     &0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1&
                     & + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(13) - 0.4D1 / 0.3D&
                     &1 * (-Q(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 &
                     &+ Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6&
                     &) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(13) -&
                     & 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(&
                     &13) / 0.2D1) * Q(12) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10&
                     &) / 0.2D1 - Q(11) * Q(13) / 0.2D1) * Q(11)) / Q(1))      
            case(2)
                sigma_vec(1) = -Q(1) * (Q(5) * Q(16) * (-0.4D1 / 0.3D1 * (-Q(5) ** 2 / 0.3&
     &D1 - Q(8) ** 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + &
     &Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10)&
     & ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(5) + 0.2D1 / 0.3D1 * (-Q(6&
     &) ** 2 / 0.3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** &
     &2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0&
     &.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(5) + 0.2D1 / 0&
     &.3D1 * (-Q(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3&
     &D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + &
     &Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(5)&
     & - 0.2D1 * (-Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q&
     &(12) / 0.2D1) * Q(6) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10&
     &) / 0.2D1 - Q(11) * Q(13) / 0.2D1) * Q(7)) / Q(1) + Q(8) * Q(16) * (&
     &-0.4D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** 2 / 0.3D1 - Q(11) *&
     &* 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 /&
     & 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6&
     &D1) * Q(8) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D1 - Q(9) ** 2 / 0.3&
     &D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + &
     &Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13&
     &) ** 2 / 0.6D1) * Q(8) + 0.2D1 / 0.3D1 * (-Q(7) ** 2 / 0.3D1 - Q(1&
     &0) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) **&
     & 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / &
     &0.6D1 + Q(12) ** 2 / 0.6D1) * Q(8) - 0.2D1 * (-Q(5) * Q(6) / 0.2D1&
     & - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.2D1) * Q(9) - 0.2D1 * (&
     &-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D1 - Q(11) * Q(13) / 0.2D&
     &1) * Q(10)) / Q(1) + Q(11) * Q(16) * (-0.4D1 / 0.3D1 * (-Q(5) ** 2 /&
     & 0.3D1 - Q(8) ** 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D&
     &1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q&
     &(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(11) + 0.2D1 / 0.3D1 * &
     &(-Q(6) ** 2 / 0.3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5&
     &) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** &
     &2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(11) + 0.2&
     &D1 / 0.3D1 * (-Q(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2&
     & / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.&
     &6D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) &
     &* Q(11) - 0.2D1 * (-Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(&
     &11) * Q(12) / 0.2D1) * Q(12) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8&
     &) * Q(10) / 0.2D1 - Q(11) * Q(13) / 0.2D1) * Q(13)) / Q(1))
      sigma_vec(2) = -Q(1) * (Q(6) * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D&
     &1 - Q(8) ** 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q&
     &(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) &
     &** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(6) - 0.4D1 / 0.3D1 * (-Q(6)&
     & ** 2 / 0.3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2&
     & / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.&
     &6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(6) + 0.2D1 / 0.&
     &3D1 * (-Q(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D&
     &1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q&
     &(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(6) &
     &- 0.2D1 * (-Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(&
     &12) / 0.2D1) * Q(5) - 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10)&
     & / 0.2D1 - Q(12) * Q(13) / 0.2D1) * Q(7)) / Q(1) + Q(9) * Q(16) * (0&
     &.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** 2 / 0.3D1 - Q(11) ** &
     &2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0&
     &.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1&
     &) * Q(9) - 0.4D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D1 - Q(9) ** 2 / 0.3D1&
     & - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(&
     &11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) &
     &** 2 / 0.6D1) * Q(9) + 0.2D1 / 0.3D1 * (-Q(7) ** 2 / 0.3D1 - Q(10)&
     & ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2&
     & / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.&
     &6D1 + Q(12) ** 2 / 0.6D1) * Q(9) - 0.2D1 * (-Q(5) * Q(6) / 0.2D1 -&
     & Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.2D1) * Q(8) - 0.2D1 * (-Q&
     &(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(13) / 0.2D1)&
     & * Q(10)) / Q(1) + Q(12) * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.&
     &3D1 - Q(8) ** 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 +&
     & Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10&
     &) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(12) - 0.4D1 / 0.3D1 * (-Q&
     &(6) ** 2 / 0.3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) *&
     &* 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 /&
     & 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(12) + 0.2D1 &
     &/ 0.3D1 * (-Q(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / &
     &0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1&
     & + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q&
     &(12) - 0.2D1 * (-Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11)&
     & * Q(12) / 0.2D1) * Q(11) - 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) *&
     & Q(10) / 0.2D1 - Q(12) * Q(13) / 0.2D1) * Q(13)) / Q(1))
      sigma_vec(3) = -Q(1) * (Q(7) * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D&
     &1 - Q(8) ** 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q&
     &(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) &
     &** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(7) + 0.2D1 / 0.3D1 * (-Q(6)&
     & ** 2 / 0.3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2&
     & / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.&
     &6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(7) - 0.4D1 / 0.&
     &3D1 * (-Q(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D&
     &1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q&
     &(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(7) &
     &- 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q&
     &(13) / 0.2D1) * Q(6) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10&
     &) / 0.2D1 - Q(11) * Q(13) / 0.2D1) * Q(5)) / Q(1) + Q(10) * Q(16) * &
     &(0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** 2 / 0.3D1 - Q(11) *&
     &* 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 /&
     & 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6&
     &D1) * Q(10) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D1 - Q(9) ** 2 / 0.&
     &3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 +&
     & Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(1&
     &3) ** 2 / 0.6D1) * Q(10) - 0.4D1 / 0.3D1 * (-Q(7) ** 2 / 0.3D1 - Q&
     &(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) &
     &** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 &
     &/ 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(10) - 0.2D1 * (-Q(6) * Q(7) / 0.&
     &2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(13) / 0.2D1) * Q(9) - 0.2D1&
     & * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D1 - Q(11) * Q(13) / &
     &0.2D1) * Q(8)) / Q(1) + Q(13) * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2&
     & / 0.3D1 - Q(8) ** 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.&
     &6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 +&
     & Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(13) + 0.2D1 / 0.3D1 &
     &* (-Q(6) ** 2 / 0.3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q&
     &(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) *&
     &* 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(13) - 0&
     &.4D1 / 0.3D1 * (-Q(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) **&
     & 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / &
     &0.6D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1&
     &) * Q(13) - 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 -&
     & Q(12) * Q(13) / 0.2D1) * Q(12) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - &
     &Q(8) * Q(10) / 0.2D1 - Q(11) * Q(13) / 0.2D1) * Q(11)) / Q(1))
      sigma_vec(4) = -Q(1) * (Q(5) * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D&
     &1 - Q(8) ** 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q&
     &(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) &
     &** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(6) - 0.4D1 / 0.3D1 * (-Q(6)&
     & ** 2 / 0.3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2&
     & / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.&
     &6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(6) + 0.2D1 / 0.&
     &3D1 * (-Q(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D&
     &1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q&
     &(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(6) &
     &- 0.2D1 * (-Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(&
     &12) / 0.2D1) * Q(5) - 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10)&
     & / 0.2D1 - Q(12) * Q(13) / 0.2D1) * Q(7)) / Q(1) + Q(8) * Q(16) * (0&
     &.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** 2 / 0.3D1 - Q(11) ** &
     &2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0&
     &.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1&
     &) * Q(9) - 0.4D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D1 - Q(9) ** 2 / 0.3D1&
     & - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(&
     &11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) &
     &** 2 / 0.6D1) * Q(9) + 0.2D1 / 0.3D1 * (-Q(7) ** 2 / 0.3D1 - Q(10)&
     & ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2&
     & / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.&
     &6D1 + Q(12) ** 2 / 0.6D1) * Q(9) - 0.2D1 * (-Q(5) * Q(6) / 0.2D1 -&
     & Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.2D1) * Q(8) - 0.2D1 * (-Q&
     &(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(13) / 0.2D1)&
     & * Q(10)) / Q(1) + Q(11) * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.&
     &3D1 - Q(8) ** 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 +&
     & Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10&
     &) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(12) - 0.4D1 / 0.3D1 * (-Q&
     &(6) ** 2 / 0.3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) *&
     &* 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 /&
     & 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(12) + 0.2D1 &
     &/ 0.3D1 * (-Q(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / &
     &0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1&
     & + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q&
     &(12) - 0.2D1 * (-Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11)&
     & * Q(12) / 0.2D1) * Q(11) - 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) *&
     & Q(10) / 0.2D1 - Q(12) * Q(13) / 0.2D1) * Q(13)) / Q(1))
      sigma_vec(5) = -Q(1) * (Q(6) * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D&
     &1 - Q(8) ** 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q&
     &(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) &
     &** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(7) + 0.2D1 / 0.3D1 * (-Q(6)&
     & ** 2 / 0.3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2&
     & / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.&
     &6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(7) - 0.4D1 / 0.&
     &3D1 * (-Q(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D&
     &1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q&
     &(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(7) &
     &- 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q&
     &(13) / 0.2D1) * Q(6) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10&
     &) / 0.2D1 - Q(11) * Q(13) / 0.2D1) * Q(5)) / Q(1) + Q(9) * Q(16) * (&
     &0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** 2 / 0.3D1 - Q(11) **&
     & 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / &
     &0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D&
     &1) * Q(10) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D1 - Q(9) ** 2 / 0.3&
     &D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + &
     &Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13&
     &) ** 2 / 0.6D1) * Q(10) - 0.4D1 / 0.3D1 * (-Q(7) ** 2 / 0.3D1 - Q(&
     &10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) *&
     &* 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 /&
     & 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(10) - 0.2D1 * (-Q(6) * Q(7) / 0.2&
     &D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(13) / 0.2D1) * Q(9) - 0.2D1 &
     &* (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D1 - Q(11) * Q(13) / 0&
     &.2D1) * Q(8)) / Q(1) + Q(12) * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 &
     &/ 0.3D1 - Q(8) ** 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6&
     &D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + &
     &Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(13) + 0.2D1 / 0.3D1 *&
     & (-Q(6) ** 2 / 0.3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(&
     &5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) **&
     & 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(13) - 0.&
     &4D1 / 0.3D1 * (-Q(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** &
     &2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0&
     &.6D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1)&
     & * Q(13) - 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - &
     &Q(12) * Q(13) / 0.2D1) * Q(12) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q&
     &(8) * Q(10) / 0.2D1 - Q(11) * Q(13) / 0.2D1) * Q(11)) / Q(1))
      sigma_vec(6) = -Q(1) * (Q(5) * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D&
     &1 - Q(8) ** 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q&
     &(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) &
     &** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(7) + 0.2D1 / 0.3D1 * (-Q(6)&
     & ** 2 / 0.3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2&
     & / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.&
     &6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(7) - 0.4D1 / 0.&
     &3D1 * (-Q(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D&
     &1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q&
     &(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(7) &
     &- 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q&
     &(13) / 0.2D1) * Q(6) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10&
     &) / 0.2D1 - Q(11) * Q(13) / 0.2D1) * Q(5)) / Q(1) + Q(8) * Q(16) * (&
     &0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** 2 / 0.3D1 - Q(11) **&
     & 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / &
     &0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D&
     &1) * Q(10) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D1 - Q(9) ** 2 / 0.3&
     &D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + &
     &Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13&
     &) ** 2 / 0.6D1) * Q(10) - 0.4D1 / 0.3D1 * (-Q(7) ** 2 / 0.3D1 - Q(&
     &10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) *&
     &* 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 /&
     & 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(10) - 0.2D1 * (-Q(6) * Q(7) / 0.2&
     &D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(13) / 0.2D1) * Q(9) - 0.2D1 &
     &* (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D1 - Q(11) * Q(13) / 0&
     &.2D1) * Q(8)) / Q(1) + Q(11) * Q(16) * (0.2D1 / 0.3D1 * (-Q(5) ** 2 &
     &/ 0.3D1 - Q(8) ** 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6&
     &D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + &
     &Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(13) + 0.2D1 / 0.3D1 *&
     & (-Q(6) ** 2 / 0.3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(&
     &5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) **&
     & 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(13) - 0.&
     &4D1 / 0.3D1 * (-Q(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** &
     &2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0&
     &.6D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1)&
     & * Q(13) - 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - &
     &Q(12) * Q(13) / 0.2D1) * Q(12) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q&
     &(8) * Q(10) / 0.2D1 - Q(11) * Q(13) / 0.2D1) * Q(11)) / Q(1))
            case default
            print *, 'GPRLEversion=',GPRLEversion, ' not implemented!'
            stop
        end select
        ! Add the pressure contribution diff_k p delta_ik
        detA = Q(5) * Q(9) * Q(13) - Q(5) * Q(10) * Q(12) - Q(8) * Q(6) * Q(13) + Q(8) * Q(7) * Q(12) + Q(11) * Q(6) * Q(10) - Q(11) * Q(7) * Q(9)

        if(addp) then
            sigma_vec(1)=+(Q(15)+2.0/3.0*Q(16))*(detA)**2*(1-detA)+sigma_vec(1)       ! A
            sigma_vec(2)=+(Q(15)+2.0/3.0*Q(16))*(detA)**2*(1-detA)+sigma_vec(2)
            sigma_vec(3)=+(Q(15)+2.0/3.0*Q(16))*(detA)**2*(1-detA)+sigma_vec(3) 
            stressnorm  = SQRT( 0.5 * ( (sigma_vec(1)-sigma_vec(2))**2 + (sigma_vec(2)-sigma_vec(3))**2 + (sigma_vec(3)-sigma_vec(1))**2 + 6.*(sigma_vec(4)**2+sigma_vec(6)**2+sigma_vec(5)**2) ) ) 
        end if
    end subroutine computeGPRLEstress
    
    subroutine ComputeAcontribution(Eshvar,Q)
        use Parameters, only : nVar

        implicit none
        real, intent(in) :: Q(nVar)
        real :: Eshvar,p, detA
        detA = Q(5) * Q(9) * Q(13) - Q(5) * Q(10) * Q(12) - Q(8) * Q(6) * Q(13) + Q(8) * Q(7) * Q(12) + Q(11) * Q(6) * Q(10) - Q(11) * Q(7) * Q(9)

        Eshvar = Q(16) / Q(1) * (Q(5) * Q(9) * Q(13) - Q(5) * Q(10) * Q(12     ) - Q(8) * Q(6) * Q(13) + Q(8) * Q(7) * Q(12) + Q(11) * Q(6) * Q(10) - Q(11) * Q(7) * Q(9)) * ((-Q(5) ** 2 / 0.3D1 - Q(8) ** 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) ** 2 + (-Q(6) ** 2 / 0.3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) ** 2 + (-Q(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) ** 2 + 0.2D1 * (-Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.2D1) ** 2 + 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(13) / 0.2D1) ** 2 + 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D1 - Q(11) * Q(13) / 0.2D1) ** 2)
        p=(Q(15)+2.0/3.0*Q(16))/Q(1)*detA*(1-detA)**2
        Eshvar=Q(1)*(Eshvar+p)
    end subroutine ComputeAcontribution

    !*********************************************************************************
    function SmoothInterface(r,ICsig,epsilon,smooth_order_in)
        implicit none
        real    :: r
        real    :: SmoothInterface,smooth_order
        real    :: eta,ICsig,xi 
        real, optional :: epsilon
        real, optional :: smooth_order_in
        
        if(.not. present(epsilon)) then
            epsilon=1.e-9    
        end if
        if(.not. present(smooth_order_in)) then
            smooth_order=4   
        end if
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
    
    function xy2t(x,ind)
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
    
    ! ********************************* GPR+DI for LE specific var and subs *********************************
    function GPRsigma2A(sigmaK,theta,lambda1,mu1,rho0,eps)
        implicit none
        real, intent(in):: sigmaK(6),theta,lambda1,mu1,rho0
        real            :: GPRsigma2A(9)
        ! Local variables
        real            :: eps,unkvarsig(6),fdf(6),fp(6),JJ(6,6),ff(6)
        integer         :: iNewton,maxNewton,ierr
        
        maxNewton=50
        !eps=1.e-8*abs(sqrt((lambda1+2.0*mu1)/rho0))
        unkvarsig=0.1
        unkvarsig(1)=1
        unkvarsig(3)=1
        unkvarsig(6)=1
        !unkvarsig=sigmaK

        do iNewton=1,maxNewton
            call A2sigmaJacobian(JJ,unkvarsig,theta,lambda1,mu1,rho0)
            call A2sigmaComputeF(fp,unkvarsig,theta,lambda1,mu1,rho0)
            ff=fp-sigmaK
            call LinSolve(6,JJ,ff,fdf)
            unkvarsig=unkvarsig-fdf
            if(sqrt(abs(sum(ff**2)))<eps) then
                goto 146    
            end if
        end do
        print *, ff
        print *, sum(ff(1:6))
        print *, 'Newton in sigma2A does not converge, err=',sqrt(abs(sum(ff(1:6)**2))),">eps=",eps
146     continue
        call A2sigmaAssignQA(GPRsigma2A,unkvarsig, theta)
    end function GPRsigma2A
    
    subroutine A2sigmaAssignQA(QA,sig, theta)
        implicit none
        real,intent(out) :: QA(9)
        real             :: Aij(3,3)
        real, intent(in) :: sig(6), theta
        integer          :: k,i,j
        Aij=0.
        Aij(1,1) = sig(1) * cos(theta)
        Aij(1,2) = -sig(1) * sin(theta)
        Aij(2,1) = sig(2) * cos(theta) + sig(3) * sin(theta)
        Aij(2,2) = -sig(2) * sin(theta) + sig(3) * cos(theta)
        Aij(3,1) = sig(4) * cos(theta) + sig(5) * sin(theta)
        Aij(3,2) = -sig(4) * sin(theta) + sig(5) * cos(theta)
        Aij(3,3) = sig(6)
        k=0
        do i=1,3
            do j=1,3
                k=k+1
                QA(k)=Aij(i,j)
            end do
        end do
        
        
    end subroutine A2sigmaAssignQA
    
    subroutine A2sigmaComputeF(fp,sig,theta,lambda1,mu1,rho0)
        implicit none
        real,intent(out) :: fp(6)
        real, intent(in) :: sig(6),theta,lambda1,mu1,rho0
        select case(GPRLEversion)
            case(0,1)
      fp(1) = (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * co&
     &s(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) *&
     & sig(1) * sin(theta) * sig(6)) * (sig(1) ** 2 * cos(theta) ** 2 + &
     &(sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + (sig(4) * cos(t&
     &heta) + sig(5) * sin(theta)) ** 2) * mu1 * (-0.2D1 / 0.3D1 * sig(1&
     &) ** 2 * cos(theta) ** 2 - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + &
     &sig(3) * sin(theta)) ** 2 - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) +&
     & sig(5) * sin(theta)) ** 2 + sig(1) ** 2 * sin(theta) ** 2 / 0.3D1&
     & + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 + (-s&
     &ig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(6) **&
     & 2 / 0.3D1) + (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + (sig(2) * &
     &cos(theta) + sig(3) * sin(theta)) * sig(1) * sin(theta) * sig(6)) &
     &** 2 * (-sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * co&
     &s(theta)) * sig(6) + 0.1D1 - (sig(2) * cos(theta) + sig(3) * sin(t&
     &heta)) * sig(1) * sin(theta) * sig(6))
      fp(2) = (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * co&
     &s(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) *&
     & sig(1) * sin(theta) * sig(6)) * (sig(1) ** 2 * sin(theta) ** 2 + &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + (-sig(4) * sin&
     &(theta) + sig(5) * cos(theta)) ** 2) * mu1 * (sig(1) ** 2 * cos(th&
     &eta) ** 2 / 0.3D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) **&
     & 2 / 0.3D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.&
     &3D1 - 0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 - 0.2D1 / 0.3D&
     &1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 - 0.2D1 / 0.&
     &3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 + sig(6) *&
     &* 2 / 0.3D1) + (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) *&
     & (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + (sig(2) *&
     & cos(theta) + sig(3) * sin(theta)) * sig(1) * sin(theta) * sig(6))&
     & ** 2 * (-sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) * sig(6) + 0.1D1 - (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) * sig(1) * sin(theta) * sig(6))
      fp(3) = (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * co&
     &s(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) *&
     & sig(1) * sin(theta) * sig(6)) * sig(6) ** 2 * mu1 * (sig(1) ** 2 &
     &* cos(theta) ** 2 / 0.3D1 + (sig(2) * cos(theta) + sig(3) * sin(th&
     &eta)) ** 2 / 0.3D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) *&
     &* 2 / 0.3D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 + (-sig(2) * s&
     &in(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 + (-sig(4) * sin(the&
     &ta) + sig(5) * cos(theta)) ** 2 / 0.3D1 - 0.2D1 / 0.3D1 * sig(6) *&
     &* 2) + (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * (-sig(2&
     &) * sin(theta) + sig(3) * cos(theta)) * sig(6) + (sig(2) * cos(the&
     &ta) + sig(3) * sin(theta)) * sig(1) * sin(theta) * sig(6)) ** 2 * &
     &(-sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta&
     &)) * sig(6) + 0.1D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) &
     &* sig(1) * sin(theta) * sig(6))
      fp(4) = 0.2D1 * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig&
     &(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(t&
     &heta)) * sig(1) * sin(theta) * sig(6)) * (-sig(1) ** 2 * cos(theta&
     &) * sin(theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-s&
     &ig(2) * sin(theta) + sig(3) * cos(theta)) + (sig(4) * cos(theta) +&
     & sig(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta&
     &)) + 0.1D1) * mu1 * (sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D1&
     & - (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(th&
     &eta) + sig(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0&
     &.2D1)
      fp(5) = -(sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) &
     &* sig(1) * sin(theta) * sig(6)) * (0.1D1 + (-sig(4) * sin(theta) +&
     & sig(5) * cos(theta)) * sig(6)) * mu1 * (-sig(4) * sin(theta) + si&
     &g(5) * cos(theta)) * sig(6)
      fp(6) = -(sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) &
     &* sig(1) * sin(theta) * sig(6)) * (0.1D1 + (sig(4) * cos(theta) + &
     &sig(5) * sin(theta)) * sig(6)) * mu1 * (sig(4) * cos(theta) + sig(&
     &5) * sin(theta)) * sig(6)
        case(2)
 fp(1) = -rho0 * (sig(1) * cos(theta) * mu1 * (-0.4D1 / 0.3D1 * (-s&
     &ig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * si&
     &n(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (&
     &-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4&
     &) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1 + sig(6) ** 2 /&
     & 0.6D1) * sig(1) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(1) ** 2 * si&
     &n(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta&
     &)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** &
     &2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(&
     &theta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) &
     &+ sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * sig(1&
     &) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + sig(1) **&
     & 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin&
     &(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)&
     &) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) &
     &* sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(&
     &theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * sig(1) * cos(theta) &
     &+ 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D1 - (sig(2)&
     & * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig&
     &(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(theta) + sig(5) * sin(th&
     &eta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.2D1) * si&
     &g(1) * sin(theta)) / rho0 + (sig(2) * cos(theta) + sig(3) * sin(th&
     &eta)) * mu1 * (-0.4D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / &
     &0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 -&
     & (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1)&
     & ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) *&
     & cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(t&
     &heta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) + 0.2D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta&
     &) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2&
     & / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3&
     &D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (sig(2) * co&
     &s(theta) + sig(3) * sin(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 / &
     &0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(thet&
     &a) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + si&
     &g(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / &
     &0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 &
     &+ (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * (si&
     &g(2) * cos(theta) + sig(3) * sin(theta)) - 0.2D1 * (sig(1) ** 2 * &
     &cos(theta) * sin(theta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * &
     &sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.2D1&
     & - (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(th&
     &eta) + sig(5) * cos(theta)) / 0.2D1) * (-sig(2) * sin(theta) + sig&
     &(3) * cos(theta))) / rho0 + (sig(4) * cos(theta) + sig(5) * sin(th&
     &eta)) * mu1 * (-0.4D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / &
     &0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 -&
     & (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1)&
     & ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) *&
     & cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(t&
     &heta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (sig(4) * cos(theta) &
     &+ sig(5) * sin(theta)) + 0.2D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta&
     &) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2&
     & / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3&
     &D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (sig(4) * co&
     &s(theta) + sig(5) * sin(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 / &
     &0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(thet&
     &a) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + si&
     &g(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / &
     &0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 &
     &+ (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * (si&
     &g(4) * cos(theta) + sig(5) * sin(theta)) - 0.2D1 * (sig(1) ** 2 * &
     &cos(theta) * sin(theta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * &
     &sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.2D1&
     & - (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(th&
     &eta) + sig(5) * cos(theta)) / 0.2D1) * (-sig(4) * sin(theta) + sig&
     &(5) * cos(theta)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * &
     &sig(6) ** 2) / rho0) + (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(&
     &theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + s&
     &ig(1) * sin(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) *&
     & sig(6)) ** 2 * (-sig(1) * cos(theta) * (-sig(2) * sin(theta) + si&
     &g(3) * cos(theta)) * sig(6) + 0.1D1 - sig(1) * sin(theta) * (sig(2&
     &) * cos(theta) + sig(3) * sin(theta)) * sig(6))
      fp(2) = -rho0 * (-sig(1) * sin(theta) * mu1 * (-0.2D1 / 0.3D1 * (-&
     &sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig&
     &(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * s&
     &in(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(&
     &4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1 + sig(6) ** 2 &
     &/ 0.6D1) * sig(1) * sin(theta) + 0.4D1 / 0.3D1 * (-sig(1) ** 2 * s&
     &in(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(thet&
     &a)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) **&
     & 2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos&
     &(theta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta)&
     & + sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * sig(&
     &1) * sin(theta) - 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + sig(1) *&
     &* 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * si&
     &n(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta&
     &)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2)&
     & * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin&
     &(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * sig(1) * sin(theta)&
     & - 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D1 - (sig(2&
     &) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta) + si&
     &g(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(theta) + sig(5) * sin(t&
     &heta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.2D1) * s&
     &ig(1) * cos(theta)) / rho0 + (-sig(2) * sin(theta) + sig(3) * cos(&
     &theta)) * mu1 * (0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 /&
     & 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 &
     &- (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1&
     &) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(&
     &theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (-sig(2) * sin(theta&
     &) + sig(3) * cos(theta)) - 0.4D1 / 0.3D1 * (-sig(1) ** 2 * sin(the&
     &ta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) **&
     & 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0&
     &.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta&
     &) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig&
     &(5) * sin(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (-sig(2) *&
     & sin(theta) + sig(3) * cos(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** 2&
     & / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(t&
     &heta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) +&
     & sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2&
     & / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6&
     &D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) - 0.2D1 * (sig(1) ** &
     &2 * cos(theta) * sin(theta) / 0.2D1 - (sig(2) * cos(theta) + sig(3&
     &) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0&
     &.2D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * si&
     &n(theta) + sig(5) * cos(theta)) / 0.2D1) * (sig(2) * cos(theta) + &
     &sig(3) * sin(theta))) / rho0 + (-sig(4) * sin(theta) + sig(5) * co&
     &s(theta)) * mu1 * (0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2&
     & / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D&
     &1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig&
     &(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3&
     &) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * co&
     &s(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (-sig(4) * sin(the&
     &ta) + sig(5) * cos(theta)) - 0.4D1 / 0.3D1 * (-sig(1) ** 2 * sin(t&
     &heta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) &
     &** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 /&
     & 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(the&
     &ta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (-sig(4)&
     & * sin(theta) + sig(5) * cos(theta)) + 0.2D1 / 0.3D1 * (-sig(6) **&
     & 2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos&
     &(theta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta)&
     & + sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) **&
     & 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0&
     &.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) &
     &* (-sig(4) * sin(theta) + sig(5) * cos(theta)) - 0.2D1 * (sig(1) *&
     &* 2 * cos(theta) * sin(theta) / 0.2D1 - (sig(2) * cos(theta) + sig&
     &(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) /&
     & 0.2D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * &
     &sin(theta) + sig(5) * cos(theta)) / 0.2D1) * (sig(4) * cos(theta) &
     &+ sig(5) * sin(theta)) + (-sig(4) * sin(theta) + sig(5) * cos(thet&
     &a)) * sig(6) ** 2) / rho0) + (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) &
     &* cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(&
     &6) + sig(1) * sin(theta) * (sig(2) * cos(theta) + sig(3) * sin(the&
     &ta)) * sig(6)) ** 2 * (-sig(1) * cos(theta) * (-sig(2) * sin(theta&
     &) + sig(3) * cos(theta)) * sig(6) + 0.1D1 - sig(1) * sin(theta) * &
     &(sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6))
      fp(3) = -sig(6) * mu1 * (0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta&
     &) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 &
     &/ 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1&
     & + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) +&
     & sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5&
     &) * cos(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * sig(6) + 0.2&
     &D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * &
     &sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(th&
     &eta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * cos(theta&
     &) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 &
     &/ 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.6D1&
     & + sig(6) ** 2 / 0.6D1) * sig(6) - 0.4D1 / 0.3D1 * (-sig(6) ** 2 /&
     & 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(the&
     &ta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 /&
     & 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1&
     & + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * si&
     &g(6) + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 * sig(6) &
     &+ (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 * sig(6)) + (lambda1&
     & + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * (-sig(2) * sin(th&
     &eta) + sig(3) * cos(theta)) * sig(6) + sig(1) * sin(theta) * (sig(&
     &2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) ** 2 * (-sig(1) *&
     & cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6&
     &) + 0.1D1 - sig(1) * sin(theta) * (sig(2) * cos(theta) + sig(3) * &
     &sin(theta)) * sig(6))
      fp(4) = -rho0 * (sig(1) * cos(theta) * mu1 * (-0.2D1 / 0.3D1 * (-s&
     &ig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * si&
     &n(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (&
     &-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4&
     &) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1 + sig(6) ** 2 /&
     & 0.6D1) * sig(1) * sin(theta) + 0.4D1 / 0.3D1 * (-sig(1) ** 2 * si&
     &n(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta&
     &)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** &
     &2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(&
     &theta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) &
     &+ sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * sig(1&
     &) * sin(theta) - 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + sig(1) **&
     & 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin&
     &(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)&
     &) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) &
     &* sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(&
     &theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * sig(1) * sin(theta) &
     &- 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D1 - (sig(2)&
     & * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig&
     &(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(theta) + sig(5) * sin(th&
     &eta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.2D1) * si&
     &g(1) * cos(theta)) / rho0 + (sig(2) * cos(theta) + sig(3) * sin(th&
     &eta)) * mu1 * (0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0&
     &.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - &
     &(sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) &
     &** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * &
     &cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(th&
     &eta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (-sig(2) * sin(theta) &
     &+ sig(3) * cos(theta)) - 0.4D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta&
     &) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2&
     & / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3&
     &D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (-sig(2) * s&
     &in(theta) + sig(3) * cos(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 /&
     & 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(the&
     &ta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 /&
     & 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1&
     & + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * (-&
     &sig(2) * sin(theta) + sig(3) * cos(theta)) - 0.2D1 * (sig(1) ** 2 &
     &* cos(theta) * sin(theta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) &
     &* sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.2&
     &D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(&
     &theta) + sig(5) * cos(theta)) / 0.2D1) * (sig(2) * cos(theta) + si&
     &g(3) * sin(theta))) / rho0 + (sig(4) * cos(theta) + sig(5) * sin(t&
     &heta)) * mu1 * (0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / &
     &0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 -&
     & (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1)&
     & ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) *&
     & cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(t&
     &heta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (-sig(4) * sin(theta)&
     & + sig(5) * cos(theta)) - 0.4D1 / 0.3D1 * (-sig(1) ** 2 * sin(thet&
     &a) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** &
     &2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.&
     &3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta)&
     & + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(&
     &5) * sin(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (-sig(4) * &
     &sin(theta) + sig(5) * cos(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 &
     &/ 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(th&
     &eta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + &
     &sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 &
     &/ 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D&
     &1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * (&
     &-sig(4) * sin(theta) + sig(5) * cos(theta)) - 0.2D1 * (sig(1) ** 2&
     & * cos(theta) * sin(theta) / 0.2D1 - (sig(2) * cos(theta) + sig(3)&
     & * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.&
     &2D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin&
     &(theta) + sig(5) * cos(theta)) / 0.2D1) * (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) + (-sig(4) * sin(theta) + sig(5) * cos(theta))&
     & * sig(6) ** 2) / rho0)
      fp(5) = -rho0 * (-sig(1) * sin(theta) * mu1 * (-(-sig(4) * sin(the&
     &ta) + sig(5) * cos(theta)) * sig(6) * sig(1) * sin(theta) + (sig(4&
     &) * cos(theta) + sig(5) * sin(theta)) * sig(6) * sig(1) * cos(thet&
     &a)) / rho0 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) * mu1 * &
     &((-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig(6) * (-sig(2) &
     &* sin(theta) + sig(3) * cos(theta)) + (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) * sig(6) * (sig(2) * cos(theta) + sig(3) * sin(the&
     &ta))) / rho0 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) * mu1 &
     &* (0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(&
     &2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * co&
     &s(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * sin(t&
     &heta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) &
     &** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 /&
     & 0.6D1 + sig(6) ** 2 / 0.6D1) * sig(6) + 0.2D1 / 0.3D1 * (-sig(1) &
     &** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * &
     &cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(th&
     &eta)) ** 2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(&
     &2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * co&
     &s(theta) + sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1&
     &) * sig(6) - 0.4D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + sig(1) ** 2 *&
     & cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(the&
     &ta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) **&
     & 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(thet&
     &a) + sig(5) * cos(theta)) ** 2 / 0.6D1) * sig(6) + (-sig(4) * sin(&
     &theta) + sig(5) * cos(theta)) ** 2 * sig(6) + (sig(4) * cos(theta)&
     & + sig(5) * sin(theta)) ** 2 * sig(6)) / rho0)
      fp(6) = -rho0 * (sig(1) * cos(theta) * mu1 * (-(-sig(4) * sin(thet&
     &a) + sig(5) * cos(theta)) * sig(6) * sig(1) * sin(theta) + (sig(4)&
     & * cos(theta) + sig(5) * sin(theta)) * sig(6) * sig(1) * cos(theta&
     &)) / rho0 + (sig(2) * cos(theta) + sig(3) * sin(theta)) * mu1 * ((&
     &-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig(6) * (-sig(2) * &
     &sin(theta) + sig(3) * cos(theta)) + (sig(4) * cos(theta) + sig(5) &
     &* sin(theta)) * sig(6) * (sig(2) * cos(theta) + sig(3) * sin(theta&
     &))) / rho0 + (sig(4) * cos(theta) + sig(5) * sin(theta)) * mu1 * (&
     &0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) &
     &* cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(t&
     &heta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * sin(thet&
     &a) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** &
     &2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.&
     &6D1 + sig(6) ** 2 / 0.6D1) * sig(6) + 0.2D1 / 0.3D1 * (-sig(1) ** &
     &2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos&
     &(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta&
     &)) ** 2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) &
     &* cos(theta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(t&
     &heta) + sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) *&
     & sig(6) - 0.4D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + sig(1) ** 2 * co&
     &s(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)&
     &) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 &
     &/ 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(t&
     &heta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) &
     &+ sig(5) * cos(theta)) ** 2 / 0.6D1) * sig(6) + (-sig(4) * sin(the&
     &ta) + sig(5) * cos(theta)) ** 2 * sig(6) + (sig(4) * cos(theta) + &
     &sig(5) * sin(theta)) ** 2 * sig(6)) / rho0)                
        end select
    end subroutine A2sigmaComputeF
    
    subroutine A2sigmaJacobian(J,sig,theta,lambda1,mu1,rho0)
        implicit none
        real :: s1,s2,s3,s4
        real,intent(out) :: J(6,6)
        real, intent(in) :: sig(6),theta,lambda1,mu1,rho0
        select case(GPRLEversion)
            case(0,1)
    J(1,1) = (cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)&
     &) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sin(the&
     &ta) * sig(6)) * (sig(1) ** 2 * cos(theta) ** 2 + (sig(2) * cos(the&
     &ta) + sig(3) * sin(theta)) ** 2 + (sig(4) * cos(theta) + sig(5) * &
     &sin(theta)) ** 2) * mu1 * (-0.2D1 / 0.3D1 * sig(1) ** 2 * cos(thet&
     &a) ** 2 - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(thet&
     &a)) ** 2 - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(the&
     &ta)) ** 2 + sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 + (-sig(2) * sin&
     &(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 + (-sig(4) * sin(theta&
     &) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(6) ** 2 / 0.3D1) + 0.2&
     &D1 * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(t&
     &heta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * si&
     &g(1) * sin(theta) * sig(6)) * sig(1) * cos(theta) ** 2 * mu1 * (-0&
     &.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 - 0.2D1 / 0.3D1 * (si&
     &g(2) * cos(theta) + sig(3) * sin(theta)) ** 2 - 0.2D1 / 0.3D1 * (s&
     &ig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 + sig(1) ** 2 * sin&
     &(theta) ** 2 / 0.3D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)&
     &) ** 2 / 0.3D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2&
     & / 0.3D1 + sig(6) ** 2 / 0.3D1) + (sig(1) * cos(theta) * (-sig(2) &
     &* sin(theta) + sig(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta&
     &) + sig(3) * sin(theta)) * sig(1) * sin(theta) * sig(6)) * (sig(1)&
     & ** 2 * cos(theta) ** 2 + (sig(2) * cos(theta) + sig(3) * sin(thet&
     &a)) ** 2 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2) * mu1&
     & * (-0.4D1 / 0.3D1 * sig(1) * cos(theta) ** 2 + 0.2D1 / 0.3D1 * si&
     &g(1) * sin(theta) ** 2) + 0.2D1 * (lambda1 + 0.2D1 / 0.3D1 * mu1) * (si&
     &g(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) *&
     & sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * s&
     &in(theta) * sig(6)) * (-sig(1) * cos(theta) * (-sig(2) * sin(theta&
     &) + sig(3) * cos(theta)) * sig(6) + 0.1D1 - (sig(2) * cos(theta) +&
     & sig(3) * sin(theta)) * sig(1) * sin(theta) * sig(6)) * (cos(theta&
     &) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + (sig(2&
     &) * cos(theta) + sig(3) * sin(theta)) * sin(theta) * sig(6)) + (lambda1&
     & + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * (-sig(2) * sin(th&
     &eta) + sig(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) * sig(1) * sin(theta) * sig(6)) ** 2 * (-cos(thet&
     &a) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) - (sig(&
     &2) * cos(theta) + sig(3) * sin(theta)) * sin(theta) * sig(6))
      J(1,2) = 0.2D1 * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + si&
     &g(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) * sig(1) * sin(theta) * sig(6)) * (sig(2) * cos(theta) + s&
     &ig(3) * sin(theta)) * cos(theta) * mu1 * (-0.2D1 / 0.3D1 * sig(1) &
     &** 2 * cos(theta) ** 2 - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) ** 2 - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) ** 2 + sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 +&
     & (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 + (-sig&
     &(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(6) ** 2&
     & / 0.3D1) + (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) * sig(1) * sin(theta) * sig(6)) * (sig(1) ** 2 * cos(theta) ** &
     &2 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + (sig(4) * c&
     &os(theta) + sig(5) * sin(theta)) ** 2) * mu1 * (-0.4D1 / 0.3D1 * (&
     &sig(2) * cos(theta) + sig(3) * sin(theta)) * cos(theta) - 0.2D1 / &
     &0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(theta))
      J(1,3) = (sig(1) * cos(theta) ** 2 * sig(6) + sin(theta) ** 2 * si&
     &g(1) * sig(6)) * (sig(1) ** 2 * cos(theta) ** 2 + (sig(2) * cos(th&
     &eta) + sig(3) * sin(theta)) ** 2 + (sig(4) * cos(theta) + sig(5) *&
     & sin(theta)) ** 2) * mu1 * (-0.2D1 / 0.3D1 * sig(1) ** 2 * cos(the&
     &ta) ** 2 - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(the&
     &ta)) ** 2 - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(th&
     &eta)) ** 2 + sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 + (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 + (-sig(4) * sin(thet&
     &a) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(6) ** 2 / 0.3D1) + 0.&
     &2D1 * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(&
     &theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * s&
     &ig(1) * sin(theta) * sig(6)) * (sig(2) * cos(theta) + sig(3) * sin&
     &(theta)) * sin(theta) * mu1 * (-0.2D1 / 0.3D1 * sig(1) ** 2 * cos(&
     &theta) ** 2 - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) ** 2 - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin&
     &(theta)) ** 2 + sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 + (-sig(2) *&
     & sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 + (-sig(4) * sin(t&
     &heta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(6) ** 2 / 0.3D1) +&
     & (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta&
     &)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1)&
     & * sin(theta) * sig(6)) * (sig(1) ** 2 * cos(theta) ** 2 + (sig(2)&
     & * cos(theta) + sig(3) * sin(theta)) ** 2 + (sig(4) * cos(theta) +&
     & sig(5) * sin(theta)) ** 2) * mu1 * (-0.4D1 / 0.3D1 * (sig(2) * co&
     &s(theta) + sig(3) * sin(theta)) * sin(theta) + 0.2D1 / 0.3D1 * (-s&
     &ig(2) * sin(theta) + sig(3) * cos(theta)) * cos(theta)) + 0.2D1 * &
     &(lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * (-sig(2) * sin&
     &(theta) + sig(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta) + s&
     &ig(3) * sin(theta)) * sig(1) * sin(theta) * sig(6)) * (-sig(1) * c&
     &os(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) &
     &+ 0.1D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * s&
     &in(theta) * sig(6)) * (sig(1) * cos(theta) ** 2 * sig(6) + sin(the&
     &ta) ** 2 * sig(1) * sig(6)) + (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1)&
     & * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig&
     &(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * sin(t&
     &heta) * sig(6)) ** 2 * (-sig(1) * cos(theta) ** 2 * sig(6) - sin(t&
     &heta) ** 2 * sig(1) * sig(6))
      J(1,4) = 0.2D1 * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + si&
     &g(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) * sig(1) * sin(theta) * sig(6)) * (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) * cos(theta) * mu1 * (-0.2D1 / 0.3D1 * sig(1) &
     &** 2 * cos(theta) ** 2 - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) ** 2 - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) ** 2 + sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 +&
     & (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 + (-sig&
     &(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(6) ** 2&
     & / 0.3D1) + (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) * sig(1) * sin(theta) * sig(6)) * (sig(1) ** 2 * cos(theta) ** &
     &2 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + (sig(4) * c&
     &os(theta) + sig(5) * sin(theta)) ** 2) * mu1 * (-0.4D1 / 0.3D1 * (&
     &sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) - 0.2D1 / &
     &0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta))
      J(1,5) = 0.2D1 * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + si&
     &g(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) * sig(1) * sin(theta) * sig(6)) * (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) * sin(theta) * mu1 * (-0.2D1 / 0.3D1 * sig(1) &
     &** 2 * cos(theta) ** 2 - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) ** 2 - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) ** 2 + sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 +&
     & (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 + (-sig&
     &(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(6) ** 2&
     & / 0.3D1) + (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) * sig(1) * sin(theta) * sig(6)) * (sig(1) ** 2 * cos(theta) ** &
     &2 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + (sig(4) * c&
     &os(theta) + sig(5) * sin(theta)) ** 2) * mu1 * (-0.4D1 / 0.3D1 * (&
     &sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) + 0.2D1 / &
     &0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * cos(theta))
      J(1,6) = (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) &
     &* sin(theta)) * (sig(1) ** 2 * cos(theta) ** 2 + (sig(2) * cos(the&
     &ta) + sig(3) * sin(theta)) ** 2 + (sig(4) * cos(theta) + sig(5) * &
     &sin(theta)) ** 2) * mu1 * (-0.2D1 / 0.3D1 * sig(1) ** 2 * cos(thet&
     &a) ** 2 - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(thet&
     &a)) ** 2 - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(the&
     &ta)) ** 2 + sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 + (-sig(2) * sin&
     &(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 + (-sig(4) * sin(theta&
     &) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(6) ** 2 / 0.3D1) + 0.2&
     &D1 / 0.3D1 * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3)&
     & * cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(thet&
     &a)) * sig(1) * sin(theta) * sig(6)) * (sig(1) ** 2 * cos(theta) **&
     & 2 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + (sig(4) * &
     &cos(theta) + sig(5) * sin(theta)) ** 2) * mu1 * sig(6) + 0.2D1 * (&
     &lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * (-sig(2) * sin(&
     &theta) + sig(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) * sig(1) * sin(theta) * sig(6)) * (-sig(1) * co&
     &s(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) +&
     & 0.1D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * si&
     &n(theta) * sig(6)) * (sig(1) * cos(theta) * (-sig(2) * sin(theta) &
     &+ sig(3) * cos(theta)) + (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) * sig(1) * sin(theta)) + (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) *&
     & cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6&
     &) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * sin(the&
     &ta) * sig(6)) ** 2 * (-sig(1) * cos(theta) * (-sig(2) * sin(theta)&
     & + sig(3) * cos(theta)) - (sig(2) * cos(theta) + sig(3) * sin(thet&
     &a)) * sig(1) * sin(theta))
      J(2,1) = (cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)&
     &) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sin(the&
     &ta) * sig(6)) * (sig(1) ** 2 * sin(theta) ** 2 + (-sig(2) * sin(th&
     &eta) + sig(3) * cos(theta)) ** 2 + (-sig(4) * sin(theta) + sig(5) &
     &* cos(theta)) ** 2) * mu1 * (sig(1) ** 2 * cos(theta) ** 2 / 0.3D1&
     & + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 + (sig&
     &(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - 0.2D1 / 0.3&
     &D1 * sig(1) ** 2 * sin(theta) ** 2 - 0.2D1 / 0.3D1 * (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta)) ** 2 - 0.2D1 / 0.3D1 * (-sig(4) * &
     &sin(theta) + sig(5) * cos(theta)) ** 2 + sig(6) ** 2 / 0.3D1) + 0.&
     &2D1 * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(&
     &theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * s&
     &ig(1) * sin(theta) * sig(6)) * sig(1) * sin(theta) ** 2 * mu1 * (s&
     &ig(1) ** 2 * cos(theta) ** 2 / 0.3D1 + (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) ** 2 / 0.3D1 + (sig(4) * cos(theta) + sig(5) * si&
     &n(theta)) ** 2 / 0.3D1 - 0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) &
     &** 2 - 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)&
     &) ** 2 - 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(thet&
     &a)) ** 2 + sig(6) ** 2 / 0.3D1) + (sig(1) * cos(theta) * (-sig(2) &
     &* sin(theta) + sig(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta&
     &) + sig(3) * sin(theta)) * sig(1) * sin(theta) * sig(6)) * (sig(1)&
     & ** 2 * sin(theta) ** 2 + (-sig(2) * sin(theta) + sig(3) * cos(the&
     &ta)) ** 2 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2) * m&
     &u1 * (0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2 - 0.4D1 / 0.3D1 * s&
     &ig(1) * sin(theta) ** 2) + 0.2D1 * (lambda1 + 0.2D1 / 0.3D1 * mu1) * (s&
     &ig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) &
     &* sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * &
     &sin(theta) * sig(6)) * (-sig(1) * cos(theta) * (-sig(2) * sin(thet&
     &a) + sig(3) * cos(theta)) * sig(6) + 0.1D1 - (sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) * sig(1) * sin(theta) * sig(6)) * (cos(thet&
     &a) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + (sig(&
     &2) * cos(theta) + sig(3) * sin(theta)) * sin(theta) * sig(6)) + (lambda1&
     & + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * (-sig(2) * sin(t&
     &heta) + sig(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig&
     &(3) * sin(theta)) * sig(1) * sin(theta) * sig(6)) ** 2 * (-cos(the&
     &ta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) - (sig&
     &(2) * cos(theta) + sig(3) * sin(theta)) * sin(theta) * sig(6))
      J(2,2) = -0.2D1 * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + s&
     &ig(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin&
     &(theta)) * sig(1) * sin(theta) * sig(6)) * (-sig(2) * sin(theta) +&
     & sig(3) * cos(theta)) * sin(theta) * mu1 * (sig(1) ** 2 * cos(thet&
     &a) ** 2 / 0.3D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2&
     & / 0.3D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D&
     &1 - 0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 - 0.2D1 / 0.3D1 &
     &* (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 - 0.2D1 / 0.3D&
     &1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 + sig(6) ** &
     &2 / 0.3D1) + (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3)&
     & * cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(thet&
     &a)) * sig(1) * sin(theta) * sig(6)) * (sig(1) ** 2 * sin(theta) **&
     & 2 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + (-sig(4) &
     &* sin(theta) + sig(5) * cos(theta)) ** 2) * mu1 * (0.2D1 / 0.3D1 *&
     & (sig(2) * cos(theta) + sig(3) * sin(theta)) * cos(theta) + 0.4D1 &
     &/ 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(theta&
     &))
      J(2,3) = (sig(1) * cos(theta) ** 2 * sig(6) + sin(theta) ** 2 * si&
     &g(1) * sig(6)) * (sig(1) ** 2 * sin(theta) ** 2 + (-sig(2) * sin(t&
     &heta) + sig(3) * cos(theta)) ** 2 + (-sig(4) * sin(theta) + sig(5)&
     & * cos(theta)) ** 2) * mu1 * (sig(1) ** 2 * cos(theta) ** 2 / 0.3D&
     &1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 + (si&
     &g(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - 0.2D1 / 0.&
     &3D1 * sig(1) ** 2 * sin(theta) ** 2 - 0.2D1 / 0.3D1 * (-sig(2) * s&
     &in(theta) + sig(3) * cos(theta)) ** 2 - 0.2D1 / 0.3D1 * (-sig(4) *&
     & sin(theta) + sig(5) * cos(theta)) ** 2 + sig(6) ** 2 / 0.3D1) + 0&
     &.2D1 * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos&
     &(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * &
     &sig(1) * sin(theta) * sig(6)) * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) * cos(theta) * mu1 * (sig(1) ** 2 * cos(theta) ** 2 / 0&
     &.3D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 + &
     &(sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - 0.2D1 /&
     & 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 - 0.2D1 / 0.3D1 * (-sig(2) &
     &* sin(theta) + sig(3) * cos(theta)) ** 2 - 0.2D1 / 0.3D1 * (-sig(4&
     &) * sin(theta) + sig(5) * cos(theta)) ** 2 + sig(6) ** 2 / 0.3D1) &
     &+ (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(thet&
     &a)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1&
     &) * sin(theta) * sig(6)) * (sig(1) ** 2 * sin(theta) ** 2 + (-sig(&
     &2) * sin(theta) + sig(3) * cos(theta)) ** 2 + (-sig(4) * sin(theta&
     &) + sig(5) * cos(theta)) ** 2) * mu1 * (0.2D1 / 0.3D1 * (sig(2) * &
     &cos(theta) + sig(3) * sin(theta)) * sin(theta) - 0.4D1 / 0.3D1 * (&
     &-sig(2) * sin(theta) + sig(3) * cos(theta)) * cos(theta)) + 0.2D1 &
     &* (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * (-sig(2) * s&
     &in(theta) + sig(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta) +&
     & sig(3) * sin(theta)) * sig(1) * sin(theta) * sig(6)) * (-sig(1) *&
     & cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6&
     &) + 0.1D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) *&
     & sin(theta) * sig(6)) * (sig(1) * cos(theta) ** 2 * sig(6) + sin(t&
     &heta) ** 2 * sig(1) * sig(6)) + (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(&
     &1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * s&
     &ig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * sin&
     &(theta) * sig(6)) ** 2 * (-sig(1) * cos(theta) ** 2 * sig(6) - sin&
     &(theta) ** 2 * sig(1) * sig(6))
      J(2,4) = -0.2D1 * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + s&
     &ig(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin&
     &(theta)) * sig(1) * sin(theta) * sig(6)) * (-sig(4) * sin(theta) +&
     & sig(5) * cos(theta)) * sin(theta) * mu1 * (sig(1) ** 2 * cos(thet&
     &a) ** 2 / 0.3D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2&
     & / 0.3D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D&
     &1 - 0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 - 0.2D1 / 0.3D1 &
     &* (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 - 0.2D1 / 0.3D&
     &1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 + sig(6) ** &
     &2 / 0.3D1) + (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3)&
     & * cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(thet&
     &a)) * sig(1) * sin(theta) * sig(6)) * (sig(1) ** 2 * sin(theta) **&
     & 2 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + (-sig(4) &
     &* sin(theta) + sig(5) * cos(theta)) ** 2) * mu1 * (0.2D1 / 0.3D1 *&
     & (sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) + 0.4D1 &
     &/ 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta&
     &))
      J(2,5) = 0.2D1 * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + si&
     &g(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) * sig(1) * sin(theta) * sig(6)) * (-sig(4) * sin(theta) + &
     &sig(5) * cos(theta)) * cos(theta) * mu1 * (sig(1) ** 2 * cos(theta&
     &) ** 2 / 0.3D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 &
     &/ 0.3D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1&
     & - 0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 - 0.2D1 / 0.3D1 *&
     & (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 - 0.2D1 / 0.3D1&
     & * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 + sig(6) ** 2&
     & / 0.3D1) + (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) * sig(1) * sin(theta) * sig(6)) * (sig(1) ** 2 * sin(theta) ** &
     &2 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + (-sig(4) *&
     & sin(theta) + sig(5) * cos(theta)) ** 2) * mu1 * (0.2D1 / 0.3D1 * &
     &(sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) - 0.4D1 /&
     & 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * cos(theta)&
     &)
      J(2,6) = (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) &
     &* sin(theta)) * (sig(1) ** 2 * sin(theta) ** 2 + (-sig(2) * sin(th&
     &eta) + sig(3) * cos(theta)) ** 2 + (-sig(4) * sin(theta) + sig(5) &
     &* cos(theta)) ** 2) * mu1 * (sig(1) ** 2 * cos(theta) ** 2 / 0.3D1&
     & + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 + (sig&
     &(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - 0.2D1 / 0.3&
     &D1 * sig(1) ** 2 * sin(theta) ** 2 - 0.2D1 / 0.3D1 * (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta)) ** 2 - 0.2D1 / 0.3D1 * (-sig(4) * &
     &sin(theta) + sig(5) * cos(theta)) ** 2 + sig(6) ** 2 / 0.3D1) + 0.&
     &2D1 / 0.3D1 * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3&
     &) * cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(the&
     &ta)) * sig(1) * sin(theta) * sig(6)) * (sig(1) ** 2 * sin(theta) *&
     &* 2 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + (-sig(4)&
     & * sin(theta) + sig(5) * cos(theta)) ** 2) * mu1 * sig(6) + 0.2D1 &
     &* (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * (-sig(2) * s&
     &in(theta) + sig(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta) +&
     & sig(3) * sin(theta)) * sig(1) * sin(theta) * sig(6)) * (-sig(1) *&
     & cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6&
     &) + 0.1D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) *&
     & sin(theta) * sig(6)) * (sig(1) * cos(theta) * (-sig(2) * sin(thet&
     &a) + sig(3) * cos(theta)) + (sig(2) * cos(theta) + sig(3) * sin(th&
     &eta)) * sig(1) * sin(theta)) + (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1&
     &) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * si&
     &g(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * sin(&
     &theta) * sig(6)) ** 2 * (-sig(1) * cos(theta) * (-sig(2) * sin(the&
     &ta) + sig(3) * cos(theta)) - (sig(2) * cos(theta) + sig(3) * sin(t&
     &heta)) * sig(1) * sin(theta))
      J(3,1) = (cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)&
     &) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sin(the&
     &ta) * sig(6)) * sig(6) ** 2 * mu1 * (sig(1) ** 2 * cos(theta) ** 2&
     & / 0.3D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D&
     &1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig&
     &(1) ** 2 * sin(theta) ** 2 / 0.3D1 + (-sig(2) * sin(theta) + sig(3&
     &) * cos(theta)) ** 2 / 0.3D1 + (-sig(4) * sin(theta) + sig(5) * co&
     &s(theta)) ** 2 / 0.3D1 - 0.2D1 / 0.3D1 * sig(6) ** 2) + (sig(1) * &
     &cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6)&
     & + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * sin(thet&
     &a) * sig(6)) * sig(6) ** 2 * mu1 * (0.2D1 / 0.3D1 * sig(1) * cos(t&
     &heta) ** 2 + 0.2D1 / 0.3D1 * sig(1) * sin(theta) ** 2) + 0.2D1 * (&
     &lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * (-sig(2) * sin(&
     &theta) + sig(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) * sig(1) * sin(theta) * sig(6)) * (-sig(1) * co&
     &s(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) +&
     & 0.1D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * si&
     &n(theta) * sig(6)) * (cos(theta) * (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) * sin(theta) * sig(6)) + (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) *&
     & cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6&
     &) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * sin(the&
     &ta) * sig(6)) ** 2 * (-cos(theta) * (-sig(2) * sin(theta) + sig(3)&
     & * cos(theta)) * sig(6) - (sig(2) * cos(theta) + sig(3) * sin(thet&
     &a)) * sin(theta) * sig(6))
      J(3,2) = (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) &
     &* sig(1) * sin(theta) * sig(6)) * sig(6) ** 2 * mu1 * (0.2D1 / 0.3&
     &D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) * cos(theta) - 0.&
     &2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(t&
     &heta))
      J(3,3) = (sig(1) * cos(theta) ** 2 * sig(6) + sin(theta) ** 2 * si&
     &g(1) * sig(6)) * sig(6) ** 2 * mu1 * (sig(1) ** 2 * cos(theta) ** &
     &2 / 0.3D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3&
     &D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + si&
     &g(1) ** 2 * sin(theta) ** 2 / 0.3D1 + (-sig(2) * sin(theta) + sig(&
     &3) * cos(theta)) ** 2 / 0.3D1 + (-sig(4) * sin(theta) + sig(5) * c&
     &os(theta)) ** 2 / 0.3D1 - 0.2D1 / 0.3D1 * sig(6) ** 2) + (sig(1) *&
     & cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6&
     &) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * sin(the&
     &ta) * sig(6)) * sig(6) ** 2 * mu1 * (0.2D1 / 0.3D1 * (sig(2) * cos&
     &(theta) + sig(3) * sin(theta)) * sin(theta) + 0.2D1 / 0.3D1 * (-si&
     &g(2) * sin(theta) + sig(3) * cos(theta)) * cos(theta)) + 0.2D1 * (&
     &lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * (-sig(2) * sin(&
     &theta) + sig(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) * sig(1) * sin(theta) * sig(6)) * (-sig(1) * co&
     &s(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) +&
     & 0.1D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * si&
     &n(theta) * sig(6)) * (sig(1) * cos(theta) ** 2 * sig(6) + sin(thet&
     &a) ** 2 * sig(1) * sig(6)) + (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) &
     &* cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(&
     &6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * sin(th&
     &eta) * sig(6)) ** 2 * (-sig(1) * cos(theta) ** 2 * sig(6) - sin(th&
     &eta) ** 2 * sig(1) * sig(6))
      J(3,4) = (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) &
     &* sig(1) * sin(theta) * sig(6)) * sig(6) ** 2 * mu1 * (0.2D1 / 0.3&
     &D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) - 0.&
     &2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(t&
     &heta))
      J(3,5) = (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) &
     &* sig(1) * sin(theta) * sig(6)) * sig(6) ** 2 * mu1 * (0.2D1 / 0.3&
     &D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) + 0.&
     &2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * cos(t&
     &heta))
      J(3,6) = (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) &
     &* sin(theta)) * sig(6) ** 2 * mu1 * (sig(1) ** 2 * cos(theta) ** 2&
     & / 0.3D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D&
     &1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig&
     &(1) ** 2 * sin(theta) ** 2 / 0.3D1 + (-sig(2) * sin(theta) + sig(3&
     &) * cos(theta)) ** 2 / 0.3D1 + (-sig(4) * sin(theta) + sig(5) * co&
     &s(theta)) ** 2 / 0.3D1 - 0.2D1 / 0.3D1 * sig(6) ** 2) + 0.2D1 * (s&
     &ig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) &
     &* sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * &
     &sin(theta) * sig(6)) * sig(6) * mu1 * (sig(1) ** 2 * cos(theta) **&
     & 2 / 0.3D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.&
     &3D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + s&
     &ig(1) ** 2 * sin(theta) ** 2 / 0.3D1 + (-sig(2) * sin(theta) + sig&
     &(3) * cos(theta)) ** 2 / 0.3D1 + (-sig(4) * sin(theta) + sig(5) * &
     &cos(theta)) ** 2 / 0.3D1 - 0.2D1 / 0.3D1 * sig(6) ** 2) - 0.4D1 / &
     &0.3D1 * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * co&
     &s(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) *&
     & sig(1) * sin(theta) * sig(6)) * sig(6) ** 3 * mu1 + 0.2D1 * (lambda1 +&
     & 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * (-sig(2) * sin(thet&
     &a) + sig(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3)&
     & * sin(theta)) * sig(1) * sin(theta) * sig(6)) * (-sig(1) * cos(th&
     &eta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + 0.1&
     &D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * sin(th&
     &eta) * sig(6)) * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + si&
     &g(3) * cos(theta)) + (sig(2) * cos(theta) + sig(3) * sin(theta)) *&
     & sig(1) * sin(theta)) + (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos&
     &(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + &
     &(sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * sin(theta) &
     &* sig(6)) ** 2 * (-sig(1) * cos(theta) * (-sig(2) * sin(theta) + s&
     &ig(3) * cos(theta)) - (sig(2) * cos(theta) + sig(3) * sin(theta)) &
     &* sig(1) * sin(theta))
      J(4,1) = 0.2D1 * (cos(theta) * (-sig(2) * sin(theta) + sig(3) * co&
     &s(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) *&
     & sin(theta) * sig(6)) * (-sig(1) ** 2 * cos(theta) * sin(theta) + &
     &(sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta&
     &) + sig(3) * cos(theta)) + (sig(4) * cos(theta) + sig(5) * sin(the&
     &ta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) + 0.1D1) * mu1&
     & * (sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D1 - (sig(2) * cos(&
     &theta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) / 0.2D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) *&
     & (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.2D1) - 0.4D1 * (&
     &sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta))&
     & * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) *&
     & sin(theta) * sig(6)) * sig(1) * cos(theta) * sin(theta) * mu1 * (&
     &sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D1 - (sig(2) * cos(thet&
     &a) + sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(t&
     &heta)) / 0.2D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-s&
     &ig(4) * sin(theta) + sig(5) * cos(theta)) / 0.2D1) + 0.2D1 * (sig(&
     &1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * s&
     &ig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * sin&
     &(theta) * sig(6)) * (-sig(1) ** 2 * cos(theta) * sin(theta) + (sig&
     &(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta) + &
     &sig(3) * cos(theta)) + (sig(4) * cos(theta) + sig(5) * sin(theta))&
     & * (-sig(4) * sin(theta) + sig(5) * cos(theta)) + 0.1D1) * mu1 * s&
     &ig(1) * cos(theta) * sin(theta)
      J(4,2) = 0.2D1 * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + si&
     &g(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) * sig(1) * sin(theta) * sig(6)) * ((-sig(2) * sin(theta) +&
     & sig(3) * cos(theta)) * cos(theta) - (sig(2) * cos(theta) + sig(3)&
     & * sin(theta)) * sin(theta)) * mu1 * (sig(1) ** 2 * cos(theta) * s&
     &in(theta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.2D1 - (sig(4) * c&
     &os(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) &
     &* cos(theta)) / 0.2D1) + 0.2D1 * (sig(1) * cos(theta) * (-sig(2) *&
     & sin(theta) + sig(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta)&
     & + sig(3) * sin(theta)) * sig(1) * sin(theta) * sig(6)) * (-sig(1)&
     & ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta) + sig(3) * &
     &sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + (sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta) + s&
     &ig(5) * cos(theta)) + 0.1D1) * mu1 * (-(-sig(2) * sin(theta) + sig&
     &(3) * cos(theta)) * cos(theta) / 0.2D1 + (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) * sin(theta) / 0.2D1)
      J(4,3) = 0.2D1 * (sig(1) * cos(theta) ** 2 * sig(6) + sin(theta) *&
     &* 2 * sig(1) * sig(6)) * (-sig(1) ** 2 * cos(theta) * sin(theta) +&
     & (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(thet&
     &a) + sig(3) * cos(theta)) + (sig(4) * cos(theta) + sig(5) * sin(th&
     &eta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) + 0.1D1) * mu&
     &1 * (sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D1 - (sig(2) * cos&
     &(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * &
     &cos(theta)) / 0.2D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) &
     &* (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.2D1) + 0.2D1 * &
     &(sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)&
     &) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) &
     &* sin(theta) * sig(6)) * ((-sig(2) * sin(theta) + sig(3) * cos(the&
     &ta)) * sin(theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * &
     &cos(theta)) * mu1 * (sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D1&
     & - (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(th&
     &eta) + sig(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0&
     &.2D1) + 0.2D1 * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig&
     &(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(t&
     &heta)) * sig(1) * sin(theta) * sig(6)) * (-sig(1) ** 2 * cos(theta&
     &) * sin(theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-s&
     &ig(2) * sin(theta) + sig(3) * cos(theta)) + (sig(4) * cos(theta) +&
     & sig(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta&
     &)) + 0.1D1) * mu1 * (-(-sig(2) * sin(theta) + sig(3) * cos(theta))&
     & * sin(theta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)&
     &) * cos(theta) / 0.2D1)
      J(4,4) = 0.2D1 * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + si&
     &g(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) * sig(1) * sin(theta) * sig(6)) * ((-sig(4) * sin(theta) +&
     & sig(5) * cos(theta)) * cos(theta) - (sig(4) * cos(theta) + sig(5)&
     & * sin(theta)) * sin(theta)) * mu1 * (sig(1) ** 2 * cos(theta) * s&
     &in(theta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.2D1 - (sig(4) * c&
     &os(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) &
     &* cos(theta)) / 0.2D1) + 0.2D1 * (sig(1) * cos(theta) * (-sig(2) *&
     & sin(theta) + sig(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta)&
     & + sig(3) * sin(theta)) * sig(1) * sin(theta) * sig(6)) * (-sig(1)&
     & ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta) + sig(3) * &
     &sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + (sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta) + s&
     &ig(5) * cos(theta)) + 0.1D1) * mu1 * (-(-sig(4) * sin(theta) + sig&
     &(5) * cos(theta)) * cos(theta) / 0.2D1 + (sig(4) * cos(theta) + si&
     &g(5) * sin(theta)) * sin(theta) / 0.2D1)
      J(4,5) = 0.2D1 * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + si&
     &g(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) * sig(1) * sin(theta) * sig(6)) * ((-sig(4) * sin(theta) +&
     & sig(5) * cos(theta)) * sin(theta) + (sig(4) * cos(theta) + sig(5)&
     & * sin(theta)) * cos(theta)) * mu1 * (sig(1) ** 2 * cos(theta) * s&
     &in(theta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.2D1 - (sig(4) * c&
     &os(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) &
     &* cos(theta)) / 0.2D1) + 0.2D1 * (sig(1) * cos(theta) * (-sig(2) *&
     & sin(theta) + sig(3) * cos(theta)) * sig(6) + (sig(2) * cos(theta)&
     & + sig(3) * sin(theta)) * sig(1) * sin(theta) * sig(6)) * (-sig(1)&
     & ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta) + sig(3) * &
     &sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + (sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta) + s&
     &ig(5) * cos(theta)) + 0.1D1) * mu1 * (-(-sig(4) * sin(theta) + sig&
     &(5) * cos(theta)) * sin(theta) / 0.2D1 - (sig(4) * cos(theta) + si&
     &g(5) * sin(theta)) * cos(theta) / 0.2D1)
      J(4,6) = 0.2D1 * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + si&
     &g(3) * cos(theta)) + (sig(2) * cos(theta) + sig(3) * sin(theta)) *&
     & sig(1) * sin(theta)) * (-sig(1) ** 2 * cos(theta) * sin(theta) + &
     &(sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta&
     &) + sig(3) * cos(theta)) + (sig(4) * cos(theta) + sig(5) * sin(the&
     &ta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) + 0.1D1) * mu1&
     & * (sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D1 - (sig(2) * cos(&
     &theta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) / 0.2D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) *&
     & (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.2D1)
      J(5,1) = -(cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta&
     &)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sin(th&
     &eta) * sig(6)) * (0.1D1 + (-sig(4) * sin(theta) + sig(5) * cos(the&
     &ta)) * sig(6)) * mu1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)&
     &) * sig(6)
      J(5,3) = -(sig(1) * cos(theta) ** 2 * sig(6) + sin(theta) ** 2 * s&
     &ig(1) * sig(6)) * (0.1D1 + (-sig(4) * sin(theta) + sig(5) * cos(th&
     &eta)) * sig(6)) * mu1 * (-sig(4) * sin(theta) + sig(5) * cos(theta&
     &)) * sig(6)
      J(5,4) = (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) &
     &* sig(1) * sin(theta) * sig(6)) * sin(theta) * sig(6) ** 2 * mu1 *&
     & (-sig(4) * sin(theta) + sig(5) * cos(theta)) + (sig(1) * cos(thet&
     &a) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + (sig(&
     &2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * sin(theta) * sig&
     &(6)) * (0.1D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig&
     &(6)) * mu1 * sin(theta) * sig(6)
      J(5,5) = -(sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * &
     &cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta))&
     & * sig(1) * sin(theta) * sig(6)) * cos(theta) * sig(6) ** 2 * mu1 &
     &* (-sig(4) * sin(theta) + sig(5) * cos(theta)) - (sig(1) * cos(the&
     &ta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + (sig&
     &(2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * sin(theta) * si&
     &g(6)) * (0.1D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) * si&
     &g(6)) * mu1 * cos(theta) * sig(6)
      J(5,6) = -(sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * &
     &cos(theta)) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1)&
     & * sin(theta)) * (0.1D1 + (-sig(4) * sin(theta) + sig(5) * cos(the&
     &ta)) * sig(6)) * mu1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)&
     &) * sig(6) - (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3)&
     & * cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(thet&
     &a)) * sig(1) * sin(theta) * sig(6)) * (-sig(4) * sin(theta) + sig(&
     &5) * cos(theta)) ** 2 * mu1 * sig(6) - (sig(1) * cos(theta) * (-si&
     &g(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + (sig(2) * cos(&
     &theta) + sig(3) * sin(theta)) * sig(1) * sin(theta) * sig(6)) * (0&
     &.1D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig(6)) * mu&
     &1 * (-sig(4) * sin(theta) + sig(5) * cos(theta))
      J(6,1) = -(cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta&
     &)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sin(th&
     &eta) * sig(6)) * (0.1D1 + (sig(4) * cos(theta) + sig(5) * sin(thet&
     &a)) * sig(6)) * mu1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) &
     &* sig(6)
      J(6,3) = -(sig(1) * cos(theta) ** 2 * sig(6) + sin(theta) ** 2 * s&
     &ig(1) * sig(6)) * (0.1D1 + (sig(4) * cos(theta) + sig(5) * sin(the&
     &ta)) * sig(6)) * mu1 * (sig(4) * cos(theta) + sig(5) * sin(theta))&
     & * sig(6)
      J(6,4) = -(sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * &
     &cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta))&
     & * sig(1) * sin(theta) * sig(6)) * cos(theta) * sig(6) ** 2 * mu1 &
     &* (sig(4) * cos(theta) + sig(5) * sin(theta)) - (sig(1) * cos(thet&
     &a) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + (sig(&
     &2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * sin(theta) * sig&
     &(6)) * (0.1D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) * sig(&
     &6)) * mu1 * cos(theta) * sig(6)
      J(6,5) = -(sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * &
     &cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta))&
     & * sig(1) * sin(theta) * sig(6)) * sin(theta) * sig(6) ** 2 * mu1 &
     &* (sig(4) * cos(theta) + sig(5) * sin(theta)) - (sig(1) * cos(thet&
     &a) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + (sig(&
     &2) * cos(theta) + sig(3) * sin(theta)) * sig(1) * sin(theta) * sig&
     &(6)) * (0.1D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) * sig(&
     &6)) * mu1 * sin(theta) * sig(6)
      J(6,6) = -(sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * &
     &cos(theta)) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(1)&
     & * sin(theta)) * (0.1D1 + (sig(4) * cos(theta) + sig(5) * sin(thet&
     &a)) * sig(6)) * mu1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) &
     &* sig(6) - (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) *&
     & cos(theta)) * sig(6) + (sig(2) * cos(theta) + sig(3) * sin(theta)&
     &) * sig(1) * sin(theta) * sig(6)) * (sig(4) * cos(theta) + sig(5) &
     &* sin(theta)) ** 2 * mu1 * sig(6) - (sig(1) * cos(theta) * (-sig(2&
     &) * sin(theta) + sig(3) * cos(theta)) * sig(6) + (sig(2) * cos(the&
     &ta) + sig(3) * sin(theta)) * sig(1) * sin(theta) * sig(6)) * (0.1D&
     &1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) * sig(6)) * mu1 * &
     &(sig(4) * cos(theta) + sig(5) * sin(theta))
            case(2)
 J(1,1) = -rho0 * (cos(theta) * mu1 * (-0.4D1 / 0.3D1 * (-sig(1) **&
     & 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin&
     &(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)&
     &) ** 2 / 0.3D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) &
     &* sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(&
     &theta) + sig(5) * cos(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) &
     &* sig(1) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta)&
     & ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 &
     &/ 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D&
     &1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) +&
     & sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5)&
     & * sin(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * sig(1) * cos(&
     &theta) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + sig(1) ** 2 * cos&
     &(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta))&
     & ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 /&
     & 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(th&
     &eta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) +&
     & sig(5) * cos(theta)) ** 2 / 0.6D1) * sig(1) * cos(theta) + 0.2D1 &
     &* (sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D1 - (sig(2) * cos(t&
     &heta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * co&
     &s(theta)) / 0.2D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) * &
     &(-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.2D1) * sig(1) * s&
     &in(theta)) / rho0 + sig(1) * cos(theta) * mu1 * (-0.4D1 / 0.3D1 * &
     &(-0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2 + sig(1) * sin(theta) *&
     &* 2 / 0.3D1) * sig(1) * cos(theta) - 0.4D1 / 0.3D1 * (-sig(1) ** 2&
     & * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(t&
     &heta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) &
     &** 2 / 0.3D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * &
     &sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(th&
     &eta) + sig(5) * cos(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * &
     &cos(theta) + 0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * sig(1) * sin(theta)&
     & ** 2 + sig(1) * cos(theta) ** 2 / 0.3D1) * sig(1) * cos(theta) + &
     &0.2D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2)&
     & * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin&
     &(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * cos(th&
     &eta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) **&
     & 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.&
     &6D1 + sig(6) ** 2 / 0.6D1) * cos(theta) + 0.2D1 / 0.3D1 * (sig(1) &
     &* cos(theta) ** 2 / 0.3D1 + sig(1) * sin(theta) ** 2 / 0.3D1) * si&
     &g(1) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + sig(1)&
     & ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * &
     &sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(the&
     &ta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(&
     &2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * s&
     &in(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * cos(theta) + 0.2D&
     &1 * sig(1) ** 2 * cos(theta) * sin(theta) ** 2 + 0.2D1 * (sig(1) *&
     &* 2 * cos(theta) * sin(theta) / 0.2D1 - (sig(2) * cos(theta) + sig&
     &(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) /&
     & 0.2D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * &
     &sin(theta) + sig(5) * cos(theta)) / 0.2D1) * sin(theta)) / rho0 + &
     &(sig(2) * cos(theta) + sig(3) * sin(theta)) * mu1 * (-0.4D1 / 0.3D&
     &1 * (-0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2 + sig(1) * sin(thet&
     &a) ** 2 / 0.3D1) * (sig(2) * cos(theta) + sig(3) * sin(theta)) + 0&
     &.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * sig(1) * sin(theta) ** 2 + sig(1)&
     & * cos(theta) ** 2 / 0.3D1) * (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) + 0.2D1 / 0.3D1 * (sig(1) * cos(theta) ** 2 / 0.3D1 + sig(&
     &1) * sin(theta) ** 2 / 0.3D1) * (sig(2) * cos(theta) + sig(3) * si&
     &n(theta)) - 0.2D1 * sig(1) * cos(theta) * sin(theta) * (-sig(2) * &
     &sin(theta) + sig(3) * cos(theta))) / rho0 + (sig(4) * cos(theta) +&
     & sig(5) * sin(theta)) * mu1 * (-0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * &
     &sig(1) * cos(theta) ** 2 + sig(1) * sin(theta) ** 2 / 0.3D1) * (si&
     &g(4) * cos(theta) + sig(5) * sin(theta)) + 0.2D1 / 0.3D1 * (-0.2D1&
     & / 0.3D1 * sig(1) * sin(theta) ** 2 + sig(1) * cos(theta) ** 2 / 0&
     &.3D1) * (sig(4) * cos(theta) + sig(5) * sin(theta)) + 0.2D1 / 0.3D&
     &1 * (sig(1) * cos(theta) ** 2 / 0.3D1 + sig(1) * sin(theta) ** 2 /&
     & 0.3D1) * (sig(4) * cos(theta) + sig(5) * sin(theta)) - 0.2D1 * si&
     &g(1) * cos(theta) * sin(theta) * (-sig(4) * sin(theta) + sig(5) * &
     &cos(theta))) / rho0) + 0.2D1 * (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1&
     &) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * si&
     &g(6) + sig(1) * sin(theta) * (sig(2) * cos(theta) + sig(3) * sin(t&
     &heta)) * sig(6)) * (-sig(1) * cos(theta) * (-sig(2) * sin(theta) +&
     & sig(3) * cos(theta)) * sig(6) + 0.1D1 - sig(1) * sin(theta) * (si&
     &g(2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) * (cos(theta) *&
     & (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + sin(theta&
     &) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) + (lambda1 + &
     &0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * (-sig(2) * sin(theta&
     &) + sig(3) * cos(theta)) * sig(6) + sig(1) * sin(theta) * (sig(2) &
     &* cos(theta) + sig(3) * sin(theta)) * sig(6)) ** 2 * (-cos(theta) &
     &* (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) - sin(thet&
     &a) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6))
      J(1,2) = -rho0 * (sig(1) * cos(theta) * mu1 * (-0.4D1 / 0.3D1 * (-&
     &0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) * cos(&
     &theta) - (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(theta)&
     & / 0.3D1) * sig(1) * cos(theta) + 0.2D1 / 0.3D1 * (0.2D1 / 0.3D1 *&
     & (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(theta) + (sig(&
     &2) * cos(theta) + sig(3) * sin(theta)) * cos(theta) / 0.3D1) * sig&
     &(1) * cos(theta) + 0.2D1 / 0.3D1 * ((sig(2) * cos(theta) + sig(3) &
     &* sin(theta)) * cos(theta) / 0.3D1 - (-sig(2) * sin(theta) + sig(3&
     &) * cos(theta)) * sin(theta) / 0.3D1) * sig(1) * cos(theta) + 0.2D&
     &1 * (-cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) / &
     &0.2D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sin(theta) /&
     & 0.2D1) * sig(1) * sin(theta)) / rho0 + cos(theta) * mu1 * (-0.4D1&
     & / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos&
     &(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta)&
     & + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * sin(theta) **&
     & 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0&
     &.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1 +&
     & sig(6) ** 2 / 0.6D1) * (sig(2) * cos(theta) + sig(3) * sin(theta)&
     &) + 0.2D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-si&
     &g(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) *&
     & sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * co&
     &s(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)&
     &) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 &
     &/ 0.6D1 + sig(6) ** 2 / 0.6D1) * (sig(2) * cos(theta) + sig(3) * s&
     &in(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + sig(1) ** 2 *&
     & cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(the&
     &ta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) **&
     & 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(thet&
     &a) + sig(5) * cos(theta)) ** 2 / 0.6D1) * (sig(2) * cos(theta) + s&
     &ig(3) * sin(theta)) - 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(thet&
     &a) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2&
     &) * sin(theta) + sig(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(thet&
     &a) + sig(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(t&
     &heta)) / 0.2D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta))) / &
     &rho0 + (sig(2) * cos(theta) + sig(3) * sin(theta)) * mu1 * (-0.4D1&
     & / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(t&
     &heta)) * cos(theta) - (-sig(2) * sin(theta) + sig(3) * cos(theta))&
     & * sin(theta) / 0.3D1) * (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) - 0.4D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (si&
     &g(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * &
     &cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * sin&
     &(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)&
     &) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2&
     & / 0.6D1 + sig(6) ** 2 / 0.6D1) * cos(theta) + 0.2D1 / 0.3D1 * (0.&
     &2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(t&
     &heta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * cos(theta) /&
     & 0.3D1) * (sig(2) * cos(theta) + sig(3) * sin(theta)) + 0.2D1 / 0.&
     &3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(the&
     &ta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + &
     &sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 &
     &/ 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.6D1&
     & + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(&
     &6) ** 2 / 0.6D1) * cos(theta) + 0.2D1 / 0.3D1 * ((sig(2) * cos(the&
     &ta) + sig(3) * sin(theta)) * cos(theta) / 0.3D1 - (-sig(2) * sin(t&
     &heta) + sig(3) * cos(theta)) * sin(theta) / 0.3D1) * (sig(2) * cos&
     &(theta) + sig(3) * sin(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0&
     &.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta&
     &) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig&
     &(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0&
     &.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 +&
     & (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * cos(&
     &theta) - 0.2D1 * (-cos(theta) * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) / 0.2D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) *&
     & sin(theta) / 0.2D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta)&
     &) + 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D1 - (sig(&
     &2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta) + s&
     &ig(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(theta) + sig(5) * sin(&
     &theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.2D1) * &
     &sin(theta)) / rho0 + (sig(4) * cos(theta) + sig(5) * sin(theta)) *&
     & mu1 * (-0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + &
     &sig(3) * sin(theta)) * cos(theta) - (-sig(2) * sin(theta) + sig(3)&
     & * cos(theta)) * sin(theta) / 0.3D1) * (sig(4) * cos(theta) + sig(&
     &5) * sin(theta)) + 0.2D1 / 0.3D1 * (0.2D1 / 0.3D1 * (-sig(2) * sin&
     &(theta) + sig(3) * cos(theta)) * sin(theta) + (sig(2) * cos(theta)&
     & + sig(3) * sin(theta)) * cos(theta) / 0.3D1) * (sig(4) * cos(thet&
     &a) + sig(5) * sin(theta)) + 0.2D1 / 0.3D1 * ((sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) * cos(theta) / 0.3D1 - (-sig(2) * sin(theta&
     &) + sig(3) * cos(theta)) * sin(theta) / 0.3D1) * (sig(4) * cos(the&
     &ta) + sig(5) * sin(theta)) - 0.2D1 * (-cos(theta) * (-sig(2) * sin&
     &(theta) + sig(3) * cos(theta)) / 0.2D1 + (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) * sin(theta) / 0.2D1) * (-sig(4) * sin(theta) +&
     & sig(5) * cos(theta))) / rho0)
      J(1,3) = -rho0 * (sig(1) * cos(theta) * mu1 * (-0.4D1 / 0.3D1 * (-&
     &0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sin(&
     &theta) + cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta))&
     & / 0.3D1) * sig(1) * cos(theta) + 0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 &
     &* cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + (sig&
     &(2) * cos(theta) + sig(3) * sin(theta)) * sin(theta) / 0.3D1) * si&
     &g(1) * cos(theta) + 0.2D1 / 0.3D1 * ((sig(2) * cos(theta) + sig(3)&
     & * sin(theta)) * sin(theta) / 0.3D1 + cos(theta) * (-sig(2) * sin(&
     &theta) + sig(3) * cos(theta)) / 0.3D1) * sig(1) * cos(theta) + 0.2&
     &D1 * (-(-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(theta) /&
     & 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * cos(theta) &
     &/ 0.2D1) * sig(1) * sin(theta)) / rho0 + sin(theta) * mu1 * (-0.4D&
     &1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * co&
     &s(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta&
     &) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * sin(theta) *&
     &* 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / &
     &0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1 &
     &+ sig(6) ** 2 / 0.6D1) * (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) + 0.2D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-s&
     &ig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) &
     &* sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * c&
     &os(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2&
     & / 0.6D1 + sig(6) ** 2 / 0.6D1) * (sig(2) * cos(theta) + sig(3) * &
     &sin(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + sig(1) ** 2 &
     &* cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(th&
     &eta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) *&
     &* 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * s&
     &in(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(the&
     &ta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * (sig(2) * cos(theta) + &
     &sig(3) * sin(theta)) - 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(the&
     &ta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(&
     &2) * sin(theta) + sig(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(the&
     &ta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(&
     &theta)) / 0.2D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta))) /&
     & rho0 + (sig(2) * cos(theta) + sig(3) * sin(theta)) * mu1 * (-0.4D&
     &1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) * sin(theta) + cos(theta) * (-sig(2) * sin(theta) + sig(3)&
     & * cos(theta)) / 0.3D1) * (sig(2) * cos(theta) + sig(3) * sin(thet&
     &a)) - 0.4D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (s&
     &ig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) *&
     & cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * si&
     &n(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta&
     &)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** &
     &2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * sin(theta) + 0.2D1 / 0.3D1 * (-&
     &0.2D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(&
     &theta)) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sin(theta)&
     & / 0.3D1) * (sig(2) * cos(theta) + sig(3) * sin(theta)) + 0.2D1 / &
     &0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(t&
     &heta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) &
     &+ sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** &
     &2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.6&
     &D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.6D1 + si&
     &g(6) ** 2 / 0.6D1) * sin(theta) + 0.2D1 / 0.3D1 * ((sig(2) * cos(t&
     &heta) + sig(3) * sin(theta)) * sin(theta) / 0.3D1 + cos(theta) * (&
     &-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.3D1) * (sig(2) * c&
     &os(theta) + sig(3) * sin(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 /&
     & 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(the&
     &ta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 /&
     & 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1&
     & + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * si&
     &n(theta) - 0.2D1 * (-(-sig(2) * sin(theta) + sig(3) * cos(theta)) &
     &* sin(theta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(theta))&
     & * cos(theta) / 0.2D1) * (-sig(2) * sin(theta) + sig(3) * cos(thet&
     &a)) - 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D1 - (si&
     &g(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta) +&
     & sig(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(theta) + sig(5) * si&
     &n(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.2D1) &
     &* cos(theta)) / rho0 + (sig(4) * cos(theta) + sig(5) * sin(theta))&
     & * mu1 * (-0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) * sin(theta) + cos(theta) * (-sig(2) * sin(&
     &theta) + sig(3) * cos(theta)) / 0.3D1) * (sig(4) * cos(theta) + si&
     &g(5) * sin(theta)) + 0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * cos(theta) &
     &* (-sig(2) * sin(theta) + sig(3) * cos(theta)) + (sig(2) * cos(the&
     &ta) + sig(3) * sin(theta)) * sin(theta) / 0.3D1) * (sig(4) * cos(t&
     &heta) + sig(5) * sin(theta)) + 0.2D1 / 0.3D1 * ((sig(2) * cos(thet&
     &a) + sig(3) * sin(theta)) * sin(theta) / 0.3D1 + cos(theta) * (-si&
     &g(2) * sin(theta) + sig(3) * cos(theta)) / 0.3D1) * (sig(4) * cos(&
     &theta) + sig(5) * sin(theta)) - 0.2D1 * (-(-sig(2) * sin(theta) + &
     &sig(3) * cos(theta)) * sin(theta) / 0.2D1 - (sig(2) * cos(theta) +&
     & sig(3) * sin(theta)) * cos(theta) / 0.2D1) * (-sig(4) * sin(theta&
     &) + sig(5) * cos(theta))) / rho0) + 0.2D1 * (lambda1 + 0.2D1 / 0.3D1 * &
     &mu1) * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos&
     &(theta)) * sig(6) + sig(1) * sin(theta) * (sig(2) * cos(theta) + s&
     &ig(3) * sin(theta)) * sig(6)) * (-sig(1) * cos(theta) * (-sig(2) *&
     & sin(theta) + sig(3) * cos(theta)) * sig(6) + 0.1D1 - sig(1) * sin&
     &(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) * &
     &(sig(1) * cos(theta) ** 2 * sig(6) + sig(1) * sin(theta) ** 2 * si&
     &g(6)) + (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * (-sig(&
     &2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + sig(1) * sin(the&
     &ta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) ** 2 *&
     & (-sig(1) * cos(theta) ** 2 * sig(6) - sig(1) * sin(theta) ** 2 * &
     &sig(6))
      J(1,4) = -rho0 * (sig(1) * cos(theta) * mu1 * (-0.4D1 / 0.3D1 * (-&
     &0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(&
     &theta) - (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta)&
     & / 0.3D1) * sig(1) * cos(theta) + 0.2D1 / 0.3D1 * (0.2D1 / 0.3D1 *&
     & (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta) + (sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) / 0.3D1) * sig&
     &(1) * cos(theta) + 0.2D1 / 0.3D1 * ((sig(4) * cos(theta) + sig(5) &
     &* sin(theta)) * cos(theta) / 0.3D1 - (-sig(4) * sin(theta) + sig(5&
     &) * cos(theta)) * sin(theta) / 0.3D1) * sig(1) * cos(theta) + 0.2D&
     &1 * (-cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / &
     &0.2D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) /&
     & 0.2D1) * sig(1) * sin(theta)) / rho0 + (sig(2) * cos(theta) + sig&
     &(3) * sin(theta)) * mu1 * (-0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig&
     &(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) - (-sig(4) * &
     &sin(theta) + sig(5) * cos(theta)) * sin(theta) / 0.3D1) * (sig(2) &
     &* cos(theta) + sig(3) * sin(theta)) + 0.2D1 / 0.3D1 * (0.2D1 / 0.3&
     &D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta) + (&
     &sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) / 0.3D1) *&
     & (sig(2) * cos(theta) + sig(3) * sin(theta)) + 0.2D1 / 0.3D1 * ((s&
     &ig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) / 0.3D1 - (&
     &-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta) / 0.3D1) &
     &* (sig(2) * cos(theta) + sig(3) * sin(theta)) - 0.2D1 * (-cos(thet&
     &a) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.2D1 + (sig(4&
     &) * cos(theta) + sig(5) * sin(theta)) * sin(theta) / 0.2D1) * (-si&
     &g(2) * sin(theta) + sig(3) * cos(theta))) / rho0 + cos(theta) * mu&
     &1 * (-0.4D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (s&
     &ig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) *&
     & cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * si&
     &n(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta&
     &)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** &
     &2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (sig(4) * cos(theta) + sig(5) *&
     & sin(theta)) + 0.2D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0&
     &.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 -&
     & (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1&
     &) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) *&
     & sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(th&
     &eta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (sig(4) * cos(theta) +&
     & sig(5) * sin(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + si&
     &g(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3&
     &) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin&
     &(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-&
     &sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4)&
     & * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * (sig(4) * cos&
     &(theta) + sig(5) * sin(theta)) - 0.2D1 * (sig(1) ** 2 * cos(theta)&
     & * sin(theta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)&
     &) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.2D1 - (sig(4)&
     & * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta) + sig&
     &(5) * cos(theta)) / 0.2D1) * (-sig(4) * sin(theta) + sig(5) * cos(&
     &theta)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * sig(6) ** &
     &2) / rho0 + (sig(4) * cos(theta) + sig(5) * sin(theta)) * mu1 * (-&
     &0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * &
     &sin(theta)) * cos(theta) - (-sig(4) * sin(theta) + sig(5) * cos(th&
     &eta)) * sin(theta) / 0.3D1) * (sig(4) * cos(theta) + sig(5) * sin(&
     &theta)) - 0.4D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 &
     &- (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 &
     &* sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(t&
     &heta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta))&
     & ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * cos(theta) + 0.2D1 / 0.3D1 &
     &* (0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * &
     &sin(theta) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(the&
     &ta) / 0.3D1) * (sig(4) * cos(theta) + sig(5) * sin(theta)) + 0.2D1&
     & / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(thet&
     &a) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * cos(theta) &
     &** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / &
     &0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.6D1 +&
     & sig(6) ** 2 / 0.6D1) * cos(theta) + 0.2D1 / 0.3D1 * ((sig(4) * co&
     &s(theta) + sig(5) * sin(theta)) * cos(theta) / 0.3D1 - (-sig(4) * &
     &sin(theta) + sig(5) * cos(theta)) * sin(theta) / 0.3D1) * (sig(4) &
     &* cos(theta) + sig(5) * sin(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** &
     &2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(&
     &theta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) &
     &+ sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** &
     &2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.&
     &6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) *&
     & cos(theta) - 0.2D1 * (-cos(theta) * (-sig(4) * sin(theta) + sig(5&
     &) * cos(theta)) / 0.2D1 + (sig(4) * cos(theta) + sig(5) * sin(thet&
     &a)) * sin(theta) / 0.2D1) * (-sig(4) * sin(theta) + sig(5) * cos(t&
     &heta)) + 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D1 - &
     &(sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta&
     &) + sig(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(theta) + sig(5) *&
     & sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.2D&
     &1) * sin(theta) + cos(theta) * sig(6) ** 2) / rho0)
      J(1,5) = -rho0 * (sig(1) * cos(theta) * mu1 * (-0.4D1 / 0.3D1 * (-&
     &0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(&
     &theta) + cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta))&
     & / 0.3D1) * sig(1) * cos(theta) + 0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 &
     &* cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) + (sig&
     &(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) / 0.3D1) * si&
     &g(1) * cos(theta) + 0.2D1 / 0.3D1 * ((sig(4) * cos(theta) + sig(5)&
     & * sin(theta)) * sin(theta) / 0.3D1 + cos(theta) * (-sig(4) * sin(&
     &theta) + sig(5) * cos(theta)) / 0.3D1) * sig(1) * cos(theta) + 0.2&
     &D1 * (-(-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta) /&
     & 0.2D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) &
     &/ 0.2D1) * sig(1) * sin(theta)) / rho0 + (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) * mu1 * (-0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (si&
     &g(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) + cos(theta)&
     & * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.3D1) * (sig(2)&
     & * cos(theta) + sig(3) * sin(theta)) + 0.2D1 / 0.3D1 * (-0.2D1 / 0&
     &.3D1 * cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) +&
     & (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) / 0.3D1)&
     & * (sig(2) * cos(theta) + sig(3) * sin(theta)) + 0.2D1 / 0.3D1 * (&
     &(sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) / 0.3D1 +&
     & cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.3D1&
     &) * (sig(2) * cos(theta) + sig(3) * sin(theta)) - 0.2D1 * (-(-sig(&
     &4) * sin(theta) + sig(5) * cos(theta)) * sin(theta) / 0.2D1 - (sig&
     &(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) / 0.2D1) * (-&
     &sig(2) * sin(theta) + sig(3) * cos(theta))) / rho0 + sin(theta) * &
     &mu1 * (-0.4D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - &
     &(sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4)&
     & * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * &
     &sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(the&
     &ta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) *&
     &* 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (sig(4) * cos(theta) + sig(5)&
     & * sin(theta)) + 0.2D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 /&
     & 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1&
     & - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig&
     &(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3)&
     & * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(&
     &theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (sig(4) * cos(theta)&
     & + sig(5) * sin(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + &
     &sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig&
     &(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * s&
     &in(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(&
     &4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * (sig(4) * c&
     &os(theta) + sig(5) * sin(theta)) - 0.2D1 * (sig(1) ** 2 * cos(thet&
     &a) * sin(theta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(thet&
     &a)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.2D1 - (sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta) + s&
     &ig(5) * cos(theta)) / 0.2D1) * (-sig(4) * sin(theta) + sig(5) * co&
     &s(theta)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * sig(6) *&
     &* 2) / rho0 + (sig(4) * cos(theta) + sig(5) * sin(theta)) * mu1 * &
     &(-0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) &
     &* sin(theta)) * sin(theta) + cos(theta) * (-sig(4) * sin(theta) + &
     &sig(5) * cos(theta)) / 0.3D1) * (sig(4) * cos(theta) + sig(5) * si&
     &n(theta)) - 0.4D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D&
     &1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (si&
     &g(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** &
     &2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos&
     &(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta&
     &)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * sin(theta) + 0.2D1 / 0.3D&
     &1 * (-0.2D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(theta) + sig(5) &
     &* cos(theta)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(&
     &theta) / 0.3D1) * (sig(4) * cos(theta) + sig(5) * sin(theta)) + 0.&
     &2D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) *&
     & sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(t&
     &heta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * cos(thet&
     &a) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2&
     & / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.6D&
     &1 + sig(6) ** 2 / 0.6D1) * sin(theta) + 0.2D1 / 0.3D1 * ((sig(4) *&
     & cos(theta) + sig(5) * sin(theta)) * sin(theta) / 0.3D1 + cos(thet&
     &a) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.3D1) * (sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) + 0.2D1 / 0.3D1 * (-sig(6) &
     &** 2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * c&
     &os(theta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(thet&
     &a) + sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) &
     &** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 /&
     & 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1&
     &) * sin(theta) - 0.2D1 * (-(-sig(4) * sin(theta) + sig(5) * cos(th&
     &eta)) * sin(theta) / 0.2D1 - (sig(4) * cos(theta) + sig(5) * sin(t&
     &heta)) * cos(theta) / 0.2D1) * (-sig(4) * sin(theta) + sig(5) * co&
     &s(theta)) - 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D1&
     & - (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(th&
     &eta) + sig(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0&
     &.2D1) * cos(theta) + sin(theta) * sig(6) ** 2) / rho0)
      J(1,6) = -rho0 * (-0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 *&
     & mu1 * sig(6) / rho0 - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) ** 2 * mu1 * sig(6) / rho0 + 0.4D1 / 0.3D1 * (sig&
     &(4) * cos(theta) + sig(5) * sin(theta)) ** 2 * mu1 * sig(6) / rho0&
     &) + 0.2D1 * (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * (-&
     &sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + sig(1) * sin&
     &(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) * &
     &(-sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta&
     &)) * sig(6) + 0.1D1 - sig(1) * sin(theta) * (sig(2) * cos(theta) +&
     & sig(3) * sin(theta)) * sig(6)) * (sig(1) * cos(theta) * (-sig(2) &
     &* sin(theta) + sig(3) * cos(theta)) + sig(1) * sin(theta) * (sig(2&
     &) * cos(theta) + sig(3) * sin(theta))) + (lambda1 + 0.2D1 / 0.3D1 * mu1&
     &) * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(th&
     &eta)) * sig(6) + sig(1) * sin(theta) * (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) * sig(6)) ** 2 * (-sig(1) * cos(theta) * (-sig(2)&
     & * sin(theta) + sig(3) * cos(theta)) - sig(1) * sin(theta) * (sig(&
     &2) * cos(theta) + sig(3) * sin(theta)))
      J(2,1) = -rho0 * (-sin(theta) * mu1 * (-0.2D1 / 0.3D1 * (-sig(1) *&
     &* 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * si&
     &n(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta&
     &)) ** 2 / 0.3D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2)&
     & * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin&
     &(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1)&
     & * sig(1) * sin(theta) + 0.4D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta&
     &) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2&
     & / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3&
     &D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * sig(1) * sin&
     &(theta) - 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + sig(1) ** 2 * co&
     &s(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)&
     &) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 &
     &/ 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(t&
     &heta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) &
     &+ sig(5) * cos(theta)) ** 2 / 0.6D1) * sig(1) * sin(theta) - 0.2D1&
     & * (sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D1 - (sig(2) * cos(&
     &theta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) / 0.2D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) *&
     & (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.2D1) * sig(1) * &
     &cos(theta)) / rho0 - sig(1) * sin(theta) * mu1 * (-0.2D1 / 0.3D1 *&
     & (-0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2 + sig(1) * sin(theta) &
     &** 2 / 0.3D1) * sig(1) * sin(theta) - 0.2D1 / 0.3D1 * (-sig(1) ** &
     &2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta))&
     & ** 2 / 0.3D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) *&
     & sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(t&
     &heta) + sig(5) * cos(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) *&
     & sin(theta) + 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * sig(1) * sin(theta&
     &) ** 2 + sig(1) * cos(theta) ** 2 / 0.3D1) * sig(1) * sin(theta) +&
     & 0.4D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2&
     &) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * si&
     &n(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * cos(t&
     &heta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) *&
     &* 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0&
     &.6D1 + sig(6) ** 2 / 0.6D1) * sin(theta) - 0.2D1 / 0.3D1 * (sig(1)&
     & * cos(theta) ** 2 / 0.3D1 + sig(1) * sin(theta) ** 2 / 0.3D1) * s&
     &ig(1) * sin(theta) - 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + sig(1&
     &) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) *&
     & sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(th&
     &eta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig&
     &(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * &
     &sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * sin(theta) - 0.2&
     &D1 * sig(1) ** 2 * cos(theta) ** 2 * sin(theta) - 0.2D1 * (sig(1) &
     &** 2 * cos(theta) * sin(theta) / 0.2D1 - (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) &
     &/ 0.2D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) *&
     & sin(theta) + sig(5) * cos(theta)) / 0.2D1) * cos(theta)) / rho0 +&
     & (-sig(2) * sin(theta) + sig(3) * cos(theta)) * mu1 * (0.2D1 / 0.3&
     &D1 * (-0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2 + sig(1) * sin(the&
     &ta) ** 2 / 0.3D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) -&
     & 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * sig(1) * sin(theta) ** 2 + sig(&
     &1) * cos(theta) ** 2 / 0.3D1) * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) + 0.2D1 / 0.3D1 * (sig(1) * cos(theta) ** 2 / 0.3D1 + s&
     &ig(1) * sin(theta) ** 2 / 0.3D1) * (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) - 0.2D1 * sig(1) * cos(theta) * sin(theta) * (sig(2)&
     & * cos(theta) + sig(3) * sin(theta))) / rho0 + (-sig(4) * sin(thet&
     &a) + sig(5) * cos(theta)) * mu1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1&
     & * sig(1) * cos(theta) ** 2 + sig(1) * sin(theta) ** 2 / 0.3D1) * &
     &(-sig(4) * sin(theta) + sig(5) * cos(theta)) - 0.4D1 / 0.3D1 * (-0&
     &.2D1 / 0.3D1 * sig(1) * sin(theta) ** 2 + sig(1) * cos(theta) ** 2&
     & / 0.3D1) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) + 0.2D1 /&
     & 0.3D1 * (sig(1) * cos(theta) ** 2 / 0.3D1 + sig(1) * sin(theta) *&
     &* 2 / 0.3D1) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) - 0.2D&
     &1 * sig(1) * cos(theta) * sin(theta) * (sig(4) * cos(theta) + sig(&
     &5) * sin(theta))) / rho0) + 0.2D1 * (lambda1 + 0.2D1 / 0.3D1 * mu1) * (&
     &sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta))&
     & * sig(6) + sig(1) * sin(theta) * (sig(2) * cos(theta) + sig(3) * &
     &sin(theta)) * sig(6)) * (-sig(1) * cos(theta) * (-sig(2) * sin(the&
     &ta) + sig(3) * cos(theta)) * sig(6) + 0.1D1 - sig(1) * sin(theta) &
     &* (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) * (cos(the&
     &ta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + sin(&
     &theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) + (&
     &lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * (-sig(2) * sin(&
     &theta) + sig(3) * cos(theta)) * sig(6) + sig(1) * sin(theta) * (si&
     &g(2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) ** 2 * (-cos(th&
     &eta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) - sin&
     &(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6))
      J(2,2) = -rho0 * (-sig(1) * sin(theta) * mu1 * (-0.2D1 / 0.3D1 * (&
     &-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) * cos&
     &(theta) - (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(theta&
     &) / 0.3D1) * sig(1) * sin(theta) + 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 &
     &* (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(theta) + (sig&
     &(2) * cos(theta) + sig(3) * sin(theta)) * cos(theta) / 0.3D1) * si&
     &g(1) * sin(theta) - 0.2D1 / 0.3D1 * ((sig(2) * cos(theta) + sig(3)&
     & * sin(theta)) * cos(theta) / 0.3D1 - (-sig(2) * sin(theta) + sig(&
     &3) * cos(theta)) * sin(theta) / 0.3D1) * sig(1) * sin(theta) - 0.2&
     &D1 * (-cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) /&
     & 0.2D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sin(theta) &
     &/ 0.2D1) * sig(1) * cos(theta)) / rho0 - sin(theta) * mu1 * (0.2D1&
     & / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos&
     &(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta)&
     & + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * sin(theta) **&
     & 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0&
     &.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1 +&
     & sig(6) ** 2 / 0.6D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta&
     &)) - 0.4D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-s&
     &ig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) &
     &* sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * c&
     &os(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2&
     & / 0.6D1 + sig(6) ** 2 / 0.6D1) * (-sig(2) * sin(theta) + sig(3) *&
     & cos(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + sig(1) ** 2&
     & * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(t&
     &heta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) &
     &** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * &
     &sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(th&
     &eta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * (-sig(2) * sin(theta) &
     &+ sig(3) * cos(theta)) - 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(t&
     &heta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-si&
     &g(2) * sin(theta) + sig(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(t&
     &heta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * co&
     &s(theta)) / 0.2D1) * (sig(2) * cos(theta) + sig(3) * sin(theta))) &
     &/ rho0 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) * mu1 * (0.2&
     &D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin&
     &(theta)) * cos(theta) - (-sig(2) * sin(theta) + sig(3) * cos(theta&
     &)) * sin(theta) / 0.3D1) * (-sig(2) * sin(theta) + sig(3) * cos(th&
     &eta)) - 0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - &
     &(sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4)&
     & * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * &
     &sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(the&
     &ta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) *&
     &* 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * sin(theta) - 0.4D1 / 0.3D1 * &
     &(0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * si&
     &n(theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * cos(theta&
     &) / 0.3D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + 0.4D1 &
     &/ 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin&
     &(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta&
     &) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * cos(theta) *&
     &* 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0&
     &.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.6D1 + &
     &sig(6) ** 2 / 0.6D1) * sin(theta) + 0.2D1 / 0.3D1 * ((sig(2) * cos&
     &(theta) + sig(3) * sin(theta)) * cos(theta) / 0.3D1 - (-sig(2) * s&
     &in(theta) + sig(3) * cos(theta)) * sin(theta) / 0.3D1) * (-sig(2) &
     &* sin(theta) + sig(3) * cos(theta)) - 0.2D1 / 0.3D1 * (-sig(6) ** &
     &2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(&
     &theta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) &
     &+ sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** &
     &2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.&
     &6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) *&
     & sin(theta) - 0.2D1 * (-cos(theta) * (-sig(2) * sin(theta) + sig(3&
     &) * cos(theta)) / 0.2D1 + (sig(2) * cos(theta) + sig(3) * sin(thet&
     &a)) * sin(theta) / 0.2D1) * (sig(2) * cos(theta) + sig(3) * sin(th&
     &eta)) - 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D1 - (&
     &sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta)&
     & + sig(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(theta) + sig(5) * &
     &sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.2D1&
     &) * cos(theta)) / rho0 + (-sig(4) * sin(theta) + sig(5) * cos(thet&
     &a)) * mu1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta&
     &) + sig(3) * sin(theta)) * cos(theta) - (-sig(2) * sin(theta) + si&
     &g(3) * cos(theta)) * sin(theta) / 0.3D1) * (-sig(4) * sin(theta) +&
     & sig(5) * cos(theta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * (-sig(2) &
     &* sin(theta) + sig(3) * cos(theta)) * sin(theta) + (sig(2) * cos(t&
     &heta) + sig(3) * sin(theta)) * cos(theta) / 0.3D1) * (-sig(4) * si&
     &n(theta) + sig(5) * cos(theta)) + 0.2D1 / 0.3D1 * ((sig(2) * cos(t&
     &heta) + sig(3) * sin(theta)) * cos(theta) / 0.3D1 - (-sig(2) * sin&
     &(theta) + sig(3) * cos(theta)) * sin(theta) / 0.3D1) * (-sig(4) * &
     &sin(theta) + sig(5) * cos(theta)) - 0.2D1 * (-cos(theta) * (-sig(2&
     &) * sin(theta) + sig(3) * cos(theta)) / 0.2D1 + (sig(2) * cos(thet&
     &a) + sig(3) * sin(theta)) * sin(theta) / 0.2D1) * (sig(4) * cos(th&
     &eta) + sig(5) * sin(theta))) / rho0)
      J(2,3) = -rho0 * (-sig(1) * sin(theta) * mu1 * (-0.2D1 / 0.3D1 * (&
     &-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sin&
     &(theta) + cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)&
     &) / 0.3D1) * sig(1) * sin(theta) + 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1&
     & * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + (si&
     &g(2) * cos(theta) + sig(3) * sin(theta)) * sin(theta) / 0.3D1) * s&
     &ig(1) * sin(theta) - 0.2D1 / 0.3D1 * ((sig(2) * cos(theta) + sig(3&
     &) * sin(theta)) * sin(theta) / 0.3D1 + cos(theta) * (-sig(2) * sin&
     &(theta) + sig(3) * cos(theta)) / 0.3D1) * sig(1) * sin(theta) - 0.&
     &2D1 * (-(-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(theta) &
     &/ 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * cos(theta)&
     & / 0.2D1) * sig(1) * cos(theta)) / rho0 + cos(theta) * mu1 * (0.2D&
     &1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * co&
     &s(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta&
     &) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * sin(theta) *&
     &* 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / &
     &0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1 &
     &+ sig(6) ** 2 / 0.6D1) * (-sig(2) * sin(theta) + sig(3) * cos(thet&
     &a)) - 0.4D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-&
     &sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4)&
     & * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * &
     &cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(thet&
     &a)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** &
     &2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + sig(1) ** &
     &2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta))&
     & ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) *&
     & sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(t&
     &heta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * (-sig(2) * sin(theta)&
     & + sig(3) * cos(theta)) - 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(&
     &theta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-s&
     &ig(2) * sin(theta) + sig(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(&
     &theta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * c&
     &os(theta)) / 0.2D1) * (sig(2) * cos(theta) + sig(3) * sin(theta)))&
     & / rho0 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) * mu1 * (0.&
     &2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * si&
     &n(theta)) * sin(theta) + cos(theta) * (-sig(2) * sin(theta) + sig(&
     &3) * cos(theta)) / 0.3D1) * (-sig(2) * sin(theta) + sig(3) * cos(t&
     &heta)) + 0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 -&
     & (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4&
     &) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 *&
     & sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(th&
     &eta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) &
     &** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * cos(theta) - 0.4D1 / 0.3D1 *&
     & (-0.2D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sin(the&
     &ta) / 0.3D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) - 0.4D&
     &1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * s&
     &in(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(the&
     &ta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * cos(theta)&
     & ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 /&
     & 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.6D1 &
     &+ sig(6) ** 2 / 0.6D1) * cos(theta) + 0.2D1 / 0.3D1 * ((sig(2) * c&
     &os(theta) + sig(3) * sin(theta)) * sin(theta) / 0.3D1 + cos(theta)&
     & * (-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.3D1) * (-sig(2&
     &) * sin(theta) + sig(3) * cos(theta)) + 0.2D1 / 0.3D1 * (-sig(6) *&
     &* 2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * co&
     &s(theta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta&
     &) + sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) *&
     &* 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / &
     &0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1)&
     & * cos(theta) - 0.2D1 * (-(-sig(2) * sin(theta) + sig(3) * cos(the&
     &ta)) * sin(theta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(th&
     &eta)) * cos(theta) / 0.2D1) * (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) - 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D1 -&
     & (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(thet&
     &a) + sig(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(theta) + sig(5) &
     &* sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.2&
     &D1) * sin(theta)) / rho0 + (-sig(4) * sin(theta) + sig(5) * cos(th&
     &eta)) * mu1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(the&
     &ta) + sig(3) * sin(theta)) * sin(theta) + cos(theta) * (-sig(2) * &
     &sin(theta) + sig(3) * cos(theta)) / 0.3D1) * (-sig(4) * sin(theta)&
     & + sig(5) * cos(theta)) - 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * cos(th&
     &eta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + (sig(2) * co&
     &s(theta) + sig(3) * sin(theta)) * sin(theta) / 0.3D1) * (-sig(4) *&
     & sin(theta) + sig(5) * cos(theta)) + 0.2D1 / 0.3D1 * ((sig(2) * co&
     &s(theta) + sig(3) * sin(theta)) * sin(theta) / 0.3D1 + cos(theta) &
     &* (-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.3D1) * (-sig(4)&
     & * sin(theta) + sig(5) * cos(theta)) - 0.2D1 * (-(-sig(2) * sin(th&
     &eta) + sig(3) * cos(theta)) * sin(theta) / 0.2D1 - (sig(2) * cos(t&
     &heta) + sig(3) * sin(theta)) * cos(theta) / 0.2D1) * (sig(4) * cos&
     &(theta) + sig(5) * sin(theta))) / rho0) + 0.2D1 * (lambda1 + 0.2D1 / 0.&
     &3D1 * mu1) * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3)&
     & * cos(theta)) * sig(6) + sig(1) * sin(theta) * (sig(2) * cos(thet&
     &a) + sig(3) * sin(theta)) * sig(6)) * (-sig(1) * cos(theta) * (-si&
     &g(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + 0.1D1 - sig(1)&
     & * sin(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(&
     &6)) * (sig(1) * cos(theta) ** 2 * sig(6) + sig(1) * sin(theta) ** &
     &2 * sig(6)) + (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + sig(1) * s&
     &in(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) &
     &** 2 * (-sig(1) * cos(theta) ** 2 * sig(6) - sig(1) * sin(theta) *&
     &* 2 * sig(6))
      J(2,4) = -rho0 * (-sig(1) * sin(theta) * mu1 * (-0.2D1 / 0.3D1 * (&
     &-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * cos&
     &(theta) - (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta&
     &) / 0.3D1) * sig(1) * sin(theta) + 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 &
     &* (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta) + (sig&
     &(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) / 0.3D1) * si&
     &g(1) * sin(theta) - 0.2D1 / 0.3D1 * ((sig(4) * cos(theta) + sig(5)&
     & * sin(theta)) * cos(theta) / 0.3D1 - (-sig(4) * sin(theta) + sig(&
     &5) * cos(theta)) * sin(theta) / 0.3D1) * sig(1) * sin(theta) - 0.2&
     &D1 * (-cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) /&
     & 0.2D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) &
     &/ 0.2D1) * sig(1) * cos(theta)) / rho0 + (-sig(2) * sin(theta) + s&
     &ig(3) * cos(theta)) * mu1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (si&
     &g(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) - (-sig(4) *&
     & sin(theta) + sig(5) * cos(theta)) * sin(theta) / 0.3D1) * (-sig(2&
     &) * sin(theta) + sig(3) * cos(theta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0&
     &.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta) +&
     & (sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) / 0.3D1)&
     & * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + 0.2D1 / 0.3D1 * &
     &((sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) / 0.3D1 &
     &- (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta) / 0.3D&
     &1) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) - 0.2D1 * (-cos(&
     &theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.2D1 + (s&
     &ig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) / 0.2D1) * &
     &(sig(2) * cos(theta) + sig(3) * sin(theta))) / rho0 - sin(theta) *&
     & mu1 * (0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - &
     &(sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4)&
     & * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * &
     &sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(the&
     &ta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) *&
     &* 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (-sig(4) * sin(theta) + sig(5&
     &) * cos(theta)) - 0.4D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 &
     &/ 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D&
     &1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + si&
     &g(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3&
     &) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin&
     &(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (-sig(4) * sin(thet&
     &a) + sig(5) * cos(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 &
     &+ sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + s&
     &ig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) *&
     & sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 &
     &+ (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-si&
     &g(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * (-sig(4) &
     &* sin(theta) + sig(5) * cos(theta)) - 0.2D1 * (sig(1) ** 2 * cos(t&
     &heta) * sin(theta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(t&
     &heta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.2D1 - (s&
     &ig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta) &
     &+ sig(5) * cos(theta)) / 0.2D1) * (sig(4) * cos(theta) + sig(5) * &
     &sin(theta)) + (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig(6&
     &) ** 2) / rho0 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) * mu&
     &1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(&
     &5) * sin(theta)) * cos(theta) - (-sig(4) * sin(theta) + sig(5) * c&
     &os(theta)) * sin(theta) / 0.3D1) * (-sig(4) * sin(theta) + sig(5) &
     &* cos(theta)) - 0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / &
     &0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 -&
     & (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1)&
     & ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) *&
     & cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(t&
     &heta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * sin(theta) - 0.4D1 / &
     &0.3D1 * (0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(thet&
     &a)) * sin(theta) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * c&
     &os(theta) / 0.3D1) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) &
     &+ 0.4D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(&
     &2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * s&
     &in(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * cos(&
     &theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) &
     &** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / &
     &0.6D1 + sig(6) ** 2 / 0.6D1) * sin(theta) + 0.2D1 / 0.3D1 * ((sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) / 0.3D1 - (-si&
     &g(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta) / 0.3D1) * (&
     &-sig(4) * sin(theta) + sig(5) * cos(theta)) - 0.2D1 / 0.3D1 * (-si&
     &g(6) ** 2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2&
     &) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos&
     &(theta) + sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(th&
     &eta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) *&
     &* 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / &
     &0.6D1) * sin(theta) - 0.2D1 * (-cos(theta) * (-sig(4) * sin(theta)&
     & + sig(5) * cos(theta)) / 0.2D1 + (sig(4) * cos(theta) + sig(5) * &
     &sin(theta)) * sin(theta) / 0.2D1) * (sig(4) * cos(theta) + sig(5) &
     &* sin(theta)) - 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(theta) / 0&
     &.2D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta))&
     & / 0.2D1) * cos(theta) - sin(theta) * sig(6) ** 2) / rho0)
      J(2,5) = -rho0 * (-sig(1) * sin(theta) * mu1 * (-0.2D1 / 0.3D1 * (&
     &-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin&
     &(theta) + cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)&
     &) / 0.3D1) * sig(1) * sin(theta) + 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1&
     & * cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) + (si&
     &g(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) / 0.3D1) * s&
     &ig(1) * sin(theta) - 0.2D1 / 0.3D1 * ((sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) * sin(theta) / 0.3D1 + cos(theta) * (-sig(4) * sin&
     &(theta) + sig(5) * cos(theta)) / 0.3D1) * sig(1) * sin(theta) - 0.&
     &2D1 * (-(-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta) &
     &/ 0.2D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta)&
     & / 0.2D1) * sig(1) * cos(theta)) / rho0 + (-sig(2) * sin(theta) + &
     &sig(3) * cos(theta)) * mu1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (s&
     &ig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) + cos(theta&
     &) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.3D1) * (-sig(&
     &2) * sin(theta) + sig(3) * cos(theta)) - 0.4D1 / 0.3D1 * (-0.2D1 /&
     & 0.3D1 * cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta))&
     & + (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) / 0.3D&
     &1) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + 0.2D1 / 0.3D1 &
     &* ((sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) / 0.3D&
     &1 + cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.&
     &3D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) - 0.2D1 * (-(-&
     &sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta) / 0.2D1 - &
     &(sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) / 0.2D1) &
     &* (sig(2) * cos(theta) + sig(3) * sin(theta))) / rho0 + cos(theta)&
     & * mu1 * (0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 &
     &- (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 &
     &* sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(t&
     &heta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta))&
     & ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (-sig(4) * sin(theta) + sig&
     &(5) * cos(theta)) - 0.4D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** &
     &2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.&
     &3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + &
     &sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig&
     &(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * s&
     &in(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (-sig(4) * sin(th&
     &eta) + sig(5) * cos(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D&
     &1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) +&
     & sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5)&
     & * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D&
     &1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-&
     &sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * (-sig(4&
     &) * sin(theta) + sig(5) * cos(theta)) - 0.2D1 * (sig(1) ** 2 * cos&
     &(theta) * sin(theta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin&
     &(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.2D1 - &
     &(sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta&
     &) + sig(5) * cos(theta)) / 0.2D1) * (sig(4) * cos(theta) + sig(5) &
     &* sin(theta)) + (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig&
     &(6) ** 2) / rho0 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) * &
     &mu1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + si&
     &g(5) * sin(theta)) * sin(theta) + cos(theta) * (-sig(4) * sin(thet&
     &a) + sig(5) * cos(theta)) / 0.3D1) * (-sig(4) * sin(theta) + sig(5&
     &) * cos(theta)) + 0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 &
     &/ 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1&
     & - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(&
     &1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3)&
     & * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos&
     &(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * cos(theta) - 0.4D1 &
     &/ 0.3D1 * (-0.2D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(theta) + s&
     &ig(5) * cos(theta)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) &
     &* sin(theta) / 0.3D1) * (-sig(4) * sin(theta) + sig(5) * cos(theta&
     &)) - 0.4D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-s&
     &ig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) &
     &* sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * c&
     &os(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2&
     & / 0.6D1 + sig(6) ** 2 / 0.6D1) * cos(theta) + 0.2D1 / 0.3D1 * ((s&
     &ig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) / 0.3D1 + c&
     &os(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.3D1) &
     &* (-sig(4) * sin(theta) + sig(5) * cos(theta)) + 0.2D1 / 0.3D1 * (&
     &-sig(6) ** 2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (si&
     &g(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * &
     &cos(theta) + sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin&
     &(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)&
     &) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2&
     & / 0.6D1) * cos(theta) - 0.2D1 * (-(-sig(4) * sin(theta) + sig(5) &
     &* cos(theta)) * sin(theta) / 0.2D1 - (sig(4) * cos(theta) + sig(5)&
     & * sin(theta)) * cos(theta) / 0.2D1) * (sig(4) * cos(theta) + sig(&
     &5) * sin(theta)) - 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(theta) &
     &/ 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) *&
     & sin(theta) + sig(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(theta) &
     &+ sig(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(thet&
     &a)) / 0.2D1) * sin(theta) + cos(theta) * sig(6) ** 2) / rho0)
      J(2,6) = -rho0 * (-0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 *&
     & mu1 * sig(6) / rho0 - 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig&
     &(3) * cos(theta)) ** 2 * mu1 * sig(6) / rho0 + 0.4D1 / 0.3D1 * (-s&
     &ig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 * mu1 * sig(6) / rh&
     &o0) + 0.2D1 * (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + sig(1) * s&
     &in(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) &
     &* (-sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(the&
     &ta)) * sig(6) + 0.1D1 - sig(1) * sin(theta) * (sig(2) * cos(theta)&
     & + sig(3) * sin(theta)) * sig(6)) * (sig(1) * cos(theta) * (-sig(2&
     &) * sin(theta) + sig(3) * cos(theta)) + sig(1) * sin(theta) * (sig&
     &(2) * cos(theta) + sig(3) * sin(theta))) + (lambda1 + 0.2D1 / 0.3D1 * m&
     &u1) * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(&
     &theta)) * sig(6) + sig(1) * sin(theta) * (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) * sig(6)) ** 2 * (-sig(1) * cos(theta) * (-sig(&
     &2) * sin(theta) + sig(3) * cos(theta)) - sig(1) * sin(theta) * (si&
     &g(2) * cos(theta) + sig(3) * sin(theta)))
      J(3,1) = -sig(6) * mu1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * sig(1)&
     & * cos(theta) ** 2 + sig(1) * sin(theta) ** 2 / 0.3D1) * sig(6) + &
     &0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * sig(1) * sin(theta) ** 2 + sig(1&
     &) * cos(theta) ** 2 / 0.3D1) * sig(6) - 0.4D1 / 0.3D1 * (sig(1) * &
     &cos(theta) ** 2 / 0.3D1 + sig(1) * sin(theta) ** 2 / 0.3D1) * sig(&
     &6)) + 0.2D1 * (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + sig(1) * s&
     &in(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) &
     &* (-sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(the&
     &ta)) * sig(6) + 0.1D1 - sig(1) * sin(theta) * (sig(2) * cos(theta)&
     & + sig(3) * sin(theta)) * sig(6)) * (cos(theta) * (-sig(2) * sin(t&
     &heta) + sig(3) * cos(theta)) * sig(6) + sin(theta) * (sig(2) * cos&
     &(theta) + sig(3) * sin(theta)) * sig(6)) + (lambda1 + 0.2D1 / 0.3D1 * m&
     &u1) * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(&
     &theta)) * sig(6) + sig(1) * sin(theta) * (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) * sig(6)) ** 2 * (-cos(theta) * (-sig(2) * sin(&
     &theta) + sig(3) * cos(theta)) * sig(6) - sin(theta) * (sig(2) * co&
     &s(theta) + sig(3) * sin(theta)) * sig(6))
      J(3,2) = -sig(6) * mu1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2&
     &) * cos(theta) + sig(3) * sin(theta)) * cos(theta) - (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta)) * sin(theta) / 0.3D1) * sig(6) + 0&
     &.2D1 / 0.3D1 * (0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) * sin(theta) + (sig(2) * cos(theta) + sig(3) * sin(thet&
     &a)) * cos(theta) / 0.3D1) * sig(6) - 0.4D1 / 0.3D1 * ((sig(2) * co&
     &s(theta) + sig(3) * sin(theta)) * cos(theta) / 0.3D1 - (-sig(2) * &
     &sin(theta) + sig(3) * cos(theta)) * sin(theta) / 0.3D1) * sig(6))
      J(3,3) = -sig(6) * mu1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2&
     &) * cos(theta) + sig(3) * sin(theta)) * sin(theta) + cos(theta) * &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.3D1) * sig(6) + 0&
     &.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(theta&
     &) + sig(3) * cos(theta)) + (sig(2) * cos(theta) + sig(3) * sin(the&
     &ta)) * sin(theta) / 0.3D1) * sig(6) - 0.4D1 / 0.3D1 * ((sig(2) * c&
     &os(theta) + sig(3) * sin(theta)) * sin(theta) / 0.3D1 + cos(theta)&
     & * (-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.3D1) * sig(6))&
     & + 0.2D1 * (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * (-s&
     &ig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + sig(1) * sin(&
     &theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) * (&
     &-sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)&
     &) * sig(6) + 0.1D1 - sig(1) * sin(theta) * (sig(2) * cos(theta) + &
     &sig(3) * sin(theta)) * sig(6)) * (sig(1) * cos(theta) ** 2 * sig(6&
     &) + sig(1) * sin(theta) ** 2 * sig(6)) + (lambda1 + 0.2D1 / 0.3D1 * mu1&
     &) * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(th&
     &eta)) * sig(6) + sig(1) * sin(theta) * (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) * sig(6)) ** 2 * (-sig(1) * cos(theta) ** 2 * sig&
     &(6) - sig(1) * sin(theta) ** 2 * sig(6))
      J(3,4) = -sig(6) * mu1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4&
     &) * cos(theta) + sig(5) * sin(theta)) * cos(theta) - (-sig(4) * si&
     &n(theta) + sig(5) * cos(theta)) * sin(theta) / 0.3D1) * sig(6) + 0&
     &.2D1 / 0.3D1 * (0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * c&
     &os(theta)) * sin(theta) + (sig(4) * cos(theta) + sig(5) * sin(thet&
     &a)) * cos(theta) / 0.3D1) * sig(6) - 0.4D1 / 0.3D1 * ((sig(4) * co&
     &s(theta) + sig(5) * sin(theta)) * cos(theta) / 0.3D1 - (-sig(4) * &
     &sin(theta) + sig(5) * cos(theta)) * sin(theta) / 0.3D1) * sig(6) -&
     & 0.2D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig(6) * s&
     &in(theta) + 0.2D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * &
     &sig(6) * cos(theta))
      J(3,5) = -sig(6) * mu1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4&
     &) * cos(theta) + sig(5) * sin(theta)) * sin(theta) + cos(theta) * &
     &(-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.3D1) * sig(6) + 0&
     &.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(theta&
     &) + sig(5) * cos(theta)) + (sig(4) * cos(theta) + sig(5) * sin(the&
     &ta)) * sin(theta) / 0.3D1) * sig(6) - 0.4D1 / 0.3D1 * ((sig(4) * c&
     &os(theta) + sig(5) * sin(theta)) * sin(theta) / 0.3D1 + cos(theta)&
     & * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.3D1) * sig(6) &
     &+ 0.2D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig(6) * &
     &cos(theta) + 0.2D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) *&
     & sig(6) * sin(theta))
      J(3,6) = -mu1 * (0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 /&
     & 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 &
     &- (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1&
     &) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(&
     &theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * sig(6) + 0.2D1 / 0.3&
     &D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(thet&
     &a) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + s&
     &ig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 /&
     & 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.6D1 &
     &+ (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(6&
     &) ** 2 / 0.6D1) * sig(6) - 0.4D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 +&
     & sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * &
     &sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 +&
     & (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig&
     &(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * sig(6) + (&
     &-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 * sig(6) + (sig(4&
     &) * cos(theta) + sig(5) * sin(theta)) ** 2 * sig(6)) - sig(6) * mu&
     &1 * (0.2D1 * sig(6) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 -&
     & (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 + 0.2D1 &
     &/ 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig(1&
     &) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) ** 2 / 0.3D1 + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta)&
     & + sig(5) * cos(theta)) ** 2) + 0.2D1 * (lambda1 + 0.2D1 / 0.3D1 * mu1)&
     & * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(the&
     &ta)) * sig(6) + sig(1) * sin(theta) * (sig(2) * cos(theta) + sig(3&
     &) * sin(theta)) * sig(6)) * (-sig(1) * cos(theta) * (-sig(2) * sin&
     &(theta) + sig(3) * cos(theta)) * sig(6) + 0.1D1 - sig(1) * sin(the&
     &ta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) * (sig&
     &(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + &
     &sig(1) * sin(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)))&
     & + (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * (-sig(2) * &
     &sin(theta) + sig(3) * cos(theta)) * sig(6) + sig(1) * sin(theta) *&
     & (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) ** 2 * (-si&
     &g(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) -&
     & sig(1) * sin(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta))&
     &)
      J(4,1) = -rho0 * (cos(theta) * mu1 * (-0.2D1 / 0.3D1 * (-sig(1) **&
     & 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin&
     &(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)&
     &) ** 2 / 0.3D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) &
     &* sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(&
     &theta) + sig(5) * cos(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) &
     &* sig(1) * sin(theta) + 0.4D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta)&
     & ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 &
     &/ 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D&
     &1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) +&
     & sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5)&
     & * sin(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * sig(1) * sin(&
     &theta) - 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + sig(1) ** 2 * cos&
     &(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta))&
     & ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 /&
     & 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(th&
     &eta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) +&
     & sig(5) * cos(theta)) ** 2 / 0.6D1) * sig(1) * sin(theta) - 0.2D1 &
     &* (sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D1 - (sig(2) * cos(t&
     &heta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * co&
     &s(theta)) / 0.2D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) * &
     &(-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.2D1) * sig(1) * c&
     &os(theta)) / rho0 + sig(1) * cos(theta) * mu1 * (-0.2D1 / 0.3D1 * &
     &(-0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2 + sig(1) * sin(theta) *&
     &* 2 / 0.3D1) * sig(1) * sin(theta) - 0.2D1 / 0.3D1 * (-sig(1) ** 2&
     & * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(t&
     &heta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) &
     &** 2 / 0.3D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * &
     &sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(th&
     &eta) + sig(5) * cos(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * &
     &sin(theta) + 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * sig(1) * sin(theta)&
     & ** 2 + sig(1) * cos(theta) ** 2 / 0.3D1) * sig(1) * sin(theta) + &
     &0.4D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2)&
     & * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin&
     &(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * cos(th&
     &eta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) **&
     & 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.&
     &6D1 + sig(6) ** 2 / 0.6D1) * sin(theta) - 0.2D1 / 0.3D1 * (sig(1) &
     &* cos(theta) ** 2 / 0.3D1 + sig(1) * sin(theta) ** 2 / 0.3D1) * si&
     &g(1) * sin(theta) - 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + sig(1)&
     & ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * &
     &sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(the&
     &ta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(&
     &2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * s&
     &in(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * sin(theta) - 0.2D&
     &1 * sig(1) ** 2 * cos(theta) ** 2 * sin(theta) - 0.2D1 * (sig(1) *&
     &* 2 * cos(theta) * sin(theta) / 0.2D1 - (sig(2) * cos(theta) + sig&
     &(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) /&
     & 0.2D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * &
     &sin(theta) + sig(5) * cos(theta)) / 0.2D1) * cos(theta)) / rho0 + &
     &(sig(2) * cos(theta) + sig(3) * sin(theta)) * mu1 * (0.2D1 / 0.3D1&
     & * (-0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2 + sig(1) * sin(theta&
     &) ** 2 / 0.3D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) - 0&
     &.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * sig(1) * sin(theta) ** 2 + sig(1)&
     & * cos(theta) ** 2 / 0.3D1) * (-sig(2) * sin(theta) + sig(3) * cos&
     &(theta)) + 0.2D1 / 0.3D1 * (sig(1) * cos(theta) ** 2 / 0.3D1 + sig&
     &(1) * sin(theta) ** 2 / 0.3D1) * (-sig(2) * sin(theta) + sig(3) * &
     &cos(theta)) - 0.2D1 * sig(1) * cos(theta) * sin(theta) * (sig(2) *&
     & cos(theta) + sig(3) * sin(theta))) / rho0 + (sig(4) * cos(theta) &
     &+ sig(5) * sin(theta)) * mu1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * &
     &sig(1) * cos(theta) ** 2 + sig(1) * sin(theta) ** 2 / 0.3D1) * (-s&
     &ig(4) * sin(theta) + sig(5) * cos(theta)) - 0.4D1 / 0.3D1 * (-0.2D&
     &1 / 0.3D1 * sig(1) * sin(theta) ** 2 + sig(1) * cos(theta) ** 2 / &
     &0.3D1) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) + 0.2D1 / 0.&
     &3D1 * (sig(1) * cos(theta) ** 2 / 0.3D1 + sig(1) * sin(theta) ** 2&
     & / 0.3D1) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) - 0.2D1 *&
     & sig(1) * cos(theta) * sin(theta) * (sig(4) * cos(theta) + sig(5) &
     &* sin(theta))) / rho0)
      J(4,2) = -rho0 * (sig(1) * cos(theta) * mu1 * (-0.2D1 / 0.3D1 * (-&
     &0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) * cos(&
     &theta) - (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(theta)&
     & / 0.3D1) * sig(1) * sin(theta) + 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 *&
     & (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(theta) + (sig(&
     &2) * cos(theta) + sig(3) * sin(theta)) * cos(theta) / 0.3D1) * sig&
     &(1) * sin(theta) - 0.2D1 / 0.3D1 * ((sig(2) * cos(theta) + sig(3) &
     &* sin(theta)) * cos(theta) / 0.3D1 - (-sig(2) * sin(theta) + sig(3&
     &) * cos(theta)) * sin(theta) / 0.3D1) * sig(1) * sin(theta) - 0.2D&
     &1 * (-cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) / &
     &0.2D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sin(theta) /&
     & 0.2D1) * sig(1) * cos(theta)) / rho0 + cos(theta) * mu1 * (0.2D1 &
     &/ 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(&
     &theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) &
     &+ sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * sin(theta) ** &
     &2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.&
     &6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1 + &
     &sig(6) ** 2 / 0.6D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta)&
     &) - 0.4D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-si&
     &g(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) *&
     & sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * co&
     &s(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)&
     &) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 &
     &/ 0.6D1 + sig(6) ** 2 / 0.6D1) * (-sig(2) * sin(theta) + sig(3) * &
     &cos(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + sig(1) ** 2 &
     &* cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(th&
     &eta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) *&
     &* 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * s&
     &in(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(the&
     &ta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * (-sig(2) * sin(theta) +&
     & sig(3) * cos(theta)) - 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(th&
     &eta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig&
     &(2) * sin(theta) + sig(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(th&
     &eta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos&
     &(theta)) / 0.2D1) * (sig(2) * cos(theta) + sig(3) * sin(theta))) /&
     & rho0 + (sig(2) * cos(theta) + sig(3) * sin(theta)) * mu1 * (0.2D1&
     & / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(t&
     &heta)) * cos(theta) - (-sig(2) * sin(theta) + sig(3) * cos(theta))&
     & * sin(theta) / 0.3D1) * (-sig(2) * sin(theta) + sig(3) * cos(thet&
     &a)) - 0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (s&
     &ig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) *&
     & cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * si&
     &n(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta&
     &)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** &
     &2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * sin(theta) - 0.4D1 / 0.3D1 * (0&
     &.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(&
     &theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * cos(theta) &
     &/ 0.3D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + 0.4D1 / &
     &0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(t&
     &heta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) &
     &+ sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** &
     &2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.6&
     &D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.6D1 + si&
     &g(6) ** 2 / 0.6D1) * sin(theta) + 0.2D1 / 0.3D1 * ((sig(2) * cos(t&
     &heta) + sig(3) * sin(theta)) * cos(theta) / 0.3D1 - (-sig(2) * sin&
     &(theta) + sig(3) * cos(theta)) * sin(theta) / 0.3D1) * (-sig(2) * &
     &sin(theta) + sig(3) * cos(theta)) - 0.2D1 / 0.3D1 * (-sig(6) ** 2 &
     &/ 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(th&
     &eta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + &
     &sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 &
     &/ 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D&
     &1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * s&
     &in(theta) - 0.2D1 * (-cos(theta) * (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) / 0.2D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)&
     &) * sin(theta) / 0.2D1) * (sig(2) * cos(theta) + sig(3) * sin(thet&
     &a)) - 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D1 - (si&
     &g(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta) +&
     & sig(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(theta) + sig(5) * si&
     &n(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.2D1) &
     &* cos(theta)) / rho0 + (sig(4) * cos(theta) + sig(5) * sin(theta))&
     & * mu1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) +&
     & sig(3) * sin(theta)) * cos(theta) - (-sig(2) * sin(theta) + sig(3&
     &) * cos(theta)) * sin(theta) / 0.3D1) * (-sig(4) * sin(theta) + si&
     &g(5) * cos(theta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * (-sig(2) * s&
     &in(theta) + sig(3) * cos(theta)) * sin(theta) + (sig(2) * cos(thet&
     &a) + sig(3) * sin(theta)) * cos(theta) / 0.3D1) * (-sig(4) * sin(t&
     &heta) + sig(5) * cos(theta)) + 0.2D1 / 0.3D1 * ((sig(2) * cos(thet&
     &a) + sig(3) * sin(theta)) * cos(theta) / 0.3D1 - (-sig(2) * sin(th&
     &eta) + sig(3) * cos(theta)) * sin(theta) / 0.3D1) * (-sig(4) * sin&
     &(theta) + sig(5) * cos(theta)) - 0.2D1 * (-cos(theta) * (-sig(2) *&
     & sin(theta) + sig(3) * cos(theta)) / 0.2D1 + (sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) * sin(theta) / 0.2D1) * (sig(4) * cos(theta&
     &) + sig(5) * sin(theta))) / rho0)
      J(4,3) = -rho0 * (sig(1) * cos(theta) * mu1 * (-0.2D1 / 0.3D1 * (-&
     &0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sin(&
     &theta) + cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta))&
     & / 0.3D1) * sig(1) * sin(theta) + 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 &
     &* cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + (sig&
     &(2) * cos(theta) + sig(3) * sin(theta)) * sin(theta) / 0.3D1) * si&
     &g(1) * sin(theta) - 0.2D1 / 0.3D1 * ((sig(2) * cos(theta) + sig(3)&
     & * sin(theta)) * sin(theta) / 0.3D1 + cos(theta) * (-sig(2) * sin(&
     &theta) + sig(3) * cos(theta)) / 0.3D1) * sig(1) * sin(theta) - 0.2&
     &D1 * (-(-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(theta) /&
     & 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * cos(theta) &
     &/ 0.2D1) * sig(1) * cos(theta)) / rho0 + sin(theta) * mu1 * (0.2D1&
     & / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos&
     &(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta)&
     & + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * sin(theta) **&
     & 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0&
     &.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1 +&
     & sig(6) ** 2 / 0.6D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta&
     &)) - 0.4D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-s&
     &ig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) &
     &* sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * c&
     &os(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2&
     & / 0.6D1 + sig(6) ** 2 / 0.6D1) * (-sig(2) * sin(theta) + sig(3) *&
     & cos(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + sig(1) ** 2&
     & * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(t&
     &heta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) &
     &** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * &
     &sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(th&
     &eta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * (-sig(2) * sin(theta) &
     &+ sig(3) * cos(theta)) - 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(t&
     &heta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-si&
     &g(2) * sin(theta) + sig(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(t&
     &heta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * co&
     &s(theta)) / 0.2D1) * (sig(2) * cos(theta) + sig(3) * sin(theta))) &
     &/ rho0 + (sig(2) * cos(theta) + sig(3) * sin(theta)) * mu1 * (0.2D&
     &1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) * sin(theta) + cos(theta) * (-sig(2) * sin(theta) + sig(3)&
     & * cos(theta)) / 0.3D1) * (-sig(2) * sin(theta) + sig(3) * cos(the&
     &ta)) + 0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (&
     &sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) &
     &* cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * s&
     &in(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(thet&
     &a)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) **&
     & 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * cos(theta) - 0.4D1 / 0.3D1 * (&
     &-0.2D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos&
     &(theta)) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * sin(theta&
     &) / 0.3D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) - 0.4D1 &
     &/ 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin&
     &(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta&
     &) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * cos(theta) *&
     &* 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0&
     &.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.6D1 + &
     &sig(6) ** 2 / 0.6D1) * cos(theta) + 0.2D1 / 0.3D1 * ((sig(2) * cos&
     &(theta) + sig(3) * sin(theta)) * sin(theta) / 0.3D1 + cos(theta) *&
     & (-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.3D1) * (-sig(2) &
     &* sin(theta) + sig(3) * cos(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** &
     &2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(&
     &theta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) &
     &+ sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** &
     &2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.&
     &6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) *&
     & cos(theta) - 0.2D1 * (-(-sig(2) * sin(theta) + sig(3) * cos(theta&
     &)) * sin(theta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(thet&
     &a)) * cos(theta) / 0.2D1) * (sig(2) * cos(theta) + sig(3) * sin(th&
     &eta)) - 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D1 - (&
     &sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta)&
     & + sig(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(theta) + sig(5) * &
     &sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.2D1&
     &) * sin(theta)) / rho0 + (sig(4) * cos(theta) + sig(5) * sin(theta&
     &)) * mu1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta)&
     & + sig(3) * sin(theta)) * sin(theta) + cos(theta) * (-sig(2) * sin&
     &(theta) + sig(3) * cos(theta)) / 0.3D1) * (-sig(4) * sin(theta) + &
     &sig(5) * cos(theta)) - 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * cos(theta&
     &) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + (sig(2) * cos(t&
     &heta) + sig(3) * sin(theta)) * sin(theta) / 0.3D1) * (-sig(4) * si&
     &n(theta) + sig(5) * cos(theta)) + 0.2D1 / 0.3D1 * ((sig(2) * cos(t&
     &heta) + sig(3) * sin(theta)) * sin(theta) / 0.3D1 + cos(theta) * (&
     &-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.3D1) * (-sig(4) * &
     &sin(theta) + sig(5) * cos(theta)) - 0.2D1 * (-(-sig(2) * sin(theta&
     &) + sig(3) * cos(theta)) * sin(theta) / 0.2D1 - (sig(2) * cos(thet&
     &a) + sig(3) * sin(theta)) * cos(theta) / 0.2D1) * (sig(4) * cos(th&
     &eta) + sig(5) * sin(theta))) / rho0)
      J(4,4) = -rho0 * (sig(1) * cos(theta) * mu1 * (-0.2D1 / 0.3D1 * (-&
     &0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(&
     &theta) - (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta)&
     & / 0.3D1) * sig(1) * sin(theta) + 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 *&
     & (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta) + (sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) / 0.3D1) * sig&
     &(1) * sin(theta) - 0.2D1 / 0.3D1 * ((sig(4) * cos(theta) + sig(5) &
     &* sin(theta)) * cos(theta) / 0.3D1 - (-sig(4) * sin(theta) + sig(5&
     &) * cos(theta)) * sin(theta) / 0.3D1) * sig(1) * sin(theta) - 0.2D&
     &1 * (-cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / &
     &0.2D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) /&
     & 0.2D1) * sig(1) * cos(theta)) / rho0 + (sig(2) * cos(theta) + sig&
     &(3) * sin(theta)) * mu1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) - (-sig(4) * s&
     &in(theta) + sig(5) * cos(theta)) * sin(theta) / 0.3D1) * (-sig(2) &
     &* sin(theta) + sig(3) * cos(theta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3&
     &D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta) + (&
     &sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) / 0.3D1) *&
     & (-sig(2) * sin(theta) + sig(3) * cos(theta)) + 0.2D1 / 0.3D1 * ((&
     &sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) / 0.3D1 - &
     &(-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta) / 0.3D1)&
     & * (-sig(2) * sin(theta) + sig(3) * cos(theta)) - 0.2D1 * (-cos(th&
     &eta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.2D1 + (sig&
     &(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) / 0.2D1) * (s&
     &ig(2) * cos(theta) + sig(3) * sin(theta))) / rho0 + cos(theta) * m&
     &u1 * (0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (s&
     &ig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) *&
     & cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * si&
     &n(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta&
     &)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** &
     &2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (-sig(4) * sin(theta) + sig(5) &
     &* cos(theta)) - 0.4D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / &
     &0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 &
     &- (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(&
     &1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) &
     &* sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(t&
     &heta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (-sig(4) * sin(theta)&
     & + sig(5) * cos(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + &
     &sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig&
     &(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * s&
     &in(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(&
     &4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * (-sig(4) * &
     &sin(theta) + sig(5) * cos(theta)) - 0.2D1 * (sig(1) ** 2 * cos(the&
     &ta) * sin(theta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(the&
     &ta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.2D1 - (sig&
     &(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta) + &
     &sig(5) * cos(theta)) / 0.2D1) * (sig(4) * cos(theta) + sig(5) * si&
     &n(theta)) + (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig(6) &
     &** 2) / rho0 + (sig(4) * cos(theta) + sig(5) * sin(theta)) * mu1 *&
     & (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) &
     &* sin(theta)) * cos(theta) - (-sig(4) * sin(theta) + sig(5) * cos(&
     &theta)) * sin(theta) / 0.3D1) * (-sig(4) * sin(theta) + sig(5) * c&
     &os(theta)) - 0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3&
     &D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (s&
     &ig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) **&
     & 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * co&
     &s(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(thet&
     &a)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * sin(theta) - 0.4D1 / 0.3&
     &D1 * (0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta))&
     & * sin(theta) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(&
     &theta) / 0.3D1) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) + 0&
     &.4D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) &
     &* sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(&
     &theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * cos(the&
     &ta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** &
     &2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.6&
     &D1 + sig(6) ** 2 / 0.6D1) * sin(theta) + 0.2D1 / 0.3D1 * ((sig(4) &
     &* cos(theta) + sig(5) * sin(theta)) * cos(theta) / 0.3D1 - (-sig(4&
     &) * sin(theta) + sig(5) * cos(theta)) * sin(theta) / 0.3D1) * (-si&
     &g(4) * sin(theta) + sig(5) * cos(theta)) - 0.2D1 / 0.3D1 * (-sig(6&
     &) ** 2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) *&
     & cos(theta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(th&
     &eta) + sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta&
     &) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2&
     & / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6&
     &D1) * sin(theta) - 0.2D1 * (-cos(theta) * (-sig(4) * sin(theta) + &
     &sig(5) * cos(theta)) / 0.2D1 + (sig(4) * cos(theta) + sig(5) * sin&
     &(theta)) * sin(theta) / 0.2D1) * (sig(4) * cos(theta) + sig(5) * s&
     &in(theta)) - 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(theta) / 0.2D&
     &1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(t&
     &heta) + sig(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(theta) + sig(&
     &5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / &
     &0.2D1) * cos(theta) - sin(theta) * sig(6) ** 2) / rho0)
      J(4,5) = -rho0 * (sig(1) * cos(theta) * mu1 * (-0.2D1 / 0.3D1 * (-&
     &0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(&
     &theta) + cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta))&
     & / 0.3D1) * sig(1) * sin(theta) + 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 &
     &* cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) + (sig&
     &(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) / 0.3D1) * si&
     &g(1) * sin(theta) - 0.2D1 / 0.3D1 * ((sig(4) * cos(theta) + sig(5)&
     & * sin(theta)) * sin(theta) / 0.3D1 + cos(theta) * (-sig(4) * sin(&
     &theta) + sig(5) * cos(theta)) / 0.3D1) * sig(1) * sin(theta) - 0.2&
     &D1 * (-(-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta) /&
     & 0.2D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) &
     &/ 0.2D1) * sig(1) * cos(theta)) / rho0 + (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) * mu1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig&
     &(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) + cos(theta) &
     &* (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.3D1) * (-sig(2)&
     & * sin(theta) + sig(3) * cos(theta)) - 0.4D1 / 0.3D1 * (-0.2D1 / 0&
     &.3D1 * cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) +&
     & (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) / 0.3D1)&
     & * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + 0.2D1 / 0.3D1 * &
     &((sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) / 0.3D1 &
     &+ cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.3D&
     &1) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) - 0.2D1 * (-(-si&
     &g(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta) / 0.2D1 - (s&
     &ig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) / 0.2D1) * &
     &(sig(2) * cos(theta) + sig(3) * sin(theta))) / rho0 + sin(theta) *&
     & mu1 * (0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - &
     &(sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4)&
     & * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * &
     &sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(the&
     &ta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) *&
     &* 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (-sig(4) * sin(theta) + sig(5&
     &) * cos(theta)) - 0.4D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 &
     &/ 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D&
     &1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + si&
     &g(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3&
     &) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin&
     &(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * (-sig(4) * sin(thet&
     &a) + sig(5) * cos(theta)) + 0.2D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 &
     &+ sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + s&
     &ig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) *&
     & sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 &
     &+ (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-si&
     &g(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * (-sig(4) &
     &* sin(theta) + sig(5) * cos(theta)) - 0.2D1 * (sig(1) ** 2 * cos(t&
     &heta) * sin(theta) / 0.2D1 - (sig(2) * cos(theta) + sig(3) * sin(t&
     &heta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.2D1 - (s&
     &ig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta) &
     &+ sig(5) * cos(theta)) / 0.2D1) * (sig(4) * cos(theta) + sig(5) * &
     &sin(theta)) + (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig(6&
     &) ** 2) / rho0 + (sig(4) * cos(theta) + sig(5) * sin(theta)) * mu1&
     & * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) * sin(theta) + cos(theta) * (-sig(4) * sin(theta) &
     &+ sig(5) * cos(theta)) / 0.3D1) * (-sig(4) * sin(theta) + sig(5) *&
     & cos(theta)) + 0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / 0&
     &.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - &
     &(sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1) &
     &** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * &
     &cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(th&
     &eta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * cos(theta) - 0.4D1 / 0&
     &.3D1 * (-0.2D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(theta) + sig(&
     &5) * cos(theta)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * s&
     &in(theta) / 0.3D1) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) &
     &- 0.4D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(&
     &2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * s&
     &in(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * cos(&
     &theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) &
     &** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / &
     &0.6D1 + sig(6) ** 2 / 0.6D1) * cos(theta) + 0.2D1 / 0.3D1 * ((sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) / 0.3D1 + cos(&
     &theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) / 0.3D1) * (&
     &-sig(4) * sin(theta) + sig(5) * cos(theta)) + 0.2D1 / 0.3D1 * (-si&
     &g(6) ** 2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2&
     &) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos&
     &(theta) + sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(th&
     &eta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) *&
     &* 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / &
     &0.6D1) * cos(theta) - 0.2D1 * (-(-sig(4) * sin(theta) + sig(5) * c&
     &os(theta)) * sin(theta) / 0.2D1 - (sig(4) * cos(theta) + sig(5) * &
     &sin(theta)) * cos(theta) / 0.2D1) * (sig(4) * cos(theta) + sig(5) &
     &* sin(theta)) - 0.2D1 * (sig(1) ** 2 * cos(theta) * sin(theta) / 0&
     &.2D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta)) / 0.2D1 - (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta))&
     & / 0.2D1) * sin(theta) + cos(theta) * sig(6) ** 2) / rho0)
      J(4,6) = -rho0 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) * mu1 *&
     & sig(6) * sin(theta) / rho0 - 0.2D1 / 0.3D1 * (sig(2) * cos(theta)&
     & + sig(3) * sin(theta)) * mu1 * sig(6) * (-sig(2) * sin(theta) + s&
     &ig(3) * cos(theta)) / rho0 + 0.4D1 / 0.3D1 * (sig(4) * cos(theta) &
     &+ sig(5) * sin(theta)) * mu1 * (-sig(4) * sin(theta) + sig(5) * co&
     &s(theta)) * sig(6) / rho0)
      J(5,1) = -rho0 * (-sin(theta) * mu1 * (-(-sig(4) * sin(theta) + si&
     &g(5) * cos(theta)) * sig(6) * sig(1) * sin(theta) + (sig(4) * cos(&
     &theta) + sig(5) * sin(theta)) * sig(6) * sig(1) * cos(theta)) / rh&
     &o0 - sig(1) * sin(theta) * mu1 * (-(-sig(4) * sin(theta) + sig(5) &
     &* cos(theta)) * sig(6) * sin(theta) + (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) * sig(6) * cos(theta)) / rho0 + (-sig(4) * sin(the&
     &ta) + sig(5) * cos(theta)) * mu1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D&
     &1 * sig(1) * cos(theta) ** 2 + sig(1) * sin(theta) ** 2 / 0.3D1) *&
     & sig(6) + 0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * sig(1) * sin(theta) **&
     & 2 + sig(1) * cos(theta) ** 2 / 0.3D1) * sig(6) - 0.4D1 / 0.3D1 * &
     &(sig(1) * cos(theta) ** 2 / 0.3D1 + sig(1) * sin(theta) ** 2 / 0.3&
     &D1) * sig(6)) / rho0)
      J(5,2) = -rho0 * (-sin(theta) * mu1 * ((-sig(4) * sin(theta) + sig&
     &(5) * cos(theta)) * sig(6) * (-sig(2) * sin(theta) + sig(3) * cos(&
     &theta)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * sig(6) * (&
     &sig(2) * cos(theta) + sig(3) * sin(theta))) / rho0 + (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta)) * mu1 * (-(-sig(4) * sin(theta) + &
     &sig(5) * cos(theta)) * sig(6) * sin(theta) + (sig(4) * cos(theta) &
     &+ sig(5) * sin(theta)) * sig(6) * cos(theta)) / rho0 + (-sig(4) * &
     &sin(theta) + sig(5) * cos(theta)) * mu1 * (0.2D1 / 0.3D1 * (-0.2D1&
     & / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) * cos(theta&
     &) - (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(theta) / 0.&
     &3D1) * sig(6) + 0.2D1 / 0.3D1 * (0.2D1 / 0.3D1 * (-sig(2) * sin(th&
     &eta) + sig(3) * cos(theta)) * sin(theta) + (sig(2) * cos(theta) + &
     &sig(3) * sin(theta)) * cos(theta) / 0.3D1) * sig(6) - 0.4D1 / 0.3D&
     &1 * ((sig(2) * cos(theta) + sig(3) * sin(theta)) * cos(theta) / 0.&
     &3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(theta) / &
     &0.3D1) * sig(6)) / rho0)
      J(5,3) = -rho0 * (cos(theta) * mu1 * ((-sig(4) * sin(theta) + sig(&
     &5) * cos(theta)) * sig(6) * (-sig(2) * sin(theta) + sig(3) * cos(t&
     &heta)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * sig(6) * (s&
     &ig(2) * cos(theta) + sig(3) * sin(theta))) / rho0 + (-sig(2) * sin&
     &(theta) + sig(3) * cos(theta)) * mu1 * ((-sig(4) * sin(theta) + si&
     &g(5) * cos(theta)) * sig(6) * cos(theta) + (sig(4) * cos(theta) + &
     &sig(5) * sin(theta)) * sig(6) * sin(theta)) / rho0 + (-sig(4) * si&
     &n(theta) + sig(5) * cos(theta)) * mu1 * (0.2D1 / 0.3D1 * (-0.2D1 /&
     & 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sin(theta) &
     &+ cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.3D&
     &1) * sig(6) + 0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * cos(theta) * (-sig&
     &(2) * sin(theta) + sig(3) * cos(theta)) + (sig(2) * cos(theta) + s&
     &ig(3) * sin(theta)) * sin(theta) / 0.3D1) * sig(6) - 0.4D1 / 0.3D1&
     & * ((sig(2) * cos(theta) + sig(3) * sin(theta)) * sin(theta) / 0.3&
     &D1 + cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0&
     &.3D1) * sig(6)) / rho0)
      J(5,4) = -rho0 * (-sig(1) * sin(theta) * mu1 * (sig(1) * cos(theta&
     &) ** 2 * sig(6) + sig(1) * sin(theta) ** 2 * sig(6)) / rho0 + (-si&
     &g(2) * sin(theta) + sig(3) * cos(theta)) * mu1 * (-sin(theta) * si&
     &g(6) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + cos(theta) *&
     & sig(6) * (sig(2) * cos(theta) + sig(3) * sin(theta))) / rho0 - si&
     &n(theta) * mu1 * (0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 &
     &/ 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1&
     & - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(&
     &1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3)&
     & * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos&
     &(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * sig(6) + 0.2D1 / 0.&
     &3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(the&
     &ta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + &
     &sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 &
     &/ 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.6D1&
     & + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(&
     &6) ** 2 / 0.6D1) * sig(6) - 0.4D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 &
     &+ sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + s&
     &ig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) *&
     & sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 &
     &+ (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-si&
     &g(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * sig(6) + &
     &(-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 * sig(6) + (sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) ** 2 * sig(6)) / rho0 + (-s&
     &ig(4) * sin(theta) + sig(5) * cos(theta)) * mu1 * (0.2D1 / 0.3D1 *&
     & (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * c&
     &os(theta) - (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(the&
     &ta) / 0.3D1) * sig(6) + 0.2D1 / 0.3D1 * (0.2D1 / 0.3D1 * (-sig(4) &
     &* sin(theta) + sig(5) * cos(theta)) * sin(theta) + (sig(4) * cos(t&
     &heta) + sig(5) * sin(theta)) * cos(theta) / 0.3D1) * sig(6) - 0.4D&
     &1 / 0.3D1 * ((sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(the&
     &ta) / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(t&
     &heta) / 0.3D1) * sig(6) - 0.2D1 * (-sig(4) * sin(theta) + sig(5) *&
     & cos(theta)) * sig(6) * sin(theta) + 0.2D1 * (sig(4) * cos(theta) &
     &+ sig(5) * sin(theta)) * sig(6) * cos(theta)) / rho0)
      J(5,5) = -rho0 * ((-sig(2) * sin(theta) + sig(3) * cos(theta)) * m&
     &u1 * (cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * &
     &sig(6) + sin(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) &
     &* sig(6)) / rho0 + cos(theta) * mu1 * (0.2D1 / 0.3D1 * (-sig(1) **&
     & 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin&
     &(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)&
     &) ** 2 / 0.3D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) &
     &* sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(&
     &theta) + sig(5) * cos(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) &
     &* sig(6) + 0.2D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1&
     & - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-s&
     &ig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) **&
     & 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin&
     &(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)&
     &) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * sig(6) - 0.4D1 / 0.3D1 * (&
     &-sig(6) ** 2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (si&
     &g(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * &
     &cos(theta) + sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin&
     &(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)&
     &) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2&
     & / 0.6D1) * sig(6) + (-sig(4) * sin(theta) + sig(5) * cos(theta)) &
     &** 2 * sig(6) + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 *&
     & sig(6)) / rho0 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) * m&
     &u1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig&
     &(5) * sin(theta)) * sin(theta) + cos(theta) * (-sig(4) * sin(theta&
     &) + sig(5) * cos(theta)) / 0.3D1) * sig(6) + 0.2D1 / 0.3D1 * (-0.2&
     &D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(the&
     &ta)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) / &
     &0.3D1) * sig(6) - 0.4D1 / 0.3D1 * ((sig(4) * cos(theta) + sig(5) *&
     & sin(theta)) * sin(theta) / 0.3D1 + cos(theta) * (-sig(4) * sin(th&
     &eta) + sig(5) * cos(theta)) / 0.3D1) * sig(6) + 0.2D1 * (-sig(4) *&
     & sin(theta) + sig(5) * cos(theta)) * sig(6) * cos(theta) + 0.2D1 *&
     & (sig(4) * cos(theta) + sig(5) * sin(theta)) * sig(6) * sin(theta)&
     &) / rho0)
      J(5,6) = -rho0 * (-sig(1) * sin(theta) * mu1 * (-(-sig(4) * sin(th&
     &eta) + sig(5) * cos(theta)) * sig(1) * sin(theta) + (sig(4) * cos(&
     &theta) + sig(5) * sin(theta)) * sig(1) * cos(theta)) / rho0 + (-si&
     &g(2) * sin(theta) + sig(3) * cos(theta)) * mu1 * ((-sig(4) * sin(t&
     &heta) + sig(5) * cos(theta)) * (-sig(2) * sin(theta) + sig(3) * co&
     &s(theta)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * (sig(2) &
     &* cos(theta) + sig(3) * sin(theta))) / rho0 + (-sig(4) * sin(theta&
     &) + sig(5) * cos(theta)) * mu1 * (0.2D1 * sig(6) ** 2 - sig(1) ** &
     &2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) ** 2 / 0.3D1 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(&
     &5) * sin(theta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-&
     &sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 + 0.2D1 / &
     &0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2) / rho0)
      J(6,1) = -rho0 * (cos(theta) * mu1 * (-(-sig(4) * sin(theta) + sig&
     &(5) * cos(theta)) * sig(6) * sig(1) * sin(theta) + (sig(4) * cos(t&
     &heta) + sig(5) * sin(theta)) * sig(6) * sig(1) * cos(theta)) / rho&
     &0 + sig(1) * cos(theta) * mu1 * (-(-sig(4) * sin(theta) + sig(5) *&
     & cos(theta)) * sig(6) * sin(theta) + (sig(4) * cos(theta) + sig(5)&
     & * sin(theta)) * sig(6) * cos(theta)) / rho0 + (sig(4) * cos(theta&
     &) + sig(5) * sin(theta)) * mu1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 &
     &* sig(1) * cos(theta) ** 2 + sig(1) * sin(theta) ** 2 / 0.3D1) * s&
     &ig(6) + 0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * sig(1) * sin(theta) ** 2&
     & + sig(1) * cos(theta) ** 2 / 0.3D1) * sig(6) - 0.4D1 / 0.3D1 * (s&
     &ig(1) * cos(theta) ** 2 / 0.3D1 + sig(1) * sin(theta) ** 2 / 0.3D1&
     &) * sig(6)) / rho0)
      J(6,2) = -rho0 * (cos(theta) * mu1 * ((-sig(4) * sin(theta) + sig(&
     &5) * cos(theta)) * sig(6) * (-sig(2) * sin(theta) + sig(3) * cos(t&
     &heta)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * sig(6) * (s&
     &ig(2) * cos(theta) + sig(3) * sin(theta))) / rho0 + (sig(2) * cos(&
     &theta) + sig(3) * sin(theta)) * mu1 * (-(-sig(4) * sin(theta) + si&
     &g(5) * cos(theta)) * sig(6) * sin(theta) + (sig(4) * cos(theta) + &
     &sig(5) * sin(theta)) * sig(6) * cos(theta)) / rho0 + (sig(4) * cos&
     &(theta) + sig(5) * sin(theta)) * mu1 * (0.2D1 / 0.3D1 * (-0.2D1 / &
     &0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) * cos(theta) -&
     & (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(theta) / 0.3D1&
     &) * sig(6) + 0.2D1 / 0.3D1 * (0.2D1 / 0.3D1 * (-sig(2) * sin(theta&
     &) + sig(3) * cos(theta)) * sin(theta) + (sig(2) * cos(theta) + sig&
     &(3) * sin(theta)) * cos(theta) / 0.3D1) * sig(6) - 0.4D1 / 0.3D1 *&
     & ((sig(2) * cos(theta) + sig(3) * sin(theta)) * cos(theta) / 0.3D1&
     & - (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(theta) / 0.3&
     &D1) * sig(6)) / rho0)
      J(6,3) = -rho0 * (sin(theta) * mu1 * ((-sig(4) * sin(theta) + sig(&
     &5) * cos(theta)) * sig(6) * (-sig(2) * sin(theta) + sig(3) * cos(t&
     &heta)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * sig(6) * (s&
     &ig(2) * cos(theta) + sig(3) * sin(theta))) / rho0 + (sig(2) * cos(&
     &theta) + sig(3) * sin(theta)) * mu1 * ((-sig(4) * sin(theta) + sig&
     &(5) * cos(theta)) * sig(6) * cos(theta) + (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) * sig(6) * sin(theta)) / rho0 + (sig(4) * cos(&
     &theta) + sig(5) * sin(theta)) * mu1 * (0.2D1 / 0.3D1 * (-0.2D1 / 0&
     &.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sin(theta) + &
     &cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.3D1)&
     & * sig(6) + 0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * cos(theta) * (-sig(2&
     &) * sin(theta) + sig(3) * cos(theta)) + (sig(2) * cos(theta) + sig&
     &(3) * sin(theta)) * sin(theta) / 0.3D1) * sig(6) - 0.4D1 / 0.3D1 *&
     & ((sig(2) * cos(theta) + sig(3) * sin(theta)) * sin(theta) / 0.3D1&
     & + cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) / 0.3&
     &D1) * sig(6)) / rho0)
      J(6,4) = -rho0 * (sig(1) * cos(theta) * mu1 * (sig(1) * cos(theta)&
     & ** 2 * sig(6) + sig(1) * sin(theta) ** 2 * sig(6)) / rho0 + (sig(&
     &2) * cos(theta) + sig(3) * sin(theta)) * mu1 * (-sin(theta) * sig(&
     &6) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + cos(theta) * s&
     &ig(6) * (sig(2) * cos(theta) + sig(3) * sin(theta))) / rho0 + cos(&
     &theta) * mu1 * (0.2D1 / 0.3D1 * (-sig(1) ** 2 * cos(theta) ** 2 / &
     &0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 -&
     & (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 + sig(1)&
     & ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) *&
     & cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(t&
     &heta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * sig(6) + 0.2D1 / 0.3D&
     &1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta&
     &) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + si&
     &g(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / &
     &0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.6D1 +&
     & (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(6)&
     & ** 2 / 0.6D1) * sig(6) - 0.4D1 / 0.3D1 * (-sig(6) ** 2 / 0.3D1 + &
     &sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig&
     &(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * s&
     &in(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(&
     &4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.6D1) * sig(6) + (-&
     &sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 * sig(6) + (sig(4)&
     & * cos(theta) + sig(5) * sin(theta)) ** 2 * sig(6)) / rho0 + (sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) * mu1 * (0.2D1 / 0.3D1 * (-&
     &0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(&
     &theta) - (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta)&
     & / 0.3D1) * sig(6) + 0.2D1 / 0.3D1 * (0.2D1 / 0.3D1 * (-sig(4) * s&
     &in(theta) + sig(5) * cos(theta)) * sin(theta) + (sig(4) * cos(thet&
     &a) + sig(5) * sin(theta)) * cos(theta) / 0.3D1) * sig(6) - 0.4D1 /&
     & 0.3D1 * ((sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta)&
     & / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(thet&
     &a) / 0.3D1) * sig(6) - 0.2D1 * (-sig(4) * sin(theta) + sig(5) * co&
     &s(theta)) * sig(6) * sin(theta) + 0.2D1 * (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) * sig(6) * cos(theta)) / rho0)
      J(6,5) = -rho0 * ((sig(2) * cos(theta) + sig(3) * sin(theta)) * mu&
     &1 * (cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * s&
     &ig(6) + sin(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) *&
     & sig(6)) / rho0 + sin(theta) * mu1 * (0.2D1 / 0.3D1 * (-sig(1) ** &
     &2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta))&
     & ** 2 / 0.3D1 + sig(1) ** 2 * sin(theta) ** 2 / 0.6D1 + (-sig(2) *&
     & sin(theta) + sig(3) * cos(theta)) ** 2 / 0.6D1 + (-sig(4) * sin(t&
     &heta) + sig(5) * cos(theta)) ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) *&
     & sig(6) + 0.2D1 / 0.3D1 * (-sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 &
     &- (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-si&
     &g(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 + sig(1) ** &
     &2 * cos(theta) ** 2 / 0.6D1 + (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) ** 2 / 0.6D1 + (sig(4) * cos(theta) + sig(5) * sin(theta))&
     & ** 2 / 0.6D1 + sig(6) ** 2 / 0.6D1) * sig(6) - 0.4D1 / 0.3D1 * (-&
     &sig(6) ** 2 / 0.3D1 + sig(1) ** 2 * cos(theta) ** 2 / 0.6D1 + (sig&
     &(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.6D1 + (sig(4) * c&
     &os(theta) + sig(5) * sin(theta)) ** 2 / 0.6D1 + sig(1) ** 2 * sin(&
     &theta) ** 2 / 0.6D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta))&
     & ** 2 / 0.6D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 &
     &/ 0.6D1) * sig(6) + (-sig(4) * sin(theta) + sig(5) * cos(theta)) *&
     &* 2 * sig(6) + (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 * &
     &sig(6)) / rho0 + (sig(4) * cos(theta) + sig(5) * sin(theta)) * mu1&
     & * (0.2D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) * sin(theta) + cos(theta) * (-sig(4) * sin(theta) &
     &+ sig(5) * cos(theta)) / 0.3D1) * sig(6) + 0.2D1 / 0.3D1 * (-0.2D1&
     & / 0.3D1 * cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta&
     &)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) / 0.&
     &3D1) * sig(6) - 0.4D1 / 0.3D1 * ((sig(4) * cos(theta) + sig(5) * s&
     &in(theta)) * sin(theta) / 0.3D1 + cos(theta) * (-sig(4) * sin(thet&
     &a) + sig(5) * cos(theta)) / 0.3D1) * sig(6) + 0.2D1 * (-sig(4) * s&
     &in(theta) + sig(5) * cos(theta)) * sig(6) * cos(theta) + 0.2D1 * (&
     &sig(4) * cos(theta) + sig(5) * sin(theta)) * sig(6) * sin(theta)) &
     &/ rho0)
      J(6,6) = -rho0 * (sig(1) * cos(theta) * mu1 * (-(-sig(4) * sin(the&
     &ta) + sig(5) * cos(theta)) * sig(1) * sin(theta) + (sig(4) * cos(t&
     &heta) + sig(5) * sin(theta)) * sig(1) * cos(theta)) / rho0 + (sig(&
     &2) * cos(theta) + sig(3) * sin(theta)) * mu1 * ((-sig(4) * sin(the&
     &ta) + sig(5) * cos(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(&
     &theta)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * (sig(2) * &
     &cos(theta) + sig(3) * sin(theta))) / rho0 + (sig(4) * cos(theta) +&
     & sig(5) * sin(theta)) * mu1 * (0.2D1 * sig(6) ** 2 - sig(1) ** 2 *&
     & cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(the&
     &ta)) ** 2 / 0.3D1 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) &
     &* sin(theta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig&
     &(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 + 0.2D1 / 0.3&
     &D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2) / rho0)               
                end select
    end subroutine A2sigmaJacobian
    ! *******************************************************************************************************
    real function LEfriction_mu(x,time,v,disp, Friction)
        implicit none
        real    :: x(3),time, v(3),disp
        !real    :: VV, LL, mu_s,mu_d
        type(tfriction) :: Friction
        ! Selfsimilar-rupture
        !VV=2000 ! m/s
        !LL=500!250 ! m
        !mu_s=0.5
        !mu_d=0.05!25
        
        LEfriction_mu=min(Friction%mu_s,max(Friction%mu_d,Friction%mu_s-(Friction%mu_s-Friction%mu_d)*(Friction%V*(time+Friction%t0)-abs(x(1)))/Friction%L))
        !LEfriction_mu=max(mu_d,mu_s-(mu_s-mu_d)*(VV*(time+0.01)-abs(x(1)))/LL)
        !if(LEfriction_mu.ge. mu_s) then
        !        LEfriction_mu=1.e+6
        !end if
        !if(VV*(time+0.5)-abs(x(1))>0) then
        !    LEfriction_mu=mu_d!0.5*(mu_s+mu_d)   
        !else
        !    LEfriction_mu=mu_s
        !end if

    end function LEfriction_mu
    
    real function mu2tauFL(mu_f, lambda, mu, rho,detA)
        implicit none
        real :: mu_f, lambda, mu, rho, detA
        real :: cs,cs_f
        cs=sqrt(mu/rho)
        cs_f=sqrt(2.0)
        mu2tauFL=mu_f*6.0/rho
        mu2tauFL=6.0/(detA*rho*cs_f**2)*mu_f
        !mu2tauFL=mu2tauFL/2.0
    end function mu2tauFL
   
    subroutine ComputeGPRLE_EA(EA,Q)
        use Parameters, only : nVar
        implicit none
        real, intent(out) :: EA(3,3)
        real, intent(in) :: Q(nVar)
        real                :: Aij(3,3),mu1,rho0
        
        ! Recontruct the matrix A
        Aij(1,:) = (/ Q( 5), Q( 6), Q( 7) /) 
        Aij(2,:) = (/ Q( 8), Q( 9), Q(10) /)
        Aij(3,:) = (/ Q(11), Q(12), Q(13) /)
        mu1=Q(16)
        rho0=Q(1)
select case(GPRLEversion)
    case(0)        
        EA(1,1) = 0.1D1 / rho0 * ((Aij(2,2) * Aij(3,3) - Aij(2,3) * Aij(3,&
     &2)) * (Aij(1,1) ** 2 + Aij(2,1) ** 2 + Aij(3,1) ** 2) * mu1 * (-0.&
     &2D1 / 0.3D1 * Aij(1,1) ** 2 - 0.2D1 / 0.3D1 * Aij(2,1) ** 2 - 0.2D&
     &1 / 0.3D1 * Aij(3,1) ** 2 + Aij(1,2) ** 2 / 0.3D1 + Aij(2,2) ** 2 &
     &/ 0.3D1 + Aij(3,2) ** 2 / 0.3D1 + Aij(1,3) ** 2 / 0.3D1 + Aij(2,3)&
     & ** 2 / 0.3D1 + Aij(3,3) ** 2 / 0.3D1) - 0.2D1 * (Aij(2,1) * Aij(3&
     &,3) - Aij(2,3) * Aij(3,1)) * (Aij(1,1) * Aij(1,2) + Aij(2,1) * Aij&
     &(2,2) + Aij(3,1) * Aij(3,2) + 0.1D1) * mu1 * (-Aij(1,1) * Aij(1,2)&
     & / 0.2D1 - Aij(2,1) * Aij(2,2) / 0.2D1 - Aij(3,1) * Aij(3,2) / 0.2&
     &D1) + 0.2D1 * (Aij(2,1) * Aij(3,2) - Aij(2,2) * Aij(3,1)) * (Aij(1&
     &,1) * Aij(1,3) + Aij(2,1) * Aij(2,3) + Aij(3,1) * Aij(3,3) + 0.1D1&
     &) * mu1 * (-Aij(1,1) * Aij(1,3) / 0.2D1 - Aij(2,1) * Aij(2,3) / 0.&
     &2D1 - Aij(3,1) * Aij(3,3) / 0.2D1))
      EA(1,2) = 0.1D1 / rho0 * (0.2D1 * (Aij(2,2) * Aij(3,3) - Aij(2,3) &
     &* Aij(3,2)) * (Aij(1,1) * Aij(1,2) + Aij(2,1) * Aij(2,2) + Aij(3,1&
     &) * Aij(3,2) + 0.1D1) * mu1 * (-Aij(1,1) * Aij(1,2) / 0.2D1 - Aij(&
     &2,1) * Aij(2,2) / 0.2D1 - Aij(3,1) * Aij(3,2) / 0.2D1) - (Aij(2,1)&
     & * Aij(3,3) - Aij(2,3) * Aij(3,1)) * (Aij(1,2) ** 2 + Aij(2,2) ** &
     &2 + Aij(3,2) ** 2) * mu1 * (Aij(1,1) ** 2 / 0.3D1 + Aij(2,1) ** 2 &
     &/ 0.3D1 + Aij(3,1) ** 2 / 0.3D1 - 0.2D1 / 0.3D1 * Aij(1,2) ** 2 - &
     &0.2D1 / 0.3D1 * Aij(2,2) ** 2 - 0.2D1 / 0.3D1 * Aij(3,2) ** 2 + Ai&
     &j(1,3) ** 2 / 0.3D1 + Aij(2,3) ** 2 / 0.3D1 + Aij(3,3) ** 2 / 0.3D&
     &1) + 0.2D1 * (Aij(2,1) * Aij(3,2) - Aij(2,2) * Aij(3,1)) * (Aij(1,&
     &2) * Aij(1,3) + Aij(2,2) * Aij(2,3) + Aij(3,2) * Aij(3,3) + 0.1D1)&
     & * mu1 * (-Aij(1,2) * Aij(1,3) / 0.2D1 - Aij(2,2) * Aij(2,3) / 0.2&
     &D1 - Aij(3,2) * Aij(3,3) / 0.2D1))
      EA(1,3) = 0.1D1 / rho0 * (0.2D1 * (Aij(2,2) * Aij(3,3) - Aij(2,3) &
     &* Aij(3,2)) * (Aij(1,1) * Aij(1,3) + Aij(2,1) * Aij(2,3) + Aij(3,1&
     &) * Aij(3,3) + 0.1D1) * mu1 * (-Aij(1,1) * Aij(1,3) / 0.2D1 - Aij(&
     &2,1) * Aij(2,3) / 0.2D1 - Aij(3,1) * Aij(3,3) / 0.2D1) - 0.2D1 * (&
     &Aij(2,1) * Aij(3,3) - Aij(2,3) * Aij(3,1)) * (Aij(1,2) * Aij(1,3) &
     &+ Aij(2,2) * Aij(2,3) + Aij(3,2) * Aij(3,3) + 0.1D1) * mu1 * (-Aij&
     &(1,2) * Aij(1,3) / 0.2D1 - Aij(2,2) * Aij(2,3) / 0.2D1 - Aij(3,2) &
     &* Aij(3,3) / 0.2D1) + (Aij(2,1) * Aij(3,2) - Aij(2,2) * Aij(3,1)) &
     &* (Aij(1,3) ** 2 + Aij(2,3) ** 2 + Aij(3,3) ** 2) * mu1 * (Aij(1,1&
     &) ** 2 / 0.3D1 + Aij(2,1) ** 2 / 0.3D1 + Aij(3,1) ** 2 / 0.3D1 + A&
     &ij(1,2) ** 2 / 0.3D1 + Aij(2,2) ** 2 / 0.3D1 + Aij(3,2) ** 2 / 0.3&
     &D1 - 0.2D1 / 0.3D1 * Aij(1,3) ** 2 - 0.2D1 / 0.3D1 * Aij(2,3) ** 2&
     & - 0.2D1 / 0.3D1 * Aij(3,3) ** 2))
      EA(2,1) = 0.1D1 / rho0 * (-(Aij(1,2) * Aij(3,3) - Aij(1,3) * Aij(3&
     &,2)) * (Aij(1,1) ** 2 + Aij(2,1) ** 2 + Aij(3,1) ** 2) * mu1 * (-0&
     &.2D1 / 0.3D1 * Aij(1,1) ** 2 - 0.2D1 / 0.3D1 * Aij(2,1) ** 2 - 0.2&
     &D1 / 0.3D1 * Aij(3,1) ** 2 + Aij(1,2) ** 2 / 0.3D1 + Aij(2,2) ** 2&
     & / 0.3D1 + Aij(3,2) ** 2 / 0.3D1 + Aij(1,3) ** 2 / 0.3D1 + Aij(2,3&
     &) ** 2 / 0.3D1 + Aij(3,3) ** 2 / 0.3D1) + 0.2D1 * (Aij(1,1) * Aij(&
     &3,3) - Aij(1,3) * Aij(3,1)) * (Aij(1,1) * Aij(1,2) + Aij(2,1) * Ai&
     &j(2,2) + Aij(3,1) * Aij(3,2) + 0.1D1) * mu1 * (-Aij(1,1) * Aij(1,2&
     &) / 0.2D1 - Aij(2,1) * Aij(2,2) / 0.2D1 - Aij(3,1) * Aij(3,2) / 0.&
     &2D1) - 0.2D1 * (Aij(1,1) * Aij(3,2) - Aij(1,2) * Aij(3,1)) * (Aij(&
     &1,1) * Aij(1,3) + Aij(2,1) * Aij(2,3) + Aij(3,1) * Aij(3,3) + 0.1D&
     &1) * mu1 * (-Aij(1,1) * Aij(1,3) / 0.2D1 - Aij(2,1) * Aij(2,3) / 0&
     &.2D1 - Aij(3,1) * Aij(3,3) / 0.2D1))
      EA(2,2) = 0.1D1 / rho0 * (-0.2D1 * (Aij(1,2) * Aij(3,3) - Aij(1,3)&
     & * Aij(3,2)) * (Aij(1,1) * Aij(1,2) + Aij(2,1) * Aij(2,2) + Aij(3,&
     &1) * Aij(3,2) + 0.1D1) * mu1 * (-Aij(1,1) * Aij(1,2) / 0.2D1 - Aij&
     &(2,1) * Aij(2,2) / 0.2D1 - Aij(3,1) * Aij(3,2) / 0.2D1) + (Aij(1,1&
     &) * Aij(3,3) - Aij(1,3) * Aij(3,1)) * (Aij(1,2) ** 2 + Aij(2,2) **&
     & 2 + Aij(3,2) ** 2) * mu1 * (Aij(1,1) ** 2 / 0.3D1 + Aij(2,1) ** 2&
     & / 0.3D1 + Aij(3,1) ** 2 / 0.3D1 - 0.2D1 / 0.3D1 * Aij(1,2) ** 2 -&
     & 0.2D1 / 0.3D1 * Aij(2,2) ** 2 - 0.2D1 / 0.3D1 * Aij(3,2) ** 2 + A&
     &ij(1,3) ** 2 / 0.3D1 + Aij(2,3) ** 2 / 0.3D1 + Aij(3,3) ** 2 / 0.3&
     &D1) - 0.2D1 * (Aij(1,1) * Aij(3,2) - Aij(1,2) * Aij(3,1)) * (Aij(1&
     &,2) * Aij(1,3) + Aij(2,2) * Aij(2,3) + Aij(3,2) * Aij(3,3) + 0.1D1&
     &) * mu1 * (-Aij(1,2) * Aij(1,3) / 0.2D1 - Aij(2,2) * Aij(2,3) / 0.&
     &2D1 - Aij(3,2) * Aij(3,3) / 0.2D1))
      EA(2,3) = 0.1D1 / rho0 * (-0.2D1 * (Aij(1,2) * Aij(3,3) - Aij(1,3)&
     & * Aij(3,2)) * (Aij(1,1) * Aij(1,3) + Aij(2,1) * Aij(2,3) + Aij(3,&
     &1) * Aij(3,3) + 0.1D1) * mu1 * (-Aij(1,1) * Aij(1,3) / 0.2D1 - Aij&
     &(2,1) * Aij(2,3) / 0.2D1 - Aij(3,1) * Aij(3,3) / 0.2D1) + 0.2D1 * &
     &(Aij(1,1) * Aij(3,3) - Aij(1,3) * Aij(3,1)) * (Aij(1,2) * Aij(1,3)&
     & + Aij(2,2) * Aij(2,3) + Aij(3,2) * Aij(3,3) + 0.1D1) * mu1 * (-Ai&
     &j(1,2) * Aij(1,3) / 0.2D1 - Aij(2,2) * Aij(2,3) / 0.2D1 - Aij(3,2)&
     & * Aij(3,3) / 0.2D1) - (Aij(1,1) * Aij(3,2) - Aij(1,2) * Aij(3,1))&
     & * (Aij(1,3) ** 2 + Aij(2,3) ** 2 + Aij(3,3) ** 2) * mu1 * (Aij(1,&
     &1) ** 2 / 0.3D1 + Aij(2,1) ** 2 / 0.3D1 + Aij(3,1) ** 2 / 0.3D1 + &
     &Aij(1,2) ** 2 / 0.3D1 + Aij(2,2) ** 2 / 0.3D1 + Aij(3,2) ** 2 / 0.&
     &3D1 - 0.2D1 / 0.3D1 * Aij(1,3) ** 2 - 0.2D1 / 0.3D1 * Aij(2,3) ** &
     &2 - 0.2D1 / 0.3D1 * Aij(3,3) ** 2))
      EA(3,1) = 0.1D1 / rho0 * ((Aij(1,2) * Aij(2,3) - Aij(1,3) * Aij(2,&
     &2)) * (Aij(1,1) ** 2 + Aij(2,1) ** 2 + Aij(3,1) ** 2) * mu1 * (-0.&
     &2D1 / 0.3D1 * Aij(1,1) ** 2 - 0.2D1 / 0.3D1 * Aij(2,1) ** 2 - 0.2D&
     &1 / 0.3D1 * Aij(3,1) ** 2 + Aij(1,2) ** 2 / 0.3D1 + Aij(2,2) ** 2 &
     &/ 0.3D1 + Aij(3,2) ** 2 / 0.3D1 + Aij(1,3) ** 2 / 0.3D1 + Aij(2,3)&
     & ** 2 / 0.3D1 + Aij(3,3) ** 2 / 0.3D1) - 0.2D1 * (Aij(1,1) * Aij(2&
     &,3) - Aij(1,3) * Aij(2,1)) * (Aij(1,1) * Aij(1,2) + Aij(2,1) * Aij&
     &(2,2) + Aij(3,1) * Aij(3,2) + 0.1D1) * mu1 * (-Aij(1,1) * Aij(1,2)&
     & / 0.2D1 - Aij(2,1) * Aij(2,2) / 0.2D1 - Aij(3,1) * Aij(3,2) / 0.2&
     &D1) + 0.2D1 * (Aij(1,1) * Aij(2,2) - Aij(1,2) * Aij(2,1)) * (Aij(1&
     &,1) * Aij(1,3) + Aij(2,1) * Aij(2,3) + Aij(3,1) * Aij(3,3) + 0.1D1&
     &) * mu1 * (-Aij(1,1) * Aij(1,3) / 0.2D1 - Aij(2,1) * Aij(2,3) / 0.&
     &2D1 - Aij(3,1) * Aij(3,3) / 0.2D1))
      EA(3,2) = 0.1D1 / rho0 * (0.2D1 * (Aij(1,2) * Aij(2,3) - Aij(1,3) &
     &* Aij(2,2)) * (Aij(1,1) * Aij(1,2) + Aij(2,1) * Aij(2,2) + Aij(3,1&
     &) * Aij(3,2) + 0.1D1) * mu1 * (-Aij(1,1) * Aij(1,2) / 0.2D1 - Aij(&
     &2,1) * Aij(2,2) / 0.2D1 - Aij(3,1) * Aij(3,2) / 0.2D1) - (Aij(1,1)&
     & * Aij(2,3) - Aij(1,3) * Aij(2,1)) * (Aij(1,2) ** 2 + Aij(2,2) ** &
     &2 + Aij(3,2) ** 2) * mu1 * (Aij(1,1) ** 2 / 0.3D1 + Aij(2,1) ** 2 &
     &/ 0.3D1 + Aij(3,1) ** 2 / 0.3D1 - 0.2D1 / 0.3D1 * Aij(1,2) ** 2 - &
     &0.2D1 / 0.3D1 * Aij(2,2) ** 2 - 0.2D1 / 0.3D1 * Aij(3,2) ** 2 + Ai&
     &j(1,3) ** 2 / 0.3D1 + Aij(2,3) ** 2 / 0.3D1 + Aij(3,3) ** 2 / 0.3D&
     &1) + 0.2D1 * (Aij(1,1) * Aij(2,2) - Aij(1,2) * Aij(2,1)) * (Aij(1,&
     &2) * Aij(1,3) + Aij(2,2) * Aij(2,3) + Aij(3,2) * Aij(3,3) + 0.1D1)&
     & * mu1 * (-Aij(1,2) * Aij(1,3) / 0.2D1 - Aij(2,2) * Aij(2,3) / 0.2&
     &D1 - Aij(3,2) * Aij(3,3) / 0.2D1))
      EA(3,3) = 0.1D1 / rho0 * (0.2D1 * (Aij(1,2) * Aij(2,3) - Aij(1,3) &
     &* Aij(2,2)) * (Aij(1,1) * Aij(1,3) + Aij(2,1) * Aij(2,3) + Aij(3,1&
     &) * Aij(3,3) + 0.1D1) * mu1 * (-Aij(1,1) * Aij(1,3) / 0.2D1 - Aij(&
     &2,1) * Aij(2,3) / 0.2D1 - Aij(3,1) * Aij(3,3) / 0.2D1) - 0.2D1 * (&
     &Aij(1,1) * Aij(2,3) - Aij(1,3) * Aij(2,1)) * (Aij(1,2) * Aij(1,3) &
     &+ Aij(2,2) * Aij(2,3) + Aij(3,2) * Aij(3,3) + 0.1D1) * mu1 * (-Aij&
     &(1,2) * Aij(1,3) / 0.2D1 - Aij(2,2) * Aij(2,3) / 0.2D1 - Aij(3,2) &
     &* Aij(3,3) / 0.2D1) + (Aij(1,1) * Aij(2,2) - Aij(1,2) * Aij(2,1)) &
     &* (Aij(1,3) ** 2 + Aij(2,3) ** 2 + Aij(3,3) ** 2) * mu1 * (Aij(1,1&
     &) ** 2 / 0.3D1 + Aij(2,1) ** 2 / 0.3D1 + Aij(3,1) ** 2 / 0.3D1 + A&
     &ij(1,2) ** 2 / 0.3D1 + Aij(2,2) ** 2 / 0.3D1 + Aij(3,2) ** 2 / 0.3&
     &D1 - 0.2D1 / 0.3D1 * Aij(1,3) ** 2 - 0.2D1 / 0.3D1 * Aij(2,3) ** 2&
     & - 0.2D1 / 0.3D1 * Aij(3,3) ** 2))
      EA=-EA
    case(1)
EA(1,1) = mu1 * (-0.4D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) **&
     & 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / &
     &0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D&
     &1 + Q(13) ** 2 / 0.6D1) * Q(5) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3&
     &D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + &
     &Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10)&
     & ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(5) + 0.2D1 / 0.3D1 * (-Q(7&
     &) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) **&
     & 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / &
     &0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(5) - 0.2D1 * (&
     &-Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.2D1&
     &) * Q(6) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D1 - &
     &Q(11) * Q(13) / 0.2D1) * Q(7)) / rho0
      EA(1,2) = mu1 * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** &
     &2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0&
     &.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1&
     & + Q(13) ** 2 / 0.6D1) * Q(6) - 0.4D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D&
     &1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q&
     &(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) &
     &** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(6) + 0.2D1 / 0.3D1 * (-Q(7)&
     & ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** &
     &2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0&
     &.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(6) - 0.2D1 * (-&
     &Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.2D1)&
     & * Q(5) - 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q&
     &(12) * Q(13) / 0.2D1) * Q(7)) / rho0
      EA(1,3) = mu1 * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** &
     &2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0&
     &.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1&
     & + Q(13) ** 2 / 0.6D1) * Q(7) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D&
     &1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q&
     &(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) &
     &** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(7) - 0.4D1 / 0.3D1 * (-Q(7)&
     & ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** &
     &2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0&
     &.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(7) - 0.2D1 * (-&
     &Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(13) / 0.2D1&
     &) * Q(6) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D1 - &
     &Q(11) * Q(13) / 0.2D1) * Q(5)) / rho0
      EA(2,1) = mu1 * (-0.4D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) **&
     & 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / &
     &0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D&
     &1 + Q(13) ** 2 / 0.6D1) * Q(8) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3&
     &D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + &
     &Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10)&
     & ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(8) + 0.2D1 / 0.3D1 * (-Q(7&
     &) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) **&
     & 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / &
     &0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(8) - 0.2D1 * (&
     &-Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.2D1&
     &) * Q(9) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D1 - &
     &Q(11) * Q(13) / 0.2D1) * Q(10)) / rho0
      EA(2,2) = mu1 * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** &
     &2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0&
     &.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1&
     & + Q(13) ** 2 / 0.6D1) * Q(9) - 0.4D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D&
     &1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q&
     &(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) &
     &** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(9) + 0.2D1 / 0.3D1 * (-Q(7)&
     & ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** &
     &2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0&
     &.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(9) - 0.2D1 * (-&
     &Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.2D1)&
     & * Q(8) - 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q&
     &(12) * Q(13) / 0.2D1) * Q(10)) / rho0
      EA(2,3) = mu1 * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** &
     &2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0&
     &.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1&
     & + Q(13) ** 2 / 0.6D1) * Q(10) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3&
     &D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + &
     &Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10)&
     & ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(10) - 0.4D1 / 0.3D1 * (-Q(&
     &7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) *&
     &* 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 /&
     & 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(10) - 0.2D1 *&
     & (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(13) / 0.&
     &2D1) * Q(9) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D1&
     & - Q(11) * Q(13) / 0.2D1) * Q(8)) / rho0
      EA(3,1) = mu1 * (-0.4D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) **&
     & 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / &
     &0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D&
     &1 + Q(13) ** 2 / 0.6D1) * Q(11) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.&
     &3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 +&
     & Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10&
     &) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(11) + 0.2D1 / 0.3D1 * (-Q&
     &(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) &
     &** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 &
     &/ 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(11) - 0.2D1 &
     &* (-Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.&
     &2D1) * Q(12) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D&
     &1 - Q(11) * Q(13) / 0.2D1) * Q(13)) / rho0
      EA(3,2) = mu1 * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** &
     &2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0&
     &.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1&
     & + Q(13) ** 2 / 0.6D1) * Q(12) - 0.4D1 / 0.3D1 * (-Q(6) ** 2 / 0.3&
     &D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + &
     &Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10)&
     & ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(12) + 0.2D1 / 0.3D1 * (-Q(&
     &7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) *&
     &* 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 /&
     & 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(12) - 0.2D1 *&
     & (-Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.2&
     &D1) * Q(11) - 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1&
     & - Q(12) * Q(13) / 0.2D1) * Q(13)) / rho0
      EA(3,3) = mu1 * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** &
     &2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0&
     &.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1&
     & + Q(13) ** 2 / 0.6D1) * Q(13) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3&
     &D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + &
     &Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10)&
     & ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(13) - 0.4D1 / 0.3D1 * (-Q(&
     &7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) *&
     &* 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 /&
     & 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(13) - 0.2D1 *&
     & (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(13) / 0.&
     &2D1) * Q(12) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D&
     &1 - Q(11) * Q(13) / 0.2D1) * Q(11)) / rho0 
    case(2)
        EA(1,1) = mu1 * (-0.4D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) **&
     & 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / &
     &0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D&
     &1 + Q(13) ** 2 / 0.6D1) * Q(5) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3&
     &D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + &
     &Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10)&
     & ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(5) + 0.2D1 / 0.3D1 * (-Q(7&
     &) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) **&
     & 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / &
     &0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(5) - 0.2D1 * (&
     &-Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.2D1&
     &) * Q(6) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D1 - &
     &Q(11) * Q(13) / 0.2D1) * Q(7)) / rho0
      EA(1,2) = mu1 * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** &
     &2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0&
     &.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1&
     & + Q(13) ** 2 / 0.6D1) * Q(6) - 0.4D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D&
     &1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q&
     &(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) &
     &** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(6) + 0.2D1 / 0.3D1 * (-Q(7)&
     & ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** &
     &2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0&
     &.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(6) - 0.2D1 * (-&
     &Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.2D1)&
     & * Q(5) - 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q&
     &(12) * Q(13) / 0.2D1) * Q(7)) / rho0
      EA(1,3) = mu1 * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** &
     &2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0&
     &.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1&
     & + Q(13) ** 2 / 0.6D1) * Q(7) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D&
     &1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q&
     &(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) &
     &** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(7) - 0.4D1 / 0.3D1 * (-Q(7)&
     & ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** &
     &2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0&
     &.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(7) - 0.2D1 * (-&
     &Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(13) / 0.2D1&
     &) * Q(6) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D1 - &
     &Q(11) * Q(13) / 0.2D1) * Q(5)) / rho0
      EA(2,1) = mu1 * (-0.4D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) **&
     & 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / &
     &0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D&
     &1 + Q(13) ** 2 / 0.6D1) * Q(8) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3&
     &D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + &
     &Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10)&
     & ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(8) + 0.2D1 / 0.3D1 * (-Q(7&
     &) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) **&
     & 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / &
     &0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(8) - 0.2D1 * (&
     &-Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.2D1&
     &) * Q(9) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D1 - &
     &Q(11) * Q(13) / 0.2D1) * Q(10)) / rho0
      EA(2,2) = mu1 * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** &
     &2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0&
     &.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1&
     & + Q(13) ** 2 / 0.6D1) * Q(9) - 0.4D1 / 0.3D1 * (-Q(6) ** 2 / 0.3D&
     &1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + Q&
     &(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) &
     &** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(9) + 0.2D1 / 0.3D1 * (-Q(7)&
     & ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) ** &
     &2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 / 0&
     &.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(9) - 0.2D1 * (-&
     &Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.2D1)&
     & * Q(8) - 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q&
     &(12) * Q(13) / 0.2D1) * Q(10)) / rho0
      EA(2,3) = mu1 * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** &
     &2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0&
     &.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1&
     & + Q(13) ** 2 / 0.6D1) * Q(10) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3&
     &D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + &
     &Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10)&
     & ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(10) - 0.4D1 / 0.3D1 * (-Q(&
     &7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) *&
     &* 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 /&
     & 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(10) - 0.2D1 *&
     & (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(13) / 0.&
     &2D1) * Q(9) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D1&
     & - Q(11) * Q(13) / 0.2D1) * Q(8)) / rho0
      EA(3,1) = mu1 * (-0.4D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) **&
     & 2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / &
     &0.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D&
     &1 + Q(13) ** 2 / 0.6D1) * Q(11) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.&
     &3D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 +&
     & Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10&
     &) ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(11) + 0.2D1 / 0.3D1 * (-Q&
     &(7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) &
     &** 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 &
     &/ 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(11) - 0.2D1 &
     &* (-Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.&
     &2D1) * Q(12) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D&
     &1 - Q(11) * Q(13) / 0.2D1) * Q(13)) / rho0
      EA(3,2) = mu1 * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** &
     &2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0&
     &.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1&
     & + Q(13) ** 2 / 0.6D1) * Q(12) - 0.4D1 / 0.3D1 * (-Q(6) ** 2 / 0.3&
     &D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + &
     &Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10)&
     & ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(12) + 0.2D1 / 0.3D1 * (-Q(&
     &7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) *&
     &* 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 /&
     & 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(12) - 0.2D1 *&
     & (-Q(5) * Q(6) / 0.2D1 - Q(8) * Q(9) / 0.2D1 - Q(11) * Q(12) / 0.2&
     &D1) * Q(11) - 0.2D1 * (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1&
     & - Q(12) * Q(13) / 0.2D1) * Q(13)) / rho0
      EA(3,3) = mu1 * (0.2D1 / 0.3D1 * (-Q(5) ** 2 / 0.3D1 - Q(8) ** &
     &2 / 0.3D1 - Q(11) ** 2 / 0.3D1 + Q(6) ** 2 / 0.6D1 + Q(9) ** 2 / 0&
     &.6D1 + Q(12) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10) ** 2 / 0.6D1&
     & + Q(13) ** 2 / 0.6D1) * Q(13) + 0.2D1 / 0.3D1 * (-Q(6) ** 2 / 0.3&
     &D1 - Q(9) ** 2 / 0.3D1 - Q(12) ** 2 / 0.3D1 + Q(5) ** 2 / 0.6D1 + &
     &Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(7) ** 2 / 0.6D1 + Q(10)&
     & ** 2 / 0.6D1 + Q(13) ** 2 / 0.6D1) * Q(13) - 0.4D1 / 0.3D1 * (-Q(&
     &7) ** 2 / 0.3D1 - Q(10) ** 2 / 0.3D1 - Q(13) ** 2 / 0.3D1 + Q(5) *&
     &* 2 / 0.6D1 + Q(8) ** 2 / 0.6D1 + Q(11) ** 2 / 0.6D1 + Q(6) ** 2 /&
     & 0.6D1 + Q(9) ** 2 / 0.6D1 + Q(12) ** 2 / 0.6D1) * Q(13) - 0.2D1 *&
     & (-Q(6) * Q(7) / 0.2D1 - Q(9) * Q(10) / 0.2D1 - Q(12) * Q(13) / 0.&
     &2D1) * Q(12) - 0.2D1 * (-Q(5) * Q(7) / 0.2D1 - Q(8) * Q(10) / 0.2D&
     &1 - Q(11) * Q(13) / 0.2D1) * Q(11)) / rho0
end select
    
    end subroutine ComputeGPRLE_EA
   
    subroutine StaticLimiterEQ99(dmpresult,lxb,ldx)
        use Parameters, only : ICType
        implicit none
        logical :: dmpresult
        real    :: lxb(3),ldx(3)
        SELECT CASE(ICType)
            CASE('SSCRACK')
                if(abs(lxb(1))<10000.0+ldx(1) .and. abs(lxb(2))<ldx(2)) then
                    dmpresult = .FALSE.
                end if
            CASE('StiffInclusion')
                ! Manually activate the limiter on the strong material jumps
                if(     ( (abs(lxb(1)-0.5)<(ldx(1)) .or. abs(lxb(1)+0.5)<(ldx(1)))  .and. (lxb(2)>-0.1-(ldx(2)) .and. lxb(2)<0.1+(ldx(2)))) ) then
                     dmpresult = .FALSE.
                end if
                if(     ( (abs(lxb(2)-0.1)<(ldx(2)) .or. abs(lxb(2)+0.1)<(ldx(2)))  .and. (lxb(1)>-0.5-(ldx(1)) .and. lxb(1)<0.5+(ldx(1)))) ) then
                     dmpresult = .FALSE.
                end if
        END SELECT
        
    end subroutine StaticLimiterEQ99
	
SUBROUTINE LinSolve(N,A,b,x)
  IMPLICIT NONE
  INTEGER       :: N
  REAL          :: A(N,N), b(N), x(N)  
  !
  INTEGER       :: i,j,ml(1) 
  REAL          :: temp(N+1),piv 
  REAL          :: C(N+1,N)  
  !
  C(1:N,:)     = TRANSPOSE(A(:,:))
  C(N+1,:)     = b 
  !    
  ! Forward elimination with column pivoting 
  ! 
  DO i = 1, N
     ! If pivot element is zero, then swap rows 
     ml = MAXLOC(ABS(C(i,i:N))) 
     j = i - 1 + ml(1) 
     temp   = C(:,j) 
     C(:,j) = C(:,i)
     C(:,i) = temp      
     IF(C(i,i).EQ.0.) THEN
        !PRINT *, 'ERROR. Matrix is singular!'
        !DO j = 1, N
        !   PRINT *, A(j,:) 
        !ENDDO
        !STOP
x = 0. 
RETURN 
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
  x = C(N+1,:)  
  !
END SUBROUTINE LinSolve
	
	
    end module SpecificVarEqn99