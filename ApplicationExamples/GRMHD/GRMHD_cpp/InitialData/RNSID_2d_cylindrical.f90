
SUBROUTINE RNSID_2d_Cartesian2Cylindrical(Q, V, xGP_loc)
	USE iso_c_binding
	IMPLICIT NONE
	!input variables
	REAL :: Q(19), V(19), xGP_loc(3)
	!
	INTENT(IN)  :: V   ! State vector in Cartesian coordinates
	INTENT(OUT) :: Q   ! State vector in Cylindrical coordinates
	INTENT(IN)  :: xGP_loc
	!
	REAL :: A_coord(3,3), iA_coord(3,3), iA_coord_contr(3,3), detA
	REAL :: rho, vel(3), p, Bmag(3), psi, lapse, shift(3), g_cov(3,3)
	
	! Gather the state
	rho        = V(1)
	vel        = V(2:4)
	p          = V(5)
	Bmag(1:3)  = V(6:8)
	psi        = V(9)
	lapse      = V(10)
	shift      = V(11:13)
	! Attention: Here we don't use the Trento Sym Matrix Linearization
	!    but instead the Pizza/SVEC one.
	!     We have (Pizza/C):  11 12 22 13 23 33
	!     We want (Trento/F): 11 12 13 22 23 33
	!
	!  Pizza Ordering         !   Trento ordering
	g_cov(1,1) = V(14)        !   g_cov(1,1) = V(14)
	g_cov(1,2) = V(15)        !   g_cov(1,2) = V(15)
	g_cov(2,2) = V(16)        !   g_cov(1,3) = V(16)
	g_cov(1,3) = V(17)        !   g_cov(2,2) = V(17)
	g_cov(2,3) = V(18)        !   g_cov(2,3) = V(18)
	g_cov(3,3) = V(19)        !   g_cov(3,3) = V(19)
	! Lower-left half:
	g_cov(2,1) = V(15)        !   g_cov(2,1) = V(15)
	g_cov(3,1) = V(17)        !   g_cov(3,1) = V(16)
	g_cov(3,2) = V(18)        !   g_cov(3,2) = V(18)
	
	! Compose the transformation matrix
        CALL RNSID_Cart2CylMatrix_cov(A_coord,xGP_loc)
        CALL MatrixInverse3x3(A_coord,iA_coord,detA)      ! now: iA_coord = Cyl2CartMatrix_cov    = Cart2CylMatrix_contr
        iA_coord_contr(:,:) = TRANSPOSE(A_coord(:,:))     ! Sph2CartMatrix_contr = TRANSPOSE( Cart2SphMatrix_cov )

        ! Note, if you get a crash because A_coord is singular,
        ! xGP_loc(1) == radius will be zero. At r=0, the coordinate
        ! system itself is singular.

        ! Transform the state: vector and matrices
        vel = MATMUL(A_coord,vel)
        Bmag = MATMUL(iA_coord,Bmag)   
        g_cov = MATMUL(MATMUL(A_coord, g_cov), iA_coord_contr)     ! HERE I NEED THE TRANSPOSE.

        ! Scatter the new state
	Q(1)     = rho
	Q(2:4)   = vel
	Q(5)     = p
	Q(6:8)   = Bmag
	Q(9)     = psi
	Q(10)    = Lapse
	Q(11:13) = shift
	!  Pizza Ordering         !   Trento ordering
	Q(14)    = g_cov(1,1)     !   Q(14)    = g_cov(1,1)
	Q(15)    = g_cov(1,2)     !   Q(15)    = g_cov(1,2)
	Q(16)    = g_cov(2,2)     !   Q(16)    = g_cov(1,3)
	Q(17)    = g_cov(1,3)     !   Q(17)    = g_cov(2,2)
	Q(18)    = g_cov(2,3)     !   Q(18)    = g_cov(2,3)
	Q(19)    = g_cov(3,3)     !   Q(19)    = g_cov(3,3)
END SUBROUTINE RNSID_2d_Cartesian2Cylindrical


SUBROUTINE RNSID_Cart2CylMatrix_cov(A,x)
	IMPLICIT NONE
	! output variables
	REAL, INTENT(IN)  :: x(3)
	REAL, INTENT(OUT) :: A(3,3)
	! local variables
	REAL :: theta,sintheta,costheta,cosphip,sinphip,phi,r,r2,ir2,datanAR,z
	!
	! x=r*cos(phi)
	! y=r*sin(phi)
	! z=z'
	r = x(1)
	z = x(2)
	phi = x(3)
	!
	cosphip = DCOS(phi)
	sinphip = DSIN(phi)
	! 
	! these are the columns!
	A(1:3,1) = (/ cosphip          ,  0.0  , - r*sinphip    /)      !:  GRAD (x)
	A(1:3,2) = (/ sinphip          ,  0.0  , + r*cosphip    /)      !:  GRAD (y)
	A(1:3,3) = (/ 0.0              ,  1.0   , 0.0           /)      !:  GRAD (z)
	!
	continue
!
END SUBROUTINE RNSID_Cart2CylMatrix_cov

 SUBROUTINE MatrixInverse3x3(M,iM,det) 
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
    ! 
    det = M(1,1)*M(2,2)*M(3,3)-M(1,1)*M(2,3)*M(3,2)-M(2,1)*M(1,2)*M(3,3)+M(2,1)*M(1,3)*M(3,2)+M(3,1)*M(1,2)*M(2,3)-M(3,1)*M(1,3)*M(2,2) !ComputeDet(M,3)
    IF(det*det.LT.1e-20) THEN
        print *, 'FATAL ERROR: det = 0'
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
    iM = iM/det
    !
    Id = MATMUL(M,iM)
    DO i=1,3
        DO j=1,3
            IF(i.eq.j) THEN
                IF((Id(i,j)-1.)**2..GT.1e-18) THEN
                    print *, 'FATAL ERROR: iM*M !=1'
                    STOP
                ENDIF
            ELSE
                IF((Id(i,j)**2).GT.1e-18) THEN
                    print *, 'FATAL ERROR: iM*M !=1'
                    STOP
                ENDIF
            ENDIF
        ENDDO
    ENDDO
    !
    CONTINUE
    !
END SUBROUTINE MatrixInverse3x3
