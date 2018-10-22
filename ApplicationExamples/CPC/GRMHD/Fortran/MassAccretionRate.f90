RECURSIVE SUBROUTINE MassAccretionRate(Q,masschange,vx,vy,vz)
  USE Parameters, ONLY: nVar, nDim
  IMPLICIT NONE
  ! Argument list declaration
  REAL :: Q(nVar), vx,vy,vz, masschange
  INTENT(IN)  :: Q
  INTENT(OUT) :: masschange
 
  

  ! Local variable declaration
  REAL :: v2,b2,e2,lf,w,ww,uem,gamma1
  REAL :: lapse, g_contr(3,3), g_cov(3,3), gp
  REAl :: shift(3),vf(3), vf_cov(3), ur_contr(3),V(nVar)
  INTEGER :: i, iErr
  REAl :: vr_contr, betar_contr, lapse_sphe, theta, phi,r
  
  ! TODO: Alejandro

  CALL PDECons2Prim(V,Q,iErr)

  !gammaij = V(14:19)
  g_cov(1,1) = V(14)
  g_cov(1,2) = V(15)
  g_cov(1,3) = V(16)
  g_cov(2,2) = V(17)
  g_cov(2,3) = V(18)
  g_cov(3,3) = V(19)
  g_cov(2,1) = V(15)
  g_cov(3,1) = V(16)
  g_cov(3,2) = V(18)

  ! Compute the inverse of the 3-metric
  CALL MatrixInverse3x3(g_cov,g_contr,gp)

  vf_cov = V(2:4)
  ! v_i --> v^i From covariant to contravariant
  vf     = MATMUL(g_contr,vf_cov)
  vx=vf(1)
  vy=vf(2)
  vz=vf(3)
  
  gp = SQRT(gp)
  lapse = V(10)
  shift = V(11:13)

  ! Computing the vector transformation

!  theta = ATAN2( sqrt(x(1)*x(1)+x(2)*x(2)), x(3) )
!  phi   = ATAN2( x(2), x(1) )


  ! Quantities in spherical simmetry
!  r= SQRT(x(1)*x(1) +  x(2)*x(2) +  x(3)*x(3))
!  vr_contr = SIN(theta)*COS(phi)*vf(1) + SIN(theta)*SIN(phi)*vf(2) + COS(phi)*vf(2)
!  betar_contr = 2.0/r
!  betar_contr = betar_contr /(1.0+2.0/r)
!  lapse_sphe  = SQRT(1.0/(1.0+2.0/r))

  ! r-component of cuadri-velocity
  ! u^i = v^i - beta^i/alpha
  
!  ur_contr = vr_contr  - betar_contr/lapse_sphe

  !  masschange = lapse_sphe *gp* Q(1) * ur_contr(1)

  masschange = Q(1) !r*r*Q(0) * ur_contr(1)

    
END SUBROUTINE MassAccretionRate
