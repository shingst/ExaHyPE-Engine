
! AdjustSolution enforcer

RECURSIVE  SUBROUTINE EnforceCCZ4Constraints(Q)
	USE typesDef, ONLY : nVar 
	IMPLICIT NONE
	! Argument list
	REAL, INTENT(INOUT)  :: Q(nVar)              ! spatial degrees of freedom 
	! Local variables
	INTEGER :: i,j,k,l,iVar,iDim, iter
	REAL    :: g_cov(3,3), det, g_contr(3,3), Aex(3,3), traceA, phi  
	
	g_cov(1,1) = Q(1)
	g_cov(1,2) = Q(2)
	g_cov(1,3) = Q(3)
	g_cov(2,1) = Q(2)
	g_cov(2,2) = Q(4)
	g_cov(2,3) = Q(5)
	g_cov(3,1) = Q(3)
	g_cov(3,2) = Q(5)
	g_cov(3,3) = Q(6)
	! This determinant should be close to unity, since we use the conformal decomposition 
	det = (Q(1)*Q(4)*Q(6)-Q(1)*Q(5)**2-Q(2)**2*Q(6)+2*Q(2)*Q(3)*Q(5)-Q(3)**2*Q(4)) 
	g_contr(1,1) = (Q(4)*Q(6)-Q(5)**2)   / det 
	g_contr(1,2) = -(Q(2)*Q(6)-Q(3)*Q(5))/ det 
	g_contr(1,3) = (Q(2)*Q(5)-Q(3)*Q(4)) / det 
	g_contr(2,1) = -(Q(2)*Q(6)-Q(3)*Q(5))/ det 
	g_contr(2,2) = (Q(1)*Q(6)-Q(3)**2)   / det 
	g_contr(2,3) = -(Q(1)*Q(5)-Q(2)*Q(3))/ det 
	g_contr(3,1) = (Q(2)*Q(5)-Q(3)*Q(4)) / det 
	g_contr(3,2) = -(Q(1)*Q(5)-Q(2)*Q(3))/ det 
	g_contr(3,3) = (Q(1)*Q(4)-Q(2)**2)   / det 
	!
	phi = det**(-1./6.) 
	g_cov = phi**2*g_cov
	det = (g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)*g_cov(3,2)-g_cov(2,1)*g_cov(1,2)*g_cov(3,3)+g_cov(2,1)*g_cov(1,3)*g_cov(3,2)+g_cov(3,1)*g_cov(1,2)*g_cov(2,3)-g_cov(3,1)*g_cov(1,3)*g_cov(2,2)) 
	g_contr(1,1) =  (g_cov(2,2)*g_cov(3,3)-g_cov(2,3)*g_cov(3,2)) / det 
	g_contr(1,2) = -(g_cov(1,2)*g_cov(3,3)-g_cov(1,3)*g_cov(3,2)) / det
	g_contr(1,3) = -(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/ det 
	g_contr(2,1) = -(g_cov(2,1)*g_cov(3,3)-g_cov(2,3)*g_cov(3,1)) / det 
	g_contr(2,2) = (g_cov(1,1)*g_cov(3,3)-g_cov(1,3)*g_cov(3,1))  / det 
	g_contr(2,3) = -(g_cov(1,1)*g_cov(2,3)-g_cov(1,3)*g_cov(2,1)) / det 
	g_contr(3,1) = -(-g_cov(2,1)*g_cov(3,2)+g_cov(2,2)*g_cov(3,1))/ det 
	g_contr(3,2) = -(g_cov(1,1)*g_cov(3,2)-g_cov(1,2)*g_cov(3,1)) / det 
	g_contr(3,3) = (g_cov(1,1)*g_cov(2,2)-g_cov(1,2)*g_cov(2,1))  / det 
	!
	Aex(1,1) = Q(7) 
	Aex(1,2) = Q(8) 
	Aex(1,3) = Q(9) 
	Aex(2,1) = Q(8) 
	Aex(2,2) = Q(10) 
	Aex(2,3) = Q(11) 
	Aex(3,1) = Q(9) 
	Aex(3,2) = Q(11) 
	Aex(3,3) = Q(12)             
	!
	traceA = SUM(g_contr*Aex) 
	!
	Aex = Aex - 1./3.*g_cov*traceA 
	!
	Q( 1) = g_cov(1,1) 
	Q( 2) = g_cov(1,2) 
	Q( 3) = g_cov(1,3) 
	Q( 4) = g_cov(2,2) 
	Q( 5) = g_cov(2,3) 
	Q( 6) = g_cov(3,3) 
	!
	Q( 7) = Aex(1,1) 
	Q( 8) = Aex(1,2) 
	Q( 9) = Aex(1,3) 
	Q(10) = Aex(2,2) 
	Q(11) = Aex(2,3) 
	Q(12) = Aex(3,3)             
	! 
END SUBROUTINE EnforceCCZ4Constraints 
