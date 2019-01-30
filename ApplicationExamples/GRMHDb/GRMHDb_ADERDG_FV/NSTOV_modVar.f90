#define RNSTOV
    
 MODULE NSTOV_modVar
    IMPLICIT NONE 
    PUBLIC
#if defined(RNSTOV) 
    REAL, DIMENSION (:), ALLOCATABLE :: t
    REAL, DIMENSION (:,:), ALLOCATABLE :: y
    !
#endif     
 END MODULE NSTOV_modVar
    
        
     