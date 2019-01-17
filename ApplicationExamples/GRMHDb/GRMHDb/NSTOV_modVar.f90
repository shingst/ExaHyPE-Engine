#define RNSTOV
    
 MODULE NSTOV_modVar
    IMPLICIT NONE 
#if defined(RNSTOV) 
    REAL, DIMENSION (:), POINTER :: t,void1,void2,void3
    REAL, DIMENSION (:,:), POINTER :: y
    REAL, DIMENSION (:,:,:), POINTER :: yDG
    !
#endif     
 END MODULE NSTOV_modVar
    
        
     