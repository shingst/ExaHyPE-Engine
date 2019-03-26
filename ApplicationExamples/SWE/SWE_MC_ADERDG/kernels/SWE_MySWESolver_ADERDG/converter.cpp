
#include "kernels/SWE_MySWESolver_ADERDG/converter.h"


void SWE::MySWESolver_ADERDG_kernels::aderdg::converter::Qhbnd_optimised2generic(const double* const opt, double* gen) {
  for(int d=0; d<4; d++) {
    for(int ij=0; ij<2; ij++) {
      #pragma omp simd aligned(gen,opt:ALIGNMENT)
      for(int n=0; n<4; n++) {
        gen[n+4*(ij+2*d)] = opt[n+4*(ij+2*d)];
      }
    }    
  }
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::converter::Qhbnd_generic2optimised(const double* const gen, double* opt) {
  for(int d=0; d<4; d++) {
    for(int ij=0; ij<2; ij++) {
      #pragma omp simd aligned(gen,opt:ALIGNMENT)
      for(int n=0; n<4; n++) {
        opt[n+4*(ij+2*d)] = gen[n+4*(ij+2*d)];
      }
    }    
  }
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::converter::Fhbnd_optimised2generic(const double* const opt, double* gen) {
  for(int d=0; d<4; d++) {
    for(int ij=0; ij<2; ij++) {
      #pragma omp simd aligned(gen,opt:ALIGNMENT)
      for(int n=0; n<4; n++) {
        gen[n+4*(ij+2*d)] = opt[n+4*(ij+2*d)];
      }
    }    
  }
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::converter::Fhbnd_generic2optimised(const double* const gen, double* opt) {
  for(int d=0; d<4; d++) {
    for(int ij=0; ij<2; ij++) {
      #pragma omp simd aligned(gen,opt:ALIGNMENT)
      for(int n=0; n<4; n++) {
        opt[n+4*(ij+2*d)] = gen[n+4*(ij+2*d)];
      }
    }    
  }
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::converter::QFace_optimised2generic(const double* const opt, double* gen) {
  for(int ij=0; ij<2; ij++) {
    #pragma omp simd aligned(opt:ALIGNMENT) //gen not always aligned if nor face0
    for(int n=0; n<4; n++) {
      gen[n+4*ij] = opt[n+4*ij];
    }
  }    
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::converter::QFace_generic2optimised(const double* const gen, double* opt) {
  for(int ij=0; ij<2; ij++) {
    #pragma omp simd aligned(opt:ALIGNMENT) //gen not always aligned if nor face0
    for(int n=0; n<4; n++) {
      opt[n+4*ij] = gen[n+4*ij];
    }
  }    
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::converter::FFace_optimised2generic(const double* const opt, double* gen) {
  for(int ij=0; ij<2; ij++) {
    #pragma omp simd aligned(opt:ALIGNMENT) //gen not always aligned if nor face0
    for(int n=0; n<4; n++) {
      gen[n+4*ij] = opt[n+4*ij];
    }
  }    
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::converter::FFace_generic2optimised(const double* const gen, double* opt) {
  for(int ij=0; ij<2; ij++) {
    #pragma omp simd aligned(opt:ALIGNMENT) //gen not always aligned if nor face0
    for(int n=0; n<4; n++) {
      opt[n+4*ij] = gen[n+4*ij];
    }
  }    
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::converter::Qhi_optimised2generic(const double* const opt, double* gen) {
  for(int j=0; j<4; j++) {
    #pragma omp simd aligned(gen,opt:ALIGNMENT)
    for(int i=0; i<4; i++) {
      gen[i+4*j] = opt[i+4*j];
    }
  }
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::converter::Qhi_generic2optimised(const double* const gen, double* opt) {
  for(int j=0; j<4; j++) {
    #pragma omp simd aligned(gen,opt:ALIGNMENT)
    for(int i=0; i<4; i++) {
      opt[i+4*j] = gen[i+4*j];
    }
  }
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::converter::Fhi_optimised2generic(const double* const opt, double* gen) {
  for(int b=0; b<3; b++) {
    for(int j=0; j<4; j++) {
      #pragma omp simd aligned(gen,opt:ALIGNMENT)
      for(int i=0; i<4; i++) {
        gen[i+4*j+b*16] = opt[i+4*j+b*16];
      }
    }
  }
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::converter::Fhi_generic2optimised(const double* const gen, double* opt) {
  for(int b=0; b<3; b++) {
    for(int j=0; j<4; j++) {
      #pragma omp simd aligned(gen,opt:ALIGNMENT)
      for(int i=0; i<4; i++) {
        opt[i+4*j+b*16] = gen[i+4*j+b*16];
      }
    }
  }
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::converter::PSi_optimised2generic(const double* const opt, double* gen) {
  for(int b=0; b<3; b++) {
    for(int j=0; j<4; j++) {
      #pragma omp simd aligned(gen,opt:ALIGNMENT)
      for(int i=0; i<4; i++) {
        gen[i+4*j+b*16] = opt[i+4*j+b*16];
      }
    }
  }
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::converter::PSi_generic2optimised(const double* const gen, double* opt) {
  for(int b=0; b<3; b++) {
    for(int j=0; j<4; j++) {
      #pragma omp simd aligned(gen,opt:ALIGNMENT)
      for(int i=0; i<4; i++) {
        opt[i+4*j+b*16] = gen[i+4*j+b*16];
      }
    }
  }
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::converter::Lduh_optimised2generic(const double* const opt, double* gen) {
  for(int zyx=0; zyx < 4; zyx++){
    #pragma omp simd aligned(gen,opt:ALIGNMENT)
    for(int n=0; n<4;n++){
      gen[zyx*4+n] = opt[zyx*4+n];
    }
  }
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::converter::Lduh_generic2optimised(const double* const gen, double* opt) {
  for(int zyx=0; zyx < 4; zyx++){
    #pragma omp simd aligned(gen,opt:ALIGNMENT)
    for(int n=0; n<4;n++){
      opt[zyx*4+n] = gen[zyx*4+n];
    }
  }
}