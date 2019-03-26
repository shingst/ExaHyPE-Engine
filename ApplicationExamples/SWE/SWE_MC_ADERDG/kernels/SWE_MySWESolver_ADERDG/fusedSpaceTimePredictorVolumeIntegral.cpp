
#include <cstring>
#include <algorithm>

#include "kernels/SWE_MySWESolver_ADERDG/Kernels.h"
#include "kernels/SWE_MySWESolver_ADERDG/DGMatrices.h"
#include "kernels/SWE_MySWESolver_ADERDG/Quadrature.h"
#include "kernels/SWE_MySWESolver_ADERDG/gemmsCPP.h"

#include "MySWESolver_ADERDG.h"


int SWE::MySWESolver_ADERDG_kernels::aderdg::fusedSpaceTimePredictorVolumeIntegral(
        SWE::MySWESolver_ADERDG& solver, 
        double* restrict lduh,
        double* restrict lQhbnd, 
        double* restrict lFhbnd,
        double* restrict lQi,
        double* restrict rhs,
        double* restrict lFi,
        double* restrict lSi,   // for NCP or Source
        double* restrict lQhi,
        double* restrict lFhi,
        double* restrict lShi,  // for NCP or Source
        double* restrict gradQ, // for NCP or Source
        const double* const restrict luh,
        const double inverseDx, //Assume dx[0] == dx[1] == dx[2]
        const double dt
) {


  //********************
  //****** Picard ******
  //********************

#ifdef __INTEL_COMPILER
  __assume_aligned(lQi, ALIGNMENT);
  __assume_aligned(rhs, ALIGNMENT);
  __assume_aligned(lFi, ALIGNMENT);
  __assume_aligned(FLCoeff, ALIGNMENT); // == F0
  __assume_aligned(Kxi, ALIGNMENT);
  __assume_aligned(iK1_T, ALIGNMENT);
  __assume_aligned(weights3, ALIGNMENT);
  __assume_aligned(iweights3, ALIGNMENT);
  __assume_aligned(luh, ALIGNMENT); //luh should be aligned, see Solver.h
  __assume_aligned(lSi, ALIGNMENT);
  __assume_aligned(gradQ, ALIGNMENT);
#endif

  // 0. Allocate local variable

  double new_lQi_slice[8] __attribute__((aligned(ALIGNMENT))); //for step 4 (computing new lQi value), doesn't update parameters
  const double dtBydx = inverseDx * dt; //Assume dx[0] == dx[1] == dx[2]
  double dudxT_by_dx[8] __attribute__((aligned(ALIGNMENT)));
  
  // 0. precompute 1/dx * dudx_T. Assume dx[0] == dx[1] == dx[2]
  #pragma omp simd aligned(dudxT_by_dx,dudx_T:ALIGNMENT)
  for(int it=0;it<8;it++) {
    dudxT_by_dx[it] = inverseDx * dudx_T[it];
  }
#if defined(USE_IPO) && ! defined(UNSAFE_IPO)
  volatile double doNotOptimizeAway_dudx_by_dt = dudxT_by_dx[0]; //used to prevent the compiler from optimizing temp array away. Needs to be volatile
#endif   

  // 1. Trivial initial guess
  for (int t = 0; t < 2; t++) {
    for (int xyz = 0; xyz < 4; xyz++) {
      std::copy_n(&luh[4*xyz], 4, &lQi[4*(xyz+4*t)]);
    }
  }
  
  // 2. Discrete Picard iterations
  constexpr int MaxIterations = 5;  
  int iter = 0;
  for (; iter < MaxIterations; iter++) {
    for (int t = 0; t < 2; t++) {  // time DOF

      { // Compute the fluxes
        double* F[2];
        for (int xyz = 0; xyz < 4; xyz++) {
          // Call PDE fluxes
          F[0] = &lFi[(t*4+xyz)*4];
          F[1] = &lFi[(t*4+xyz)*4+32];
          #ifdef USE_IPO
              #pragma forceinline recursive
          #endif
          solver.SWE::MySWESolver_ADERDG::flux(lQi+(t*4+xyz)*4, F);
        }
      }

      // Compute the contribution of the initial condition uh to the right-hand side (rhs)
      for (int xyz = 0; xyz < 4; xyz++) {
        const double weight = weights3[xyz] * FLCoeff[t];
        #pragma omp simd aligned(rhs,luh:ALIGNMENT)
        for (int n = 0; n < 4; n++) {
          rhs[n+4*(xyz+4*t)] = weight * luh[n+4*xyz];
        }
      }
      
      //set gradQ to 0
      std::memset(gradQ, 0, 32 * sizeof(double));
      
      // Compute the "derivatives" (contributions of the stiffness matrix)      
      // x direction (independent from the y and z derivatives)
      for (int z = 0; z < 1; z++) {
        for (int y = 0; y < 2; y++) {
          double coeffRhsX[8] __attribute__((aligned(ALIGNMENT)));
          #pragma omp simd aligned(coeffRhsX,Kxi:ALIGNMENT)
          for (int it = 0; it < 8; it++) {
            coeffRhsX[it] = - weights3[t*2+z*2+y] * dtBydx * Kxi[it];
          }
          #if defined(USE_IPO) && !defined(UNSAFE_IPO)
          volatile double doNotOptimizeAway_coeffRhsX = coeffRhsX[0]; //used to prevent the compiler from optimizing temp array away. Needs to be volatile
          #endif
          #ifdef USE_IPO
          #pragma forceinline
          #endif
          gemm_4_2_2_rhs_x(lFi+((t*1+z)*2+y)*8, coeffRhsX, rhs+((t*1+z)*2+y)*8);
          #ifdef USE_IPO
          #pragma forceinline
          #endif
          gemm_4_2_2_gradQ_x(lQi+((t*1+z)*2+y)*8, dudxT_by_dx, gradQ+(z*2+y)*8);
        }
      }
      
      // y direction (independent from the x and z derivatives)
      for (int z = 0; z < 1; z++) {
        for (int x = 0; x < 2; x++) {
          double coeffRhsY[8] __attribute__((aligned(ALIGNMENT)));
          #pragma omp simd aligned(coeffRhsY,Kxi:ALIGNMENT)
          for (int it = 0; it < 8; it++) {
            coeffRhsY[it] = - weights3[t*2+z*2+x] * dtBydx * Kxi[it];
          }
          #if defined(USE_IPO) && !defined(UNSAFE_IPO)
          volatile double doNotOptimizeAway_coeffRhsY = coeffRhsY[0]; //used to prevent the compiler from optimizing temp array away. Needs to be volatile
          #endif
          #ifdef USE_IPO
          #pragma forceinline
          #endif
          gemm_4_2_2_rhs_y(lFi+((t*1+z)*4+x)*4+32, coeffRhsY, rhs+((t*1+z)*4+x)*4);
          #ifdef USE_IPO
          #pragma forceinline
          #endif
          gemm_4_2_2_gradQ_y(lQi+((t*1+z)*4+x)*4, dudxT_by_dx, gradQ+(z*4+x)*4+16);
        }
      }
       

      // Compute the Nonconservative part NCP + Source
      double gradQ_PDE[8] __attribute__((aligned(ALIGNMENT)));
      double tmp_ncp_output[4] __attribute__((aligned(ALIGNMENT))) = {0.}; //initialize for padding
      for(int xyz = 0; xyz < 4; xyz++) {
        // Remove padding to use the same user function as generic kernel
        double gradQ_PDE[8]; 
        std::copy_n(gradQ+4*xyz, 4, gradQ_PDE); //x
        std::copy_n(gradQ+4*xyz+16, 4, gradQ_PDE+4); //y
        
        // NCP
        #ifdef USE_IPO
          #pragma forceinline recursive
        #endif
        solver.SWE::MySWESolver_ADERDG::nonConservativeProduct(lQi+(t*4+xyz)*4, gradQ_PDE, tmp_ncp_output);
        #pragma omp simd aligned(lSi,tmp_ncp_output:ALIGNMENT)
        for(int n = 0; n<4; n++) {
          lSi[n+(t*4+xyz)*4] = - tmp_ncp_output[n];  }

        // Update rhs
        const double updateSize = weights1[t] * weights3[xyz] * dt;
        #pragma omp simd aligned(rhs,lSi:ALIGNMENT)
        for (int n = 0; n < 4; n++) {
          rhs[n+(t*4+xyz)*4] += updateSize * lSi[n+(t*4+xyz)*4];
        }
      }

    }  // end time dof

    // 3. Multiply with (K1)^(-1) to get the discrete time integral of the
    // discrete Picard iteration
    double sq_res = 0.0;
    for (int xyz = 0; xyz < 4; xyz++) {
      double s_m_QSlice[8] __attribute__((aligned(ALIGNMENT)));
      #pragma omp simd aligned(s_m_QSlice,iK1_T:ALIGNMENT)
      for (int it = 0; it < 8; it++) {
        s_m_QSlice[it] = iweights3[xyz] * iK1_T[it];
      }
      #if defined(USE_IPO) && !defined(UNSAFE_IPO)
      volatile double doNotOptimizeAway_s_m_QSlice = s_m_QSlice[0]; //used to prevent the compiler from optimizing temp array away. Needs to be volatile
      #endif
      #ifdef USE_IPO
      #pragma forceinline
      #endif
      gemm_4_2_2_lqi(rhs+4*xyz, s_m_QSlice, new_lQi_slice);
      for(int t = 0; t < 2; t++) {
        for(int n=0; n<4; n++) { //only copy and change the variables, skip parameters
          sq_res += (new_lQi_slice[n+4*t] - lQi[n+4*(xyz+4*t)]) * (new_lQi_slice[n+4*t] - lQi[n+4*(xyz+4*t)]);
          lQi[n+4*(xyz+4*t)] = new_lQi_slice[n+4*t];
        }
      }
    }

    // 4. Exit condition
    constexpr double tol2 = 1e-7 * 1e-7;
    if (sq_res < tol2) {
      break;
    }
  }  // end iter

  //***********************
  //****** Predictor ******
  //***********************
  
#ifdef __INTEL_COMPILER
  __assume_aligned(lQi, ALIGNMENT);
  __assume_aligned(lQhi, ALIGNMENT);
  __assume_aligned(lFi, ALIGNMENT);
  __assume_aligned(lFhi, ALIGNMENT);
  __assume_aligned(weights1, ALIGNMENT);
  __assume_aligned(lSi, ALIGNMENT);
  __assume_aligned(lShi, ALIGNMENT);
#endif  

  std::memset(lQhi, 0, 16 * sizeof(double));
  std::memset(lFhi, 0, 32 * sizeof(double));
  std::memset(lShi, 0, 16 * sizeof(double));

  for (int z=0; z<1; z++) {
    for (int y=0; y<2; y++) {
      for (int x=0; x<2; x++) {
        
        // Matrix-Vector Products
        for (int t=0; t<2; t++) {
          #pragma omp simd aligned(lQhi,lQi:ALIGNMENT)
          for (int n=0; n<4; n++) {
            // Fortran: lQhi(:,x,y,z) = lQi(:,:,x,y,z) * wGPN(:)
            lQhi[((z*2+y)*2+x)*4+n] += weights1[t] *
                lQi[(((t*1+z)*2+y)*2+x)*4+n];
          }
          #pragma omp simd aligned(lFhi,lFi:ALIGNMENT)
          for (int n=0; n<4; n++) {
            // Fortran: lFhi_x(:,x,y,z) = lFh(:,1,x,y,z,:) * wGPN(:)
            lFhi[((z*2+y)*2+x)*4+n+0] += weights1[t] *
                lFi[(((t*1+z)*2+y)*2+x)*4+n+0];
          }  
          #pragma omp simd aligned(lFhi,lFi:ALIGNMENT)
          for (int n=0; n<4; n++) {
            // Fortran: lFhi_y(:,y,x,z) = lFh(:,2,:x,y,z,:) * wGPN(:)
            lFhi[((z*2+x)*2+y)*4+n+16] += weights1[t] *
                lFi[(((t*1+z)*2+y)*2+x)*4+n+32];
          }  
            
          #pragma omp simd aligned(lShi,lSi:ALIGNMENT)
          for (int n=0; n<4; n++) {
            // Fortran: lFhi_S(:,x,y,z) = lSh(:,x,y,z,:) * wGPN(:)
            lShi[((z*2+y)*2+x)*4+n] += weights1[t] *
              lSi[(((t*1+z)*2+y)*2+x)*4+n];
          }
        }
      
      }
    }
  }
  
  //**************************
  //****** Extrapolator ******
  //**************************
  
#ifdef __INTEL_COMPILER
  __assume_aligned(lQhi, ALIGNMENT);
  __assume_aligned(lQhbnd, ALIGNMENT);
  __assume_aligned(lFhi, ALIGNMENT);
  __assume_aligned(lFhbnd, ALIGNMENT);
  __assume_aligned(FRCoeff, ALIGNMENT);
  __assume_aligned(FLCoeff, ALIGNMENT);
#endif
  
  std::memset(lQhbnd, 0, 32 * sizeof(double));
  std::memset(lFhbnd, 0, 32 * sizeof(double));

  // x-direction: face 1 (left) and face 2 (right)
  for (int yz = 0; yz < 2; yz++) {
    // Matrix-Vector Products
    for (int x = 0; x < 2; x++) {
      #pragma omp simd aligned(lQhbnd,lQhi:ALIGNMENT)
      for (int n = 0; n < 4; n++) {    
        // Fortran: lQhbnd(:,j,i,1) = lQhi(:,:,j,i) * FLCoeff(:)
        lQhbnd[n+4*yz+0] +=
            lQhi[n+4*(x+2*yz)] * FLCoeff[x];

        // Fortran: lQhbnd(:,j,i,2) = lQhi(:,:,j,i) * FRCoeff(:)
        lQhbnd[n+4*yz+8] +=
            lQhi[n+4*(x+2*yz)] * FRCoeff[x];
        // Fortran: lFhbnd(:,j,i,1) = lFhi_x(:,:,j,i) * FLCoeff(:)
        lFhbnd[n+4*yz+0] +=
            lFhi[n+4*(x+2*yz)] * FLCoeff[x];

        // Fortran: lFhbnd(:,j,i,2) = lFhi_x(:,:,j,i) * FRCoeff(:)
        lFhbnd[n+4*yz+8] +=
            lFhi[n+4*(x+2*yz)] * FRCoeff[x];
      }
    }
  }

  // y-direction: face 3 (left) and face 4 (right)
  for (int xz = 0; xz < 2; xz++) {  
    // Matrix-Vector Products
    for (int y = 0; y < 2; y++) {
      const int z = 0;
      const int x = xz;
      #pragma omp simd aligned(lQhbnd,lQhi:ALIGNMENT)
      for (int n = 0; n < 4; n++) {
        // Fortran: lQhbnd(:,j,i,3) = lQhi(:,j,:,i) * FLCoeff(:)
        lQhbnd[n+4*xz+16] +=
            lQhi[n+4*(x+2*(y+1*z))] * FLCoeff[y];

        // Fortran: lQhbnd(:,j,i,4) = lQhi(:,j,:,i) * FRCoeff(:)
        lQhbnd[n+4*xz+24] +=
            lQhi[n+4*(x+2*(y+1*z))] * FRCoeff[y];
        // Fortran: lFhbnd(:,j,i,3) = lFhi_y(:,:,j,i) * FLCoeff(:)
        lFhbnd[n+4*xz+16] +=
            lFhi[n+4*(y+2*xz)+16] * FLCoeff[y];

        // Fortran: lFhbnd(:,j,i,4) = lFhi_y(:,:,j,i) * FRCoeff(:)
        lFhbnd[n+4*xz+24] +=
            lFhi[n+4*(y+2*xz)+16] * FRCoeff[y];
      }
    }
  }

  

  
  //*****************************
  //****** Volume Integral ******
  //*****************************
  

  memset(lduh, 0, 16*sizeof(double));

#ifdef __INTEL_COMPILER
  __assume_aligned(lFhi,     ALIGNMENT);
  __assume_aligned(Kxi_T,    ALIGNMENT);
  __assume_aligned(weights2, ALIGNMENT);
  __assume_aligned(lduh,     ALIGNMENT); //lduh should be aligned, see Solver.h
  __assume_aligned(weights3, ALIGNMENT);
  __assume_aligned(lShi,     ALIGNMENT);
#endif
  
  // Assume equispaced mesh, dx[0] == dx[1] == dx[2]
  for (int j=0; j<1; j++) {
    for (int i=0; i<2; i++) {
      
      //x, also define coefficient matrix coeffVolume
      double coeffVolume[8] __attribute__((aligned(ALIGNMENT)));
      #pragma omp simd aligned(coeffVolume,Kxi_T:ALIGNMENT)
      for (int it = 0; it < 8; it++) {
        coeffVolume[it] = weights2[i+j*2] * inverseDx * Kxi_T[it];
      }
      #if defined(USE_IPO) && !defined(UNSAFE_IPO)
      volatile double doNotOptimizeAway_coeffVolume = coeffVolume[0]; //used to prevent the compiler from optimizing temp array away. Needs to be volatile
      #endif
      #ifdef USE_IPO
      #pragma forceinline
      #endif
      gemm_4_2_2_lduh_x(lFhi+(j*2+i)*8+0, coeffVolume, lduh+(j*2+i)*8);

      //y, reuse coeffVolume
      #ifdef USE_IPO
      #pragma forceinline
      #endif
      gemm_4_2_2_lduh_y(lFhi+(j*2+i)*8+16, coeffVolume, lduh+(j*4+i)*4);

    }
  }
  // source
  for (int xyz = 0; xyz < 4; xyz++) {
    // Fortran: lduh(:,k,j,i) += w * lShi(:,k,j,i)
    #pragma omp simd aligned(lduh,lShi:ALIGNMENT)
    for (int n = 0; n < 4; n++) {
      lduh[xyz*4+n] += weights3[xyz] * lShi[xyz*4+n];
    }
  }

  return std::min(iter+1, MaxIterations); //return number of Picard iterations, min to avoid doing a +1 if the loop wasn't exited early
}