void SWE::MySWESolver_ADERDG_kernels::aderdg::gemm_4_2_2_rhs_x(const double* A, const double* B, double* C) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r12\n\t"
                       "movq $0, %%r13\n\t"
                       "movq $0, %%r14\n\t"
                       "33:\n\t"
                       "addq $2, %%r13\n\t"
                       "movq $0, %%r12\n\t"
                       "34:\n\t"
                       "addq $4, %%r12\n\t"
                       "vmovapd 0(%%rdx), %%ymm14\n\t"
                       "vmovapd 32(%%rdx), %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm2\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm2\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovapd %%ymm14, 0(%%rdx)\n\t"
                       "vmovapd %%ymm15, 32(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $32, %%rdi\n\t"
                       "cmpq $4, %%r12\n\t"
                       "jl 34b\n\t"
                       "addq $32, %%rdx\n\t"
                       "addq $64, %%rsi\n\t"
                       "subq $32, %%rdi\n\t"
                       "cmpq $2, %%r13\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 32;
#endif
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::gemm_4_2_2_rhs_y(const double* A, const double* B, double* C) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r12\n\t"
                       "movq $0, %%r13\n\t"
                       "movq $0, %%r14\n\t"
                       "33:\n\t"
                       "addq $2, %%r13\n\t"
                       "movq $0, %%r12\n\t"
                       "34:\n\t"
                       "addq $4, %%r12\n\t"
                       "vmovapd 0(%%rdx), %%ymm14\n\t"
                       "vmovapd 64(%%rdx), %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm2\n\t"
                       "addq $64, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm2\n\t"
                       "addq $64, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovapd %%ymm14, 0(%%rdx)\n\t"
                       "vmovapd %%ymm15, 64(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $96, %%rdi\n\t"
                       "cmpq $4, %%r12\n\t"
                       "jl 34b\n\t"
                       "addq $96, %%rdx\n\t"
                       "addq $64, %%rsi\n\t"
                       "subq $32, %%rdi\n\t"
                       "cmpq $2, %%r13\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 32;
#endif
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::gemm_4_2_2_lduh_x(const double* A, const double* B, double* C) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r12\n\t"
                       "movq $0, %%r13\n\t"
                       "movq $0, %%r14\n\t"
                       "33:\n\t"
                       "addq $2, %%r13\n\t"
                       "movq $0, %%r12\n\t"
                       "34:\n\t"
                       "addq $4, %%r12\n\t"
                       "vmovupd 0(%%rdx), %%ymm14\n\t"
                       "vmovupd 32(%%rdx), %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm2\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm2\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovupd %%ymm14, 0(%%rdx)\n\t"
                       "vmovupd %%ymm15, 32(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $32, %%rdi\n\t"
                       "cmpq $4, %%r12\n\t"
                       "jl 34b\n\t"
                       "addq $32, %%rdx\n\t"
                       "addq $64, %%rsi\n\t"
                       "subq $32, %%rdi\n\t"
                       "cmpq $2, %%r13\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 32;
#endif
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::gemm_4_2_2_lduh_y(const double* A, const double* B, double* C) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r12\n\t"
                       "movq $0, %%r13\n\t"
                       "movq $0, %%r14\n\t"
                       "33:\n\t"
                       "addq $2, %%r13\n\t"
                       "movq $0, %%r12\n\t"
                       "34:\n\t"
                       "addq $4, %%r12\n\t"
                       "vmovupd 0(%%rdx), %%ymm14\n\t"
                       "vmovupd 64(%%rdx), %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm2\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm2\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovupd %%ymm14, 0(%%rdx)\n\t"
                       "vmovupd %%ymm15, 64(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $32, %%rdi\n\t"
                       "cmpq $4, %%r12\n\t"
                       "jl 34b\n\t"
                       "addq $96, %%rdx\n\t"
                       "addq $64, %%rsi\n\t"
                       "subq $32, %%rdi\n\t"
                       "cmpq $2, %%r13\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 32;
#endif
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::gemm_4_2_2_gradQ_x(const double* A, const double* B, double* C) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r12\n\t"
                       "movq $0, %%r13\n\t"
                       "movq $0, %%r14\n\t"
                       "33:\n\t"
                       "addq $2, %%r13\n\t"
                       "movq $0, %%r12\n\t"
                       "34:\n\t"
                       "addq $4, %%r12\n\t"
                       "vmovapd 0(%%rdx), %%ymm14\n\t"
                       "vmovapd 32(%%rdx), %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm2\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm2\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovapd %%ymm14, 0(%%rdx)\n\t"
                       "vmovapd %%ymm15, 32(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $32, %%rdi\n\t"
                       "cmpq $4, %%r12\n\t"
                       "jl 34b\n\t"
                       "addq $32, %%rdx\n\t"
                       "addq $64, %%rsi\n\t"
                       "subq $32, %%rdi\n\t"
                       "cmpq $2, %%r13\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 32;
#endif
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::gemm_4_2_2_gradQ_y(const double* A, const double* B, double* C) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r12\n\t"
                       "movq $0, %%r13\n\t"
                       "movq $0, %%r14\n\t"
                       "33:\n\t"
                       "addq $2, %%r13\n\t"
                       "movq $0, %%r12\n\t"
                       "34:\n\t"
                       "addq $4, %%r12\n\t"
                       "vmovapd 0(%%rdx), %%ymm14\n\t"
                       "vmovapd 64(%%rdx), %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm2\n\t"
                       "addq $64, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm2\n\t"
                       "addq $64, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovapd %%ymm14, 0(%%rdx)\n\t"
                       "vmovapd %%ymm15, 64(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $96, %%rdi\n\t"
                       "cmpq $4, %%r12\n\t"
                       "jl 34b\n\t"
                       "addq $96, %%rdx\n\t"
                       "addq $64, %%rsi\n\t"
                       "subq $32, %%rdi\n\t"
                       "cmpq $2, %%r13\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 32;
#endif
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::gemm_4_2_2_lqi(const double* A, const double* B, double* C) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r12\n\t"
                       "movq $0, %%r13\n\t"
                       "movq $0, %%r14\n\t"
                       "33:\n\t"
                       "addq $2, %%r13\n\t"
                       "movq $0, %%r12\n\t"
                       "34:\n\t"
                       "addq $4, %%r12\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm2\n\t"
                       "addq $128, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm2\n\t"
                       "addq $128, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovapd %%ymm14, 0(%%rdx)\n\t"
                       "vmovapd %%ymm15, 32(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $224, %%rdi\n\t"
                       "cmpq $4, %%r12\n\t"
                       "jl 34b\n\t"
                       "addq $32, %%rdx\n\t"
                       "addq $64, %%rsi\n\t"
                       "subq $32, %%rdi\n\t"
                       "cmpq $2, %%r13\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 32;
#endif
}

