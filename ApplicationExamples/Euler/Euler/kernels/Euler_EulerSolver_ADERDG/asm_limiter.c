void Euler::EulerSolver_ADERDG_kernels::aderdg::gemm_5_3_2_dg2fv_x(const double* A, const double* B, double* C) {
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
                       "addq $3, %%r13\n\t"
                       "movq $0, %%r12\n\t"
                       "34:\n\t"
                       "addq $4, %%r12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm3\n\t"
                       "addq $40, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm3\n\t"
                       "addq $40, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm13, 0(%%rdx)\n\t"
                       "vmovapd %%ymm14, 64(%%rdx)\n\t"
                       "vmovapd %%ymm15, 128(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $48, %%rdi\n\t"
                       "cmpq $4, %%r12\n\t"
                       "jl 34b\n\t"
                       "34:\n\t"
                       "addq $1, %%r12\n\t"
                       "vxorpd %%xmm13, %%xmm13, %%xmm13\n\t"
                       "vxorpd %%xmm14, %%xmm14, %%xmm14\n\t"
                       "vxorpd %%xmm15, %%xmm15, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm3\n\t"
                       "addq $40, %%rdi\n\t"
                       "vmovsd 0(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm3, %%xmm0, %%xmm13\n\t"
                       "vmovsd 32(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm3, %%xmm1, %%xmm14\n\t"
                       "vmovsd 64(%%rsi), %%xmm2\n\t"
                       "vfmadd231sd %%xmm3, %%xmm2, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm3\n\t"
                       "addq $40, %%rdi\n\t"
                       "vmovsd 8(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm3, %%xmm0, %%xmm13\n\t"
                       "vmovsd 40(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm3, %%xmm1, %%xmm14\n\t"
                       "vmovsd 72(%%rsi), %%xmm2\n\t"
                       "vfmadd231sd %%xmm3, %%xmm2, %%xmm15\n\t"
                       "vmovsd %%xmm13, 0(%%rdx)\n\t"
                       "vmovsd %%xmm14, 64(%%rdx)\n\t"
                       "vmovsd %%xmm15, 128(%%rdx)\n\t"
                       "addq $8, %%rdx\n\t"
                       "subq $72, %%rdi\n\t"
                       "cmpq $5, %%r12\n\t"
                       "jl 34b\n\t"
                       "addq $152, %%rdx\n\t"
                       "addq $96, %%rsi\n\t"
                       "subq $40, %%rdi\n\t"
                       "cmpq $3, %%r13\n\t"
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
libxsmm_num_total_flops += 60;
#endif
}

void Euler::EulerSolver_ADERDG_kernels::aderdg::gemm_8_3_2_dg2fv_y(const double* A, const double* B, double* C) {
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
                       "addq $3, %%r13\n\t"
                       "movq $0, %%r12\n\t"
                       "34:\n\t"
                       "addq $8, %%r12\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm10, 0(%%rdx)\n\t"
                       "vmovapd %%ymm11, 32(%%rdx)\n\t"
                       "vmovapd %%ymm12, 192(%%rdx)\n\t"
                       "vmovapd %%ymm13, 224(%%rdx)\n\t"
                       "vmovapd %%ymm14, 384(%%rdx)\n\t"
                       "vmovapd %%ymm15, 416(%%rdx)\n\t"
                       "addq $64, %%rdx\n\t"
                       "subq $320, %%rdi\n\t"
                       "cmpq $8, %%r12\n\t"
                       "jl 34b\n\t"
                       "addq $512, %%rdx\n\t"
                       "addq $96, %%rsi\n\t"
                       "subq $64, %%rdi\n\t"
                       "cmpq $3, %%r13\n\t"
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
libxsmm_num_total_flops += 96;
#endif
}

void Euler::EulerSolver_ADERDG_kernels::aderdg::gemm_5_3_2_dg2fv_z(const double* A, const double* B, double* C) {
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
                       "addq $3, %%r13\n\t"
                       "movq $0, %%r12\n\t"
                       "34:\n\t"
                       "addq $4, %%r12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $576, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $576, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovupd %%ymm13, 0(%%rdx)\n\t"
                       "vmovupd %%ymm14, 360(%%rdx)\n\t"
                       "vmovupd %%ymm15, 720(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $1120, %%rdi\n\t"
                       "cmpq $4, %%r12\n\t"
                       "jl 34b\n\t"
                       "34:\n\t"
                       "addq $1, %%r12\n\t"
                       "vxorpd %%xmm13, %%xmm13, %%xmm13\n\t"
                       "vxorpd %%xmm14, %%xmm14, %%xmm14\n\t"
                       "vxorpd %%xmm15, %%xmm15, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm3\n\t"
                       "addq $576, %%rdi\n\t"
                       "vmovsd 0(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm3, %%xmm0, %%xmm13\n\t"
                       "vmovsd 32(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm3, %%xmm1, %%xmm14\n\t"
                       "vmovsd 64(%%rsi), %%xmm2\n\t"
                       "vfmadd231sd %%xmm3, %%xmm2, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm3\n\t"
                       "addq $576, %%rdi\n\t"
                       "vmovsd 8(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm3, %%xmm0, %%xmm13\n\t"
                       "vmovsd 40(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm3, %%xmm1, %%xmm14\n\t"
                       "vmovsd 72(%%rsi), %%xmm2\n\t"
                       "vfmadd231sd %%xmm3, %%xmm2, %%xmm15\n\t"
                       "vmovsd %%xmm13, 0(%%rdx)\n\t"
                       "vmovsd %%xmm14, 360(%%rdx)\n\t"
                       "vmovsd %%xmm15, 720(%%rdx)\n\t"
                       "addq $8, %%rdx\n\t"
                       "subq $1144, %%rdi\n\t"
                       "cmpq $5, %%r12\n\t"
                       "jl 34b\n\t"
                       "addq $1040, %%rdx\n\t"
                       "addq $96, %%rsi\n\t"
                       "subq $40, %%rdi\n\t"
                       "cmpq $3, %%r13\n\t"
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
libxsmm_num_total_flops += 60;
#endif
}

void Euler::EulerSolver_ADERDG_kernels::aderdg::gemm_5_2_3_fv2dg_x(const double* A, const double* B, double* C) {
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
                       "vmovupd 0(%%rdi), %%ymm2\n\t"
                       "addq $40, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm2\n\t"
                       "addq $40, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm2\n\t"
                       "addq $40, %%rdi\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovapd %%ymm14, 0(%%rdx)\n\t"
                       "vmovapd %%ymm15, 64(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $88, %%rdi\n\t"
                       "cmpq $4, %%r12\n\t"
                       "jl 34b\n\t"
                       "34:\n\t"
                       "addq $1, %%r12\n\t"
                       "vxorpd %%xmm14, %%xmm14, %%xmm14\n\t"
                       "vxorpd %%xmm15, %%xmm15, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm2\n\t"
                       "addq $40, %%rdi\n\t"
                       "vmovsd 0(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm2, %%xmm0, %%xmm14\n\t"
                       "vmovsd 32(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm2, %%xmm1, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm2\n\t"
                       "addq $40, %%rdi\n\t"
                       "vmovsd 8(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm2, %%xmm0, %%xmm14\n\t"
                       "vmovsd 40(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm2, %%xmm1, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm2\n\t"
                       "addq $40, %%rdi\n\t"
                       "vmovsd 16(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm2, %%xmm0, %%xmm14\n\t"
                       "vmovsd 48(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm2, %%xmm1, %%xmm15\n\t"
                       "vmovsd %%xmm14, 0(%%rdx)\n\t"
                       "vmovsd %%xmm15, 64(%%rdx)\n\t"
                       "addq $8, %%rdx\n\t"
                       "subq $112, %%rdi\n\t"
                       "cmpq $5, %%r12\n\t"
                       "jl 34b\n\t"
                       "addq $88, %%rdx\n\t"
                       "addq $64, %%rsi\n\t"
                       "subq $40, %%rdi\n\t"
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
libxsmm_num_total_flops += 60;
#endif
}

void Euler::EulerSolver_ADERDG_kernels::aderdg::gemm_8_2_3_fv2dg_y(const double* A, const double* B, double* C) {
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
                       "addq $8, %%r12\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm1\n\t"
                       "vmovapd 0(%%rdi), %%ymm2\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm12\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm14\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm1\n\t"
                       "vmovapd 0(%%rdi), %%ymm2\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm12\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm14\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm1\n\t"
                       "vmovapd 0(%%rdi), %%ymm2\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm12\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm14\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm15\n\t"
                       "vmovapd %%ymm12, 0(%%rdx)\n\t"
                       "vmovapd %%ymm13, 32(%%rdx)\n\t"
                       "vmovapd %%ymm14, 128(%%rdx)\n\t"
                       "vmovapd %%ymm15, 160(%%rdx)\n\t"
                       "addq $64, %%rdx\n\t"
                       "subq $320, %%rdi\n\t"
                       "cmpq $8, %%r12\n\t"
                       "jl 34b\n\t"
                       "addq $192, %%rdx\n\t"
                       "addq $64, %%rsi\n\t"
                       "subq $64, %%rdi\n\t"
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
libxsmm_num_total_flops += 96;
#endif
}

void Euler::EulerSolver_ADERDG_kernels::aderdg::gemm_5_2_3_fv2dg_z(const double* A, const double* B, double* C) {
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
                       "addq $256, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm2\n\t"
                       "addq $256, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm2\n\t"
                       "addq $256, %%rdi\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovupd %%ymm14, 0(%%rdx)\n\t"
                       "vmovupd %%ymm15, 160(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $736, %%rdi\n\t"
                       "cmpq $4, %%r12\n\t"
                       "jl 34b\n\t"
                       "34:\n\t"
                       "addq $1, %%r12\n\t"
                       "vxorpd %%xmm14, %%xmm14, %%xmm14\n\t"
                       "vxorpd %%xmm15, %%xmm15, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm2\n\t"
                       "addq $256, %%rdi\n\t"
                       "vmovsd 0(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm2, %%xmm0, %%xmm14\n\t"
                       "vmovsd 32(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm2, %%xmm1, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm2\n\t"
                       "addq $256, %%rdi\n\t"
                       "vmovsd 8(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm2, %%xmm0, %%xmm14\n\t"
                       "vmovsd 40(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm2, %%xmm1, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm2\n\t"
                       "addq $256, %%rdi\n\t"
                       "vmovsd 16(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm2, %%xmm0, %%xmm14\n\t"
                       "vmovsd 48(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm2, %%xmm1, %%xmm15\n\t"
                       "vmovsd %%xmm14, 0(%%rdx)\n\t"
                       "vmovsd %%xmm15, 160(%%rdx)\n\t"
                       "addq $8, %%rdx\n\t"
                       "subq $760, %%rdi\n\t"
                       "cmpq $5, %%r12\n\t"
                       "jl 34b\n\t"
                       "addq $280, %%rdx\n\t"
                       "addq $64, %%rsi\n\t"
                       "subq $40, %%rdi\n\t"
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
libxsmm_num_total_flops += 60;
#endif
}

void Euler::EulerSolver_ADERDG_kernels::aderdg::gemm_5_2_2_uh2lob_x(const double* A, const double* B, double* C) {
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
                       "vmovupd 0(%%rdi), %%ymm2\n\t"
                       "addq $40, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm2\n\t"
                       "addq $40, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm15\n\t"
                       "vmovapd %%ymm14, 0(%%rdx)\n\t"
                       "vmovapd %%ymm15, 64(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $48, %%rdi\n\t"
                       "cmpq $4, %%r12\n\t"
                       "jl 34b\n\t"
                       "34:\n\t"
                       "addq $1, %%r12\n\t"
                       "vxorpd %%xmm14, %%xmm14, %%xmm14\n\t"
                       "vxorpd %%xmm15, %%xmm15, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm2\n\t"
                       "addq $40, %%rdi\n\t"
                       "vmovsd 0(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm2, %%xmm0, %%xmm14\n\t"
                       "vmovsd 32(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm2, %%xmm1, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm2\n\t"
                       "addq $40, %%rdi\n\t"
                       "vmovsd 8(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm2, %%xmm0, %%xmm14\n\t"
                       "vmovsd 40(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm2, %%xmm1, %%xmm15\n\t"
                       "vmovsd %%xmm14, 0(%%rdx)\n\t"
                       "vmovsd %%xmm15, 64(%%rdx)\n\t"
                       "addq $8, %%rdx\n\t"
                       "subq $72, %%rdi\n\t"
                       "cmpq $5, %%r12\n\t"
                       "jl 34b\n\t"
                       "addq $88, %%rdx\n\t"
                       "addq $64, %%rsi\n\t"
                       "subq $40, %%rdi\n\t"
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
libxsmm_num_total_flops += 40;
#endif
}

void Euler::EulerSolver_ADERDG_kernels::aderdg::gemm_8_2_2_uh2lob_y(const double* A, const double* B, double* C) {
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
                       "addq $8, %%r12\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm1\n\t"
                       "vmovapd 0(%%rdi), %%ymm2\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm12\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm14\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm1\n\t"
                       "vmovapd 0(%%rdi), %%ymm2\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm12\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm14\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm15\n\t"
                       "vmovapd %%ymm12, 0(%%rdx)\n\t"
                       "vmovapd %%ymm13, 32(%%rdx)\n\t"
                       "vmovapd %%ymm14, 128(%%rdx)\n\t"
                       "vmovapd %%ymm15, 160(%%rdx)\n\t"
                       "addq $64, %%rdx\n\t"
                       "subq $192, %%rdi\n\t"
                       "cmpq $8, %%r12\n\t"
                       "jl 34b\n\t"
                       "addq $192, %%rdx\n\t"
                       "addq $64, %%rsi\n\t"
                       "subq $64, %%rdi\n\t"
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
libxsmm_num_total_flops += 64;
#endif
}

void Euler::EulerSolver_ADERDG_kernels::aderdg::gemm_8_2_2_uh2lob_z_slice(const double* A, const double* B, double* C) {
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
                       "addq $8, %%r12\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm1\n\t"
                       "vmovapd 0(%%rdi), %%ymm2\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm12\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm14\n\t"
                       "addq $256, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm1\n\t"
                       "vmovapd 0(%%rdi), %%ymm2\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm12\n\t"
                       "vfmadd231pd %%ymm2, %%ymm1, %%ymm14\n\t"
                       "addq $256, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm15\n\t"
                       "vmovapd %%ymm12, 0(%%rdx)\n\t"
                       "vmovapd %%ymm13, 32(%%rdx)\n\t"
                       "vmovapd %%ymm14, 64(%%rdx)\n\t"
                       "vmovapd %%ymm15, 96(%%rdx)\n\t"
                       "addq $64, %%rdx\n\t"
                       "subq $448, %%rdi\n\t"
                       "cmpq $8, %%r12\n\t"
                       "jl 34b\n\t"
                       "addq $64, %%rdx\n\t"
                       "addq $64, %%rsi\n\t"
                       "subq $64, %%rdi\n\t"
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
libxsmm_num_total_flops += 64;
#endif
}

void Euler::EulerSolver_ADERDG_kernels::aderdg::gemm_8_3_2_dg2fv_z_slice(const double* A, const double* B, double* C) {
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
                       "addq $3, %%r13\n\t"
                       "movq $0, %%r12\n\t"
                       "34:\n\t"
                       "addq $8, %%r12\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $576, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $576, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm10, 0(%%rdx)\n\t"
                       "vmovapd %%ymm11, 32(%%rdx)\n\t"
                       "vmovapd %%ymm12, 64(%%rdx)\n\t"
                       "vmovapd %%ymm13, 96(%%rdx)\n\t"
                       "vmovapd %%ymm14, 128(%%rdx)\n\t"
                       "vmovapd %%ymm15, 160(%%rdx)\n\t"
                       "addq $64, %%rdx\n\t"
                       "subq $1088, %%rdi\n\t"
                       "cmpq $8, %%r12\n\t"
                       "jl 34b\n\t"
                       "addq $128, %%rdx\n\t"
                       "addq $96, %%rsi\n\t"
                       "subq $64, %%rdi\n\t"
                       "cmpq $3, %%r13\n\t"
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
libxsmm_num_total_flops += 96;
#endif
}

