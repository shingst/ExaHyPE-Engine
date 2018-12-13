/* Simple AVX test code for very small matrix-matrix multiplication needed for the space-time predictor of ADER-DG P3 */
/* The vectorlength 4 is for AVX256 and the nVar is only for FO-CCZ4 with AVX256 padding */
/* written by Michael Dumbser, University of Trento */

#define ND3D 64   // number of 3d degrees of freedom
#define ND2D 16   // number of 2d degrees of freedom
#define ND1D 4    // number of 1d degrees of freedom
#define SOD  8 

#define GRMHD

#ifdef AVX512 

#define VECTORLENGTH 8 
#ifdef EULER
#define nVar 8  
#endif 
#ifdef CCZ4EINSTEIN
#define nVar 64 
#endif
#ifdef GRMHD 
#define nVar 24   
#endif 
#ifdef CCZ4GRMHD
#define nVar 72  
#endif
#ifdef GRGPR
#define nVar 32  
#endif 
#ifdef CCZ4GRGPR 
#define nVar 80  
#endif

#else

#define VECTORLENGTH 4 
#ifdef EULER
#define nVar 8  
#endif 
#ifdef CCZ4EINSTEIN
#define nVar 64 
#endif
#ifdef GRMHD 
#define nVar 20   
#endif 
#ifdef CCZ4GRMHD
#define nVar 68  
#endif
#ifdef GRGPR
#define nVar 32  
#endif 
#ifdef CCZ4GRGPR 
#define nVar 80  
#endif 

#endif

#include<immintrin.h>


// compute all matmuls in x direction on a 4D data set 
#ifdef AVX512

extern "C" void __cdecl vsmm4x4cs4d_x(double *cptr, const double *aptr, double b[4][4], double s[1])
{
	__assume_aligned(aptr, 64);
	__assume_aligned(cptr, 64);
	__assume_aligned(b, 64);

	// matrix b zmm registers (our small stiffness matrix). We can use the fact that the matrix has anti-symmetric entries, as detailed below. 
	// zmm0 zmm4  -zmm7 -zmm3  
	// zmm1 zmm5  -zmm6 -zmm2 
	// zmm2 zmm6  -zmm5 -zmm1 
	// zmm3 zmm7  -zmm4 -zmm0  
	// 
	// matrix a registers 
	// 
	// zmm28 zmm29 zmm30 zmm31   
	// 
	// result registers 
	// first half of the operations 
	// zmm8  zmm9  zmm10 zmm11 
	// zmm12 zmm13 zmm14 zmm15 
	// second half of the operations 
	// zmm16 zmm17 zmm18 zmm19  
	// zmm20 zmm21 zmm22 zmm23 
	 
	__asm { 

		mov rax, aptr                                       // pointer to a 
		mov rbx, b											// pointer to b 
		mov rcx, cptr										// pointer to c 
		mov r13, s											// pointer to s  
		// mov rdx, dptr									// pointer to d  ( only for debugging ) 
		 
		vbroadcastsd zmm31, QWORD PTR[r13]					// load s in vector register 
		
		prefetcht0[rax + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8 ]
		prefetcht0[rax + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 0 * nVar * 8 ]

		vbroadcastsd zmm20, QWORD PTR[rbx + 0 * 8 ]			// load b[0][0] 
		vbroadcastsd zmm21, QWORD PTR[rbx + 1 * 8 ]			// load b[0][1] 
		vbroadcastsd zmm22, QWORD PTR[rbx + 2 * 8 ]			// load b[0][2] 
		vbroadcastsd zmm23, QWORD PTR[rbx + 3 * 8 ]			// load b[0][3] 
		
		vmulpd zmm0,  zmm20, zmm31                          // multiply column of b with scalar s 
		vmulpd zmm1,  zmm21, zmm31
		vmulpd zmm2,  zmm22, zmm31
		vmulpd zmm3,  zmm23, zmm31

		prefetcht0[rax + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 0 * nVar * 8 ]
		prefetcht0[rax + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 0 * nVar * 8 ]

		vbroadcastsd zmm20, QWORD PTR[rbx + 4  * 8 ]		// load b[1][0] 
		vbroadcastsd zmm21, QWORD PTR[rbx + 5  * 8 ]		// load b[1][1] 
		vbroadcastsd zmm22, QWORD PTR[rbx + 6  * 8 ]		// load b[1][2] 
		vbroadcastsd zmm23, QWORD PTR[rbx + 7  * 8 ]		// load b[1][3] 
		
		vmulpd zmm4,  zmm20, zmm31                          // multiply column of b with scalar s 
		vmulpd zmm5,  zmm21, zmm31
		vmulpd zmm6,  zmm22, zmm31
		vmulpd zmm7,  zmm23, zmm31
		
		mov r8, 0
loopt:
		mov r10, 0
loopz:
		mov r15, 0   
loopvar: 		
			//  
			// compute address of the pointers 
			mov r14,  0  
			add r14,  r8         // time loop index 
			add r14,  r10        // z loop index  
			add r14,  r15        // variable loop index 
			// 
			mov r12, rax	     // compute a base pointer 
			add r12, r14        
			mov r13, rcx         // compute c base pointer 
			add r13, r14
			
			//
			// only one part   
			// 
			vmovapd		zmm28,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8 ]		
			vmovapd		zmm29,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 0 * nVar * 8 ]			
			vmulpd		zmm8,   zmm0,  zmm28
			vmulpd      zmm12,  zmm0,  zmm29 			
			vmovapd		zmm30,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 0 * nVar * 8 ]		
			vmovapd		zmm31,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 0 * nVar * 8 ]			
			vmulpd		zmm16,  zmm0,  zmm30 
			vmulpd      zmm20,  zmm0,  zmm31 			
			//prefetcht0[r12 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8 ]
			//prefetcht0[r12 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 1 * nVar * 8 ]			
			vmulpd		zmm9,   zmm4,  zmm28
			vmulpd      zmm13,  zmm4,  zmm29 
			//prefetcht0[r12 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 1 * nVar * 8 ]
			//prefetcht0[r12 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 1 * nVar * 8 ]			
			vmulpd		zmm17,  zmm4,  zmm30 
			vmulpd      zmm21,  zmm4,  zmm31 			
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8 ]
			vxorpd      zmm10,  zmm10, zmm10  
			vfnmadd231pd	zmm10,  zmm7,  zmm28
			vxorpd      zmm14,  zmm14, zmm14  
			vfnmadd231pd	zmm14,  zmm7,  zmm29
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 0 * nVar * 8 ]
			vxorpd      zmm18,  zmm18, zmm18  
			vfnmadd231pd	zmm18,  zmm7,  zmm30
			vxorpd      zmm22,  zmm22, zmm22  
			vfnmadd231pd	zmm22,  zmm7,  zmm31
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 0 * nVar * 8 ]
			vxorpd      zmm11,  zmm11, zmm11  
			vfnmadd231pd	zmm11,  zmm3,  zmm28
			vxorpd      zmm15,  zmm15, zmm15  
			vfnmadd231pd	zmm15,  zmm3,  zmm29
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 0 * nVar * 8 ]
			vxorpd      zmm19,  zmm19, zmm19  
			vfnmadd231pd	zmm19,  zmm3,  zmm30
			vxorpd      zmm23,  zmm23, zmm23  
			vfnmadd231pd	zmm23,  zmm3,  zmm31
			
			vmovapd		zmm28,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8 ]		
			vmovapd		zmm29,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 1 * nVar * 8 ]			
			vfmadd231pd		zmm8,   zmm1,  zmm28
			vfmadd231pd 	zmm12,  zmm1,  zmm29 
			vmovapd		zmm30,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 1 * nVar * 8 ]		
			vmovapd		zmm31,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 1 * nVar * 8 ]			
			vfmadd231pd		zmm16,  zmm1,  zmm30 
			vfmadd231pd     zmm20,  zmm1,  zmm31 
			//prefetcht0[r12 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8 ]
			//prefetcht0[r12 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 2 * nVar * 8 ]			
			vfmadd231pd		zmm9,   zmm5,  zmm28
			vfmadd231pd     zmm13,  zmm5,  zmm29 
			//prefetcht0[r12 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 2 * nVar * 8 ]
			//prefetcht0[r12 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 2 * nVar * 8 ]			
			vfmadd231pd		zmm17,  zmm5,  zmm30 
			vfmadd231pd     zmm21,  zmm5,  zmm31 			
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8 ]
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8 ]
			vfnmadd231pd	zmm10,  zmm6,  zmm28
			vfnmadd231pd	zmm14,  zmm6,  zmm29
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 1 * nVar * 8 ]
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 2 * nVar * 8 ]
			vfnmadd231pd	zmm18,  zmm6,  zmm30
			vfnmadd231pd	zmm22,  zmm6,  zmm31
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 1 * nVar * 8 ]
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 2 * nVar * 8 ]
			vfnmadd231pd	zmm11,  zmm2,  zmm28
			vfnmadd231pd	zmm15,  zmm2,  zmm29
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 1 * nVar * 8 ]
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 2 * nVar * 8 ]
			vfnmadd231pd	zmm19,  zmm2,  zmm30
			vfnmadd231pd	zmm23,  zmm2,  zmm31
			
			vmovapd		zmm28,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8 ]		
			vmovapd		zmm29,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 2 * nVar * 8 ]			
			vfmadd231pd		zmm8,   zmm2,  zmm28
			vfmadd231pd 	zmm12,  zmm2,  zmm29 
			vmovapd		zmm30,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 2 * nVar * 8 ]		
			vmovapd		zmm31,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 2 * nVar * 8 ]			
			vfmadd231pd		zmm16,  zmm2,  zmm30 
			vfmadd231pd     zmm20,  zmm2,  zmm31 			
			//prefetcht0[r12 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8 ]
			//prefetcht0[r12 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 3 * nVar * 8 ]			
			vfmadd231pd		zmm9,   zmm6,  zmm28
			vfmadd231pd     zmm13,  zmm6,  zmm29 
			//prefetcht0[r12 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 3 * nVar * 8 ]
			//prefetcht0[r12 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 3 * nVar * 8 ]			
			vfmadd231pd		zmm17,  zmm6,  zmm30 
			vfmadd231pd     zmm21,  zmm6,  zmm31 			
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8 ]
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 3 * nVar * 8 ]
			vfnmadd231pd	zmm10,  zmm5,  zmm28
			vfnmadd231pd	zmm14,  zmm5,  zmm29
			vfnmadd231pd	zmm18,  zmm5,  zmm30
			vfnmadd231pd	zmm22,  zmm5,  zmm31
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 3 * nVar * 8 ]
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 3 * nVar * 8 ]
			vfnmadd231pd	zmm11,  zmm1,  zmm28
			vfnmadd231pd	zmm15,  zmm1,  zmm29
			vfnmadd231pd	zmm19,  zmm1,  zmm30
			vfnmadd231pd	zmm23,  zmm1,  zmm31

			// 
			// store the result (very slow)  
			// 	
			vmovapd		zmm28,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8 ]		
			vmovapd		zmm29,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 3 * nVar * 8 ]			
			vmovapd		zmm30,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 3 * nVar * 8 ]		
			vmovapd		zmm31,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 3 * nVar * 8 ]			
			vfmadd231pd		zmm8,   zmm3,  zmm28
			vfmadd231pd 	zmm12,  zmm3,  zmm29 
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8], zmm8
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 0 * nVar * 8], zmm12
			vfmadd231pd		zmm16,  zmm3,  zmm30 
			vfmadd231pd     zmm20,  zmm3,  zmm31 			
			//prefetcht0[r12 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 4 * nVar * 8 ]
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 0 * nVar * 8], zmm16
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 0 * nVar * 8], zmm20
			vfmadd231pd		zmm9,   zmm7,  zmm28
			vfmadd231pd     zmm13,  zmm7,  zmm29 
			//prefetcht0[r12 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 4 * nVar * 8 ]			
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8], zmm9
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 1 * nVar * 8], zmm13
			vfmadd231pd		zmm17,  zmm7,  zmm30 
			vfmadd231pd     zmm21,  zmm7,  zmm31 			
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 1 * nVar * 8], zmm17
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 1 * nVar * 8], zmm21
			vfnmadd231pd	zmm10,  zmm4,  zmm28
			vfnmadd231pd	zmm14,  zmm4,  zmm29
			//prefetcht0[r12 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 4 * nVar * 8 ]
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8], zmm10
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 2 * nVar * 8], zmm14
			vfnmadd231pd	zmm18,  zmm4,  zmm30
			vfnmadd231pd	zmm22,  zmm4,  zmm31
			//prefetcht0[r12 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 4 * nVar * 8 ]			
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 2 * nVar * 8], zmm18
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 2 * nVar * 8], zmm22
			vfnmadd231pd	zmm11,  zmm0,  zmm28
			vfnmadd231pd	zmm15,  zmm0,  zmm29
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8], zmm11
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 3 * nVar * 8], zmm15
			vfnmadd231pd	zmm19,  zmm0,  zmm30
			vfnmadd231pd	zmm23,  zmm0,  zmm31
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 3 * nVar * 8], zmm19
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 3 * nVar * 8], zmm23			
			
		add r15, VECTORLENGTH*8 
		cmp r15, nVar*8 
		jl loopvar

		add r10, ND2D*nVar*8   
		cmp r10, 4*ND2D*nVar*8 
		jl loopz

		add r8, ND3D*nVar*8
		cmp r8, 4*ND3D*nVar*8 
		jl loopt 

	}
	
#else
extern "C" void __cdecl vsmm4x4cs4d_x(double c[4][4][4][4][nVar], double a[4][4][4][4][nVar], double b[4][4], double s[1])
{ 
	__assume_aligned(c, 32);
	__assume_aligned(a, 32);
	__assume_aligned(b, 32);
	int i, j, k, iy, iz, it;
	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	__m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;

	// optimized AVX code for matmul on all tensor slices 
	// load the scalar into register 14 
	ymm14 = _mm256_broadcast_sd(&s[0]);
	// load first column of the small matrix b and scale with the scalar s 
	ymm4 = _mm256_broadcast_sd(&b[0][0]);
	ymm5 = _mm256_broadcast_sd(&b[0][1]);
	ymm6 = _mm256_broadcast_sd(&b[0][2]);
	ymm7 = _mm256_broadcast_sd(&b[0][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load second column of the small matrix b and scale with the scalar s  
	ymm4 = _mm256_broadcast_sd(&b[1][0]);
	ymm5 = _mm256_broadcast_sd(&b[1][1]);
	ymm6 = _mm256_broadcast_sd(&b[1][2]);
	ymm7 = _mm256_broadcast_sd(&b[1][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);
	for (it = 0; it < 4; it++)
	{
		for (iz = 0; iz < 4; iz++)
		{
			for (iy = 0; iy < 4; iy++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of matrix a on position i 
					ymm4 = _mm256_load_pd(&a[it][iz][iy][0][i]);
					ymm5 = _mm256_load_pd(&a[it][iz][iy][1][i]);
					ymm6 = _mm256_load_pd(&a[it][iz][iy][2][i]);
					ymm7 = _mm256_load_pd(&a[it][iz][iy][3][i]);
					// set the result register 1 to zero 
					ymm9 = _mm256_setzero_pd();
					// set the result register 2 to zero
					ymm15 = _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the first and second column of b and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
					// multiply the row of a with the first and second column of b and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 
					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result in c 
					_mm256_store_pd(&c[it][iz][iy][0][i], ymm9);
					_mm256_store_pd(&c[it][iz][iy][1][i], ymm15);
				}
			}
		}
	}
	// load third column of the small matrix b and scale with the scalar s 
	ymm4 = _mm256_broadcast_sd(&b[2][0]);
	ymm5 = _mm256_broadcast_sd(&b[2][1]);
	ymm6 = _mm256_broadcast_sd(&b[2][2]);
	ymm7 = _mm256_broadcast_sd(&b[2][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load fourth column of the small matrix b and scale with the scalar s  
	ymm4 = _mm256_broadcast_sd(&b[3][0]);
	ymm5 = _mm256_broadcast_sd(&b[3][1]);
	ymm6 = _mm256_broadcast_sd(&b[3][2]);
	ymm7 = _mm256_broadcast_sd(&b[3][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);
	for (it = 0; it < 4; it++)
	{
		for (iz = 0; iz < 4; iz++)
		{
			for (iy = 0; iy < 4; iy++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of a on position i 
					ymm4 = _mm256_load_pd(&a[it][iz][iy][0][i]);
					ymm5 = _mm256_load_pd(&a[it][iz][iy][1][i]);
					ymm6 = _mm256_load_pd(&a[it][iz][iy][2][i]);
					ymm7 = _mm256_load_pd(&a[it][iz][iy][3][i]);
					// set the result register 1 to zero 
					ymm9 = _mm256_setzero_pd();
					// set the result register 2 to zero
					ymm15 = _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the first and second column of b and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
					// multiply the row of a with the first and second column of b and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 

					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result
					_mm256_store_pd(&c[it][iz][iy][2][i], ymm9);
					_mm256_store_pd(&c[it][iz][iy][3][i], ymm15);
				}
			}
		}
	}
#endif 
}

// compute all matmuls in y direction on a 4D data set 
#ifdef AVX512
extern "C" void __cdecl vsmm4x4cs4d_y(double *cptr, const double *aptr, double b[4][4], double s[1])
{
	__assume_aligned(cptr, 64);
	__assume_aligned(aptr, 64);
	__assume_aligned(b, 64);

	// matrix b zmm registers (our small stiffness matrix). We can use the fact that the matrix has anti-symmetric entries, as detailed below. 
	// zmm0 zmm4  -zmm7 -zmm3  
	// zmm1 zmm5  -zmm6 -zmm2 
	// zmm2 zmm6  -zmm5 -zmm1 
	// zmm3 zmm7  -zmm4 -zmm0  
	// 
	// matrix a registers 
	// 
	// zmm28 zmm29 zmm30 zmm31   
	// 
	// result registers 
	// first half of the operations 
	// zmm8  zmm9  zmm10 zmm11 
	// zmm12 zmm13 zmm14 zmm15 
	// second half of the operations 
	// zmm16 zmm17 zmm18 zmm19  
	// zmm20 zmm21 zmm22 zmm23 
	 
	__asm { 

		mov rax, aptr                                       // pointer to a 
		mov rbx, b											// pointer to b 
		mov rcx, cptr										// pointer to c 
		mov r13, s											// pointer to s  
		// mov rdx, dptr									// pointer to d  ( only for debugging ) 
		 
		vbroadcastsd zmm31, QWORD PTR[r13]					// load s in vector register 
		
		prefetcht0[rax + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8 ]
		prefetcht0[rax + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8 ]

		vbroadcastsd zmm20, QWORD PTR[rbx + 0 * 8 ]			// load b[0][0] 
		vbroadcastsd zmm21, QWORD PTR[rbx + 1 * 8 ]			// load b[0][1] 
		vbroadcastsd zmm22, QWORD PTR[rbx + 2 * 8 ]			// load b[0][2] 
		vbroadcastsd zmm23, QWORD PTR[rbx + 3 * 8 ]			// load b[0][3] 
		
		vmulpd zmm0,  zmm20, zmm31                          // multiply column of b with scalar s 
		vmulpd zmm1,  zmm21, zmm31
		vmulpd zmm2,  zmm22, zmm31
		vmulpd zmm3,  zmm23, zmm31

		prefetcht0[rax + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8 ]
		prefetcht0[rax + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8 ]

		vbroadcastsd zmm20, QWORD PTR[rbx + 4  * 8 ]		// load b[1][0] 
		vbroadcastsd zmm21, QWORD PTR[rbx + 5  * 8 ]		// load b[1][1] 
		vbroadcastsd zmm22, QWORD PTR[rbx + 6  * 8 ]		// load b[1][2] 
		vbroadcastsd zmm23, QWORD PTR[rbx + 7  * 8 ]		// load b[1][3] 
		
		vmulpd zmm4,  zmm20, zmm31                          // multiply column of b with scalar s 
		vmulpd zmm5,  zmm21, zmm31
		vmulpd zmm6,  zmm22, zmm31
		vmulpd zmm7,  zmm23, zmm31
		
		mov r8, 0
loopt:
		mov r10, 0
loopz:
		mov r15, 0   
loopvar: 		
			//  
			// compute address of the pointers 
			mov r14,  0  
			add r14,  r8         // time loop index 
			add r14,  r10        // z loop index  
			add r14,  r15        // variable loop index 
			// 
			mov r12, rax	     // compute a base pointer 
			add r12, r14        
			mov r13, rcx         // compute c base pointer 
			add r13, r14
			
			//
			// only one part   
			// 
			vmovapd		zmm28,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8 ]		
			vmovapd		zmm29,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8 ]			
			vmulpd		zmm8,   zmm0,  zmm28
			vmulpd      zmm12,  zmm0,  zmm29 			
			vmovapd		zmm30,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8 ]		
			vmovapd		zmm31,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8 ]			
			vmulpd		zmm16,  zmm0,  zmm30 
			vmulpd      zmm20,  zmm0,  zmm31 			
			vmulpd		zmm9,   zmm4,  zmm28
			vmulpd      zmm13,  zmm4,  zmm29 
			vmulpd		zmm17,  zmm4,  zmm30 
			vmulpd      zmm21,  zmm4,  zmm31 			
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8 ]
			vxorpd      zmm10,  zmm10, zmm10  
			vfnmadd231pd	zmm10,  zmm7,  zmm28
			vxorpd      zmm14,  zmm14, zmm14  
			vfnmadd231pd	zmm14,  zmm7,  zmm29
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8 ]
			vxorpd      zmm18,  zmm18, zmm18  
			vfnmadd231pd	zmm18,  zmm7,  zmm30
			vxorpd      zmm22,  zmm22, zmm22  
			vfnmadd231pd	zmm22,  zmm7,  zmm31
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8 ]
			vxorpd      zmm11,  zmm11, zmm11  
			vfnmadd231pd	zmm11,  zmm3,  zmm28
			vxorpd      zmm15,  zmm15, zmm15  
			vfnmadd231pd	zmm15,  zmm3,  zmm29
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8 ]
			vxorpd      zmm19,  zmm19, zmm19  
			vfnmadd231pd	zmm19,  zmm3,  zmm30
			vxorpd      zmm23,  zmm23, zmm23  
			vfnmadd231pd	zmm23,  zmm3,  zmm31
			
			vmovapd		zmm28,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 0 * nVar * 8 ]		
			vmovapd		zmm29,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 1 * nVar * 8 ]			
			vfmadd231pd		zmm8,   zmm1,  zmm28
			vfmadd231pd 	zmm12,  zmm1,  zmm29 
			vmovapd		zmm30,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 2 * nVar * 8 ]		
			vmovapd		zmm31,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 3 * nVar * 8 ]			
			vfmadd231pd		zmm16,  zmm1,  zmm30 
			vfmadd231pd     zmm20,  zmm1,  zmm31 
			vfmadd231pd		zmm9,   zmm5,  zmm28
			vfmadd231pd     zmm13,  zmm5,  zmm29 
			vfmadd231pd		zmm17,  zmm5,  zmm30 
			vfmadd231pd     zmm21,  zmm5,  zmm31 			
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 0 * nVar * 8 ]
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 0 * nVar * 8 ]
			vfnmadd231pd	zmm10,  zmm6,  zmm28
			vfnmadd231pd	zmm14,  zmm6,  zmm29
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 1 * nVar * 8 ]
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 1 * nVar * 8 ]
			vfnmadd231pd	zmm18,  zmm6,  zmm30
			vfnmadd231pd	zmm22,  zmm6,  zmm31
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 2 * nVar * 8 ]
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 2 * nVar * 8 ]
			vfnmadd231pd	zmm11,  zmm2,  zmm28
			vfnmadd231pd	zmm15,  zmm2,  zmm29
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 3 * nVar * 8 ]
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 3 * nVar * 8 ]
			vfnmadd231pd	zmm19,  zmm2,  zmm30
			vfnmadd231pd	zmm23,  zmm2,  zmm31
			
			vmovapd		zmm28,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 0 * nVar * 8 ]		
			vmovapd		zmm29,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 1 * nVar * 8 ]			
			vfmadd231pd		zmm8,   zmm2,  zmm28
			vfmadd231pd 	zmm12,  zmm2,  zmm29 
			vmovapd		zmm30,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 2 * nVar * 8 ]		
			vmovapd		zmm31,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 3 * nVar * 8 ]			
			vfmadd231pd		zmm16,  zmm2,  zmm30 
			vfmadd231pd     zmm20,  zmm2,  zmm31 			
			vfmadd231pd		zmm9,   zmm6,  zmm28
			vfmadd231pd     zmm13,  zmm6,  zmm29 
			vfmadd231pd		zmm17,  zmm6,  zmm30 
			vfmadd231pd     zmm21,  zmm6,  zmm31 			
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 0 * nVar * 8 ]
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 1 * nVar * 8 ]
			vfnmadd231pd	zmm10,  zmm5,  zmm28
			vfnmadd231pd	zmm14,  zmm5,  zmm29
			vfnmadd231pd	zmm18,  zmm5,  zmm30
			vfnmadd231pd	zmm22,  zmm5,  zmm31
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 2 * nVar * 8 ]
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 3 * nVar * 8 ]
			vfnmadd231pd	zmm11,  zmm1,  zmm28
			vfnmadd231pd	zmm15,  zmm1,  zmm29
			vfnmadd231pd	zmm19,  zmm1,  zmm30
			vfnmadd231pd	zmm23,  zmm1,  zmm31

			// 
			// store the result (very slow)  
			// 	
			vmovapd		zmm28,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 0 * nVar * 8 ]		
			vmovapd		zmm29,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 1 * nVar * 8 ]			
			vmovapd		zmm30,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 2 * nVar * 8 ]		
			vmovapd		zmm31,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 3 * nVar * 8 ]			
			vfmadd231pd		zmm8,   zmm3,  zmm28
			vfmadd231pd 	zmm12,  zmm3,  zmm29 
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8], zmm8
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8], zmm12
			vfmadd231pd		zmm16,  zmm3,  zmm30 
			vfmadd231pd     zmm20,  zmm3,  zmm31 			
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8], zmm16
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8], zmm20
			vfmadd231pd		zmm9,   zmm7,  zmm28
			vfmadd231pd     zmm13,  zmm7,  zmm29 
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 0 * nVar * 8], zmm9
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 1 * nVar * 8], zmm13
			vfmadd231pd		zmm17,  zmm7,  zmm30 
			vfmadd231pd     zmm21,  zmm7,  zmm31 			
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 2 * nVar * 8], zmm17
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 1 * ND1D*nVar * 8 + 3 * nVar * 8], zmm21
			vfnmadd231pd	zmm10,  zmm4,  zmm28
			vfnmadd231pd	zmm14,  zmm4,  zmm29
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 0 * nVar * 8], zmm10
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 1 * nVar * 8], zmm14
			vfnmadd231pd	zmm18,  zmm4,  zmm30
			vfnmadd231pd	zmm22,  zmm4,  zmm31
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 2 * nVar * 8], zmm18
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 2 * ND1D*nVar * 8 + 3 * nVar * 8], zmm22
			vfnmadd231pd	zmm11,  zmm0,  zmm28
			vfnmadd231pd	zmm15,  zmm0,  zmm29
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 0 * nVar * 8], zmm11
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 1 * nVar * 8], zmm15
			vfnmadd231pd	zmm19,  zmm0,  zmm30
			vfnmadd231pd	zmm23,  zmm0,  zmm31
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 2 * nVar * 8], zmm19
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 3 * ND1D*nVar * 8 + 3 * nVar * 8], zmm23			
			
		add r15, VECTORLENGTH*8 
		cmp r15, nVar*8 
		jl loopvar

		add r10, ND2D*nVar*8   
		cmp r10, 4*ND2D*nVar*8 
		jl loopz

		add r8, ND3D*nVar*8
		cmp r8, 4*ND3D*nVar*8 
		jl loopt 

	}

}


#else 
extern "C" void __cdecl vsmm4x4cs4d_y(double c[4][4][4][4][nVar], const double a[4][4][4][4][nVar], const double b[4][4], const double s[1])
{
	__assume_aligned(c, 32);
	__assume_aligned(a, 32);
	__assume_aligned(b, 32);
	int i, j, k, ix, iz, it;
	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	__m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;

	// optimized AVX code for matmul on all tensor slices 
	// load the scalar into register 14 
	ymm14 = _mm256_broadcast_sd(&s[0]);
	// load first column of the small matrix b and scale with scalar s 
	ymm4 = _mm256_broadcast_sd(&b[0][0]);
	ymm5 = _mm256_broadcast_sd(&b[0][1]);
	ymm6 = _mm256_broadcast_sd(&b[0][2]);
	ymm7 = _mm256_broadcast_sd(&b[0][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load second column of the small matrix b and scale with scalar s  
	ymm4 = _mm256_broadcast_sd(&b[1][0]);
	ymm5 = _mm256_broadcast_sd(&b[1][1]);
	ymm6 = _mm256_broadcast_sd(&b[1][2]);
	ymm7 = _mm256_broadcast_sd(&b[1][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);
	for (it = 0; it < 4; it++)
	{
		for (iz = 0; iz < 4; iz++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of matrix a on position i 
					ymm4 = _mm256_load_pd(&a[it][iz][0][ix][i]);
					ymm5 = _mm256_load_pd(&a[it][iz][1][ix][i]);
					ymm6 = _mm256_load_pd(&a[it][iz][2][ix][i]);
					ymm7 = _mm256_load_pd(&a[it][iz][3][ix][i]);
					// set the result register 1 to zero 
					ymm9 = _mm256_setzero_pd();
					// set the result register 2 to zero
					ymm15 = _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the first and second column of b and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
					// multiply the row of a with the first and second column of b and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 
					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result in c 
					_mm256_store_pd(&c[it][iz][0][ix][i], ymm9);
					_mm256_store_pd(&c[it][iz][1][ix][i], ymm15);
				}
			}
		}
	}
	// load third column of the small matrix b 
	ymm4 = _mm256_broadcast_sd(&b[2][0]);
	ymm5 = _mm256_broadcast_sd(&b[2][1]);
	ymm6 = _mm256_broadcast_sd(&b[2][2]);
	ymm7 = _mm256_broadcast_sd(&b[2][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load fourth column of the small matrix b 
	ymm4 = _mm256_broadcast_sd(&b[3][0]);
	ymm5 = _mm256_broadcast_sd(&b[3][1]);
	ymm6 = _mm256_broadcast_sd(&b[3][2]);
	ymm7 = _mm256_broadcast_sd(&b[3][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);
	for (it = 0; it < 4; it++)
	{
		for (iz = 0; iz < 4; iz++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of a on position i 
					ymm4 = _mm256_load_pd(&a[it][iz][0][ix][i]);
					ymm5 = _mm256_load_pd(&a[it][iz][1][ix][i]);
					ymm6 = _mm256_load_pd(&a[it][iz][2][ix][i]);
					ymm7 = _mm256_load_pd(&a[it][iz][3][ix][i]);
					// set the result register 1 to zero 
					ymm9 = _mm256_setzero_pd();
					// set the result register 2 to zero
					ymm15 = _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the first and second column of b and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
					// multiply the row of a with the first and second column of b and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 

					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result
					_mm256_store_pd(&c[it][iz][2][ix][i], ymm9);
					_mm256_store_pd(&c[it][iz][3][ix][i], ymm15);
				} 
			}
		}
	}

}
#endif 

#ifdef AVX512
extern "C" void __cdecl vsmm4x4cs4d_z(double *cptr, const double *aptr, double b[4][4], double s[1])
{
	__assume_aligned(cptr, 64);
	__assume_aligned(aptr, 64);
	__assume_aligned(b, 64);

	// matrix b zmm registers (our small stiffness matrix). We can use the fact that the matrix has anti-symmetric entries, as detailed below. 
	// zmm0 zmm4  -zmm7 -zmm3  
	// zmm1 zmm5  -zmm6 -zmm2 
	// zmm2 zmm6  -zmm5 -zmm1 
	// zmm3 zmm7  -zmm4 -zmm0  
	// 
	// matrix a registers 
	// 
	// zmm28 zmm29 zmm30 zmm31   
	// 
	// result registers 
	// first half of the operations 
	// zmm8  zmm9  zmm10 zmm11 
	// zmm12 zmm13 zmm14 zmm15 
	// second half of the operations 
	// zmm16 zmm17 zmm18 zmm19  
	// zmm20 zmm21 zmm22 zmm23 
	 
	__asm { 

		mov rax, aptr                                       // pointer to a 
		mov rbx, b											// pointer to b 
		mov rcx, cptr										// pointer to c 
		mov r13, s											// pointer to s  
		// mov rdx, dptr									// pointer to d  ( only for debugging ) 
		 
		vbroadcastsd zmm31, QWORD PTR[r13]					// load s in vector register 
		
		prefetcht0[rax + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8 ]
		prefetcht0[rax + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8 ]

		vbroadcastsd zmm20, QWORD PTR[rbx + 0 * 8 ]			// load b[0][0] 
		vbroadcastsd zmm21, QWORD PTR[rbx + 1 * 8 ]			// load b[0][1] 
		vbroadcastsd zmm22, QWORD PTR[rbx + 2 * 8 ]			// load b[0][2] 
		vbroadcastsd zmm23, QWORD PTR[rbx + 3 * 8 ]			// load b[0][3] 
		
		vmulpd zmm0,  zmm20, zmm31                          // multiply column of b with scalar s 
		vmulpd zmm1,  zmm21, zmm31
		vmulpd zmm2,  zmm22, zmm31
		vmulpd zmm3,  zmm23, zmm31

		prefetcht0[rax + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8 ]
		prefetcht0[rax + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8 ]

		vbroadcastsd zmm20, QWORD PTR[rbx + 4  * 8 ]		// load b[1][0] 
		vbroadcastsd zmm21, QWORD PTR[rbx + 5  * 8 ]		// load b[1][1] 
		vbroadcastsd zmm22, QWORD PTR[rbx + 6  * 8 ]		// load b[1][2] 
		vbroadcastsd zmm23, QWORD PTR[rbx + 7  * 8 ]		// load b[1][3] 
		
		vmulpd zmm4,  zmm20, zmm31                          // multiply column of b with scalar s 
		vmulpd zmm5,  zmm21, zmm31
		vmulpd zmm6,  zmm22, zmm31
		vmulpd zmm7,  zmm23, zmm31
		
		mov r8, 0
loopt:
		mov r10, 0
loopy:
		mov r15, 0   
loopvar: 		
			//  
			// compute address of the pointers 
			mov r14,  0  
			add r14,  r8         // time loop index 
			add r14,  r10        // y loop index  
			add r14,  r15        // variable loop index 
			// 
			mov r12, rax	     // compute a base pointer 
			add r12, r14        
			mov r13, rcx         // compute c base pointer 
			add r13, r14
			
			//
			// only one part   
			// 
			vmovapd		zmm28,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8 ]		
			vmovapd		zmm29,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8 ]			
			vmulpd		zmm8,   zmm0,  zmm28
			vmulpd      zmm12,  zmm0,  zmm29 			
			vmovapd		zmm30,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8 ]		
			vmovapd		zmm31,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8 ]			
			vmulpd		zmm16,  zmm0,  zmm30 
			vmulpd      zmm20,  zmm0,  zmm31 			
			vmulpd		zmm9,   zmm4,  zmm28
			vmulpd      zmm13,  zmm4,  zmm29 
			vmulpd		zmm17,  zmm4,  zmm30 
			vmulpd      zmm21,  zmm4,  zmm31 			
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8 ]
			vxorpd      zmm10,  zmm10, zmm10  
			vfnmadd231pd	zmm10,  zmm7,  zmm28
			vxorpd      zmm14,  zmm14, zmm14  
			vfnmadd231pd	zmm14,  zmm7,  zmm29
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8 ]
			vxorpd      zmm18,  zmm18, zmm18  
			vfnmadd231pd	zmm18,  zmm7,  zmm30
			vxorpd      zmm22,  zmm22, zmm22  
			vfnmadd231pd	zmm22,  zmm7,  zmm31
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8 ]
			vxorpd      zmm11,  zmm11, zmm11  
			vfnmadd231pd	zmm11,  zmm3,  zmm28
			vxorpd      zmm15,  zmm15, zmm15  
			vfnmadd231pd	zmm15,  zmm3,  zmm29
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8 ]
			vxorpd      zmm19,  zmm19, zmm19  
			vfnmadd231pd	zmm19,  zmm3,  zmm30
			vxorpd      zmm23,  zmm23, zmm23  
			vfnmadd231pd	zmm23,  zmm3,  zmm31
			
			vmovapd		zmm28,  ZMMWORD PTR[r12 + 1 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8 ]		
			vmovapd		zmm29,  ZMMWORD PTR[r12 + 1 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8 ]			
			vfmadd231pd		zmm8,   zmm1,  zmm28
			vfmadd231pd 	zmm12,  zmm1,  zmm29 
			vmovapd		zmm30,  ZMMWORD PTR[r12 + 1 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8 ]		
			vmovapd		zmm31,  ZMMWORD PTR[r12 + 1 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8 ]			
			vfmadd231pd		zmm16,  zmm1,  zmm30 
			vfmadd231pd     zmm20,  zmm1,  zmm31 
			vfmadd231pd		zmm9,   zmm5,  zmm28
			vfmadd231pd     zmm13,  zmm5,  zmm29 
			vfmadd231pd		zmm17,  zmm5,  zmm30 
			vfmadd231pd     zmm21,  zmm5,  zmm31 			
			prefetchw[r13 + 1 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8 ]
			prefetchw[r13 + 2 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8 ]
			vfnmadd231pd	zmm10,  zmm6,  zmm28
			vfnmadd231pd	zmm14,  zmm6,  zmm29
			prefetchw[r13 + 1 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8 ]
			prefetchw[r13 + 2 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8 ]
			vfnmadd231pd	zmm18,  zmm6,  zmm30
			vfnmadd231pd	zmm22,  zmm6,  zmm31
			prefetchw[r13 + 1 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8 ]
			prefetchw[r13 + 2 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8 ]
			vfnmadd231pd	zmm11,  zmm2,  zmm28
			vfnmadd231pd	zmm15,  zmm2,  zmm29
			prefetchw[r13 + 1 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8 ]
			prefetchw[r13 + 2 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8 ]
			vfnmadd231pd	zmm19,  zmm2,  zmm30
			vfnmadd231pd	zmm23,  zmm2,  zmm31
			
			vmovapd		zmm28,  ZMMWORD PTR[r12 + 2 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8 ]		
			vmovapd		zmm29,  ZMMWORD PTR[r12 + 2 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8 ]			
			vfmadd231pd		zmm8,   zmm2,  zmm28
			vfmadd231pd 	zmm12,  zmm2,  zmm29 
			vmovapd		zmm30,  ZMMWORD PTR[r12 + 2 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8 ]		
			vmovapd		zmm31,  ZMMWORD PTR[r12 + 2 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8 ]			
			vfmadd231pd		zmm16,  zmm2,  zmm30 
			vfmadd231pd     zmm20,  zmm2,  zmm31 			
			vfmadd231pd		zmm9,   zmm6,  zmm28
			vfmadd231pd     zmm13,  zmm6,  zmm29 
			vfmadd231pd		zmm17,  zmm6,  zmm30 
			vfmadd231pd     zmm21,  zmm6,  zmm31 			
			prefetchw[r13 + 3 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8 ]
			prefetchw[r13 + 3 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8 ]
			vfnmadd231pd	zmm10,  zmm5,  zmm28
			vfnmadd231pd	zmm14,  zmm5,  zmm29
			vfnmadd231pd	zmm18,  zmm5,  zmm30
			vfnmadd231pd	zmm22,  zmm5,  zmm31
			prefetchw[r13 + 3 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8 ]
			prefetchw[r13 + 3 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8 ]
			vfnmadd231pd	zmm11,  zmm1,  zmm28
			vfnmadd231pd	zmm15,  zmm1,  zmm29
			vfnmadd231pd	zmm19,  zmm1,  zmm30
			vfnmadd231pd	zmm23,  zmm1,  zmm31

			// 
			// store the result (very slow)  
			// 	
			vmovapd		zmm28,  ZMMWORD PTR[r12 + 3 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8 ]		
			vmovapd		zmm29,  ZMMWORD PTR[r12 + 3 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8 ]			
			vmovapd		zmm30,  ZMMWORD PTR[r12 + 3 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8 ]		
			vmovapd		zmm31,  ZMMWORD PTR[r12 + 3 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8 ]			
			vfmadd231pd		zmm8,   zmm3,  zmm28
			vfmadd231pd 	zmm12,  zmm3,  zmm29 
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8], zmm8
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8], zmm12
			vfmadd231pd		zmm16,  zmm3,  zmm30 
			vfmadd231pd     zmm20,  zmm3,  zmm31 			
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8], zmm16
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8], zmm20
			vfmadd231pd		zmm9,   zmm7,  zmm28
			vfmadd231pd     zmm13,  zmm7,  zmm29 
			vmovapd		 ZMMWORD PTR[r13 + 1 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8], zmm9
			vmovapd		 ZMMWORD PTR[r13 + 1 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8], zmm13
			vfmadd231pd		zmm17,  zmm7,  zmm30 
			vfmadd231pd     zmm21,  zmm7,  zmm31 			
			vmovapd		 ZMMWORD PTR[r13 + 1 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8], zmm17
			vmovapd		 ZMMWORD PTR[r13 + 1 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8], zmm21
			vfnmadd231pd	zmm10,  zmm4,  zmm28
			vfnmadd231pd	zmm14,  zmm4,  zmm29
			vmovapd		 ZMMWORD PTR[r13 + 2 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8], zmm10
			vmovapd		 ZMMWORD PTR[r13 + 2 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8], zmm14
			vfnmadd231pd	zmm18,  zmm4,  zmm30
			vfnmadd231pd	zmm22,  zmm4,  zmm31
			vmovapd		 ZMMWORD PTR[r13 + 2 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8], zmm18
			vmovapd		 ZMMWORD PTR[r13 + 2 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8], zmm22
			vfnmadd231pd	zmm11,  zmm0,  zmm28
			vfnmadd231pd	zmm15,  zmm0,  zmm29
			vmovapd		 ZMMWORD PTR[r13 + 3 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8], zmm11
			vmovapd		 ZMMWORD PTR[r13 + 3 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8], zmm15
			vfnmadd231pd	zmm19,  zmm0,  zmm30
			vfnmadd231pd	zmm23,  zmm0,  zmm31
			vmovapd		 ZMMWORD PTR[r13 + 3 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8], zmm19
			vmovapd		 ZMMWORD PTR[r13 + 3 * ND2D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8], zmm23			
			
		add r15, VECTORLENGTH*8 
		cmp r15, nVar*8 
		jl loopvar

		add r10, ND1D*nVar*8   
		cmp r10, 4*ND1D*nVar*8 
		jl loopy

		add r8, ND3D*nVar*8
		cmp r8, 4*ND3D*nVar*8 
		jl loopt 

	}

}

#else 

// compute all matmuls in z direction on a 4D data set 
extern "C" void __cdecl vsmm4x4cs4d_z(double c[4][4][4][4][nVar], const double a[4][4][4][4][nVar], const double b[4][4], const double s[1])
{
	__assume_aligned(c, 32);
	__assume_aligned(a, 32);
	__assume_aligned(b, 32);
	int i, j, k, ix, iy, it;
	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	__m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;

	// optimized AVX code for matmul on all tensor slices 
	// load the scalar into register 14 
	ymm14 = _mm256_broadcast_sd(&s[0]);
	// load first column of the small matrix b and scale with scalar s 
	ymm4 = _mm256_broadcast_sd(&b[0][0]);
	ymm5 = _mm256_broadcast_sd(&b[0][1]);
	ymm6 = _mm256_broadcast_sd(&b[0][2]);
	ymm7 = _mm256_broadcast_sd(&b[0][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load second column of the small matrix b and scale with scalar s  
	ymm4 = _mm256_broadcast_sd(&b[1][0]);
	ymm5 = _mm256_broadcast_sd(&b[1][1]);
	ymm6 = _mm256_broadcast_sd(&b[1][2]);
	ymm7 = _mm256_broadcast_sd(&b[1][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);
	for (it = 0; it < 4; it++)
	{
		for (iy = 0; iy < 4; iy++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of matrix a on position i 
					ymm4 = _mm256_load_pd(&a[it][0][iy][ix][i]);
					ymm5 = _mm256_load_pd(&a[it][1][iy][ix][i]);
					ymm6 = _mm256_load_pd(&a[it][2][iy][ix][i]);
					ymm7 = _mm256_load_pd(&a[it][3][iy][ix][i]);
					// set the result register 1 to zero 
					ymm9 = _mm256_setzero_pd();
					// set the result register 2 to zero
					ymm15 = _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the first and second column of b and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
					// multiply the row of a with the first and second column of b and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 
					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result in c 
					_mm256_store_pd(&c[it][0][iy][ix][i], ymm9);
					_mm256_store_pd(&c[it][1][iy][ix][i], ymm15);
				}
			}
		}
	}
	// load third column of the small matrix b 
	ymm4 = _mm256_broadcast_sd(&b[2][0]);
	ymm5 = _mm256_broadcast_sd(&b[2][1]);
	ymm6 = _mm256_broadcast_sd(&b[2][2]);
	ymm7 = _mm256_broadcast_sd(&b[2][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load fourth column of the small matrix b 
	ymm4 = _mm256_broadcast_sd(&b[3][0]);
	ymm5 = _mm256_broadcast_sd(&b[3][1]);
	ymm6 = _mm256_broadcast_sd(&b[3][2]);
	ymm7 = _mm256_broadcast_sd(&b[3][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);
	for (it = 0; it < 4; it++)
	{
		for (iy = 0; iy < 4; iy++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of a on position i 
					ymm4 = _mm256_load_pd(&a[it][0][iy][ix][i]);
					ymm5 = _mm256_load_pd(&a[it][1][iy][ix][i]);
					ymm6 = _mm256_load_pd(&a[it][2][iy][ix][i]);
					ymm7 = _mm256_load_pd(&a[it][3][iy][ix][i]);
					// set the result register 1 to zero 
					ymm9 = _mm256_setzero_pd();
					// set the result register 2 to zero
					ymm15 = _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the first and second column of b and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
					// multiply the row of a with the first and second column of b and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 

					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result
					_mm256_store_pd(&c[it][2][iy][ix][i], ymm9);
					_mm256_store_pd(&c[it][3][iy][ix][i], ymm15);
				}
			}
		}
	}

}
#endif 



// -----------------------------------------------------------------------------------------------------------

// compute all matmuls in x direction on a 4D data set with vector weights (very useful for stiffness matrix multiplications) 
extern "C" void __cdecl vsmm4x4c4d_kx(double c[4][4][4][4][nVar], double a[4][4][4][4][nVar], double b[4][4], double s[1] )
{
#ifdef AVX512

	/* 
	__assume_aligned(c, 64);
	__assume_aligned(a, 64);
	__assume_aligned(b, 64);

	

	__asm {
		mov rax, a											// pointer to a 
		mov rbx, b											// pointer to b 
		mov rcx, c											// pointer to c 

		mov r8, 0
loopt:
		mov r10, 0
loopz:
		mov r15, 0  

		mov r13, s											// pointer to s  
		mov r11, wGP										// pointer to d  

		vbroadcastsd zmm14, QWORD PTR[r13]					// load s in vector register 

		vbroadcastsd zmm4, QWORD PTR[rbx + 0 * 8]			// load b[0][0] 
		vbroadcastsd zmm5, QWORD PTR[rbx + 1 * 8]			// load b[0][1] 
		vbroadcastsd zmm6, QWORD PTR[rbx + 2 * 8]			// load b[0][2] 
		vbroadcastsd zmm7, QWORD PTR[rbx + 3 * 8]			// load b[0][3] 

		vmulpd zmm0, zmm4, zmm14                            // multiply column of b with scalar s 
		vmulpd zmm1, zmm5, zmm14
		vmulpd zmm2, zmm6, zmm14
		vmulpd zmm3, zmm7, zmm14

		vbroadcastsd zmm4, QWORD PTR[rbx + 4 * 8]			// load b[1][0] 
		vbroadcastsd zmm5, QWORD PTR[rbx + 5 * 8]			// load b[1][1] 
		vbroadcastsd zmm6, QWORD PTR[rbx + 6 * 8]			// load b[1][2] 
		vbroadcastsd zmm7, QWORD PTR[rbx + 7 * 8]			// load b[1][3] 

		vmulpd zmm10, zmm4, zmm14							// multiply column of b with scalar s  
		vmulpd zmm11, zmm5, zmm14
		vmulpd zmm12, zmm6, zmm14
		vmulpd zmm13, zmm7, zmm14

		vbroadcastsd zmm4, QWORD PTR[rbx + 8 * 8]			// load b[1][0] 
		vbroadcastsd zmm5, QWORD PTR[rbx + 9 * 8]			// load b[1][1] 
		vbroadcastsd zmm6, QWORD PTR[rbx + 10 * 8]			// load b[1][2] 
		vbroadcastsd zmm7, QWORD PTR[rbx + 11 * 8]			// load b[1][3] 

		vmulpd zmm20, zmm4, zmm14							// multiply column of b with scalar s  
		vmulpd zmm21, zmm5, zmm14
		vmulpd zmm22, zmm6, zmm14
		vmulpd zmm23, zmm7, zmm14

		vbroadcastsd zmm4, QWORD PTR[rbx + 12 * 8]		    // load b[1][0] 
		vbroadcastsd zmm5, QWORD PTR[rbx + 13 * 8]		    // load b[1][1] 
		vbroadcastsd zmm6, QWORD PTR[rbx + 14 * 8]		    // load b[1][2] 
		vbroadcastsd zmm7, QWORD PTR[rbx + 15 * 8]		    // load b[1][3] 

		vmulpd zmm24, zmm4, zmm14							// multiply column of b with scalar s  
		vmulpd zmm25, zmm5, zmm14
		vmulpd zmm26, zmm6, zmm14
		vmulpd zmm27, zmm7, zmm14

		vxorpd zmm4, zmm4, zmm4
		vxorpd zmm5, zmm5, zmm5
		vxorpd zmm6, zmm6, zmm6
		vxorpd zmm7, zmm7, zmm7
		vxorpd zmm8, zmm8, zmm8
		vxorpd zmm14, zmm14, zmm14
		vxorpd zmm18, zmm18, zmm18
		vxorpd zmm19, zmm19, zmm19

loopvar: 

			// compute address of the pointers 
			mov r11, rax		// save rax 
			mov r14,  0
			mov rax,  r8 
			mov r9, ND3D*nVar * 8 
			mul r9 
			add r14,  rax        // time loop index 
			mov rax, r10 
			mov r9, ND2D*nVar * 8 
			mul r9 
			add r14,  rax        // z loop index 
			add r14,  r15        // variable loop index 
			mov rax, r11 

			mov r12, rax 	     // compute a base pointer 
			add r12, r14        
			mov r13, rcx         // compute c base pointer 
			add r13, r14

			// first half 
			vmovapd		zmm4,  ZMMWORD PTR[r12 + 0 * ND1D*nVar*8 + 0 * nVar*8]		    // zmm4 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 0 * nVar + i * VECTORLENGTH);  // &a[it][iz][0][0][i] 
			vmovapd		zmm5,  ZMMWORD PTR[r12 + 0 * ND1D*nVar*8 + 1 * nVar*8]		    // zmm5 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 1 * nVar + i * VECTORLENGTH);  // &a[it][iz][0][1][i]  
			vmovapd		zmm6,  ZMMWORD PTR[r12 + 0 * ND1D*nVar*8 + 2 * nVar*8]		    // zmm6 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 2 * nVar + i * VECTORLENGTH);  // &a[it][iz][0][2][i] 
			vmovapd		zmm7,  ZMMWORD PTR[r12 + 0 * ND1D*nVar*8 + 3 * nVar*8]		    // zmm7 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 3 * nVar + i * VECTORLENGTH);  // &a[it][iz][0][3][i] 

			vmovapd		zmm9, ZMMWORD PTR[r13 + 0 * ND1D*nVar * 8 + 0 * nVar * 8]		// zmm4 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 0 * nVar + i * VECTORLENGTH);  // &a[it][iz][0][0][i] 
			vmovapd		zmm15, ZMMWORD PTR[r13 + 0 * ND1D*nVar * 8 + 1 * nVar * 8]		// zmm5 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 1 * nVar + i * VECTORLENGTH);  // &a[it][iz][0][1][i]  
			vmovapd		zmm16, ZMMWORD PTR[r13 + 0 * ND1D*nVar * 8 + 2 * nVar * 8]		// zmm6 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 2 * nVar + i * VECTORLENGTH);  // &a[it][iz][0][2][i] 
			vmovapd		zmm17, ZMMWORD PTR[r13 + 0 * ND1D*nVar * 8 + 3 * nVar * 8]		// zmm7 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 3 * nVar + i * VECTORLENGTH);  // &a[it][iz][0][3][i] 

			vfmadd231pd   zmm9, zmm4,  zmm0									// zmm9 = _mm512_fmadd_pd(zmm4, zmm0, zmm9);  // 09-4 		  				  		 
			vfmadd231pd   zmm15, zmm5, zmm11								// zmm15 = _mm512_fmadd_pd(zmm5, zmm11, zmm15); // 15-5			  				  
			vfnmadd231pd  zmm16, zmm6, zmm11								// zmm16 = _mm512_fmadd_pd(zmm6, zmm22, zmm16); // 16-6
			vfnmadd231pd  zmm17, zmm7, zmm0 								// zmm17 = _mm512_fmadd_pd(zmm7, zmm27, zmm17);   // 17-7 

			vmovapd		zmm8,  ZMMWORD PTR[r12 + 1 * ND1D*nVar*8 + 0 * nVar*8]     	    // zmm8 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 0 * nVar + i * VECTORLENGTH); // &a[it][iz][1][0][i] 
			vmovapd		zmm14, ZMMWORD PTR[r12 + 1 * ND1D*nVar*8 + 1 * nVar*8]		    // zmm14 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 1 * nVar + i * VECTORLENGTH); // & a[it][iz][1][1][i]
			vmovapd		zmm18, ZMMWORD PTR[r12 + 1 * ND1D*nVar*8 + 2 * nVar*8]     	    // zmm18 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 2 * nVar + i * VECTORLENGTH); // &a[it][iz][1][2][i] 
			vmovapd		zmm19, ZMMWORD PTR[r12 + 1 * ND1D*nVar*8 + 3 * nVar*8]		    // zmm19 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 3 * nVar + i * VECTORLENGTH);  // &a[it][iz][1][3][i] 

			vmovapd		zmm28, ZMMWORD PTR[r13 + 1 * ND1D*nVar * 8 + 0 * nVar * 8]     	// zmm8 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 0 * nVar + i * VECTORLENGTH); // &a[it][iz][1][0][i] 
			vmovapd		zmm29, ZMMWORD PTR[r13 + 1 * ND1D*nVar * 8 + 1 * nVar * 8]		// zmm14 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 1 * nVar + i * VECTORLENGTH); // & a[it][iz][1][1][i]
			vmovapd		zmm30, ZMMWORD PTR[r13 + 1 * ND1D*nVar * 8 + 2 * nVar * 8]     	// zmm18 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 2 * nVar + i * VECTORLENGTH); // &a[it][iz][1][2][i] 
			vmovapd		zmm31, ZMMWORD PTR[r13 + 1 * ND1D*nVar * 8 + 3 * nVar * 8]		// zmm19 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 3 * nVar + i * VECTORLENGTH);  // &a[it][iz][1][3][i] 

			vfmadd231pd   zmm28, zmm8,  zmm0							// zmm28 = _mm512_fmadd_pd(zmm8, zmm0, zmm28);     // 28-4 			  		 
			vfmadd231pd   zmm29, zmm14, zmm11							// zmm29 = _mm512_fmadd_pd(zmm14, zmm11, zmm29);   // 29-5		 		  				  
			vfnmadd231pd  zmm30, zmm18, zmm11							// zmm30 = _mm512_fmadd_pd(zmm18, zmm22, zmm30);   // 30-6	 
			vfnmadd231pd  zmm31, zmm19, zmm0							// zmm31 = _mm512_fmadd_pd(zmm19, zmm27, zmm31);   // 31-7 

			prefetchw[r13 + 0 * ND1D*nVar * 8 + 0 * nVar * 8]
			prefetchw[r13 + 0 * ND1D*nVar * 8 + 1 * nVar * 8]
			prefetchw[r13 + 0 * ND1D*nVar * 8 + 2 * nVar * 8]
			prefetchw[r13 + 0 * ND1D*nVar * 8 + 3 * nVar * 8]
			prefetchw[r13 + 1 * ND1D*nVar * 8 + 0 * nVar * 8]
			prefetchw[r13 + 1 * ND1D*nVar * 8 + 1 * nVar * 8]
			prefetchw[r13 + 1 * ND1D*nVar * 8 + 2 * nVar * 8]
			prefetchw[r13 + 1 * ND1D*nVar * 8 + 3 * nVar * 8]

			vfmadd231pd  zmm9,  zmm5,  zmm1							// zmm9 = _mm512_fmadd_pd(zmm5, zmm1, zmm9);      // 09-5 		  				  		 
			vfmadd231pd  zmm15, zmm6,  zmm12						// zmm15 = _mm512_fmadd_pd(zmm6, zmm12, zmm15);    // 15-6	  				  
			vfnmadd231pd zmm16, zmm7,  zmm10						// zmm16 = _mm512_fmadd_pd(zmm7, zmm23, zmm16);    // 16-7
			vfnmadd231pd zmm17, zmm4,  zmm3 						// zmm17 = _mm512_fmadd_pd(zmm4, zmm24, zmm17);    // 17-4
			vfmadd231pd  zmm28, zmm14, zmm1							// zmm28 = _mm512_fmadd_pd(zmm14, zmm1, zmm28);    // 28-5		  		 
			vfmadd231pd  zmm29, zmm18, zmm12						// zmm29 = _mm512_fmadd_pd(zmm18, zmm12, zmm29);   // 29-6	 		  				  
			vfnmadd231pd zmm30, zmm19, zmm10						// zmm30 = _mm512_fmadd_pd(zmm19, zmm23, zmm30);   // 30-7
			vfnmadd231pd zmm31, zmm8,  zmm3							// zmm31 = _mm512_fmadd_pd(zmm8, zmm24, zmm31);    // 31-4

			vfmadd231pd  zmm9,  zmm6,  zmm2							// zmm9 = _mm512_fmadd_pd(zmm6, zmm2, zmm9);	    // 09-6		  				  		 
			vfmadd231pd  zmm15, zmm7,  zmm13						// zmm15 = _mm512_fmadd_pd(zmm7, zmm13, zmm15);    // 15-7	  				  
			vfnmadd231pd zmm16, zmm4,  zmm13						// zmm16 = _mm512_fmadd_pd(zmm4, zmm20, zmm16);    // 16-4
			vfnmadd231pd zmm17, zmm5,  zmm2							// zmm17 = _mm512_fmadd_pd(zmm5, zmm25, zmm17);    // 17-5
			vfmadd231pd  zmm28, zmm18, zmm2							// zmm28 = _mm512_fmadd_pd(zmm18, zmm2, zmm28);    // 28-6		  		 
			vfmadd231pd  zmm29, zmm19, zmm13						// zmm29 = _mm512_fmadd_pd(zmm19, zmm13, zmm29);   // 29-7	 		  				  
			vfnmadd231pd zmm30, zmm8,  zmm13						// zmm30 = _mm512_fmadd_pd(zmm8, zmm20, zmm30);    // 30-4
			vfnmadd231pd zmm31, zmm14, zmm2							// zmm31 = _mm512_fmadd_pd(zmm14, zmm25, zmm31);   // 31-5

			vfmadd231pd  zmm9, zmm7,  zmm3							// zmm9 = _mm512_fmadd_pd(zmm7, zmm3, zmm9);				  		 
			vfmadd231pd  zmm15, zmm4, zmm10							// zmm15 = _mm512_fmadd_pd(zmm4, zmm10, zmm15);				  
			vfnmadd231pd zmm16, zmm5, zmm12							// zmm16 = _mm512_fmadd_pd(zmm5, zmm21, zmm16);
			vfnmadd231pd zmm17, zmm6, zmm1 							// zmm17 = _mm512_fmadd_pd(zmm6, zmm26, zmm17);

			vmovapd		ZMMWORD PTR[r13 + 0 * ND1D*nVar * 8 + 0 * nVar * 8], zmm9     			        // _mm512_store_pd(cptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 0 * nVar + i * VECTORLENGTH, zmm9);     // &c[it][iz][0][0][i]
			vmovapd		ZMMWORD PTR[r13 + 0 * ND1D*nVar * 8 + 1 * nVar * 8], zmm15						// _mm512_store_pd(cptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 1 * nVar + i * VECTORLENGTH, zmm15);    // &c[it][iz][0][1][i]
			vmovapd		ZMMWORD PTR[r13 + 0 * ND1D*nVar * 8 + 2 * nVar * 8], zmm16    			        // _mm512_store_pd(cptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 2 * nVar + i * VECTORLENGTH, zmm16);    // &c[it][iz][0][2][i]
			vmovapd		ZMMWORD PTR[r13 + 0 * ND1D*nVar * 8 + 3 * nVar * 8], zmm17						// _mm512_store_pd(cptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 3 * nVar + i * VECTORLENGTH, zmm17);    // &c[it][iz][0][3][i]

			vfmadd231pd  zmm28, zmm19, zmm3							// zmm28 = _mm512_fmadd_pd(zmm19, zmm3, zmm28);	  		 
			vfmadd231pd  zmm29, zmm8, zmm10							// zmm29 = _mm512_fmadd_pd(zmm8, zmm10, zmm29);  		  				  
			vfnmadd231pd zmm30, zmm14, zmm12						// zmm30 = _mm512_fmadd_pd(zmm14, zmm21, zmm30);
			vfnmadd231pd zmm31, zmm18, zmm1							// zmm31 = _mm512_fmadd_pd(zmm18, zmm26, zmm31); 

			vmovapd		ZMMWORD PTR[r13 + 1 * ND1D*nVar * 8 + 0 * nVar * 8], zmm28    			        // _mm512_store_pd(cptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 0 * nVar + i * VECTORLENGTH, zmm28);    // &c[it][iz][1][0][i] 
			vmovapd		ZMMWORD PTR[r13 + 1 * ND1D*nVar * 8 + 1 * nVar * 8], zmm29						// _mm512_store_pd(cptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 1 * nVar + i * VECTORLENGTH, zmm29);    // &c[it][iz][1][1][i]
			vmovapd		ZMMWORD PTR[r13 + 1 * ND1D*nVar * 8 + 2 * nVar * 8], zmm30    			        // _mm512_store_pd(cptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 2 * nVar + i * VECTORLENGTH, zmm30);    // &c[it][iz][1][2][i]
			vmovapd		ZMMWORD PTR[r13 + 1 * ND1D*nVar * 8 + 3 * nVar * 8], zmm31						// _mm512_store_pd(cptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 3 * nVar + i * VECTORLENGTH, zmm31);    // &c[it][iz][1][3][i]

			// second half 

			vmovapd		zmm4,  ZMMWORD PTR[r12 + 2 * ND1D*nVar*8 + 0 * nVar*8]     					// zmm4 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 0 * nVar + i * VECTORLENGTH);  // &a[it][iz][0][0][i] 
			vmovapd		zmm5,  ZMMWORD PTR[r12 + 2 * ND1D*nVar*8 + 1 * nVar*8]						// zmm5 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 1 * nVar + i * VECTORLENGTH);  // &a[it][iz][0][1][i]  
			vmovapd		zmm6,  ZMMWORD PTR[r12 + 2 * ND1D*nVar*8 + 2 * nVar*8]     					// zmm6 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 2 * nVar + i * VECTORLENGTH);  // &a[it][iz][0][2][i] 
			vmovapd		zmm7,  ZMMWORD PTR[r12 + 2 * ND1D*nVar*8 + 3 * nVar*8]						// zmm7 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 3 * nVar + i * VECTORLENGTH);  // &a[it][iz][0][3][i] 
			
			vmovapd		zmm9, ZMMWORD PTR[r13 + 2 * ND1D*nVar * 8 + 0 * nVar * 8]		// zmm4 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 0 * nVar + i * VECTORLENGTH);  // &a[it][iz][0][0][i] 
			vmovapd		zmm15, ZMMWORD PTR[r13 + 2 * ND1D*nVar * 8 + 1 * nVar * 8]		// zmm5 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 1 * nVar + i * VECTORLENGTH);  // &a[it][iz][0][1][i]  
			vmovapd		zmm16, ZMMWORD PTR[r13 + 2 * ND1D*nVar * 8 + 2 * nVar * 8]		// zmm6 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 2 * nVar + i * VECTORLENGTH);  // &a[it][iz][0][2][i] 
			vmovapd		zmm17, ZMMWORD PTR[r13 + 2 * ND1D*nVar * 8 + 3 * nVar * 8]		// zmm7 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 3 * nVar + i * VECTORLENGTH);  // &a[it][iz][0][3][i] 

			vfmadd231pd   zmm9, zmm4,  zmm0							    // zmm9 = _mm512_fmadd_pd(zmm4, zmm0, zmm9);  // 09-4 		  				  		 
			vfmadd231pd   zmm15, zmm5, zmm11						// zmm15 = _mm512_fmadd_pd(zmm5, zmm11, zmm15); // 15-5			  				  
			vfnmadd231pd  zmm16, zmm6, zmm11						// zmm16 = _mm512_fmadd_pd(zmm6, zmm22, zmm16); // 16-6
			vfnmadd231pd  zmm17, zmm7, zmm0 						// zmm17 = _mm512_fmadd_pd(zmm7, zmm27, zmm17);   // 17-7 
				
			vmovapd		zmm8,  ZMMWORD PTR[r12 + 3 * ND1D*nVar*8 + 0 * nVar*8]     					// zmm8 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 0 * nVar + i * VECTORLENGTH); // &a[it][iz][1][0][i] 
			vmovapd		zmm14, ZMMWORD PTR[r12 + 3 * ND1D*nVar*8 + 1 * nVar*8]						// zmm14 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 1 * nVar + i * VECTORLENGTH); // & a[it][iz][1][1][i]
			vmovapd		zmm18, ZMMWORD PTR[r12 + 3 * ND1D*nVar*8 + 2 * nVar*8]     					// zmm18 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 2 * nVar + i * VECTORLENGTH); // &a[it][iz][1][2][i] 
			vmovapd		zmm19, ZMMWORD PTR[r12 + 3 * ND1D*nVar*8 + 3 * nVar*8]						// zmm19 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 3 * nVar + i * VECTORLENGTH);  // &a[it][iz][1][3][i] 
			
			vmovapd		zmm28, ZMMWORD PTR[r13 + 3 * ND1D*nVar * 8 + 0 * nVar * 8]     	// zmm8 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 0 * nVar + i * VECTORLENGTH); // &a[it][iz][1][0][i] 
			vmovapd		zmm29, ZMMWORD PTR[r13 + 3 * ND1D*nVar * 8 + 1 * nVar * 8]		// zmm14 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 1 * nVar + i * VECTORLENGTH); // & a[it][iz][1][1][i]
			vmovapd		zmm30, ZMMWORD PTR[r13 + 3 * ND1D*nVar * 8 + 2 * nVar * 8]     	// zmm18 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 2 * nVar + i * VECTORLENGTH); // &a[it][iz][1][2][i] 
			vmovapd		zmm31, ZMMWORD PTR[r13 + 3 * ND1D*nVar * 8 + 3 * nVar * 8]		// zmm19 = _mm512_load_pd(aptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 3 * nVar + i * VECTORLENGTH);  // &a[it][iz][1][3][i] 

			vfmadd231pd   zmm28, zmm8,  zmm0						// zmm28 = _mm512_fmadd_pd(zmm8, zmm0, zmm28);     // 28-4 			  		 
			vfmadd231pd   zmm29, zmm14, zmm11						// zmm29 = _mm512_fmadd_pd(zmm14, zmm11, zmm29);   // 29-5		 		  				  
			vfnmadd231pd  zmm30, zmm18, zmm11						// zmm30 = _mm512_fmadd_pd(zmm18, zmm22, zmm30);   // 30-6	 
			vfnmadd231pd  zmm31, zmm19, zmm0						// zmm31 = _mm512_fmadd_pd(zmm19, zmm27, zmm31);   // 31-7 

			prefetchw[r13 + 2 * ND1D*nVar * 8 + 0 * nVar * 8]
			prefetchw[r13 + 2 * ND1D*nVar * 8 + 1 * nVar * 8]
			prefetchw[r13 + 2 * ND1D*nVar * 8 + 2 * nVar * 8]
			prefetchw[r13 + 2 * ND1D*nVar * 8 + 3 * nVar * 8]
			prefetchw[r13 + 3 * ND1D*nVar * 8 + 0 * nVar * 8]
			prefetchw[r13 + 3 * ND1D*nVar * 8 + 1 * nVar * 8]
			prefetchw[r13 + 3 * ND1D*nVar * 8 + 2 * nVar * 8]
			prefetchw[r13 + 3 * ND1D*nVar * 8 + 3 * nVar * 8]

			vfmadd231pd  zmm9,  zmm5,  zmm1							// zmm9 = _mm512_fmadd_pd(zmm5, zmm1, zmm9);      // 09-5 		  				  		 
			vfmadd231pd  zmm15, zmm6,  zmm12						// zmm15 = _mm512_fmadd_pd(zmm6, zmm12, zmm15);    // 15-6	  				  
			vfnmadd231pd zmm16, zmm7,  zmm10						// zmm16 = _mm512_fmadd_pd(zmm7, zmm23, zmm16);    // 16-7
			vfnmadd231pd zmm17, zmm4,  zmm3 						// zmm17 = _mm512_fmadd_pd(zmm4, zmm24, zmm17);    // 17-4
			vfmadd231pd  zmm28, zmm14, zmm1							// zmm28 = _mm512_fmadd_pd(zmm14, zmm1, zmm28);    // 28-5		  		 
			vfmadd231pd  zmm29, zmm18, zmm12						// zmm29 = _mm512_fmadd_pd(zmm18, zmm12, zmm29);   // 29-6	 		  				  
			vfnmadd231pd zmm30, zmm19, zmm10						// zmm30 = _mm512_fmadd_pd(zmm19, zmm23, zmm30);   // 30-7
			vfnmadd231pd zmm31, zmm8,  zmm3							// zmm31 = _mm512_fmadd_pd(zmm8, zmm24, zmm31);    // 31-4

			vfmadd231pd  zmm9,  zmm6,  zmm2							// zmm9 = _mm512_fmadd_pd(zmm6, zmm2, zmm9);	    // 09-6		  				  		 
			vfmadd231pd  zmm15, zmm7,  zmm13						// zmm15 = _mm512_fmadd_pd(zmm7, zmm13, zmm15);    // 15-7	  				  
			vfnmadd231pd zmm16, zmm4,  zmm13						// zmm16 = _mm512_fmadd_pd(zmm4, zmm20, zmm16);    // 16-4
			vfnmadd231pd zmm17, zmm5,  zmm2							// zmm17 = _mm512_fmadd_pd(zmm5, zmm25, zmm17);    // 17-5
			vfmadd231pd  zmm28, zmm18, zmm2							// zmm28 = _mm512_fmadd_pd(zmm18, zmm2, zmm28);    // 28-6		  		 
			vfmadd231pd  zmm29, zmm19, zmm13						// zmm29 = _mm512_fmadd_pd(zmm19, zmm13, zmm29);   // 29-7	 		  				  
			vfnmadd231pd zmm30, zmm8,  zmm13						// zmm30 = _mm512_fmadd_pd(zmm8, zmm20, zmm30);    // 30-4
			vfnmadd231pd zmm31, zmm14, zmm2							// zmm31 = _mm512_fmadd_pd(zmm14, zmm25, zmm31);   // 31-5

			vfmadd231pd  zmm9, zmm7,  zmm3							// zmm9 = _mm512_fmadd_pd(zmm7, zmm3, zmm9);				  		 
			vfmadd231pd  zmm15, zmm4, zmm10							// zmm15 = _mm512_fmadd_pd(zmm4, zmm10, zmm15);				  
			vfnmadd231pd zmm16, zmm5, zmm12							// zmm16 = _mm512_fmadd_pd(zmm5, zmm21, zmm16);
			vfnmadd231pd zmm17, zmm6, zmm1 							// zmm17 = _mm512_fmadd_pd(zmm6, zmm26, zmm17);

			vmovapd		ZMMWORD PTR[r13 + 2 * ND1D*nVar*8 + 0 * nVar*8], zmm9     			        // _mm512_store_pd(cptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 0 * nVar + i * VECTORLENGTH, zmm9);     // &c[it][iz][0][0][i]
			vmovapd		ZMMWORD PTR[r13 + 2 * ND1D*nVar*8 + 1 * nVar*8], zmm15						// _mm512_store_pd(cptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 1 * nVar + i * VECTORLENGTH, zmm15);    // &c[it][iz][0][1][i]
			vmovapd		ZMMWORD PTR[r13 + 2 * ND1D*nVar*8 + 2 * nVar*8], zmm16    			        // _mm512_store_pd(cptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 2 * nVar + i * VECTORLENGTH, zmm16);    // &c[it][iz][0][2][i]
			vmovapd		ZMMWORD PTR[r13 + 2 * ND1D*nVar*8 + 3 * nVar*8], zmm17						// _mm512_store_pd(cptr + it * ND3D*nVar + iz * ND2D*nVar + 0 * ND1D*nVar + 3 * nVar + i * VECTORLENGTH, zmm17);    // &c[it][iz][0][3][i]

			vfmadd231pd  zmm28, zmm19, zmm3							// zmm28 = _mm512_fmadd_pd(zmm19, zmm3, zmm28);	  		 
			vfmadd231pd  zmm29, zmm8, zmm10							// zmm29 = _mm512_fmadd_pd(zmm8, zmm10, zmm29);  		  				  
			vfnmadd231pd zmm30, zmm14, zmm12						// zmm30 = _mm512_fmadd_pd(zmm14, zmm21, zmm30);
			vfnmadd231pd zmm31, zmm18, zmm1							// zmm31 = _mm512_fmadd_pd(zmm18, zmm26, zmm31); 
		
			vmovapd		ZMMWORD PTR[r13 + 3 * ND1D*nVar*8 + 0 * nVar*8], zmm28    			        // _mm512_store_pd(cptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 0 * nVar + i * VECTORLENGTH, zmm28);    // &c[it][iz][1][0][i] 
			vmovapd		ZMMWORD PTR[r13 + 3 * ND1D*nVar*8 + 1 * nVar*8], zmm29						// _mm512_store_pd(cptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 1 * nVar + i * VECTORLENGTH, zmm29);    // &c[it][iz][1][1][i]
			vmovapd		ZMMWORD PTR[r13 + 3 * ND1D*nVar*8 + 2 * nVar*8], zmm30    			        // _mm512_store_pd(cptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 2 * nVar + i * VECTORLENGTH, zmm30);    // &c[it][iz][1][2][i]
			vmovapd		ZMMWORD PTR[r13 + 3 * ND1D*nVar*8 + 3 * nVar*8], zmm31						// _mm512_store_pd(cptr + it * ND3D*nVar + iz * ND2D*nVar + 1 * ND1D*nVar + 3 * nVar + i * VECTORLENGTH, zmm31);    // &c[it][iz][1][3][i]

		add r15, VECTORLENGTH*8 
		cmp r15, nVar*8  
		jl loopvar

		add r10, 1 // ND2D*nVar*8   
		cmp r10, 4 // 4*ND2D*nVar*8 
		jl loopz

		add r8, 1  // ND3D*nVar*8
		cmp r8, 4  // 4*ND3D*nVar*8 
		jl loopt 



	}
	
	*/ 

	int i, j, k, it, iy, iz;
	double weight[1];
	__m512d zmm0, zmm1, zmm2, zmm3, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15;
	__m512d zmm20, zmm21, zmm22, zmm23, zmm24, zmm25, zmm26, zmm27; 
	__m512d zmm4, zmm5, zmm6, zmm7, zmm8, zmm9, zmm16, zmm17, zmm18, zmm19, zmm28, zmm29, zmm30, zmm31;
	__m128d xmm;

	// load the scalar into register 14 
	xmm = _mm_load1_pd(&s[0]);
	zmm14 = _mm512_broadcastsd_pd(xmm);
	// load first column of the small matrix b and scale with the scalar s 
	xmm = _mm_load1_pd(&b[0][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm0 = _mm512_mul_pd(zmm4, zmm14);
	zmm1 = _mm512_mul_pd(zmm5, zmm14);
	zmm2 = _mm512_mul_pd(zmm6, zmm14);
	zmm3 = _mm512_mul_pd(zmm7, zmm14);
	// load second column of the small matrix b and scale with the scalar s  
	xmm = _mm_load1_pd(&b[1][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm10 = _mm512_mul_pd(zmm4, zmm14);
	zmm11 = _mm512_mul_pd(zmm5, zmm14);
	zmm12 = _mm512_mul_pd(zmm6, zmm14);
	zmm13 = _mm512_mul_pd(zmm7, zmm14);
	// load third column of the small matrix b and scale with the scalar s 
	xmm = _mm_load1_pd(&b[2][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm20 = _mm512_mul_pd(zmm4, zmm14);
	zmm21 = _mm512_mul_pd(zmm5, zmm14);
	zmm22 = _mm512_mul_pd(zmm6, zmm14);
	zmm23 = _mm512_mul_pd(zmm7, zmm14);
	// load fourth column of the small matrix b and scale with the scalar s  
	xmm = _mm_load1_pd(&b[3][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm24 = _mm512_mul_pd(zmm4, zmm14);
	zmm25 = _mm512_mul_pd(zmm5, zmm14);
	zmm26 = _mm512_mul_pd(zmm6, zmm14);
	zmm27 = _mm512_mul_pd(zmm7, zmm14);

	for (it = 0; it < 4; it++)
	{
		for (iz = 0; iz < 4; iz++)
		{
			for (iy = 0; iy < 4; iy++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of matrix a on position i 
					// load c into the result registers 
					// multiply the row of a with the first and second column of b and add the result to c using FMA 
					// multiply the row of a with the third and fourth column of b and add the result to c using FMA 
					zmm4 = _mm512_load_pd(&a[it][iz][iy][0][i]);
					zmm9 = _mm512_load_pd(&c[it][iz][iy][0][i]);
					zmm9 = _mm512_fmadd_pd(zmm4, zmm0, zmm9);
					zmm5 = _mm512_load_pd(&a[it][iz][iy][1][i]);
					zmm15 = _mm512_load_pd(&c[it][iz][iy][1][i]);
					zmm15 = _mm512_fmadd_pd(zmm5, zmm11, zmm15);
					zmm6 = _mm512_load_pd(&a[it][iz][iy][2][i]);
					zmm16 = _mm512_load_pd(&c[it][iz][iy][2][i]);
					zmm16 = _mm512_fmadd_pd(zmm6, zmm22, zmm16);
					zmm7 = _mm512_load_pd(&a[it][iz][iy][3][i]);
					zmm17 = _mm512_load_pd(&c[it][iz][iy][3][i]);
					zmm17 = _mm512_fmadd_pd(zmm7, zmm27, zmm17);
					// 
					zmm9 = _mm512_fmadd_pd(zmm5, zmm1, zmm9);
					zmm15 = _mm512_fmadd_pd(zmm6, zmm12, zmm15);
					zmm16 = _mm512_fmadd_pd(zmm7, zmm23, zmm16);
					zmm17 = _mm512_fmadd_pd(zmm4, zmm24, zmm17);
					// 
					zmm9 = _mm512_fmadd_pd(zmm6, zmm2, zmm9);
					zmm15 = _mm512_fmadd_pd(zmm7, zmm13, zmm15);
					zmm16 = _mm512_fmadd_pd(zmm4, zmm20, zmm16);
					zmm17 = _mm512_fmadd_pd(zmm5, zmm25, zmm17);
					// store the result in c 
					zmm9 = _mm512_fmadd_pd(zmm7, zmm3, zmm9);
					_mm512_store_pd(&c[it][iz][iy][0][i], zmm9);
  					zmm15 = _mm512_fmadd_pd(zmm4, zmm10, zmm15);
					_mm512_store_pd(&c[it][iz][iy][1][i], zmm15);
					zmm16 = _mm512_fmadd_pd(zmm5, zmm21, zmm16);
					_mm512_store_pd(&c[it][iz][iy][2][i], zmm16);
					zmm17 = _mm512_fmadd_pd(zmm6, zmm26, zmm17);
					_mm512_store_pd(&c[it][iz][iy][3][i], zmm17);
				}
			}
		}
	}

	 

#else 

	__assume_aligned(c, 32);
	__assume_aligned(a, 32);
	__assume_aligned(b, 32);
	int i, j, k, iy, iz, it;
	double weight[1];
	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	__m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;

	// load the scalar into register 14 
	ymm14 = _mm256_broadcast_sd(&s[0]);
	// load first column of the small matrix b and scale with scalar s 
	ymm4 = _mm256_broadcast_sd(&b[0][0]);
	ymm5 = _mm256_broadcast_sd(&b[0][1]);
	ymm6 = _mm256_broadcast_sd(&b[0][2]);
	ymm7 = _mm256_broadcast_sd(&b[0][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load second column of the small matrix b and scale with scalar s  
	ymm4 = _mm256_broadcast_sd(&b[1][0]);
	ymm5 = _mm256_broadcast_sd(&b[1][1]);
	ymm6 = _mm256_broadcast_sd(&b[1][2]);
	ymm7 = _mm256_broadcast_sd(&b[1][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);

	for (it = 0; it < 4; it++) 
	{
		for (iz = 0; iz < 4; iz++)
		{
			for (iy = 0; iy < 4; iy++)
			{
				// optimized AVX code for matmul on all tensor slices 
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of matrix a on position i 
					ymm4 = _mm256_load_pd(&a[it][iz][iy][0][i]);
					ymm5 = _mm256_load_pd(&a[it][iz][iy][1][i]);
					ymm6 = _mm256_load_pd(&a[it][iz][iy][2][i]);
					ymm7 = _mm256_load_pd(&a[it][iz][iy][3][i]);
					// set the result register 1 to zero 
					ymm9 = _mm256_load_pd(&c[it][iz][iy][0][i]); //  _mm256_setzero_pd();
					// set the result register 2 to zero
					ymm15 = _mm256_load_pd(&c[it][iz][iy][1][i]);  // _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the first and second column of b and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
					// multiply the row of a with the first and second column of b and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 
					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result in c 
					_mm256_store_pd(&c[it][iz][iy][0][i], ymm9);
					_mm256_store_pd(&c[it][iz][iy][1][i], ymm15);
				}
			}
		}
	}

	ymm14 = _mm256_broadcast_sd(&s[0]);
	// load third column of the small matrix b 
	ymm4 = _mm256_broadcast_sd(&b[2][0]);
	ymm5 = _mm256_broadcast_sd(&b[2][1]);
	ymm6 = _mm256_broadcast_sd(&b[2][2]);
	ymm7 = _mm256_broadcast_sd(&b[2][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load fourth column of the small matrix b 
	ymm4 = _mm256_broadcast_sd(&b[3][0]);
	ymm5 = _mm256_broadcast_sd(&b[3][1]);
	ymm6 = _mm256_broadcast_sd(&b[3][2]);
	ymm7 = _mm256_broadcast_sd(&b[3][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);


	for (it = 0; it < 4; it++)
	{
		for (iz = 0; iz < 4; iz++)
		{
			for (iy = 0; iy < 4; iy++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of a on position i 
					ymm4 = _mm256_load_pd(&a[it][iz][iy][0][i]);
					ymm5 = _mm256_load_pd(&a[it][iz][iy][1][i]);
					ymm6 = _mm256_load_pd(&a[it][iz][iy][2][i]);
					ymm7 = _mm256_load_pd(&a[it][iz][iy][3][i]);
					// set the result register 1 to zero 
					ymm9 = _mm256_load_pd(&c[it][iz][iy][2][i]);   // _mm256_setzero_pd();
					// set the result register 2 to zero
					ymm15 = _mm256_load_pd(&c[it][iz][iy][3][i]);  // _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the first and second column of b and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
					// multiply the row of a with the first and second column of b and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 

					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result
					_mm256_store_pd(&c[it][iz][iy][2][i], ymm9);
					_mm256_store_pd(&c[it][iz][iy][3][i], ymm15);
				}
			}
		}
	}
#endif 
}

// compute all matmuls in y direction on a 4D data set with vector weights (very useful for stiffness matrix multiplications) 
extern "C" void __cdecl vsmm4x4c4d_ky(double c[4][4][4][4][nVar], double a[4][4][4][4][nVar], double b[4][4], double s[1] )
{
#ifdef AVX512
	__assume_aligned(c, 64);
	__assume_aligned(a, 64);
	__assume_aligned(b, 64);
	int i, j, k, it, ix, iz;
	double weight[1];
	__m512d zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7;
	__m512d zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15;
	__m512d zmm16, zmm17;
	__m512d zmm20, zmm21, zmm22, zmm23, zmm24, zmm25, zmm26, zmm27;
	__m128d xmm;

	// load the scalar into register 14 
	xmm = _mm_load1_pd(&s[0]);
	zmm14 = _mm512_broadcastsd_pd(xmm);
	// load first column of the small matrix b and scale with the scalar s 
	xmm = _mm_load1_pd(&b[0][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm0 = _mm512_mul_pd(zmm4, zmm14);
	zmm1 = _mm512_mul_pd(zmm5, zmm14);
	zmm2 = _mm512_mul_pd(zmm6, zmm14);
	zmm3 = _mm512_mul_pd(zmm7, zmm14);
	// load second column of the small matrix b and scale with the scalar s  
	xmm = _mm_load1_pd(&b[1][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm10 = _mm512_mul_pd(zmm4, zmm14);
	zmm11 = _mm512_mul_pd(zmm5, zmm14);
	zmm12 = _mm512_mul_pd(zmm6, zmm14);
	zmm13 = _mm512_mul_pd(zmm7, zmm14);
	// load third column of the small matrix b and scale with the scalar s 
	xmm = _mm_load1_pd(&b[2][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm20 = _mm512_mul_pd(zmm4, zmm14);
	zmm21 = _mm512_mul_pd(zmm5, zmm14);
	zmm22 = _mm512_mul_pd(zmm6, zmm14);
	zmm23 = _mm512_mul_pd(zmm7, zmm14);
	// load fourth column of the small matrix b and scale with the scalar s  
	xmm = _mm_load1_pd(&b[3][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm24 = _mm512_mul_pd(zmm4, zmm14);
	zmm25 = _mm512_mul_pd(zmm5, zmm14);
	zmm26 = _mm512_mul_pd(zmm6, zmm14);
	zmm27 = _mm512_mul_pd(zmm7, zmm14);

	for (it = 0; it < 4; it++)
	{
		for (iz = 0; iz < 4; iz++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of matrix a on position i 
					// load c into the result registers 
					// multiply the row of a with the first and second column of b and add the result to c using FMA 
					// multiply the row of a with the third and fourth column of b and add the result to c using FMA 
					zmm4 = _mm512_load_pd(&a[it][iz][0][ix][i]);
					zmm9 = _mm512_load_pd(&c[it][iz][0][ix][i]);
					zmm9 = _mm512_fmadd_pd(zmm4, zmm0, zmm9);
					zmm5 = _mm512_load_pd(&a[it][iz][1][ix][i]);
					zmm15 = _mm512_load_pd(&c[it][iz][1][ix][i]);
					zmm15 = _mm512_fmadd_pd(zmm5, zmm11, zmm15);
					zmm6 = _mm512_load_pd(&a[it][iz][2][ix][i]);
					zmm16 = _mm512_load_pd(&c[it][iz][2][ix][i]);
					zmm16 = _mm512_fmadd_pd(zmm6, zmm22, zmm16);
					zmm7 = _mm512_load_pd(&a[it][iz][3][ix][i]);
					zmm17 = _mm512_load_pd(&c[it][iz][3][ix][i]);
					zmm17 = _mm512_fmadd_pd(zmm7, zmm27, zmm17);
					// 
					zmm9 = _mm512_fmadd_pd(zmm5, zmm1, zmm9);
					zmm15 = _mm512_fmadd_pd(zmm6, zmm12, zmm15);
					zmm16 = _mm512_fmadd_pd(zmm7, zmm23, zmm16);
					zmm17 = _mm512_fmadd_pd(zmm4, zmm24, zmm17);
					// 
					zmm9 = _mm512_fmadd_pd(zmm6, zmm2, zmm9);
					zmm15 = _mm512_fmadd_pd(zmm7, zmm13, zmm15);
					zmm16 = _mm512_fmadd_pd(zmm4, zmm20, zmm16);
					zmm17 = _mm512_fmadd_pd(zmm5, zmm25, zmm17);
					// store the result in c 
					zmm9 = _mm512_fmadd_pd(zmm7, zmm3, zmm9);
					_mm512_store_pd(&c[it][iz][0][ix][i], zmm9);
					zmm15 = _mm512_fmadd_pd(zmm4, zmm10, zmm15);
					_mm512_store_pd(&c[it][iz][1][ix][i], zmm15);
					zmm16 = _mm512_fmadd_pd(zmm5, zmm21, zmm16);
					_mm512_store_pd(&c[it][iz][2][ix][i], zmm16);
					zmm17 = _mm512_fmadd_pd(zmm6, zmm26, zmm17);
					_mm512_store_pd(&c[it][iz][3][ix][i], zmm17);
				}
			}
		}
	}

#else 
	__assume_aligned(c, 32);
	__assume_aligned(a, 32);
	__assume_aligned(b, 32);
	int i, j, k, ix, iz, it;
	double weight[1];
	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	__m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;

	// load the scalar into register 14 
	ymm14 = _mm256_broadcast_sd(&s[0]);
	// load first column of the small matrix b and scale with scalar s 
	ymm4 = _mm256_broadcast_sd(&b[0][0]);
	ymm5 = _mm256_broadcast_sd(&b[0][1]);
	ymm6 = _mm256_broadcast_sd(&b[0][2]);
	ymm7 = _mm256_broadcast_sd(&b[0][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load second column of the small matrix b and scale with scalar s  
	ymm4 = _mm256_broadcast_sd(&b[1][0]);
	ymm5 = _mm256_broadcast_sd(&b[1][1]);
	ymm6 = _mm256_broadcast_sd(&b[1][2]);
	ymm7 = _mm256_broadcast_sd(&b[1][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);

	for (it = 0; it < 4; it++)
	{
		for (iz = 0; iz < 4; iz++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				// optimized AVX code for matmul on all tensor slices 
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of matrix a on position i 
					ymm4 = _mm256_load_pd(&a[it][iz][0][ix][i]);
					ymm5 = _mm256_load_pd(&a[it][iz][1][ix][i]);
					ymm6 = _mm256_load_pd(&a[it][iz][2][ix][i]);
					ymm7 = _mm256_load_pd(&a[it][iz][3][ix][i]);
					// load current value of the the result into result register 1  
					ymm9 = _mm256_load_pd(&c[it][iz][0][ix][i]);    //  _mm256_setzero_pd();
					// set the result register 2 to zero
					ymm15 = _mm256_load_pd(&c[it][iz][1][ix][i]);   // _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the first and second column of b and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
					// multiply the row of a with the first and second column of b and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 
					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result in c 
					_mm256_store_pd(&c[it][iz][0][ix][i], ymm9);
					_mm256_store_pd(&c[it][iz][1][ix][i], ymm15);
				}
			}
		}
	}

	ymm14 = _mm256_broadcast_sd(&s[0]);
	// load third column of the small matrix b 
	ymm4 = _mm256_broadcast_sd(&b[2][0]);
	ymm5 = _mm256_broadcast_sd(&b[2][1]);
	ymm6 = _mm256_broadcast_sd(&b[2][2]);
	ymm7 = _mm256_broadcast_sd(&b[2][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load fourth column of the small matrix b 
	ymm4 = _mm256_broadcast_sd(&b[3][0]);
	ymm5 = _mm256_broadcast_sd(&b[3][1]);
	ymm6 = _mm256_broadcast_sd(&b[3][2]);
	ymm7 = _mm256_broadcast_sd(&b[3][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);

	for (it = 0; it < 4; it++)
	{
		for (iz = 0; iz < 4; iz++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of a on position i 
					ymm4 = _mm256_load_pd(&a[it][iz][0][ix][i]);
					ymm5 = _mm256_load_pd(&a[it][iz][1][ix][i]);
					ymm6 = _mm256_load_pd(&a[it][iz][2][ix][i]);
					ymm7 = _mm256_load_pd(&a[it][iz][3][ix][i]);
					// load the result register 1  
					ymm9 = _mm256_load_pd(&c[it][iz][2][ix][i]);  // _mm256_setzero_pd();
					// load the result register 2  
					ymm15 = _mm256_load_pd(&c[it][iz][3][ix][i]);  // _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the first and second column of b and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
					// multiply the row of a with the first and second column of b and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 

					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result
					_mm256_store_pd(&c[it][iz][2][ix][i], ymm9);
					_mm256_store_pd(&c[it][iz][3][ix][i], ymm15);
				}
			}
		}
	}
#endif 
}


// compute all matmuls in z direction on a 4D data set with vector weights (very useful for stiffness matrix multiplications) 
extern "C" void __cdecl vsmm4x4c4d_kz(double c[4][4][4][4][nVar], double a[4][4][4][4][nVar], double b[4][4], double s[1], double wGP[4])
{
#ifdef AVX512
	__assume_aligned(c, 64);
	__assume_aligned(a, 64);
	__assume_aligned(b, 64);
	int i, j, k, ix, iy, it;
	double weight[1];
	__m512d zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7;
	__m512d zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15;
	__m512d zmm16, zmm17;
	__m512d zmm20, zmm21, zmm22, zmm23, zmm24, zmm25, zmm26, zmm27;
	__m128d xmm;

	// load the scalar into register 14 
	xmm = _mm_load1_pd(&s[0]);
	zmm14 = _mm512_broadcastsd_pd(xmm);
	// load first column of the small matrix b and scale with the scalar s 
	xmm = _mm_load1_pd(&b[0][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm0 = _mm512_mul_pd(zmm4, zmm14);
	zmm1 = _mm512_mul_pd(zmm5, zmm14);
	zmm2 = _mm512_mul_pd(zmm6, zmm14);
	zmm3 = _mm512_mul_pd(zmm7, zmm14);
	// load second column of the small matrix b and scale with the scalar s  
	xmm = _mm_load1_pd(&b[1][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm10 = _mm512_mul_pd(zmm4, zmm14);
	zmm11 = _mm512_mul_pd(zmm5, zmm14);
	zmm12 = _mm512_mul_pd(zmm6, zmm14);
	zmm13 = _mm512_mul_pd(zmm7, zmm14);
	// load third column of the small matrix b and scale with the scalar s 
	xmm = _mm_load1_pd(&b[2][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm20 = _mm512_mul_pd(zmm4, zmm14);
	zmm21 = _mm512_mul_pd(zmm5, zmm14);
	zmm22 = _mm512_mul_pd(zmm6, zmm14);
	zmm23 = _mm512_mul_pd(zmm7, zmm14);
	// load fourth column of the small matrix b and scale with the scalar s  
	xmm = _mm_load1_pd(&b[3][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm24 = _mm512_mul_pd(zmm4, zmm14);
	zmm25 = _mm512_mul_pd(zmm5, zmm14);
	zmm26 = _mm512_mul_pd(zmm6, zmm14);
	zmm27 = _mm512_mul_pd(zmm7, zmm14);

	for (it = 0; it < 4; it++)
	{
		for (iy = 0; iy < 4; iy++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				// optimized AVX code for matmul on all tensor slices 
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of matrix a on position i 
					// load c into the result registers 
					// multiply the row of a with the first and second column of b and add the result to c using FMA 
					// multiply the row of a with the third and fourth column of b and add the result to c using FMA 
					zmm4 = _mm512_load_pd(&a[it][0][iy][ix][i]);
					zmm9 = _mm512_load_pd(&c[it][0][iy][ix][i]);
					zmm9 = _mm512_fmadd_pd(zmm4, zmm0, zmm9);
					zmm5 = _mm512_load_pd(&a[it][1][iy][ix][i]);
					zmm15 = _mm512_load_pd(&c[it][1][iy][ix][i]);
					zmm15 = _mm512_fmadd_pd(zmm5, zmm11, zmm15);
					zmm6 = _mm512_load_pd(&a[it][2][iy][ix][i]);
					zmm16 = _mm512_load_pd(&c[it][2][iy][ix][i]);
					zmm16 = _mm512_fmadd_pd(zmm6, zmm22, zmm16);
					zmm7 = _mm512_load_pd(&a[it][3][iy][ix][i]);
					zmm17 = _mm512_load_pd(&c[it][3][iy][ix][i]);
					zmm17 = _mm512_fmadd_pd(zmm7, zmm27, zmm17);
					// 
					zmm9 = _mm512_fmadd_pd(zmm5, zmm1, zmm9);
					zmm15 = _mm512_fmadd_pd(zmm6, zmm12, zmm15);
					zmm16 = _mm512_fmadd_pd(zmm7, zmm23, zmm16);
					zmm17 = _mm512_fmadd_pd(zmm4, zmm24, zmm17);
					// 
					zmm9 = _mm512_fmadd_pd(zmm6, zmm2, zmm9);
					zmm15 = _mm512_fmadd_pd(zmm7, zmm13, zmm15);
					zmm16 = _mm512_fmadd_pd(zmm4, zmm20, zmm16);
					zmm17 = _mm512_fmadd_pd(zmm5, zmm25, zmm17);
					// store the result in c 
					zmm9 = _mm512_fmadd_pd(zmm7, zmm3, zmm9);
					_mm512_store_pd(&c[it][0][iy][ix][i], zmm9);
					zmm15 = _mm512_fmadd_pd(zmm4, zmm10, zmm15);
					_mm512_store_pd(&c[it][1][iy][ix][i], zmm15);
					zmm16 = _mm512_fmadd_pd(zmm5, zmm21, zmm16);
					_mm512_store_pd(&c[it][2][iy][ix][i], zmm16);
					zmm17 = _mm512_fmadd_pd(zmm6, zmm26, zmm17);
					_mm512_store_pd(&c[it][3][iy][ix][i], zmm17);
				}
			}
		}
	}

#else 
	__assume_aligned(c, 32);
	__assume_aligned(a, 32);
	__assume_aligned(b, 32);
	int i, j, k, ix, iy, it;
	double weight[1];
	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	__m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;

	// load the scalar into register 14 
	ymm14 = _mm256_broadcast_sd(&s[0]);
	// load first column of the small matrix b and scale with scalar s 
	ymm4 = _mm256_broadcast_sd(&b[0][0]);
	ymm5 = _mm256_broadcast_sd(&b[0][1]);
	ymm6 = _mm256_broadcast_sd(&b[0][2]);
	ymm7 = _mm256_broadcast_sd(&b[0][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load second column of the small matrix b and scale with scalar s  
	ymm4 = _mm256_broadcast_sd(&b[1][0]);
	ymm5 = _mm256_broadcast_sd(&b[1][1]);
	ymm6 = _mm256_broadcast_sd(&b[1][2]);
	ymm7 = _mm256_broadcast_sd(&b[1][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);

	for (it = 0; it < 4; it++)
	{
		for (iy = 0; iy < 4; iy++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				// optimized AVX code for matmul on all tensor slices 

				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of matrix a on position i 
					ymm4 = _mm256_load_pd(&a[it][0][iy][ix][i]);
					ymm5 = _mm256_load_pd(&a[it][1][iy][ix][i]);
					ymm6 = _mm256_load_pd(&a[it][2][iy][ix][i]);
					ymm7 = _mm256_load_pd(&a[it][3][iy][ix][i]);
					// set the result register 1 to zero 
					ymm9 = _mm256_load_pd(&c[it][0][iy][ix][i]);  // _mm256_setzero_pd();
					// set the result register 2 to zero
					ymm15 = _mm256_load_pd(&c[it][1][iy][ix][i]);  // _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the first and second column of b and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
					// multiply the row of a with the first and second column of b and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 
					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result in c 
					_mm256_store_pd(&c[it][0][iy][ix][i], ymm9);
					_mm256_store_pd(&c[it][1][iy][ix][i], ymm15);
				}
			}
		}
	}

	ymm14 = _mm256_broadcast_sd(&s[0]);
	// load third column of the small matrix b 
	ymm4 = _mm256_broadcast_sd(&b[2][0]);
	ymm5 = _mm256_broadcast_sd(&b[2][1]);
	ymm6 = _mm256_broadcast_sd(&b[2][2]);
	ymm7 = _mm256_broadcast_sd(&b[2][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load fourth column of the small matrix b 
	ymm4 = _mm256_broadcast_sd(&b[3][0]);
	ymm5 = _mm256_broadcast_sd(&b[3][1]);
	ymm6 = _mm256_broadcast_sd(&b[3][2]);
	ymm7 = _mm256_broadcast_sd(&b[3][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);

	for (it = 0; it < 4; it++)
	{
		for (iy = 0; iy < 4; iy++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of a on position i 
					ymm4 = _mm256_load_pd(&a[it][0][iy][ix][i]);
					ymm5 = _mm256_load_pd(&a[it][1][iy][ix][i]);
					ymm6 = _mm256_load_pd(&a[it][2][iy][ix][i]);
					ymm7 = _mm256_load_pd(&a[it][3][iy][ix][i]);
					// set the result register 1 to zero 
					ymm9 = _mm256_load_pd(&c[it][2][iy][ix][i]);  // _mm256_setzero_pd();
					// set the result register 2 to zero
					ymm15 = _mm256_load_pd(&c[it][3][iy][ix][i]); // _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the first and second column of b and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
					// multiply the row of a with the first and second column of b and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 

					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result
					_mm256_store_pd(&c[it][2][iy][ix][i], ymm9);
					_mm256_store_pd(&c[it][3][iy][ix][i], ymm15);
				}
			}
		}
	}
#endif 
}


// compute all matmuls in t direction on a 4D data set with vector weights (very useful for stiffness matrix multiplications) 
extern "C" void __cdecl vsmm4x4c4d_kt(double c[4][4][4][4][nVar], double a[4][4][4][4][nVar], double b[4][4])
{

#ifdef AVX512

	__assume_aligned(c, 64);
	__assume_aligned(a, 64);
	__assume_aligned(b, 64);
	int i, j, k, ix, iy, iz;
	double weight[1]; 
	__m512d zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7;
	__m512d zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15;
	__m512d zmm16, zmm17;
	__m512d zmm20, zmm21, zmm22, zmm23, zmm24, zmm25, zmm26, zmm27;
	__m128d xmm;

	// optimized AVX code for matmul on all tensor slices 
	// load first column of the small matrix b and scale with the scalar s 
	xmm = _mm_load1_pd(&b[0][0]);
	zmm0 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][1]);
	zmm1 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][2]);
	zmm2 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][3]);
	zmm3 = _mm512_broadcastsd_pd(xmm);
	// load second column of the small matrix b and scale with the scalar s  
	xmm = _mm_load1_pd(&b[1][0]);
	zmm10 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][1]);
	zmm11 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][2]);
	zmm12 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][3]);
	zmm13 = _mm512_broadcastsd_pd(xmm);
	// load third column of the small matrix b and scale with the scalar s 
	xmm = _mm_load1_pd(&b[2][0]);
	zmm20 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][1]);
	zmm21 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][2]);
	zmm22 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][3]);
	zmm23 = _mm512_broadcastsd_pd(xmm);
	// load fourth column of the small matrix b and scale with the scalar s  
	xmm = _mm_load1_pd(&b[3][0]);
	zmm24 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][1]);
	zmm25 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][2]);
	zmm26 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][3]);
	zmm27 = _mm512_broadcastsd_pd(xmm);


	for (iz = 0; iz < 4; iz++) 
	{
		for (iy = 0; iy < 4; iy++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// set the result registers to zero 
					// load VECTORLENGTH rows of matrix a on position i 
					// multiply the row of a with the first and second column of b and add the result to c using FMA 
					// multiply the row of a with the third and fourth column of b and add the result to c using FMA 
					zmm9 = _mm512_setzero_pd();
					zmm4 = _mm512_load_pd(&a[0][iz][iy][ix][i]);
					zmm9 = _mm512_fmadd_pd(zmm4, zmm0, zmm9);      // 09-4 
					zmm15 = _mm512_setzero_pd();
					zmm5 = _mm512_load_pd(&a[1][iz][iy][ix][i]);
					zmm15 = _mm512_fmadd_pd(zmm5, zmm11, zmm15);   // 15-5
					zmm16 = _mm512_setzero_pd();
					zmm6 = _mm512_load_pd(&a[2][iz][iy][ix][i]);
					zmm16 = _mm512_fmadd_pd(zmm6, zmm22, zmm16);   // 16-6
					zmm17 = _mm512_setzero_pd();
					zmm7 = _mm512_load_pd(&a[3][iz][iy][ix][i]);
					zmm17 = _mm512_fmadd_pd(zmm7, zmm27, zmm17);   // 17-7
					// 
					zmm9 = _mm512_fmadd_pd(zmm5, zmm1, zmm9);      // 09-5
					zmm15 = _mm512_fmadd_pd(zmm6, zmm12, zmm15);   // 15-6
					zmm16 = _mm512_fmadd_pd(zmm7, zmm23, zmm16);   // 16-7
					zmm17 = _mm512_fmadd_pd(zmm4, zmm24, zmm17);   // 17-4 
					// 
					zmm9 = _mm512_fmadd_pd(zmm6, zmm2, zmm9);	   // 09-6
					zmm15 = _mm512_fmadd_pd(zmm7, zmm13, zmm15);   // 15-7
					zmm16 = _mm512_fmadd_pd(zmm4, zmm20, zmm16);   // 16-4
					zmm17 = _mm512_fmadd_pd(zmm5, zmm25, zmm17);   // 17-5
					// 
					// store the result in c 
					zmm9 = _mm512_fmadd_pd(zmm7, zmm3, zmm9);	   // 09-7 
					_mm512_store_pd(&c[0][iz][iy][ix][i], zmm9);
					zmm15 = _mm512_fmadd_pd(zmm4, zmm10, zmm15);   // 15-4
					_mm512_store_pd(&c[1][iz][iy][ix][i], zmm15);
					zmm16 = _mm512_fmadd_pd(zmm5, zmm21, zmm16);   // 16-5
					_mm512_store_pd(&c[2][iz][iy][ix][i], zmm16);
					zmm17 = _mm512_fmadd_pd(zmm6, zmm26, zmm17);   // 17-6
					_mm512_store_pd(&c[3][iz][iy][ix][i], zmm17);
				}
			}
		}
	}


#else 

	__assume_aligned(c, 32);
	__assume_aligned(a, 32);
	__assume_aligned(b, 32);
	int i, j, k, ix, iy, iz;
	double weight[1]; 
	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	__m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;

	// load first column of the small matrix b and scale with scalar s 
	ymm0 = _mm256_broadcast_sd(&b[0][0]);
	ymm1 = _mm256_broadcast_sd(&b[0][1]);
	ymm2 = _mm256_broadcast_sd(&b[0][2]);
	ymm3 = _mm256_broadcast_sd(&b[0][3]);
	// load second column of the small matrix b and scale with scalar s  
	ymm10 = _mm256_broadcast_sd(&b[1][0]);
	ymm11 = _mm256_broadcast_sd(&b[1][1]);
	ymm12 = _mm256_broadcast_sd(&b[1][2]);
	ymm13 = _mm256_broadcast_sd(&b[1][3]);

	for (iz = 0; iz < 4; iz++)
	{
		for (iy = 0; iy < 4; iy++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				// optimized AVX code for matmul on all tensor slices 
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of matrix a on position i 
					ymm4 = _mm256_load_pd(&a[0][iz][iy][ix][i]);
					ymm5 = _mm256_load_pd(&a[1][iz][iy][ix][i]);
					ymm6 = _mm256_load_pd(&a[2][iz][iy][ix][i]);
					ymm7 = _mm256_load_pd(&a[3][iz][iy][ix][i]);
					// set the result register 1 to zero 
					ymm9 = _mm256_setzero_pd();
					// set the result register 2 to zero
					ymm15 = _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the first and second column of b and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
					// multiply the row of a with the first and second column of b and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 
					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result in c 
					_mm256_store_pd(&c[0][iz][iy][ix][i], ymm9);
					_mm256_store_pd(&c[1][iz][iy][ix][i], ymm15);
				}
			}
		}
	}

	// load third column of the small matrix b and scale with scalar s 
	ymm0 = _mm256_broadcast_sd(&b[2][0]);
	ymm1 = _mm256_broadcast_sd(&b[2][1]);
	ymm2 = _mm256_broadcast_sd(&b[2][2]);
	ymm3 = _mm256_broadcast_sd(&b[2][3]);
	// load fourth column of the small matrix b and scale with scalar s  
	ymm10 = _mm256_broadcast_sd(&b[3][0]);
	ymm11 = _mm256_broadcast_sd(&b[3][1]);
	ymm12 = _mm256_broadcast_sd(&b[3][2]);
	ymm13 = _mm256_broadcast_sd(&b[3][3]);

	for (iz = 0; iz < 4; iz++)
	{
		for (iy = 0; iy < 4; iy++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of a on position i 
					ymm4 = _mm256_load_pd(&a[0][iz][iy][ix][i]);
					ymm5 = _mm256_load_pd(&a[1][iz][iy][ix][i]);
					ymm6 = _mm256_load_pd(&a[2][iz][iy][ix][i]);
					ymm7 = _mm256_load_pd(&a[3][iz][iy][ix][i]);
					// set the result register 1 to zero 
					ymm9 = _mm256_setzero_pd();
					// set the result register 2 to zero
					ymm15 = _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the first and second column of b and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
					// multiply the row of a with the first and second column of b and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 

					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result
					_mm256_store_pd(&c[2][iz][iy][ix][i], ymm9);
					_mm256_store_pd(&c[3][iz][iy][ix][i], ymm15);
				}
			}
		}
	}

#endif 

}



// ------------------------------------------------------------------------------------------------------------- 

// compute all matmuls in x direction on a 3D data set 
#ifdef AVX512
extern "C" void __cdecl vsmm4x4cs3d_x(double *cptr, const double *aptr, double b[4][4], double s[1])
{
	__assume_aligned(aptr, 64);
	__assume_aligned(cptr, 64);
	__assume_aligned(b, 64);

	__asm { 
		mov rax, aptr                                       // pointer to a 
		mov rbx, b											// pointer to b 
		mov rcx, cptr										// pointer to c 
		mov r13, s											// pointer to s  
		// mov rdx, dptr									// pointer to d  ( only for debugging ) 
		 
		vbroadcastsd zmm14, QWORD PTR[r13]					// load s in vector register 

		vbroadcastsd zmm4, QWORD PTR[rbx + 0 * 8 ]			// load b[0][0] 
		vbroadcastsd zmm5, QWORD PTR[rbx + 1 * 8 ]			// load b[0][1] 
		vbroadcastsd zmm6, QWORD PTR[rbx + 2 * 8 ]			// load b[0][2] 
		vbroadcastsd zmm7, QWORD PTR[rbx + 3 * 8 ]			// load b[0][3] 

		vmulpd zmm0, zmm4, zmm14                            // multiply column of b with scalar s 
		vmulpd zmm1, zmm5, zmm14
		vmulpd zmm2, zmm6, zmm14
		vmulpd zmm3, zmm7, zmm14

		vbroadcastsd zmm4, QWORD PTR[rbx + 4 * 8]			// load b[1][0] 
		vbroadcastsd zmm5, QWORD PTR[rbx + 5 * 8]			// load b[1][1] 
		vbroadcastsd zmm6, QWORD PTR[rbx + 6 * 8]			// load b[1][2] 
		vbroadcastsd zmm7, QWORD PTR[rbx + 7 * 8]			// load b[1][3] 

		vmulpd zmm10, zmm4, zmm14							// multiply column of b with scalar s  
		vmulpd zmm11, zmm5, zmm14
		vmulpd zmm12, zmm6, zmm14
		vmulpd zmm13, zmm7, zmm14

		mov r10, 0
loopz:
		mov r15, 0  
loopvar: 

			// compute address of the pointers 
			mov r14,  0
			add r14,  r10        // z loop index 
			add r14,  r15        // variable loop index 

			mov r12, rax	     // compute a base pointer 
			add r12, r14        
			mov r13, rcx         // compute c base pointer 
			add r13, r14

			// first half 
			vmovapd		zmm4,  ZMMWORD PTR[r12 + 0 * ND1D*nVar*8 + 0 * nVar*8]		
			vmovapd		zmm5,  ZMMWORD PTR[r12 + 0 * ND1D*nVar*8 + 1 * nVar*8]		
			vmulpd			zmm9, zmm4, zmm0
			vmulpd			zmm15, zmm5, zmm11
			vmovapd		zmm6,  ZMMWORD PTR[r12 + 0 * ND1D*nVar*8 + 2 * nVar*8]		
			vmovapd		zmm7,  ZMMWORD PTR[r12 + 0 * ND1D*nVar*8 + 3 * nVar*8]		
			vxorpd			zmm16, zmm16, zmm16
			vfnmadd231pd	zmm16, zmm6, zmm11
			vxorpd			zmm17, zmm17, zmm17
			vfnmadd231pd	zmm17, zmm7, zmm0
			vmovapd		zmm8,  ZMMWORD PTR[r12 + 1 * ND1D*nVar*8 + 0 * nVar*8]     	
			vmovapd		zmm14, ZMMWORD PTR[r12 + 1 * ND1D*nVar*8 + 1 * nVar*8]	
			vmulpd			zmm28, zmm8, zmm0
			vmulpd			zmm29, zmm14, zmm11
			vmovapd		zmm18, ZMMWORD PTR[r12 + 1 * ND1D*nVar*8 + 2 * nVar*8]     	
			vmovapd		zmm19, ZMMWORD PTR[r12 + 1 * ND1D*nVar*8 + 3 * nVar*8]		

			vxorpd			zmm30, zmm30, zmm30											 
			vfnmadd231pd	zmm30, zmm18, zmm11											 
			vxorpd			zmm31, zmm31, zmm31											 
			vfnmadd231pd	zmm31, zmm19, zmm0											 

			prefetchw[r13 + 0 * ND1D*nVar * 8 + 0 * nVar * 8]
			prefetchw[r13 + 0 * ND1D*nVar * 8 + 1 * nVar * 8]
			prefetchw[r13 + 0 * ND1D*nVar * 8 + 2 * nVar * 8]
			prefetchw[r13 + 0 * ND1D*nVar * 8 + 3 * nVar * 8]
			prefetchw[r13 + 1 * ND1D*nVar * 8 + 0 * nVar * 8]
			prefetchw[r13 + 1 * ND1D*nVar * 8 + 1 * nVar * 8]
			prefetchw[r13 + 1 * ND1D*nVar * 8 + 2 * nVar * 8]
			prefetchw[r13 + 1 * ND1D*nVar * 8 + 3 * nVar * 8]

			vfmadd231pd  zmm9,  zmm5,  zmm1							 
			vfmadd231pd  zmm15, zmm6,  zmm12						 
			vfnmadd231pd zmm16, zmm7,  zmm10						 
			vfnmadd231pd zmm17, zmm4,  zmm3
			vfmadd231pd  zmm28, zmm14, zmm1							 
			vfmadd231pd  zmm29, zmm18, zmm12						  
			vfnmadd231pd zmm30, zmm19, zmm10
			vfnmadd231pd zmm31, zmm8,  zmm3

			vfmadd231pd  zmm9, zmm6,   zmm2							   		 
			vfmadd231pd  zmm15, zmm7,  zmm13							 
			vfnmadd231pd zmm16, zmm4,  zmm13
			vfnmadd231pd zmm17, zmm5,  zmm2
			vfmadd231pd  zmm28, zmm18, zmm2							 
			vfmadd231pd  zmm29, zmm19, zmm13							 	  
			vfnmadd231pd zmm30, zmm8,  zmm13
			vfnmadd231pd zmm31, zmm14, zmm2
																	 
			vfmadd231pd  zmm9,  zmm7, zmm3							 
			vfmadd231pd  zmm15, zmm4, zmm10							 
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND1D*nVar * 8 + 0 * nVar * 8], zmm9
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND1D*nVar * 8 + 1 * nVar * 8], zmm15
			vfnmadd231pd zmm16, zmm5, zmm12
			vfnmadd231pd zmm17, zmm6, zmm1							 								    
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND1D*nVar * 8 + 2 * nVar * 8], zmm16
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND1D*nVar * 8 + 3 * nVar * 8], zmm17
			vfmadd231pd  zmm28, zmm19, zmm3
			vfmadd231pd  zmm29, zmm8,  zmm10							 		  				  		    
			vmovapd		 ZMMWORD PTR[r13 + 1 * ND1D*nVar * 8 + 0 * nVar * 8], zmm28
			vmovapd		 ZMMWORD PTR[r13 + 1 * ND1D*nVar * 8 + 1 * nVar * 8], zmm29
			vfnmadd231pd zmm30, zmm14, zmm12
			vfnmadd231pd zmm31, zmm18, zmm1 		
			vmovapd		 ZMMWORD PTR[r13 + 1 * ND1D*nVar * 8 + 2 * nVar * 8], zmm30
			vmovapd		 ZMMWORD PTR[r13 + 1 * ND1D*nVar * 8 + 3 * nVar * 8], zmm31

			// second half 

			vmovapd		zmm4,  ZMMWORD PTR[r12 + 2 * ND1D*nVar*8 + 0 * nVar*8]     				 
			vmovapd		zmm5,  ZMMWORD PTR[r12 + 2 * ND1D*nVar*8 + 1 * nVar*8]					 
			vmulpd			zmm20, zmm4, zmm0
			vmulpd			zmm21, zmm5, zmm11
			vmovapd		zmm6,  ZMMWORD PTR[r12 + 2 * ND1D*nVar*8 + 2 * nVar*8]
			vmovapd		zmm7,  ZMMWORD PTR[r12 + 2 * ND1D*nVar*8 + 3 * nVar*8]					 
			vxorpd			zmm22, zmm22, zmm22
			vfnmadd231pd	zmm22, zmm6, zmm11
			vxorpd			zmm23, zmm23, zmm23
			vfnmadd231pd	zmm23, zmm7, zmm0

			vmovapd		zmm8,  ZMMWORD PTR[r12 + 3 * ND1D*nVar * 8 + 0 * nVar * 8]
			vmovapd		zmm14, ZMMWORD PTR[r12 + 3 * ND1D*nVar * 8 + 1 * nVar * 8]
			vmulpd			zmm24, zmm8, zmm0
			vmulpd			zmm25, zmm14, zmm11
			vmovapd		zmm18, ZMMWORD PTR[r12 + 3 * ND1D*nVar * 8 + 2 * nVar * 8]
			vmovapd		zmm19, ZMMWORD PTR[r12 + 3 * ND1D*nVar * 8 + 3 * nVar * 8]

			vxorpd			zmm26, zmm26, zmm26
			vfnmadd231pd	zmm26, zmm18, zmm11
			vxorpd			zmm27, zmm27, zmm27 
			vfnmadd231pd	zmm27, zmm19, zmm0

			prefetchw[r13 + 2 * ND1D*nVar * 8 + 0 * nVar * 8]
			prefetchw[r13 + 2 * ND1D*nVar * 8 + 1 * nVar * 8]
			prefetchw[r13 + 2 * ND1D*nVar * 8 + 2 * nVar * 8]
			prefetchw[r13 + 2 * ND1D*nVar * 8 + 3 * nVar * 8]
			prefetchw[r13 + 3 * ND1D*nVar * 8 + 0 * nVar * 8]
			prefetchw[r13 + 3 * ND1D*nVar * 8 + 1 * nVar * 8]
			prefetchw[r13 + 3 * ND1D*nVar * 8 + 2 * nVar * 8]
			prefetchw[r13 + 3 * ND1D*nVar * 8 + 3 * nVar * 8]

			vfmadd231pd  zmm20, zmm5,  zmm1							 
			vfmadd231pd  zmm21, zmm6,  zmm12						 
			vfnmadd231pd zmm22, zmm7,  zmm10						 
			vfnmadd231pd zmm23, zmm4,  zmm3
			vfmadd231pd  zmm24, zmm14, zmm1							 
			vfmadd231pd  zmm25, zmm18, zmm12						  
			vfnmadd231pd zmm26, zmm19, zmm10
			vfnmadd231pd zmm27, zmm8,  zmm3

			vfmadd231pd  zmm20, zmm6,   zmm2	 						   		 
			vfmadd231pd  zmm21, zmm7,  zmm13							 
			vfnmadd231pd zmm22, zmm4,  zmm13
			vfnmadd231pd zmm23, zmm5,  zmm2
			vfmadd231pd  zmm24, zmm18, zmm2							 
			vfmadd231pd  zmm25, zmm19, zmm13							 	  
			vfnmadd231pd zmm26, zmm8,  zmm13
			vfnmadd231pd zmm27, zmm14, zmm2

			vfmadd231pd zmm20, zmm7, zmm3							 		  		 
			vfmadd231pd zmm21, zmm4, zmm10							 		  
			vmovapd		ZMMWORD PTR[r13 + 2 * ND1D*nVar*8 + 0 * nVar*8], zmm20     			      
			vmovapd		ZMMWORD PTR[r13 + 2 * ND1D*nVar*8 + 1 * nVar*8], zmm21					 
			vfnmadd231pd zmm22, zmm5, zmm12							 							  
			vfnmadd231pd zmm23, zmm6, zmm1							 							  
			vmovapd		ZMMWORD PTR[r13 + 2 * ND1D*nVar*8 + 2 * nVar*8], zmm22    			      
			vmovapd		ZMMWORD PTR[r13 + 2 * ND1D*nVar*8 + 3 * nVar*8], zmm23					 
			vfmadd231pd zmm24, zmm19, zmm3							  							  
			vfmadd231pd zmm25, zmm8, zmm10							   				  			  
			vmovapd		ZMMWORD PTR[r13 + 3 * ND1D*nVar*8 + 0 * nVar*8], zmm24    			      
			vmovapd		ZMMWORD PTR[r13 + 3 * ND1D*nVar*8 + 1 * nVar*8], zmm25					 
			vfnmadd231pd zmm26, zmm14, zmm12							 							  
			vfnmadd231pd zmm27, zmm18, zmm1							 							  
			vmovapd		ZMMWORD PTR[r13 + 3 * ND1D*nVar*8 + 2 * nVar*8], zmm26    			      
			vmovapd		ZMMWORD PTR[r13 + 3 * ND1D*nVar*8 + 3 * nVar*8], zmm27					 

		add r15, VECTORLENGTH*8 
		cmp r15, nVar*8  
		jl loopvar

		add r10, ND2D*nVar*8   
		cmp r10, 4*ND2D*nVar*8 
		jl loopz

	}

	/* 

	__assume_aligned(c, 64);
	__assume_aligned(a, 64);
	__assume_aligned(b, 64);
	int i, j, k, iy, iz, it;
	__m512d zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7;
	__m512d zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15;
	__m512d zmm16, zmm17;
	__m512d zmm20, zmm21, zmm22, zmm23, zmm24, zmm25, zmm26, zmm27;
	__m128d xmm;


	// optimized AVX code for matmul on all tensor slices 
	// load the scalar into register 14 
	xmm = _mm_load1_pd(&s[0]);
	zmm14 = _mm512_broadcastsd_pd(xmm);
	// load first column of the small matrix b and scale with the scalar s 
	xmm = _mm_load1_pd(&b[0][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm0 = _mm512_mul_pd(zmm4, zmm14);
	zmm1 = _mm512_mul_pd(zmm5, zmm14);
	zmm2 = _mm512_mul_pd(zmm6, zmm14);
	zmm3 = _mm512_mul_pd(zmm7, zmm14);
	// load second column of the small matrix b and scale with the scalar s  
	xmm = _mm_load1_pd(&b[1][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm10 = _mm512_mul_pd(zmm4, zmm14);
	zmm11 = _mm512_mul_pd(zmm5, zmm14);
	zmm12 = _mm512_mul_pd(zmm6, zmm14);
	zmm13 = _mm512_mul_pd(zmm7, zmm14);
	// load third column of the small matrix b and scale with the scalar s 
	xmm = _mm_load1_pd(&b[2][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm20 = _mm512_mul_pd(zmm4, zmm14);
	zmm21 = _mm512_mul_pd(zmm5, zmm14);
	zmm22 = _mm512_mul_pd(zmm6, zmm14);
	zmm23 = _mm512_mul_pd(zmm7, zmm14);
	// load fourth column of the small matrix b and scale with the scalar s  
	xmm = _mm_load1_pd(&b[3][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm24 = _mm512_mul_pd(zmm4, zmm14);
	zmm25 = _mm512_mul_pd(zmm5, zmm14);
	zmm26 = _mm512_mul_pd(zmm6, zmm14);
	zmm27 = _mm512_mul_pd(zmm7, zmm14);

	for (iz = 0; iz < 4; iz++)
	{
		for (iy = 0; iy < 4; iy++)
		{
			for (i = 0; i < nVar; i = i + VECTORLENGTH)
			{
				// load VECTORLENGTH rows of matrix a on position i 
				// multiply the row of a with the first and second column of b and add the result to c using FMA 
				// multiply the row of a with the third and fourth column of b and add the result to c using FMA 
				// set the result registers 1-4 to zero 
				zmm4 = _mm512_load_pd(&a[iz][iy][0][i]);
				zmm9 = _mm512_setzero_pd();
				zmm9 = _mm512_fmadd_pd(zmm4, zmm0, zmm9);
				zmm5 = _mm512_load_pd(&a[iz][iy][1][i]);
				zmm15 = _mm512_setzero_pd();
				zmm15 = _mm512_fmadd_pd(zmm5, zmm11, zmm15);
				zmm6 = _mm512_load_pd(&a[iz][iy][2][i]);
				zmm16 = _mm512_setzero_pd();
				zmm16 = _mm512_fmadd_pd(zmm6, zmm22, zmm16);
				zmm7 = _mm512_load_pd(&a[iz][iy][3][i]);
				zmm17 = _mm512_setzero_pd();
				zmm17 = _mm512_fmadd_pd(zmm7, zmm27, zmm17);
				// 
				zmm9 = _mm512_fmadd_pd(zmm5, zmm1, zmm9);
				zmm15 = _mm512_fmadd_pd(zmm6, zmm12, zmm15);
				zmm16 = _mm512_fmadd_pd(zmm7, zmm23, zmm16);
				zmm17 = _mm512_fmadd_pd(zmm4, zmm24, zmm17);
				// 
				zmm9 = _mm512_fmadd_pd(zmm6, zmm2, zmm9);
				zmm15 = _mm512_fmadd_pd(zmm7, zmm13, zmm15);
				zmm16 = _mm512_fmadd_pd(zmm4, zmm20, zmm16);
				zmm17 = _mm512_fmadd_pd(zmm5, zmm25, zmm17);
				// store the result in c 
				zmm9 = _mm512_fmadd_pd(zmm7, zmm3, zmm9);
				_mm512_store_pd(&c[iz][iy][0][i], zmm9);
				zmm15 = _mm512_fmadd_pd(zmm4, zmm10, zmm15);
				_mm512_store_pd(&c[iz][iy][1][i], zmm15);
				zmm16 = _mm512_fmadd_pd(zmm5, zmm21, zmm16);
				_mm512_store_pd(&c[iz][iy][2][i], zmm16);
				zmm17 = _mm512_fmadd_pd(zmm6, zmm26, zmm17);
				_mm512_store_pd(&c[iz][iy][3][i], zmm17);
			}
		}
	}
	
	*/

#else 
extern "C" void __cdecl vsmm4x4cs3d_x(double c[4][4][4][nVar], double a[4][4][4][nVar], double b[4][4], double s[1])
{
	__assume_aligned(c, 32);
	__assume_aligned(a, 32);
	__assume_aligned(b, 32);
	int i, j, k, iy, iz, it;
	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	__m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;

	// optimized AVX code for matmul on all tensor slices 
	// load the scalar into register 14 
	ymm14 = _mm256_broadcast_sd(&s[0]);
	// load first column of the small matrix b and scale with the scalar s 
	ymm4 = _mm256_broadcast_sd(&b[0][0]);
	//ymm4 = _mm256_set_pd(b[0][0], b[0][1], b[0][2], b[0][3]);
	ymm5 = _mm256_broadcast_sd(&b[0][1]);
	ymm6 = _mm256_broadcast_sd(&b[0][2]);
	ymm7 = _mm256_broadcast_sd(&b[0][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load second column of the small matrix b and scale with the scalar s  
	ymm4 = _mm256_broadcast_sd(&b[1][0]);
	ymm5 = _mm256_broadcast_sd(&b[1][1]);
	ymm6 = _mm256_broadcast_sd(&b[1][2]);
	ymm7 = _mm256_broadcast_sd(&b[1][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);
		for (iz = 0; iz < 4; iz++)
		{
			for (iy = 0; iy < 4; iy++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of matrix a on position i 
					ymm4 = _mm256_load_pd(&a[iz][iy][0][i]);
					ymm5 = _mm256_load_pd(&a[iz][iy][1][i]);
					ymm6 = _mm256_load_pd(&a[iz][iy][2][i]);
					ymm7 = _mm256_load_pd(&a[iz][iy][3][i]);
					// set the result register 1 to zero 
					ymm9 = _mm256_setzero_pd();
					// set the result register 2 to zero
					ymm15 = _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the first and second column of b and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
					// multiply the row of a with the first and second column of b and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 
					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result in c 
					_mm256_store_pd(&c[iz][iy][0][i], ymm9);
					_mm256_store_pd(&c[iz][iy][1][i], ymm15);
				}
			}
		}

	// load third column of the small matrix b and scale with the scalar s 
	ymm4 = _mm256_broadcast_sd(&b[2][0]);
	ymm5 = _mm256_broadcast_sd(&b[2][1]);
	ymm6 = _mm256_broadcast_sd(&b[2][2]);
	ymm7 = _mm256_broadcast_sd(&b[2][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load fourth column of the small matrix b and scale with the scalar s  
	ymm4 = _mm256_broadcast_sd(&b[3][0]);
	ymm5 = _mm256_broadcast_sd(&b[3][1]);
	ymm6 = _mm256_broadcast_sd(&b[3][2]);
	ymm7 = _mm256_broadcast_sd(&b[3][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);

		for (iz = 0; iz < 4; iz++)
		{
			for (iy = 0; iy < 4; iy++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of a on position i 
					ymm4 = _mm256_load_pd(&a[iz][iy][0][i]);
					ymm5 = _mm256_load_pd(&a[iz][iy][1][i]);
					ymm6 = _mm256_load_pd(&a[iz][iy][2][i]);
					ymm7 = _mm256_load_pd(&a[iz][iy][3][i]);
					// set the result register 1 to zero 
					ymm9 = _mm256_setzero_pd();
					// set the result register 2 to zero
					ymm15 = _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the first and second column of b and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
					// multiply the row of a with the first and second column of b and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 

					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result
					_mm256_store_pd(&c[iz][iy][2][i], ymm9);
					_mm256_store_pd(&c[iz][iy][3][i], ymm15);
				}
			}
		}
#endif 	

}

// compute all matmuls in y direction on a 3D data set 
#ifdef AVX512
extern "C" void __cdecl vsmm4x4cs3d_y(double *cptr, const double *aptr, double b[4][4], double s[1])
{
	__assume_aligned(aptr, 64);
	__assume_aligned(cptr, 64);
	__assume_aligned(b, 64);

	__asm { 
		mov rax, aptr                                       // pointer to a 
		mov rbx, b											// pointer to b 
		mov rcx, cptr										// pointer to c 
		mov r13, s											// pointer to s  
		// mov rdx, dptr									// pointer to d  ( only for debugging ) 
		 
		vbroadcastsd zmm14, QWORD PTR[r13]					// load s in vector register 

		vbroadcastsd zmm4, QWORD PTR[rbx + 0 * 8 ]			// load b[0][0] 
		vbroadcastsd zmm5, QWORD PTR[rbx + 1 * 8 ]			// load b[0][1] 
		vbroadcastsd zmm6, QWORD PTR[rbx + 2 * 8 ]			// load b[0][2] 
		vbroadcastsd zmm7, QWORD PTR[rbx + 3 * 8 ]			// load b[0][3] 

		vmulpd zmm0, zmm4, zmm14                            // multiply column of b with scalar s 
		vmulpd zmm1, zmm5, zmm14
		vmulpd zmm2, zmm6, zmm14
		vmulpd zmm3, zmm7, zmm14

		vbroadcastsd zmm4, QWORD PTR[rbx + 4 * 8]			// load b[1][0] 
		vbroadcastsd zmm5, QWORD PTR[rbx + 5 * 8]			// load b[1][1] 
		vbroadcastsd zmm6, QWORD PTR[rbx + 6 * 8]			// load b[1][2] 
		vbroadcastsd zmm7, QWORD PTR[rbx + 7 * 8]			// load b[1][3] 

		vmulpd zmm10, zmm4, zmm14							// multiply column of b with scalar s  
		vmulpd zmm11, zmm5, zmm14
		vmulpd zmm12, zmm6, zmm14
		vmulpd zmm13, zmm7, zmm14

		mov r10, 0
loopz:
		mov r15, 0  
loopvar: 

			// compute address of the pointers 
			mov r14,  0
			add r14,  r10        // x loop index 
			add r14,  r15        // variable loop index 

			mov r12, rax	     // compute a base pointer 
			add r12, r14        
			mov r13, rcx         // compute c base pointer 
			add r13, r14

			// first half 
			vmovapd		zmm4,  ZMMWORD PTR[r12 + 0 * ND1D*nVar*8 + 0 * nVar*8]		
			vmovapd		zmm5,  ZMMWORD PTR[r12 + 1 * ND1D*nVar*8 + 0 * nVar*8]		
			vmulpd			zmm9, zmm4, zmm0
			vmulpd			zmm15, zmm5, zmm11
			vmovapd		zmm6,  ZMMWORD PTR[r12 + 2 * ND1D*nVar*8 + 0 * nVar*8]		
			vmovapd		zmm7,  ZMMWORD PTR[r12 + 3 * ND1D*nVar*8 + 0 * nVar*8]		
			vxorpd			zmm16, zmm16, zmm16
			vfnmadd231pd	zmm16, zmm6, zmm11
			vxorpd			zmm17, zmm17, zmm17
			vfnmadd231pd	zmm17, zmm7, zmm0
			vmovapd		zmm8,  ZMMWORD PTR[r12 + 0 * ND1D*nVar*8 + 1 * nVar*8]     	
			vmovapd		zmm14, ZMMWORD PTR[r12 + 1 * ND1D*nVar*8 + 1 * nVar*8]	
			vmulpd			zmm28, zmm8, zmm0
			vmulpd			zmm29, zmm14, zmm11
			vmovapd		zmm18, ZMMWORD PTR[r12 + 2 * ND1D*nVar*8 + 1 * nVar*8]     	
			vmovapd		zmm19, ZMMWORD PTR[r12 + 3 * ND1D*nVar*8 + 1 * nVar*8]		

			vxorpd			zmm30, zmm30, zmm30											 
			vfnmadd231pd	zmm30, zmm18, zmm11											 
			vxorpd			zmm31, zmm31, zmm31											 
			vfnmadd231pd	zmm31, zmm19, zmm0											 

			prefetchw[r13 + 0 * ND1D*nVar * 8 + 0 * nVar * 8]
			prefetchw[r13 + 1 * ND1D*nVar * 8 + 0 * nVar * 8]
			prefetchw[r13 + 2 * ND1D*nVar * 8 + 0 * nVar * 8]
			prefetchw[r13 + 3 * ND1D*nVar * 8 + 0 * nVar * 8]
			prefetchw[r13 + 0 * ND1D*nVar * 8 + 1 * nVar * 8]
			prefetchw[r13 + 1 * ND1D*nVar * 8 + 1 * nVar * 8]
			prefetchw[r13 + 2 * ND1D*nVar * 8 + 1 * nVar * 8]
			prefetchw[r13 + 3 * ND1D*nVar * 8 + 1 * nVar * 8]

			vfmadd231pd  zmm9,  zmm5,  zmm1							 
			vfmadd231pd  zmm15, zmm6,  zmm12						 
			vfnmadd231pd zmm16, zmm7,  zmm10						 
			vfnmadd231pd zmm17, zmm4,  zmm3
			vfmadd231pd  zmm28, zmm14, zmm1							 
			vfmadd231pd  zmm29, zmm18, zmm12						  
			vfnmadd231pd zmm30, zmm19, zmm10
			vfnmadd231pd zmm31, zmm8,  zmm3

			vfmadd231pd  zmm9, zmm6,   zmm2							   		 
			vfmadd231pd  zmm15, zmm7,  zmm13							 
			vfnmadd231pd zmm16, zmm4,  zmm13
			vfnmadd231pd zmm17, zmm5,  zmm2
			vfmadd231pd  zmm28, zmm18, zmm2							 
			vfmadd231pd  zmm29, zmm19, zmm13							 	  
			vfnmadd231pd zmm30, zmm8,  zmm13
			vfnmadd231pd zmm31, zmm14, zmm2
																	 
			vfmadd231pd  zmm9,  zmm7, zmm3							 
			vfmadd231pd  zmm15, zmm4, zmm10							 
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND1D*nVar * 8 + 0 * nVar * 8], zmm9
			vmovapd		 ZMMWORD PTR[r13 + 1 * ND1D*nVar * 8 + 0 * nVar * 8], zmm15
			vfnmadd231pd zmm16, zmm5, zmm12
			vfnmadd231pd zmm17, zmm6, zmm1							 								    
			vmovapd		 ZMMWORD PTR[r13 + 2 * ND1D*nVar * 8 + 0 * nVar * 8], zmm16
			vmovapd		 ZMMWORD PTR[r13 + 3 * ND1D*nVar * 8 + 0 * nVar * 8], zmm17
			vfmadd231pd  zmm28, zmm19, zmm3
			vfmadd231pd  zmm29, zmm8,  zmm10							 		  				  		    
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND1D*nVar * 8 + 1 * nVar * 8], zmm28
			vmovapd		 ZMMWORD PTR[r13 + 1 * ND1D*nVar * 8 + 1 * nVar * 8], zmm29
			vfnmadd231pd zmm30, zmm14, zmm12
			vfnmadd231pd zmm31, zmm18, zmm1 		
			vmovapd		 ZMMWORD PTR[r13 + 2 * ND1D*nVar * 8 + 1 * nVar * 8], zmm30
			vmovapd		 ZMMWORD PTR[r13 + 3 * ND1D*nVar * 8 + 1 * nVar * 8], zmm31

			// second half 

			vmovapd		zmm4,  ZMMWORD PTR[r12 + 0 * ND1D*nVar*8 + 2 * nVar*8]     				 
			vmovapd		zmm5,  ZMMWORD PTR[r12 + 1 * ND1D*nVar*8 + 2 * nVar*8]					 
			vmulpd			zmm20, zmm4, zmm0
			vmulpd			zmm21, zmm5, zmm11
			vmovapd		zmm6,  ZMMWORD PTR[r12 + 2 * ND1D*nVar*8 + 2 * nVar*8]
			vmovapd		zmm7,  ZMMWORD PTR[r12 + 3 * ND1D*nVar*8 + 2 * nVar*8]					 
			vxorpd			zmm22, zmm22, zmm22
			vfnmadd231pd	zmm22, zmm6, zmm11
			vxorpd			zmm23, zmm23, zmm23
			vfnmadd231pd	zmm23, zmm7, zmm0

			vmovapd		zmm8,  ZMMWORD PTR[r12 + 0 * ND1D*nVar * 8 + 3 * nVar * 8]
			vmovapd		zmm14, ZMMWORD PTR[r12 + 1 * ND1D*nVar * 8 + 3 * nVar * 8]
			vmulpd			zmm24, zmm8, zmm0
			vmulpd			zmm25, zmm14, zmm11
			vmovapd		zmm18, ZMMWORD PTR[r12 + 2 * ND1D*nVar * 8 + 3 * nVar * 8]
			vmovapd		zmm19, ZMMWORD PTR[r12 + 3 * ND1D*nVar * 8 + 3 * nVar * 8]

			vxorpd			zmm26, zmm26, zmm26
			vfnmadd231pd	zmm26, zmm18, zmm11
			vxorpd			zmm27, zmm27, zmm27 
			vfnmadd231pd	zmm27, zmm19, zmm0

			prefetchw[r13 + 0 * ND1D*nVar * 8 + 2 * nVar * 8]
			prefetchw[r13 + 1 * ND1D*nVar * 8 + 2 * nVar * 8]
			prefetchw[r13 + 2 * ND1D*nVar * 8 + 2 * nVar * 8]
			prefetchw[r13 + 3 * ND1D*nVar * 8 + 2 * nVar * 8]
			prefetchw[r13 + 0 * ND1D*nVar * 8 + 3 * nVar * 8]
			prefetchw[r13 + 1 * ND1D*nVar * 8 + 3 * nVar * 8]
			prefetchw[r13 + 2 * ND1D*nVar * 8 + 3 * nVar * 8]
			prefetchw[r13 + 3 * ND1D*nVar * 8 + 3 * nVar * 8]

			vfmadd231pd  zmm20, zmm5,  zmm1							 
			vfmadd231pd  zmm21, zmm6,  zmm12						 
			vfnmadd231pd zmm22, zmm7,  zmm10						 
			vfnmadd231pd zmm23, zmm4,  zmm3
			vfmadd231pd  zmm24, zmm14, zmm1							 
			vfmadd231pd  zmm25, zmm18, zmm12						  
			vfnmadd231pd zmm26, zmm19, zmm10
			vfnmadd231pd zmm27, zmm8,  zmm3

			vfmadd231pd  zmm20, zmm6,   zmm2	 						   		 
			vfmadd231pd  zmm21, zmm7,  zmm13							 
			vfnmadd231pd zmm22, zmm4,  zmm13
			vfnmadd231pd zmm23, zmm5,  zmm2
			vfmadd231pd  zmm24, zmm18, zmm2							 
			vfmadd231pd  zmm25, zmm19, zmm13							 	  
			vfnmadd231pd zmm26, zmm8,  zmm13
			vfnmadd231pd zmm27, zmm14, zmm2

			vfmadd231pd zmm20, zmm7, zmm3							 		  		 
			vfmadd231pd zmm21, zmm4, zmm10							 		  
			vmovapd		ZMMWORD PTR[r13 + 0 * ND1D*nVar*8 + 2 * nVar*8], zmm20     			      
			vmovapd		ZMMWORD PTR[r13 + 1 * ND1D*nVar*8 + 2 * nVar*8], zmm21					 
			vfnmadd231pd zmm22, zmm5, zmm12							 							  
			vfnmadd231pd zmm23, zmm6, zmm1							 							  
			vmovapd		ZMMWORD PTR[r13 + 2 * ND1D*nVar*8 + 2 * nVar*8], zmm22    			      
			vmovapd		ZMMWORD PTR[r13 + 3 * ND1D*nVar*8 + 2 * nVar*8], zmm23					 
			vfmadd231pd zmm24, zmm19, zmm3							  							  
			vfmadd231pd zmm25, zmm8, zmm10							   				  			  
			vmovapd		ZMMWORD PTR[r13 + 0 * ND1D*nVar*8 + 3 * nVar*8], zmm24    			      
			vmovapd		ZMMWORD PTR[r13 + 1 * ND1D*nVar*8 + 3 * nVar*8], zmm25					 
			vfnmadd231pd zmm26, zmm14, zmm12							 							  
			vfnmadd231pd zmm27, zmm18, zmm1							 							  
			vmovapd		ZMMWORD PTR[r13 + 2 * ND1D*nVar*8 + 3 * nVar*8], zmm26    			      
			vmovapd		ZMMWORD PTR[r13 + 3 * ND1D*nVar*8 + 3 * nVar*8], zmm27					 

		add r15, VECTORLENGTH*8 
		cmp r15, nVar*8  
		jl loopvar

		add r10, ND2D*nVar*8   
		cmp r10, 4*ND2D*nVar*8 
		jl loopz

	}

	/*

	__assume_aligned(c, 64);
	__assume_aligned(a, 64);
	__assume_aligned(b, 64);
	int i, j, k, ix, iz, it;
	__m512d zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7;
	__m512d zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15;
	__m512d zmm16, zmm17;
	__m512d zmm20, zmm21, zmm22, zmm23, zmm24, zmm25, zmm26, zmm27;
	__m128d xmm;


	// optimized AVX code for matmul on all tensor slices 
	// load the scalar into register 14 
	xmm = _mm_load1_pd(&s[0]);
	zmm14 = _mm512_broadcastsd_pd(xmm);
	// load first column of the small matrix b and scale with the scalar s 
	xmm = _mm_load1_pd(&b[0][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm0 = _mm512_mul_pd(zmm4, zmm14);
	zmm1 = _mm512_mul_pd(zmm5, zmm14);
	zmm2 = _mm512_mul_pd(zmm6, zmm14);
	zmm3 = _mm512_mul_pd(zmm7, zmm14);
	// load second column of the small matrix b and scale with the scalar s  
	xmm = _mm_load1_pd(&b[1][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm10 = _mm512_mul_pd(zmm4, zmm14);
	zmm11 = _mm512_mul_pd(zmm5, zmm14);
	zmm12 = _mm512_mul_pd(zmm6, zmm14);
	zmm13 = _mm512_mul_pd(zmm7, zmm14);
	// load third column of the small matrix b and scale with the scalar s 
	xmm = _mm_load1_pd(&b[2][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm20 = _mm512_mul_pd(zmm4, zmm14);
	zmm21 = _mm512_mul_pd(zmm5, zmm14);
	zmm22 = _mm512_mul_pd(zmm6, zmm14);
	zmm23 = _mm512_mul_pd(zmm7, zmm14);
	// load fourth column of the small matrix b and scale with the scalar s  
	xmm = _mm_load1_pd(&b[3][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm24 = _mm512_mul_pd(zmm4, zmm14);
	zmm25 = _mm512_mul_pd(zmm5, zmm14);
	zmm26 = _mm512_mul_pd(zmm6, zmm14);
	zmm27 = _mm512_mul_pd(zmm7, zmm14);

	for (iz = 0; iz < 4; iz++)
	{
		for (ix = 0; ix < 4; ix++)
		{
			for (i = 0; i < nVar; i = i + VECTORLENGTH)
			{
				// load VECTORLENGTH rows of matrix a on position i 
				// multiply the row of a with the first and second column of b and add the result to c using FMA 
				// multiply the row of a with the third and fourth column of b and add the result to c using FMA 
				// set the result registers 1-4 to zero 
				zmm4 = _mm512_load_pd(&a[iz][0][ix][i]);
				zmm9 = _mm512_setzero_pd();
				zmm9 = _mm512_fmadd_pd(zmm4, zmm0, zmm9);
				zmm5 = _mm512_load_pd(&a[iz][1][ix][i]);
				zmm15 = _mm512_setzero_pd();
				zmm15 = _mm512_fmadd_pd(zmm5, zmm11, zmm15);
				zmm6 = _mm512_load_pd(&a[iz][2][ix][i]);
				zmm16 = _mm512_setzero_pd();
				zmm16 = _mm512_fmadd_pd(zmm6, zmm22, zmm16);
				zmm7 = _mm512_load_pd(&a[iz][3][ix][i]);
				zmm17 = _mm512_setzero_pd();
				zmm17 = _mm512_fmadd_pd(zmm7, zmm27, zmm17);
				// 
				zmm9 = _mm512_fmadd_pd(zmm5, zmm1, zmm9);
				zmm15 = _mm512_fmadd_pd(zmm6, zmm12, zmm15);
				zmm16 = _mm512_fmadd_pd(zmm7, zmm23, zmm16);
				zmm17 = _mm512_fmadd_pd(zmm4, zmm24, zmm17);
				// 
				zmm9 = _mm512_fmadd_pd(zmm6, zmm2, zmm9);
				zmm15 = _mm512_fmadd_pd(zmm7, zmm13, zmm15);
				zmm16 = _mm512_fmadd_pd(zmm4, zmm20, zmm16);
				zmm17 = _mm512_fmadd_pd(zmm5, zmm25, zmm17);
				// store the result in c 
				zmm9 = _mm512_fmadd_pd(zmm7, zmm3, zmm9);
				_mm512_store_pd(&c[iz][0][ix][i], zmm9);
				zmm15 = _mm512_fmadd_pd(zmm4, zmm10, zmm15);
				_mm512_store_pd(&c[iz][1][ix][i], zmm15);
				zmm16 = _mm512_fmadd_pd(zmm5, zmm21, zmm16);
				_mm512_store_pd(&c[iz][2][ix][i], zmm16);
				zmm17 = _mm512_fmadd_pd(zmm6, zmm26, zmm17);
				_mm512_store_pd(&c[iz][3][ix][i], zmm17);
			}
		}
	}

	*/

#else 
extern "C" void __cdecl vsmm4x4cs3d_y(double c[4][4][4][nVar], double a[4][4][4][nVar], double b[4][4], double s[1])
{
	__assume_aligned(c, 32);
	__assume_aligned(a, 32);
	__assume_aligned(b, 32);
	int i, j, k, ix, iz, it;
	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	__m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;

	// optimized AVX code for matmul on all tensor slices 
	// load the scalar into register 14 
	ymm14 = _mm256_broadcast_sd(&s[0]);
	// load first column of the small matrix b and scale with scalar s 
	ymm4 = _mm256_broadcast_sd(&b[0][0]);
	ymm5 = _mm256_broadcast_sd(&b[0][1]);
	ymm6 = _mm256_broadcast_sd(&b[0][2]);
	ymm7 = _mm256_broadcast_sd(&b[0][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load second column of the small matrix b and scale with scalar s  
	ymm4 = _mm256_broadcast_sd(&b[1][0]);
	ymm5 = _mm256_broadcast_sd(&b[1][1]);
	ymm6 = _mm256_broadcast_sd(&b[1][2]);
	ymm7 = _mm256_broadcast_sd(&b[1][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);
		for (iz = 0; iz < 4; iz++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of matrix a on position i 
					ymm4 = _mm256_load_pd(&a[iz][0][ix][i]);
					ymm5 = _mm256_load_pd(&a[iz][1][ix][i]);
					ymm6 = _mm256_load_pd(&a[iz][2][ix][i]);
					ymm7 = _mm256_load_pd(&a[iz][3][ix][i]);
					// set the result register 1 to zero 
					ymm9 = _mm256_setzero_pd();
					// set the result register 2 to zero
					ymm15 = _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the first and second column of b and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
					// multiply the row of a with the first and second column of b and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 
					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result in c 
					_mm256_store_pd(&c[iz][0][ix][i], ymm9);
					_mm256_store_pd(&c[iz][1][ix][i], ymm15);
				}
			}
		}
	
	// load third column of the small matrix b 
	ymm4 = _mm256_broadcast_sd(&b[2][0]);
	ymm5 = _mm256_broadcast_sd(&b[2][1]);
	ymm6 = _mm256_broadcast_sd(&b[2][2]);
	ymm7 = _mm256_broadcast_sd(&b[2][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load fourth column of the small matrix b 
	ymm4 = _mm256_broadcast_sd(&b[3][0]);
	ymm5 = _mm256_broadcast_sd(&b[3][1]);
	ymm6 = _mm256_broadcast_sd(&b[3][2]);
	ymm7 = _mm256_broadcast_sd(&b[3][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);

		for (iz = 0; iz < 4; iz++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of a on position i 
					ymm4 = _mm256_load_pd(&a[iz][0][ix][i]);
					ymm5 = _mm256_load_pd(&a[iz][1][ix][i]);
					ymm6 = _mm256_load_pd(&a[iz][2][ix][i]);
					ymm7 = _mm256_load_pd(&a[iz][3][ix][i]);
					// set the result register 1 to zero 
					ymm9 = _mm256_setzero_pd();
					// set the result register 2 to zero
					ymm15 = _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the first and second column of b and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
					// multiply the row of a with the first and second column of b and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 

					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result
					_mm256_store_pd(&c[iz][2][ix][i], ymm9);
					_mm256_store_pd(&c[iz][3][ix][i], ymm15);
				}
			}
		}
#endif 	
}

#ifdef AVX512
extern "C" void __cdecl vsmm4x4cs3d_z(double *cptr, const double *aptr, double b[4][4], double s[1])
{
	__assume_aligned(aptr, 64);
	__assume_aligned(cptr, 64);
	__assume_aligned(b, 64);

	__asm { 
		mov rax, aptr                                       // pointer to a 
		mov rbx, b											// pointer to b 
		mov rcx, cptr										// pointer to c 
		mov r13, s											// pointer to s  
		// mov rdx, dptr									// pointer to d  ( only for debugging ) 
		 
		vbroadcastsd zmm14, QWORD PTR[r13]					// load s in vector register 

		vbroadcastsd zmm4, QWORD PTR[rbx + 0 * 8 ]			// load b[0][0] 
		vbroadcastsd zmm5, QWORD PTR[rbx + 1 * 8 ]			// load b[0][1] 
		vbroadcastsd zmm6, QWORD PTR[rbx + 2 * 8 ]			// load b[0][2] 
		vbroadcastsd zmm7, QWORD PTR[rbx + 3 * 8 ]			// load b[0][3] 

		vmulpd zmm0, zmm4, zmm14                            // multiply column of b with scalar s 
		vmulpd zmm1, zmm5, zmm14
		vmulpd zmm2, zmm6, zmm14
		vmulpd zmm3, zmm7, zmm14

		vbroadcastsd zmm4, QWORD PTR[rbx + 4 * 8]			// load b[1][0] 
		vbroadcastsd zmm5, QWORD PTR[rbx + 5 * 8]			// load b[1][1] 
		vbroadcastsd zmm6, QWORD PTR[rbx + 6 * 8]			// load b[1][2] 
		vbroadcastsd zmm7, QWORD PTR[rbx + 7 * 8]			// load b[1][3] 

		vmulpd zmm10, zmm4, zmm14							// multiply column of b with scalar s  
		vmulpd zmm11, zmm5, zmm14
		vmulpd zmm12, zmm6, zmm14
		vmulpd zmm13, zmm7, zmm14

		mov r10, 0
loopy:
		mov r15, 0  
loopvar: 

			// compute address of the pointers 
			mov r14,  0
			add r14,  r10        // x loop index 
			add r14,  r15        // variable loop index 

			mov r12, rax	     // compute a base pointer 
			add r12, r14        
			mov r13, rcx         // compute c base pointer 
			add r13, r14

			// first half 
			vmovapd		zmm4,  ZMMWORD PTR[r12 + 0 * ND2D*nVar*8 + 0 * nVar*8]		
			vmovapd		zmm5,  ZMMWORD PTR[r12 + 1 * ND2D*nVar*8 + 0 * nVar*8]		
			vmulpd			zmm9, zmm4, zmm0
			vmulpd			zmm15, zmm5, zmm11
			vmovapd		zmm6,  ZMMWORD PTR[r12 + 2 * ND2D*nVar*8 + 0 * nVar*8]		
			vmovapd		zmm7,  ZMMWORD PTR[r12 + 3 * ND2D*nVar*8 + 0 * nVar*8]		
			vxorpd			zmm16, zmm16, zmm16
			vfnmadd231pd	zmm16, zmm6, zmm11
			vxorpd			zmm17, zmm17, zmm17
			vfnmadd231pd	zmm17, zmm7, zmm0
			vmovapd		zmm8,  ZMMWORD PTR[r12 + 0 * ND2D*nVar*8 + 1 * nVar*8]     	
			vmovapd		zmm14, ZMMWORD PTR[r12 + 1 * ND2D*nVar*8 + 1 * nVar*8]	
			vmulpd			zmm28, zmm8, zmm0
			vmulpd			zmm29, zmm14, zmm11
			vmovapd		zmm18, ZMMWORD PTR[r12 + 2 * ND2D*nVar*8 + 1 * nVar*8]     	
			vmovapd		zmm19, ZMMWORD PTR[r12 + 3 * ND2D*nVar*8 + 1 * nVar*8]		

			vxorpd			zmm30, zmm30, zmm30											 
			vfnmadd231pd	zmm30, zmm18, zmm11											 
			vxorpd			zmm31, zmm31, zmm31											 
			vfnmadd231pd	zmm31, zmm19, zmm0											 

			prefetchw[r13 + 0 * ND2D*nVar * 8 + 0 * nVar * 8]
			prefetchw[r13 + 1 * ND2D*nVar * 8 + 0 * nVar * 8]
			prefetchw[r13 + 2 * ND2D*nVar * 8 + 0 * nVar * 8]
			prefetchw[r13 + 3 * ND2D*nVar * 8 + 0 * nVar * 8]
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 1 * nVar * 8]
			prefetchw[r13 + 1 * ND2D*nVar * 8 + 1 * nVar * 8]
			prefetchw[r13 + 2 * ND2D*nVar * 8 + 1 * nVar * 8]
			prefetchw[r13 + 3 * ND2D*nVar * 8 + 1 * nVar * 8]

			vfmadd231pd  zmm9,  zmm5,  zmm1							 
			vfmadd231pd  zmm15, zmm6,  zmm12						 
			vfnmadd231pd zmm16, zmm7,  zmm10						 
			vfnmadd231pd zmm17, zmm4,  zmm3
			vfmadd231pd  zmm28, zmm14, zmm1							 
			vfmadd231pd  zmm29, zmm18, zmm12						  
			vfnmadd231pd zmm30, zmm19, zmm10
			vfnmadd231pd zmm31, zmm8,  zmm3

			vfmadd231pd  zmm9, zmm6,   zmm2							   		 
			vfmadd231pd  zmm15, zmm7,  zmm13							 
			vfnmadd231pd zmm16, zmm4,  zmm13
			vfnmadd231pd zmm17, zmm5,  zmm2
			vfmadd231pd  zmm28, zmm18, zmm2							 
			vfmadd231pd  zmm29, zmm19, zmm13							 	  
			vfnmadd231pd zmm30, zmm8,  zmm13
			vfnmadd231pd zmm31, zmm14, zmm2
																	 
			vfmadd231pd  zmm9,  zmm7, zmm3							 
			vfmadd231pd  zmm15, zmm4, zmm10							 
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 0 * nVar * 8], zmm9
			vmovapd		 ZMMWORD PTR[r13 + 1 * ND2D*nVar * 8 + 0 * nVar * 8], zmm15
			vfnmadd231pd zmm16, zmm5, zmm12
			vfnmadd231pd zmm17, zmm6, zmm1							 								    
			vmovapd		 ZMMWORD PTR[r13 + 2 * ND2D*nVar * 8 + 0 * nVar * 8], zmm16
			vmovapd		 ZMMWORD PTR[r13 + 3 * ND2D*nVar * 8 + 0 * nVar * 8], zmm17
			vfmadd231pd  zmm28, zmm19, zmm3
			vfmadd231pd  zmm29, zmm8,  zmm10							 		  				  		    
			vmovapd		 ZMMWORD PTR[r13 + 0 * ND2D*nVar * 8 + 1 * nVar * 8], zmm28
			vmovapd		 ZMMWORD PTR[r13 + 1 * ND2D*nVar * 8 + 1 * nVar * 8], zmm29
			vfnmadd231pd zmm30, zmm14, zmm12
			vfnmadd231pd zmm31, zmm18, zmm1 		
			vmovapd		 ZMMWORD PTR[r13 + 2 * ND2D*nVar * 8 + 1 * nVar * 8], zmm30
			vmovapd		 ZMMWORD PTR[r13 + 3 * ND2D*nVar * 8 + 1 * nVar * 8], zmm31

			// second half 

			vmovapd		zmm4,  ZMMWORD PTR[r12 + 0 * ND2D*nVar*8 + 2 * nVar*8]     				 
			vmovapd		zmm5,  ZMMWORD PTR[r12 + 1 * ND2D*nVar*8 + 2 * nVar*8]					 
			vmulpd			zmm20, zmm4, zmm0
			vmulpd			zmm21, zmm5, zmm11
			vmovapd		zmm6,  ZMMWORD PTR[r12 + 2 * ND2D*nVar*8 + 2 * nVar*8]
			vmovapd		zmm7,  ZMMWORD PTR[r12 + 3 * ND2D*nVar*8 + 2 * nVar*8]					 
			vxorpd			zmm22, zmm22, zmm22
			vfnmadd231pd	zmm22, zmm6, zmm11
			vxorpd			zmm23, zmm23, zmm23
			vfnmadd231pd	zmm23, zmm7, zmm0

			vmovapd		zmm8,  ZMMWORD PTR[r12 + 0 * ND2D*nVar * 8 + 3 * nVar * 8]
			vmovapd		zmm14, ZMMWORD PTR[r12 + 1 * ND2D*nVar * 8 + 3 * nVar * 8]
			vmulpd			zmm24, zmm8, zmm0
			vmulpd			zmm25, zmm14, zmm11
			vmovapd		zmm18, ZMMWORD PTR[r12 + 2 * ND2D*nVar * 8 + 3 * nVar * 8]
			vmovapd		zmm19, ZMMWORD PTR[r12 + 3 * ND2D*nVar * 8 + 3 * nVar * 8]

			vxorpd			zmm26, zmm26, zmm26
			vfnmadd231pd	zmm26, zmm18, zmm11
			vxorpd			zmm27, zmm27, zmm27 
			vfnmadd231pd	zmm27, zmm19, zmm0

			prefetchw[r13 + 0 * ND2D*nVar * 8 + 2 * nVar * 8]
			prefetchw[r13 + 1 * ND2D*nVar * 8 + 2 * nVar * 8]
			prefetchw[r13 + 2 * ND2D*nVar * 8 + 2 * nVar * 8]
			prefetchw[r13 + 3 * ND2D*nVar * 8 + 2 * nVar * 8]
			prefetchw[r13 + 0 * ND2D*nVar * 8 + 3 * nVar * 8]
			prefetchw[r13 + 1 * ND2D*nVar * 8 + 3 * nVar * 8]
			prefetchw[r13 + 2 * ND2D*nVar * 8 + 3 * nVar * 8]
			prefetchw[r13 + 3 * ND2D*nVar * 8 + 3 * nVar * 8]

			vfmadd231pd  zmm20, zmm5,  zmm1							 
			vfmadd231pd  zmm21, zmm6,  zmm12						 
			vfnmadd231pd zmm22, zmm7,  zmm10						 
			vfnmadd231pd zmm23, zmm4,  zmm3
			vfmadd231pd  zmm24, zmm14, zmm1							 
			vfmadd231pd  zmm25, zmm18, zmm12						  
			vfnmadd231pd zmm26, zmm19, zmm10
			vfnmadd231pd zmm27, zmm8,  zmm3

			vfmadd231pd  zmm20, zmm6,   zmm2	 						   		 
			vfmadd231pd  zmm21, zmm7,  zmm13							 
			vfnmadd231pd zmm22, zmm4,  zmm13
			vfnmadd231pd zmm23, zmm5,  zmm2
			vfmadd231pd  zmm24, zmm18, zmm2							 
			vfmadd231pd  zmm25, zmm19, zmm13							 	  
			vfnmadd231pd zmm26, zmm8,  zmm13
			vfnmadd231pd zmm27, zmm14, zmm2

			vfmadd231pd zmm20, zmm7, zmm3							 		  		 
			vfmadd231pd zmm21, zmm4, zmm10							 		  
			vmovapd		ZMMWORD PTR[r13 + 0 * ND2D*nVar*8 + 2 * nVar*8], zmm20     			      
			vmovapd		ZMMWORD PTR[r13 + 1 * ND2D*nVar*8 + 2 * nVar*8], zmm21					 
			vfnmadd231pd zmm22, zmm5, zmm12							 							  
			vfnmadd231pd zmm23, zmm6, zmm1							 							  
			vmovapd		ZMMWORD PTR[r13 + 2 * ND2D*nVar*8 + 2 * nVar*8], zmm22    			      
			vmovapd		ZMMWORD PTR[r13 + 3 * ND2D*nVar*8 + 2 * nVar*8], zmm23					 
			vfmadd231pd zmm24, zmm19, zmm3							  							  
			vfmadd231pd zmm25, zmm8, zmm10							   				  			  
			vmovapd		ZMMWORD PTR[r13 + 0 * ND2D*nVar*8 + 3 * nVar*8], zmm24    			      
			vmovapd		ZMMWORD PTR[r13 + 1 * ND2D*nVar*8 + 3 * nVar*8], zmm25					 
			vfnmadd231pd zmm26, zmm14, zmm12							 							  
			vfnmadd231pd zmm27, zmm18, zmm1							 							  
			vmovapd		ZMMWORD PTR[r13 + 2 * ND2D*nVar*8 + 3 * nVar*8], zmm26    			      
			vmovapd		ZMMWORD PTR[r13 + 3 * ND2D*nVar*8 + 3 * nVar*8], zmm27					 

		add r15, VECTORLENGTH*8 
		cmp r15, nVar*8  
		jl loopvar

		add r10, ND1D*nVar*8   
		cmp r10, 4*ND1D*nVar*8 
		jl loopy

	}


	/*
	__assume_aligned(c, 64);
	__assume_aligned(a, 64);
	__assume_aligned(b, 64);
	int i, j, k, ix, iy, it;
	__m512d zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7;
	__m512d zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15;
	__m512d zmm16, zmm17;
	__m512d zmm20, zmm21, zmm22, zmm23, zmm24, zmm25, zmm26, zmm27;
	__m128d xmm;


	// optimized AVX code for matmul on all tensor slices 
	// load the scalar into register 14 
	xmm = _mm_load1_pd(&s[0]);
	zmm14 = _mm512_broadcastsd_pd(xmm);
	// load first column of the small matrix b and scale with the scalar s 
	xmm = _mm_load1_pd(&b[0][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm0 = _mm512_mul_pd(zmm4, zmm14);
	zmm1 = _mm512_mul_pd(zmm5, zmm14);
	zmm2 = _mm512_mul_pd(zmm6, zmm14);
	zmm3 = _mm512_mul_pd(zmm7, zmm14);
	// load second column of the small matrix b and scale with the scalar s  
	xmm = _mm_load1_pd(&b[1][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm10 = _mm512_mul_pd(zmm4, zmm14);
	zmm11 = _mm512_mul_pd(zmm5, zmm14);
	zmm12 = _mm512_mul_pd(zmm6, zmm14);
	zmm13 = _mm512_mul_pd(zmm7, zmm14);
	// load third column of the small matrix b and scale with the scalar s 
	xmm = _mm_load1_pd(&b[2][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm20 = _mm512_mul_pd(zmm4, zmm14);
	zmm21 = _mm512_mul_pd(zmm5, zmm14);
	zmm22 = _mm512_mul_pd(zmm6, zmm14);
	zmm23 = _mm512_mul_pd(zmm7, zmm14);
	// load fourth column of the small matrix b and scale with the scalar s  
	xmm = _mm_load1_pd(&b[3][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm24 = _mm512_mul_pd(zmm4, zmm14);
	zmm25 = _mm512_mul_pd(zmm5, zmm14);
	zmm26 = _mm512_mul_pd(zmm6, zmm14);
	zmm27 = _mm512_mul_pd(zmm7, zmm14);

	for (iy = 0; iy < 4; iy++)
	{
		for (ix = 0; ix < 4; ix++)
		{
			for (i = 0; i < nVar; i = i + VECTORLENGTH)
			{
				// load VECTORLENGTH rows of matrix a on position i 
				// multiply the row of a with the first and second column of b and add the result to c using FMA 
				// multiply the row of a with the third and fourth column of b and add the result to c using FMA 
				// set the result registers 1-4 to zero 
				zmm4 = _mm512_load_pd(&a[0][iy][ix][i]);
				zmm9 = _mm512_setzero_pd();
				zmm9 = _mm512_fmadd_pd(zmm4, zmm0, zmm9);
				zmm5 = _mm512_load_pd(&a[1][iy][ix][i]);
				zmm15 = _mm512_setzero_pd();
				zmm15 = _mm512_fmadd_pd(zmm5, zmm11, zmm15);
				zmm6 = _mm512_load_pd(&a[2][iy][ix][i]);
				zmm16 = _mm512_setzero_pd();
				zmm16 = _mm512_fmadd_pd(zmm6, zmm22, zmm16);
				zmm7 = _mm512_load_pd(&a[3][iy][ix][i]);
				zmm17 = _mm512_setzero_pd();
				zmm17 = _mm512_fmadd_pd(zmm7, zmm27, zmm17);
				// 
				zmm9 = _mm512_fmadd_pd(zmm5, zmm1, zmm9);
				zmm15 = _mm512_fmadd_pd(zmm6, zmm12, zmm15);
				zmm16 = _mm512_fmadd_pd(zmm7, zmm23, zmm16);
				zmm17 = _mm512_fmadd_pd(zmm4, zmm24, zmm17);
				// 
				zmm9 = _mm512_fmadd_pd(zmm6, zmm2, zmm9);
				zmm15 = _mm512_fmadd_pd(zmm7, zmm13, zmm15);
				zmm16 = _mm512_fmadd_pd(zmm4, zmm20, zmm16);
				zmm17 = _mm512_fmadd_pd(zmm5, zmm25, zmm17);
				// store the result in c 
				zmm9 = _mm512_fmadd_pd(zmm7, zmm3, zmm9);
				_mm512_store_pd(&c[0][iy][ix][i], zmm9);
				zmm15 = _mm512_fmadd_pd(zmm4, zmm10, zmm15);
				_mm512_store_pd(&c[1][iy][ix][i], zmm15);
				zmm16 = _mm512_fmadd_pd(zmm5, zmm21, zmm16);
				_mm512_store_pd(&c[2][iy][ix][i], zmm16);
				zmm17 = _mm512_fmadd_pd(zmm6, zmm26, zmm17);
				_mm512_store_pd(&c[3][iy][ix][i], zmm17);
			}
		}
	}
	*/

#else 
// compute all matmuls in z direction on a 3D data set 
extern "C" void __cdecl vsmm4x4cs3d_z(double c[4][4][4][nVar], double a[4][4][4][nVar], double b[4][4], double s[1])
{
	__assume_aligned(c, 32);
	__assume_aligned(a, 32);
	__assume_aligned(b, 32);
	int i, j, k, ix, iy, it;
	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	__m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;

	// optimized AVX code for matmul on all tensor slices 
	// load the scalar into register 14 
	ymm14 = _mm256_broadcast_sd(&s[0]);
	// load first column of the small matrix b and scale with scalar s 
	ymm4 = _mm256_broadcast_sd(&b[0][0]);
	ymm5 = _mm256_broadcast_sd(&b[0][1]);
	ymm6 = _mm256_broadcast_sd(&b[0][2]);
	ymm7 = _mm256_broadcast_sd(&b[0][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load second column of the small matrix b and scale with scalar s  
	ymm4 = _mm256_broadcast_sd(&b[1][0]);
	ymm5 = _mm256_broadcast_sd(&b[1][1]);
	ymm6 = _mm256_broadcast_sd(&b[1][2]);
	ymm7 = _mm256_broadcast_sd(&b[1][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);

		for (iy = 0; iy < 4; iy++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of matrix a on position i 
					ymm4 = _mm256_load_pd(&a[0][iy][ix][i]);
					ymm5 = _mm256_load_pd(&a[1][iy][ix][i]);
					ymm6 = _mm256_load_pd(&a[2][iy][ix][i]);
					ymm7 = _mm256_load_pd(&a[3][iy][ix][i]);
					// set the result register 1 to zero 
					ymm9 = _mm256_setzero_pd();
					// set the result register 2 to zero
					ymm15 = _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the first and second column of b and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
					// multiply the row of a with the first and second column of b and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 
					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result in c 
					_mm256_store_pd(&c[0][iy][ix][i], ymm9);
					_mm256_store_pd(&c[1][iy][ix][i], ymm15);
				}
			}
		}
	
	// load third column of the small matrix b 
	ymm4 = _mm256_broadcast_sd(&b[2][0]);
	ymm5 = _mm256_broadcast_sd(&b[2][1]);
	ymm6 = _mm256_broadcast_sd(&b[2][2]);
	ymm7 = _mm256_broadcast_sd(&b[2][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load fourth column of the small matrix b 
	ymm4 = _mm256_broadcast_sd(&b[3][0]);
	ymm5 = _mm256_broadcast_sd(&b[3][1]);
	ymm6 = _mm256_broadcast_sd(&b[3][2]);
	ymm7 = _mm256_broadcast_sd(&b[3][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);

		for (iy = 0; iy < 4; iy++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of a on position i 
					ymm4 = _mm256_load_pd(&a[0][iy][ix][i]);
					ymm5 = _mm256_load_pd(&a[1][iy][ix][i]);
					ymm6 = _mm256_load_pd(&a[2][iy][ix][i]);
					ymm7 = _mm256_load_pd(&a[3][iy][ix][i]);
					// set the result register 1 to zero 
					ymm9 = _mm256_setzero_pd();
					// set the result register 2 to zero
					ymm15 = _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the first and second column of b and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
					// multiply the row of a with the first and second column of b and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 

					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result
					_mm256_store_pd(&c[2][iy][ix][i], ymm9);
					_mm256_store_pd(&c[3][iy][ix][i], ymm15);
				}
			}
		}
#endif 

}


// compute all matmuls in x direction on a 3D data set with vector weights (very useful for stiffness matrix multiplications) 
extern "C" void __cdecl vsmm4x4c3d_kx(double c[4][4][4][nVar], double a[4][4][4][nVar], double b[4][4], double s[1])
{
#ifdef AVX512

	int i, j, k, it, iy, iz;
	double weight[1];
	__m512d zmm0, zmm1, zmm2, zmm3, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15;
	__m512d zmm20, zmm21, zmm22, zmm23, zmm24, zmm25, zmm26, zmm27;
	__m512d zmm4, zmm5, zmm6, zmm7, zmm8, zmm9, zmm16, zmm17, zmm18, zmm19, zmm28, zmm29, zmm30, zmm31;
	__m128d xmm;

	// load the scalar into register 14 
	xmm = _mm_load1_pd(&s[0]);
	zmm14 = _mm512_broadcastsd_pd(xmm);
	// load first column of the small matrix b and scale with the scalar s 
	xmm = _mm_load1_pd(&b[0][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm0 = _mm512_mul_pd(zmm4, zmm14);
	zmm1 = _mm512_mul_pd(zmm5, zmm14);
	zmm2 = _mm512_mul_pd(zmm6, zmm14);
	zmm3 = _mm512_mul_pd(zmm7, zmm14);
	// load second column of the small matrix b and scale with the scalar s  
	xmm = _mm_load1_pd(&b[1][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm10 = _mm512_mul_pd(zmm4, zmm14);
	zmm11 = _mm512_mul_pd(zmm5, zmm14);
	zmm12 = _mm512_mul_pd(zmm6, zmm14);
	zmm13 = _mm512_mul_pd(zmm7, zmm14);
	// load third column of the small matrix b and scale with the scalar s 
	xmm = _mm_load1_pd(&b[2][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm20 = _mm512_mul_pd(zmm4, zmm14);
	zmm21 = _mm512_mul_pd(zmm5, zmm14);
	zmm22 = _mm512_mul_pd(zmm6, zmm14);
	zmm23 = _mm512_mul_pd(zmm7, zmm14);
	// load fourth column of the small matrix b and scale with the scalar s  
	xmm = _mm_load1_pd(&b[3][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm24 = _mm512_mul_pd(zmm4, zmm14);
	zmm25 = _mm512_mul_pd(zmm5, zmm14);
	zmm26 = _mm512_mul_pd(zmm6, zmm14);
	zmm27 = _mm512_mul_pd(zmm7, zmm14);

	for (iz = 0; iz < 4; iz++)
	{
		for (iy = 0; iy < 4; iy++)
		{
			for (i = 0; i < nVar; i = i + VECTORLENGTH)
			{
				// load VECTORLENGTH rows of matrix a on position i 
				// load c into the result registers 
				// multiply the row of a with the first and second column of b and add the result to c using FMA 
				// multiply the row of a with the third and fourth column of b and add the result to c using FMA 
				zmm4 = _mm512_load_pd(&a[iz][iy][0][i]);
				zmm9 = _mm512_load_pd(&c[iz][iy][0][i]);
				zmm9 = _mm512_fmadd_pd(zmm4, zmm0, zmm9);
				zmm5 = _mm512_load_pd(&a[iz][iy][1][i]);
				zmm15 = _mm512_load_pd(&c[iz][iy][1][i]);
				zmm15 = _mm512_fmadd_pd(zmm5, zmm11, zmm15);
				zmm6 = _mm512_load_pd(&a[iz][iy][2][i]);
				zmm16 = _mm512_load_pd(&c[iz][iy][2][i]);
				zmm16 = _mm512_fmadd_pd(zmm6, zmm22, zmm16);
				zmm7 = _mm512_load_pd(&a[iz][iy][3][i]);
				zmm17 = _mm512_load_pd(&c[iz][iy][3][i]);
				zmm17 = _mm512_fmadd_pd(zmm7, zmm27, zmm17);
				// 
				zmm9 = _mm512_fmadd_pd(zmm5, zmm1, zmm9);
				zmm15 = _mm512_fmadd_pd(zmm6, zmm12, zmm15);
				zmm16 = _mm512_fmadd_pd(zmm7, zmm23, zmm16);
				zmm17 = _mm512_fmadd_pd(zmm4, zmm24, zmm17);
				// 
				zmm9 = _mm512_fmadd_pd(zmm6, zmm2, zmm9);
				zmm15 = _mm512_fmadd_pd(zmm7, zmm13, zmm15);
				zmm16 = _mm512_fmadd_pd(zmm4, zmm20, zmm16);
				zmm17 = _mm512_fmadd_pd(zmm5, zmm25, zmm17);
				// store the result in c 
				zmm9 = _mm512_fmadd_pd(zmm7, zmm3, zmm9);
				_mm512_store_pd(&c[iz][iy][0][i], zmm9);
				zmm15 = _mm512_fmadd_pd(zmm4, zmm10, zmm15);
				_mm512_store_pd(&c[iz][iy][1][i], zmm15);
				zmm16 = _mm512_fmadd_pd(zmm5, zmm21, zmm16);
				_mm512_store_pd(&c[iz][iy][2][i], zmm16);
				zmm17 = _mm512_fmadd_pd(zmm6, zmm26, zmm17);
				_mm512_store_pd(&c[iz][iy][3][i], zmm17);
			}
		}
	}
	



#else 

	__assume_aligned(c, 32);
	__assume_aligned(a, 32);
	__assume_aligned(b, 32);
	int i, j, k, iy, iz, it;
	double weight[1];
	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	__m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;

	// load the scalar into register 14 
	ymm14 = _mm256_broadcast_sd(&s[0]);
	// load first column of the small matrix b and scale with scalar s 
	ymm4 = _mm256_broadcast_sd(&b[0][0]);
	ymm5 = _mm256_broadcast_sd(&b[0][1]);
	ymm6 = _mm256_broadcast_sd(&b[0][2]);
	ymm7 = _mm256_broadcast_sd(&b[0][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load second column of the small matrix b and scale with scalar s  
	ymm4 = _mm256_broadcast_sd(&b[1][0]);
	ymm5 = _mm256_broadcast_sd(&b[1][1]);
	ymm6 = _mm256_broadcast_sd(&b[1][2]);
	ymm7 = _mm256_broadcast_sd(&b[1][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);

	for (iz = 0; iz < 4; iz++)
	{
		for (iy = 0; iy < 4; iy++)
		{
			// optimized AVX code for matmul on all tensor slices 
			for (i = 0; i < nVar; i = i + VECTORLENGTH)
			{
				// load VECTORLENGTH rows of matrix a on position i 
				ymm4 = _mm256_load_pd(&a[iz][iy][0][i]);
				ymm5 = _mm256_load_pd(&a[iz][iy][1][i]);
				ymm6 = _mm256_load_pd(&a[iz][iy][2][i]);
				ymm7 = _mm256_load_pd(&a[iz][iy][3][i]);
				// set the result register 1 to zero 
				ymm9 = _mm256_load_pd(&c[iz][iy][0][i]); //  _mm256_setzero_pd();
															 // set the result register 2 to zero
				ymm15 = _mm256_load_pd(&c[iz][iy][1][i]);  // _mm256_setzero_pd();
#ifdef AVX2 
																   // multiply the row of a with the first and second column of b and add the result to c  
																   // using FMA 
				ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
				ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
				ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
				ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
				ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
				ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
				ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
				ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
																   // multiply the row of a with the first and second column of b and add the result to c 
																   // using vectorized multiply and add for CPU where FMA is not available 
				ymm8 = _mm256_mul_pd(ymm4, ymm0);
				ymm9 = _mm256_add_pd(ymm8, ymm9);

				ymm8 = _mm256_mul_pd(ymm7, ymm13);
				ymm15 = _mm256_add_pd(ymm8, ymm15);

				ymm8 = _mm256_mul_pd(ymm5, ymm1);
				ymm9 = _mm256_add_pd(ymm8, ymm9);

				ymm8 = _mm256_mul_pd(ymm6, ymm12);
				ymm15 = _mm256_add_pd(ymm8, ymm15);

				ymm8 = _mm256_mul_pd(ymm6, ymm2);
				ymm9 = _mm256_add_pd(ymm8, ymm9);

				ymm8 = _mm256_mul_pd(ymm5, ymm11);
				ymm15 = _mm256_add_pd(ymm8, ymm15);

				ymm8 = _mm256_mul_pd(ymm7, ymm3);
				ymm9 = _mm256_add_pd(ymm8, ymm9);

				ymm8 = _mm256_mul_pd(ymm4, ymm10);
				ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
				// store the result in c 
				_mm256_store_pd(&c[iz][iy][0][i], ymm9);
				_mm256_store_pd(&c[iz][iy][1][i], ymm15);
			}
		}
	}


	ymm14 = _mm256_broadcast_sd(&s[0]);
	// load third column of the small matrix b 
	ymm4 = _mm256_broadcast_sd(&b[2][0]);
	ymm5 = _mm256_broadcast_sd(&b[2][1]);
	ymm6 = _mm256_broadcast_sd(&b[2][2]);
	ymm7 = _mm256_broadcast_sd(&b[2][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load fourth column of the small matrix b 
	ymm4 = _mm256_broadcast_sd(&b[3][0]);
	ymm5 = _mm256_broadcast_sd(&b[3][1]);
	ymm6 = _mm256_broadcast_sd(&b[3][2]);
	ymm7 = _mm256_broadcast_sd(&b[3][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);


	for (iz = 0; iz < 4; iz++)
	{
		for (iy = 0; iy < 4; iy++)
		{
			for (i = 0; i < nVar; i = i + VECTORLENGTH)
			{
				// load VECTORLENGTH rows of a on position i 
				ymm4 = _mm256_load_pd(&a[iz][iy][0][i]);
				ymm5 = _mm256_load_pd(&a[iz][iy][1][i]);
				ymm6 = _mm256_load_pd(&a[iz][iy][2][i]);
				ymm7 = _mm256_load_pd(&a[iz][iy][3][i]);
				// set the result register 1 to zero 
				ymm9 = _mm256_load_pd(&c[iz][iy][2][i]);   // _mm256_setzero_pd();
															   // set the result register 2 to zero
				ymm15 = _mm256_load_pd(&c[iz][iy][3][i]);  // _mm256_setzero_pd();
#ifdef AVX2 
																   // multiply the row of a with the first and second column of b and add the result to c  
																   // using FMA 
				ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
				ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
				ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
				ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
				ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
				ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
				ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
				ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
																   // multiply the row of a with the first and second column of b and add the result to c 
																   // using vectorized multiply and add for CPU where FMA is not available 

				ymm8 = _mm256_mul_pd(ymm4, ymm0);
				ymm9 = _mm256_add_pd(ymm8, ymm9);

				ymm8 = _mm256_mul_pd(ymm7, ymm13);
				ymm15 = _mm256_add_pd(ymm8, ymm15);

				ymm8 = _mm256_mul_pd(ymm5, ymm1);
				ymm9 = _mm256_add_pd(ymm8, ymm9);

				ymm8 = _mm256_mul_pd(ymm6, ymm12);
				ymm15 = _mm256_add_pd(ymm8, ymm15);

				ymm8 = _mm256_mul_pd(ymm6, ymm2);
				ymm9 = _mm256_add_pd(ymm8, ymm9);

				ymm8 = _mm256_mul_pd(ymm5, ymm11);
				ymm15 = _mm256_add_pd(ymm8, ymm15);

				ymm8 = _mm256_mul_pd(ymm7, ymm3);
				ymm9 = _mm256_add_pd(ymm8, ymm9);

				ymm8 = _mm256_mul_pd(ymm4, ymm10);
				ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
				// store the result
				_mm256_store_pd(&c[iz][iy][2][i], ymm9);
				_mm256_store_pd(&c[iz][iy][3][i], ymm15);
			}
		}
	}
	
#endif 
}

// compute all matmuls in y direction on a 3D data set with vector weights (very useful for stiffness matrix multiplications) 
extern "C" void __cdecl vsmm4x4c3d_ky(double c[4][4][4][nVar], double a[4][4][4][nVar], double b[4][4], double s[1])
{
#ifdef AVX512
	__assume_aligned(c, 64);
	__assume_aligned(a, 64);
	__assume_aligned(b, 64);
	int i, j, k, it, ix, iz;
	double weight[1];
	__m512d zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7;
	__m512d zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15;
	__m512d zmm16, zmm17;
	__m512d zmm20, zmm21, zmm22, zmm23, zmm24, zmm25, zmm26, zmm27;
	__m128d xmm;

	// load the scalar into register 14 
	xmm = _mm_load1_pd(&s[0]);
	zmm14 = _mm512_broadcastsd_pd(xmm);
	// load first column of the small matrix b and scale with the scalar s 
	xmm = _mm_load1_pd(&b[0][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm0 = _mm512_mul_pd(zmm4, zmm14);
	zmm1 = _mm512_mul_pd(zmm5, zmm14);
	zmm2 = _mm512_mul_pd(zmm6, zmm14);
	zmm3 = _mm512_mul_pd(zmm7, zmm14);
	// load second column of the small matrix b and scale with the scalar s  
	xmm = _mm_load1_pd(&b[1][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm10 = _mm512_mul_pd(zmm4, zmm14);
	zmm11 = _mm512_mul_pd(zmm5, zmm14);
	zmm12 = _mm512_mul_pd(zmm6, zmm14);
	zmm13 = _mm512_mul_pd(zmm7, zmm14);
	// load third column of the small matrix b and scale with the scalar s 
	xmm = _mm_load1_pd(&b[2][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm20 = _mm512_mul_pd(zmm4, zmm14);
	zmm21 = _mm512_mul_pd(zmm5, zmm14);
	zmm22 = _mm512_mul_pd(zmm6, zmm14);
	zmm23 = _mm512_mul_pd(zmm7, zmm14);
	// load fourth column of the small matrix b and scale with the scalar s  
	xmm = _mm_load1_pd(&b[3][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm24 = _mm512_mul_pd(zmm4, zmm14);
	zmm25 = _mm512_mul_pd(zmm5, zmm14);
	zmm26 = _mm512_mul_pd(zmm6, zmm14);
	zmm27 = _mm512_mul_pd(zmm7, zmm14);

		for (iz = 0; iz < 4; iz++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of matrix a on position i 
					// load c into the result registers 
					// multiply the row of a with the first and second column of b and add the result to c using FMA 
					// multiply the row of a with the third and fourth column of b and add the result to c using FMA 
					zmm4 = _mm512_load_pd(&a[iz][0][ix][i]);
					zmm9 = _mm512_load_pd(&c[iz][0][ix][i]);
					zmm9 = _mm512_fmadd_pd(zmm4, zmm0, zmm9);
					zmm5 = _mm512_load_pd(&a[iz][1][ix][i]);
					zmm15 = _mm512_load_pd(&c[iz][1][ix][i]);
					zmm15 = _mm512_fmadd_pd(zmm5, zmm11, zmm15);
					zmm6 = _mm512_load_pd(&a[iz][2][ix][i]);
					zmm16 = _mm512_load_pd(&c[iz][2][ix][i]);
					zmm16 = _mm512_fmadd_pd(zmm6, zmm22, zmm16);
					zmm7 = _mm512_load_pd(&a[iz][3][ix][i]);
					zmm17 = _mm512_load_pd(&c[iz][3][ix][i]);
					zmm17 = _mm512_fmadd_pd(zmm7, zmm27, zmm17);
					// 
					zmm9 = _mm512_fmadd_pd(zmm5, zmm1, zmm9);
					zmm15 = _mm512_fmadd_pd(zmm6, zmm12, zmm15);
					zmm16 = _mm512_fmadd_pd(zmm7, zmm23, zmm16);
					zmm17 = _mm512_fmadd_pd(zmm4, zmm24, zmm17);
					// 
					zmm9 = _mm512_fmadd_pd(zmm6, zmm2, zmm9);
					zmm15 = _mm512_fmadd_pd(zmm7, zmm13, zmm15);
					zmm16 = _mm512_fmadd_pd(zmm4, zmm20, zmm16);
					zmm17 = _mm512_fmadd_pd(zmm5, zmm25, zmm17);
					// store the result in c 
					zmm9 = _mm512_fmadd_pd(zmm7, zmm3, zmm9);
					_mm512_store_pd(&c[iz][0][ix][i], zmm9);
					zmm15 = _mm512_fmadd_pd(zmm4, zmm10, zmm15);
					_mm512_store_pd(&c[iz][1][ix][i], zmm15);
					zmm16 = _mm512_fmadd_pd(zmm5, zmm21, zmm16);
					_mm512_store_pd(&c[iz][2][ix][i], zmm16);
					zmm17 = _mm512_fmadd_pd(zmm6, zmm26, zmm17);
					_mm512_store_pd(&c[iz][3][ix][i], zmm17);
				}
			}
		}

#else 
	__assume_aligned(c, 32);
	__assume_aligned(a, 32);
	__assume_aligned(b, 32);
	int i, j, k, ix, iz, it;
	double weight[1];
	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	__m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;

	// load the scalar into register 14 
	ymm14 = _mm256_broadcast_sd(&s[0]);
	// load first column of the small matrix b and scale with scalar s 
	ymm4 = _mm256_broadcast_sd(&b[0][0]);
	ymm5 = _mm256_broadcast_sd(&b[0][1]);
	ymm6 = _mm256_broadcast_sd(&b[0][2]);
	ymm7 = _mm256_broadcast_sd(&b[0][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load second column of the small matrix b and scale with scalar s  
	ymm4 = _mm256_broadcast_sd(&b[1][0]);
	ymm5 = _mm256_broadcast_sd(&b[1][1]);
	ymm6 = _mm256_broadcast_sd(&b[1][2]);
	ymm7 = _mm256_broadcast_sd(&b[1][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);

		for (iz = 0; iz < 4; iz++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				// optimized AVX code for matmul on all tensor slices 
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of matrix a on position i 
					ymm4 = _mm256_load_pd(&a[iz][0][ix][i]);
					ymm5 = _mm256_load_pd(&a[iz][1][ix][i]);
					ymm6 = _mm256_load_pd(&a[iz][2][ix][i]);
					ymm7 = _mm256_load_pd(&a[iz][3][ix][i]);
					// load current value of the the result into result register 1  
					ymm9 = _mm256_load_pd(&c[iz][0][ix][i]);    //  _mm256_setzero_pd();
																	// set the result register 2 to zero
					ymm15 = _mm256_load_pd(&c[iz][1][ix][i]);   // _mm256_setzero_pd();
#ifdef AVX2 
																	// multiply the row of a with the first and second column of b and add the result to c  
																	// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
																	// multiply the row of a with the first and second column of b and add the result to c 
																	// using vectorized multiply and add for CPU where FMA is not available 
					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result in c 
					_mm256_store_pd(&c[iz][0][ix][i], ymm9);
					_mm256_store_pd(&c[iz][1][ix][i], ymm15);
				}
			}
		}

	ymm14 = _mm256_broadcast_sd(&s[0]);
	// load third column of the small matrix b 
	ymm4 = _mm256_broadcast_sd(&b[2][0]);
	ymm5 = _mm256_broadcast_sd(&b[2][1]);
	ymm6 = _mm256_broadcast_sd(&b[2][2]);
	ymm7 = _mm256_broadcast_sd(&b[2][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load fourth column of the small matrix b 
	ymm4 = _mm256_broadcast_sd(&b[3][0]);
	ymm5 = _mm256_broadcast_sd(&b[3][1]);
	ymm6 = _mm256_broadcast_sd(&b[3][2]);
	ymm7 = _mm256_broadcast_sd(&b[3][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);

		for (iz = 0; iz < 4; iz++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of a on position i 
					ymm4 = _mm256_load_pd(&a[iz][0][ix][i]);
					ymm5 = _mm256_load_pd(&a[iz][1][ix][i]);
					ymm6 = _mm256_load_pd(&a[iz][2][ix][i]);
					ymm7 = _mm256_load_pd(&a[iz][3][ix][i]);
					// load the result register 1  
					ymm9 = _mm256_load_pd(&c[iz][2][ix][i]);  // _mm256_setzero_pd();
																  // load the result register 2  
					ymm15 = _mm256_load_pd(&c[iz][3][ix][i]);  // _mm256_setzero_pd();
#ifdef AVX2 
																   // multiply the row of a with the first and second column of b and add the result to c  
																   // using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
																   // multiply the row of a with the first and second column of b and add the result to c 
																   // using vectorized multiply and add for CPU where FMA is not available 

					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result
					_mm256_store_pd(&c[iz][2][ix][i], ymm9);
					_mm256_store_pd(&c[iz][3][ix][i], ymm15);
				}
			}
		}
#endif 
}


// compute all matmuls in z direction on a 4D data set with vector weights (very useful for stiffness matrix multiplications) 
extern "C" void __cdecl vsmm4x4c3d_kz(double c[4][4][4][nVar], double a[4][4][4][nVar], double b[4][4], double s[1], double wGP[4])
{
#ifdef AVX512
	__assume_aligned(c, 64);
	__assume_aligned(a, 64);
	__assume_aligned(b, 64);
	int i, j, k, ix, iy, it;
	double weight[1];
	__m512d zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7;
	__m512d zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15;
	__m512d zmm16, zmm17;
	__m512d zmm20, zmm21, zmm22, zmm23, zmm24, zmm25, zmm26, zmm27;
	__m128d xmm;

	// load the scalar into register 14 
	xmm = _mm_load1_pd(&s[0]);
	zmm14 = _mm512_broadcastsd_pd(xmm);
	// load first column of the small matrix b and scale with the scalar s 
	xmm = _mm_load1_pd(&b[0][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[0][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm0 = _mm512_mul_pd(zmm4, zmm14);
	zmm1 = _mm512_mul_pd(zmm5, zmm14);
	zmm2 = _mm512_mul_pd(zmm6, zmm14);
	zmm3 = _mm512_mul_pd(zmm7, zmm14);
	// load second column of the small matrix b and scale with the scalar s  
	xmm = _mm_load1_pd(&b[1][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[1][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm10 = _mm512_mul_pd(zmm4, zmm14);
	zmm11 = _mm512_mul_pd(zmm5, zmm14);
	zmm12 = _mm512_mul_pd(zmm6, zmm14);
	zmm13 = _mm512_mul_pd(zmm7, zmm14);
	// load third column of the small matrix b and scale with the scalar s 
	xmm = _mm_load1_pd(&b[2][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[2][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm20 = _mm512_mul_pd(zmm4, zmm14);
	zmm21 = _mm512_mul_pd(zmm5, zmm14);
	zmm22 = _mm512_mul_pd(zmm6, zmm14);
	zmm23 = _mm512_mul_pd(zmm7, zmm14);
	// load fourth column of the small matrix b and scale with the scalar s  
	xmm = _mm_load1_pd(&b[3][0]);
	zmm4 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][1]);
	zmm5 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][2]);
	zmm6 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&b[3][3]);
	zmm7 = _mm512_broadcastsd_pd(xmm);
	zmm24 = _mm512_mul_pd(zmm4, zmm14);
	zmm25 = _mm512_mul_pd(zmm5, zmm14);
	zmm26 = _mm512_mul_pd(zmm6, zmm14);
	zmm27 = _mm512_mul_pd(zmm7, zmm14);

		for (iy = 0; iy < 4; iy++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				// optimized AVX code for matmul on all tensor slices 
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of matrix a on position i 
					// load c into the result registers 
					// multiply the row of a with the first and second column of b and add the result to c using FMA 
					// multiply the row of a with the third and fourth column of b and add the result to c using FMA 
					zmm4 = _mm512_load_pd(&a[0][iy][ix][i]);
					zmm9 = _mm512_load_pd(&c[0][iy][ix][i]);
					zmm9 = _mm512_fmadd_pd(zmm4, zmm0, zmm9);
					zmm5 = _mm512_load_pd(&a[1][iy][ix][i]);
					zmm15 = _mm512_load_pd(&c[1][iy][ix][i]);
					zmm15 = _mm512_fmadd_pd(zmm5, zmm11, zmm15);
					zmm6 = _mm512_load_pd(&a[2][iy][ix][i]);
					zmm16 = _mm512_load_pd(&c[2][iy][ix][i]);
					zmm16 = _mm512_fmadd_pd(zmm6, zmm22, zmm16);
					zmm7 = _mm512_load_pd(&a[3][iy][ix][i]);
					zmm17 = _mm512_load_pd(&c[3][iy][ix][i]);
					zmm17 = _mm512_fmadd_pd(zmm7, zmm27, zmm17);
					// 
					zmm9 = _mm512_fmadd_pd(zmm5, zmm1, zmm9);
					zmm15 = _mm512_fmadd_pd(zmm6, zmm12, zmm15);
					zmm16 = _mm512_fmadd_pd(zmm7, zmm23, zmm16);
					zmm17 = _mm512_fmadd_pd(zmm4, zmm24, zmm17);
					// 
					zmm9 = _mm512_fmadd_pd(zmm6, zmm2, zmm9);
					zmm15 = _mm512_fmadd_pd(zmm7, zmm13, zmm15);
					zmm16 = _mm512_fmadd_pd(zmm4, zmm20, zmm16);
					zmm17 = _mm512_fmadd_pd(zmm5, zmm25, zmm17);
					// store the result in c 
					zmm9 = _mm512_fmadd_pd(zmm7, zmm3, zmm9);
					_mm512_store_pd(&c[0][iy][ix][i], zmm9);
					zmm15 = _mm512_fmadd_pd(zmm4, zmm10, zmm15);
					_mm512_store_pd(&c[1][iy][ix][i], zmm15);
					zmm16 = _mm512_fmadd_pd(zmm5, zmm21, zmm16);
					_mm512_store_pd(&c[2][iy][ix][i], zmm16);
					zmm17 = _mm512_fmadd_pd(zmm6, zmm26, zmm17);
					_mm512_store_pd(&c[3][iy][ix][i], zmm17);
				}
			}
		}

#else 
	__assume_aligned(c, 32);
	__assume_aligned(a, 32);
	__assume_aligned(b, 32);
	int i, j, k, ix, iy, it;
	double weight[1];
	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	__m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;

	// load the scalar into register 14 
	ymm14 = _mm256_broadcast_sd(&s[0]);
	// load first column of the small matrix b and scale with scalar s 
	ymm4 = _mm256_broadcast_sd(&b[0][0]);
	ymm5 = _mm256_broadcast_sd(&b[0][1]);
	ymm6 = _mm256_broadcast_sd(&b[0][2]);
	ymm7 = _mm256_broadcast_sd(&b[0][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load second column of the small matrix b and scale with scalar s  
	ymm4 = _mm256_broadcast_sd(&b[1][0]);
	ymm5 = _mm256_broadcast_sd(&b[1][1]);
	ymm6 = _mm256_broadcast_sd(&b[1][2]);
	ymm7 = _mm256_broadcast_sd(&b[1][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);

		for (iy = 0; iy < 4; iy++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				// optimized AVX code for matmul on all tensor slices 

				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of matrix a on position i 
					ymm4 = _mm256_load_pd(&a[0][iy][ix][i]);
					ymm5 = _mm256_load_pd(&a[1][iy][ix][i]);
					ymm6 = _mm256_load_pd(&a[2][iy][ix][i]);
					ymm7 = _mm256_load_pd(&a[3][iy][ix][i]);
					// set the result register 1 to zero 
					ymm9 = _mm256_load_pd(&c[0][iy][ix][i]);  // _mm256_setzero_pd();
																  // set the result register 2 to zero
					ymm15 = _mm256_load_pd(&c[1][iy][ix][i]);  // _mm256_setzero_pd();
#ifdef AVX2 
																   // multiply the row of a with the first and second column of b and add the result to c  
																   // using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
																   // multiply the row of a with the first and second column of b and add the result to c 
																   // using vectorized multiply and add for CPU where FMA is not available 
					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result in c 
					_mm256_store_pd(&c[0][iy][ix][i], ymm9);
					_mm256_store_pd(&c[1][iy][ix][i], ymm15);
				}
			}
		}

	ymm14 = _mm256_broadcast_sd(&s[0]);
	// load third column of the small matrix b 
	ymm4 = _mm256_broadcast_sd(&b[2][0]);
	ymm5 = _mm256_broadcast_sd(&b[2][1]);
	ymm6 = _mm256_broadcast_sd(&b[2][2]);
	ymm7 = _mm256_broadcast_sd(&b[2][3]);
	ymm0 = _mm256_mul_pd(ymm4, ymm14);
	ymm1 = _mm256_mul_pd(ymm5, ymm14);
	ymm2 = _mm256_mul_pd(ymm6, ymm14);
	ymm3 = _mm256_mul_pd(ymm7, ymm14);
	// load fourth column of the small matrix b 
	ymm4 = _mm256_broadcast_sd(&b[3][0]);
	ymm5 = _mm256_broadcast_sd(&b[3][1]);
	ymm6 = _mm256_broadcast_sd(&b[3][2]);
	ymm7 = _mm256_broadcast_sd(&b[3][3]);
	ymm10 = _mm256_mul_pd(ymm4, ymm14);
	ymm11 = _mm256_mul_pd(ymm5, ymm14);
	ymm12 = _mm256_mul_pd(ymm6, ymm14);
	ymm13 = _mm256_mul_pd(ymm7, ymm14);

		for (iy = 0; iy < 4; iy++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of a on position i 
					ymm4 = _mm256_load_pd(&a[0][iy][ix][i]);
					ymm5 = _mm256_load_pd(&a[1][iy][ix][i]);
					ymm6 = _mm256_load_pd(&a[2][iy][ix][i]);
					ymm7 = _mm256_load_pd(&a[3][iy][ix][i]);
					// set the result register 1 to zero 
					ymm9 = _mm256_load_pd(&c[2][iy][ix][i]);  // _mm256_setzero_pd();
																  // set the result register 2 to zero
					ymm15 = _mm256_load_pd(&c[3][iy][ix][i]); // _mm256_setzero_pd();
#ifdef AVX2 
																  // multiply the row of a with the first and second column of b and add the result to c  
																  // using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm7, ymm13, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm6, ymm12, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm5, ymm11, ymm15);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
					ymm15 = _mm256_fmadd_pd(ymm4, ymm10, ymm15);
#else
																  // multiply the row of a with the first and second column of b and add the result to c 
																  // using vectorized multiply and add for CPU where FMA is not available 

					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm7, ymm13);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm6, ymm12);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm5, ymm11);
					ymm15 = _mm256_add_pd(ymm8, ymm15);

					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);

					ymm8 = _mm256_mul_pd(ymm4, ymm10);
					ymm15 = _mm256_add_pd(ymm8, ymm15);
#endif 
					// store the result
					_mm256_store_pd(&c[2][iy][ix][i], ymm9);
					_mm256_store_pd(&c[3][iy][ix][i], ymm15);
				}
			}
		}
#endif 
}



// ****************************************** dot-product reductions 

// compute dot-product reduction in x direction on a 3D data set 
extern "C" void __cdecl vsmm4x4ra_x(double c[4][4][nVar], double a[4][4][4][nVar], double wGP[4])
{
#ifdef AVX512 
	__assume_aligned(c, 64);
	__assume_aligned(a, 64);
	int i, j, k, iy, iz;
	double weight[1];
	__m512d zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7;
	__m512d zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15;
	__m128d xmm; 

	// load weights into the vector registers 
	xmm  = _mm_load1_pd(&wGP[0]);
	zmm0 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&wGP[1]);
	zmm1 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&wGP[2]);
	zmm2 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&wGP[3]);
	zmm3 = _mm512_broadcastsd_pd(xmm);

	for (iz = 0; iz < 4; iz++)
	{
		for (iy = 0; iy < 4; iy++)
		{
			// optimized AVX code for matmul on all tensor slices 

			for (i = 0; i < nVar; i = i + VECTORLENGTH)
			{
				// load VECTORLENGTH rows of matrix a on position i 
				// multiply the row of a with the weights wGP and add the result to c  
				// using FMA 
				zmm4 = _mm512_load_pd(&a[iz][iy][0][i]);
				zmm9 = _mm512_mul_pd(zmm4, zmm0);
				zmm5 = _mm512_load_pd(&a[iz][iy][1][i]);
				zmm9 = _mm512_fmadd_pd(zmm5, zmm1, zmm9);
				zmm6 = _mm512_load_pd(&a[iz][iy][2][i]);
				zmm9 = _mm512_fmadd_pd(zmm6, zmm2, zmm9);
				zmm7 = _mm512_load_pd(&a[iz][iy][3][i]);
				zmm9 = _mm512_fmadd_pd(zmm7, zmm3, zmm9);
				// store the result in c 
				_mm512_store_pd(&c[iz][iy][i], zmm9);
			}
		}
	}
#else 
	__assume_aligned(c, 32);
	__assume_aligned(a, 32);
	int i, j, k, iy, iz;
	double weight[1];
	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	__m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;

	// load weights into the vector registers 
	ymm0 = _mm256_broadcast_sd(&wGP[0]);
	ymm1 = _mm256_broadcast_sd(&wGP[1]);
	ymm2 = _mm256_broadcast_sd(&wGP[2]);
	ymm3 = _mm256_broadcast_sd(&wGP[3]);

		for (iz = 0; iz < 4; iz++)
		{
			for (iy = 0; iy < 4; iy++)
			{
				// optimized AVX code for matmul on all tensor slices 

				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of matrix a on position i 
					ymm4 = _mm256_load_pd(&a[iz][iy][0][i]);
					ymm5 = _mm256_load_pd(&a[iz][iy][1][i]);
					ymm6 = _mm256_load_pd(&a[iz][iy][2][i]);
					ymm7 = _mm256_load_pd(&a[iz][iy][3][i]);
					// set the result register to zero 
					ymm9 = _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the weights wGP and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
#else
					// multiply the row of a with the weights wGP and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 
					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);
					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);
					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);
					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);
#endif 
					// store the result in c 
					_mm256_store_pd(&c[iz][iy][i], ymm9);
				}
			}
		}
#endif 
}

// compute dot-product reduction in y direction on a 3D data set 
extern "C" void __cdecl vsmm4x4ra_y(double c[4][4][nVar], double a[4][4][4][nVar], double wGP[4])
{
#ifdef AVX512 
	__assume_aligned(c, 64);
	__assume_aligned(a, 64);
	int i, j, k, ix, iz;
	double weight[1];
	__m512d zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7;
	__m512d zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15;
	__m128d xmm;

	// load weights into the vector registers 
	xmm = _mm_load1_pd(&wGP[0]);
	zmm0 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&wGP[1]);
	zmm1 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&wGP[2]);
	zmm2 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&wGP[3]);
	zmm3 = _mm512_broadcastsd_pd(xmm);

	for (iz = 0; iz < 4; iz++)
	{
		for (ix = 0; ix < 4; ix++)
		{
			// optimized AVX code for matmul on all tensor slices 

			for (i = 0; i < nVar; i = i + VECTORLENGTH)
			{
				// load VECTORLENGTH rows of matrix a on position i 
				// multiply the row of a with the weights wGP and add the result to c  
				// using FMA 
				zmm4 = _mm512_load_pd(&a[iz][0][ix][i]);
				zmm9 = _mm512_mul_pd(zmm4, zmm0);
				zmm5 = _mm512_load_pd(&a[iz][1][ix][i]);
				zmm9 = _mm512_fmadd_pd(zmm5, zmm1, zmm9);
				zmm6 = _mm512_load_pd(&a[iz][2][ix][i]);
				zmm9 = _mm512_fmadd_pd(zmm6, zmm2, zmm9);
				zmm7 = _mm512_load_pd(&a[iz][3][ix][i]);
				zmm9 = _mm512_fmadd_pd(zmm7, zmm3, zmm9);
				// store the result in c 
				_mm512_store_pd(&c[iz][ix][i], zmm9);
			}
		}
	}
#else
	__assume_aligned(c, 32);
	__assume_aligned(a, 32);
	int i, j, k, ix, iz;
	double weight[1];
	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	__m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;

	// load weights into the vector registers 
	ymm0 = _mm256_broadcast_sd(&wGP[0]);
	ymm1 = _mm256_broadcast_sd(&wGP[1]);
	ymm2 = _mm256_broadcast_sd(&wGP[2]);
	ymm3 = _mm256_broadcast_sd(&wGP[3]);

	for (iz = 0; iz < 4; iz++)
	{
		for (ix = 0; ix < 4; ix++)
		{
			// optimized AVX code for matmul on all tensor slices 

			for (i = 0; i < nVar; i = i + VECTORLENGTH)
			{
				// load VECTORLENGTH rows of matrix a on position i 
				ymm4 = _mm256_load_pd(&a[iz][0][ix][i]);
				ymm5 = _mm256_load_pd(&a[iz][1][ix][i]);
				ymm6 = _mm256_load_pd(&a[iz][2][ix][i]);
				ymm7 = _mm256_load_pd(&a[iz][3][ix][i]);
				// set the result register to zero 
				ymm9 = _mm256_setzero_pd();
#ifdef AVX2 
				// multiply the row of a with the weights wGP and add the result to c  
				// using FMA 
				ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
				ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
				ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
				ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
#else
				// multiply the row of a with the weights wGP and add the result to c 
				// using vectorized multiply and add for CPU where FMA is not available 
				ymm8 = _mm256_mul_pd(ymm4, ymm0);
				ymm9 = _mm256_add_pd(ymm8, ymm9);
				ymm8 = _mm256_mul_pd(ymm5, ymm1);
				ymm9 = _mm256_add_pd(ymm8, ymm9);
				ymm8 = _mm256_mul_pd(ymm6, ymm2);
				ymm9 = _mm256_add_pd(ymm8, ymm9);
				ymm8 = _mm256_mul_pd(ymm7, ymm3);
				ymm9 = _mm256_add_pd(ymm8, ymm9);
#endif 
				// store the result in c 
				_mm256_store_pd(&c[iz][ix][i], ymm9);
			}
		}
	}
#endif 
}
// compute dot-product reduction in z direction on a 3D data set 
extern "C" void __cdecl vsmm4x4ra_z(double c[4][4][nVar], double a[4][4][4][nVar], double wGP[4])
{
#ifdef AVX512 
	__assume_aligned(c, 64);
	__assume_aligned(a, 64);
	int i, j, k, ix, iy;
	double weight[1];
	__m512d zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7;
	__m512d zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15;
	__m128d xmm;

	// load weights into the vector registers 
	xmm = _mm_load1_pd(&wGP[0]);
	zmm0 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&wGP[1]);
	zmm1 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&wGP[2]);
	zmm2 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&wGP[3]);
	zmm3 = _mm512_broadcastsd_pd(xmm);

	for (iy = 0; iy < 4; iy++)
	{
		for (ix = 0; ix < 4; ix++)
		{
			// optimized AVX code for matmul on all tensor slices 

			for (i = 0; i < nVar; i = i + VECTORLENGTH) 
			{
				// load VECTORLENGTH rows of matrix a on position i 
				// multiply the row of a with the weights wGP and add the result to c  
				// using FMA 
				zmm4 = _mm512_load_pd(&a[0][iy][ix][i]);
				zmm9 = _mm512_mul_pd(zmm4, zmm0);
				zmm5 = _mm512_load_pd(&a[1][iy][ix][i]);
				zmm9 = _mm512_fmadd_pd(zmm5, zmm1, zmm9);
				zmm6 = _mm512_load_pd(&a[2][iy][ix][i]);
				zmm9 = _mm512_fmadd_pd(zmm6, zmm2, zmm9);
				zmm7 = _mm512_load_pd(&a[3][iy][ix][i]);
				zmm9 = _mm512_fmadd_pd(zmm7, zmm3, zmm9);
				// store the result in c 
				_mm512_store_pd(&c[iy][ix][i], zmm9);
			}
		}
	}
#else
	__assume_aligned(c, 32);
	__assume_aligned(a, 32);
	int i, j, k, ix, iy;
	double weight[1];
	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	__m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;

	// load weights into the vector registers 
	ymm0 = _mm256_broadcast_sd(&wGP[0]);
	ymm1 = _mm256_broadcast_sd(&wGP[1]);
	ymm2 = _mm256_broadcast_sd(&wGP[2]);
	ymm3 = _mm256_broadcast_sd(&wGP[3]);

	for (iy = 0; iy < 4; iy++)
	{
		for (ix = 0; ix < 4; ix++)
		{
			// optimized AVX code for matmul on all tensor slices 

			for (i = 0; i < nVar; i = i + VECTORLENGTH)
			{
				// load VECTORLENGTH rows of matrix a on position i 
				ymm4 = _mm256_load_pd(&a[0][iy][ix][i]);
				ymm5 = _mm256_load_pd(&a[1][iy][ix][i]);
				ymm6 = _mm256_load_pd(&a[2][iy][ix][i]);
				ymm7 = _mm256_load_pd(&a[3][iy][ix][i]);
				// set the result register to zero 
				ymm9 = _mm256_setzero_pd();
#ifdef AVX2 
				// multiply the row of a with the weights wGP and add the result to c  
				// using FMA 
				ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
				ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
				ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
				ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
#else
				// multiply the row of a with the weights wGP and add the result to c 
				// using vectorized multiply and add for CPU where FMA is not available 
				ymm8 = _mm256_mul_pd(ymm4, ymm0);
				ymm9 = _mm256_add_pd(ymm8, ymm9);
				ymm8 = _mm256_mul_pd(ymm5, ymm1);
				ymm9 = _mm256_add_pd(ymm8, ymm9);
				ymm8 = _mm256_mul_pd(ymm6, ymm2);
				ymm9 = _mm256_add_pd(ymm8, ymm9);
				ymm8 = _mm256_mul_pd(ymm7, ymm3);
				ymm9 = _mm256_add_pd(ymm8, ymm9);
#endif 
				// store the result in c 
				_mm256_store_pd(&c[iy][ix][i], ymm9);
			}
		}
	}
#endif 
}

// compute dot-product reduction in t direction on a 4D data set 
extern "C" void __cdecl vsmm4x4ra_t(double c[4][4][4][nVar], double a[4][4][4][4][nVar], double wGP[4])
{
#ifdef AVX512 
	__assume_aligned(c, 64);
	__assume_aligned(a, 64);


	__asm { 
		mov rax, a	                                        // pointer to a 
		mov rbx, wGP										// pointer to wGP  
		mov rcx, c											// pointer to c 
		 
		vbroadcastsd zmm0, QWORD PTR[rbx + 0 * 8 ]			// load wGP[0] 
		vbroadcastsd zmm1, QWORD PTR[rbx + 1 * 8 ]			// load wGP[1] 
		vbroadcastsd zmm2, QWORD PTR[rbx + 2 * 8 ]			// load wGP[2] 
		vbroadcastsd zmm3, QWORD PTR[rbx + 3 * 8 ]			// load wGP[3] 

		mov r10, 0
loopz:
		mov r15, 0  
loopvar: 

			// compute address of the pointers 
			mov r14,  0
			add r14,  r10        // z loop index 
			add r14,  r15        // variable loop index 

			mov r12, rax	     // compute a base pointer 
			add r12, r14        
			mov r13, rcx         // compute c base pointer 
			add r13, r14

			// first part  

			vmovapd		zmm4,  ZMMWORD PTR[r12 + 0 * ND3D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8]		
			vmovapd		zmm8,  ZMMWORD PTR[r12 + 0 * ND3D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8]
			vmovapd		zmm12, ZMMWORD PTR[r12 + 0 * ND3D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8]
			vmovapd		zmm16, ZMMWORD PTR[r12 + 0 * ND3D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8]
			vmovapd		zmm5,  ZMMWORD PTR[r12 + 1 * ND3D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8]		
			vmovapd		zmm9,  ZMMWORD PTR[r12 + 1 * ND3D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8]
			vmovapd		zmm13, ZMMWORD PTR[r12 + 1 * ND3D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8]
			vmovapd		zmm17, ZMMWORD PTR[r12 + 1 * ND3D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8]
			vmovapd		zmm6,  ZMMWORD PTR[r12 + 2 * ND3D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8]		
			vmovapd		zmm10, ZMMWORD PTR[r12 + 2 * ND3D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8]
			vmovapd		zmm14, ZMMWORD PTR[r12 + 2 * ND3D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8]
			vmovapd		zmm18, ZMMWORD PTR[r12 + 2 * ND3D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8]
			vmovapd		zmm7,  ZMMWORD PTR[r12 + 3 * ND3D*nVar * 8 + 0 * ND1D*nVar * 8 + 0 * nVar * 8]		
			vmovapd		zmm11, ZMMWORD PTR[r12 + 3 * ND3D*nVar * 8 + 0 * ND1D*nVar * 8 + 1 * nVar * 8]
			vmovapd		zmm15, ZMMWORD PTR[r12 + 3 * ND3D*nVar * 8 + 0 * ND1D*nVar * 8 + 2 * nVar * 8]
			vmovapd		zmm19, ZMMWORD PTR[r12 + 3 * ND3D*nVar * 8 + 0 * ND1D*nVar * 8 + 3 * nVar * 8]
			vmovapd		zmm28, ZMMWORD PTR[r12 + 0 * ND3D*nVar * 8 + 1 * ND1D*nVar * 8 + 3 * nVar * 8]
			vmovapd		zmm29, ZMMWORD PTR[r12 + 1 * ND3D*nVar * 8 + 1 * ND1D*nVar * 8 + 3 * nVar * 8]
			vmovapd		zmm30, ZMMWORD PTR[r12 + 2 * ND3D*nVar * 8 + 1 * ND1D*nVar * 8 + 3 * nVar * 8]
			vmovapd		zmm31, ZMMWORD PTR[r12 + 3 * ND3D*nVar * 8 + 1 * ND1D*nVar * 8 + 3 * nVar * 8]

			vmulpd		zmm20, zmm0, zmm4
			vmulpd		zmm21, zmm0, zmm8
			vmulpd		zmm22, zmm0, zmm12
			vmulpd		zmm23, zmm0, zmm16
			vmulpd      zmm24, zmm0, ZMMWORD PTR[r12 + 0 * ND3D*nVar * 8 + 1 * ND1D*nVar * 8 + 0 * nVar * 8]
			vmulpd      zmm25, zmm0, ZMMWORD PTR[r12 + 0 * ND3D*nVar * 8 + 1 * ND1D*nVar * 8 + 1 * nVar * 8]
			vmulpd      zmm26, zmm0, ZMMWORD PTR[r12 + 0 * ND3D*nVar * 8 + 1 * ND1D*nVar * 8 + 2 * nVar * 8]
			vmulpd      zmm27, zmm0, zmm28 

			prefetchw [r13 + 0 * ND1D*nVar * 8 + 0 * nVar * 8]
			prefetchw [r13 + 0 * ND1D*nVar * 8 + 1 * nVar * 8]
			prefetchw [r13 + 0 * ND1D*nVar * 8 + 2 * nVar * 8]
			prefetchw [r13 + 0 * ND1D*nVar * 8 + 3 * nVar * 8]
			prefetchw [r13 + 1 * ND1D*nVar * 8 + 0 * nVar * 8]
			prefetchw [r13 + 1 * ND1D*nVar * 8 + 1 * nVar * 8]
			prefetchw [r13 + 1 * ND1D*nVar * 8 + 2 * nVar * 8]
			prefetchw [r13 + 1 * ND1D*nVar * 8 + 3 * nVar * 8]

			vfmadd231pd zmm20, zmm1, zmm5
			vfmadd231pd zmm21, zmm1, zmm9
			vfmadd231pd zmm22, zmm1, zmm13
			vfmadd231pd zmm23, zmm1, zmm17
			vfmadd231pd zmm24, zmm1, ZMMWORD PTR[r12 + 1 * ND3D*nVar * 8 + 1 * ND1D*nVar * 8 + 0 * nVar * 8]
			vfmadd231pd zmm25, zmm1, ZMMWORD PTR[r12 + 1 * ND3D*nVar * 8 + 1 * ND1D*nVar * 8 + 1 * nVar * 8]
			vfmadd231pd zmm26, zmm1, ZMMWORD PTR[r12 + 1 * ND3D*nVar * 8 + 1 * ND1D*nVar * 8 + 2 * nVar * 8]
			vfmadd231pd zmm27, zmm1, zmm29 

			vfmadd231pd zmm20, zmm2, zmm6
			vfmadd231pd zmm21, zmm2, zmm10
			vfmadd231pd zmm22, zmm2, zmm14
			vfmadd231pd zmm23, zmm2, zmm18
			vfmadd231pd zmm24, zmm2, ZMMWORD PTR[r12 + 2 * ND3D*nVar * 8 + 1 * ND1D*nVar * 8 + 0 * nVar * 8]
			vfmadd231pd zmm25, zmm2, ZMMWORD PTR[r12 + 2 * ND3D*nVar * 8 + 1 * ND1D*nVar * 8 + 1 * nVar * 8]
			vfmadd231pd zmm26, zmm2, ZMMWORD PTR[r12 + 2 * ND3D*nVar * 8 + 1 * ND1D*nVar * 8 + 2 * nVar * 8]
			vfmadd231pd zmm27, zmm2, zmm30 

			vfmadd231pd zmm20, zmm3, zmm7
			vfmadd231pd zmm21, zmm3, zmm11
			vfmadd231pd zmm22, zmm3, zmm15
			vfmadd231pd zmm23, zmm3, zmm19
			vfmadd231pd zmm24, zmm3, ZMMWORD PTR[r12 + 3 * ND3D*nVar * 8 + 1 * ND1D*nVar * 8 + 0 * nVar * 8]
			vfmadd231pd zmm25, zmm3, ZMMWORD PTR[r12 + 3 * ND3D*nVar * 8 + 1 * ND1D*nVar * 8 + 1 * nVar * 8]
			vfmadd231pd zmm26, zmm3, ZMMWORD PTR[r12 + 3 * ND3D*nVar * 8 + 1 * ND1D*nVar * 8 + 2 * nVar * 8]
			vfmadd231pd zmm27, zmm3, zmm31 

			vmovapd		ZMMWORD PTR[r13 + 0 * ND1D*nVar * 8 + 0 * nVar * 8], zmm20 
			vmovapd		ZMMWORD PTR[r13 + 0 * ND1D*nVar * 8 + 1 * nVar * 8], zmm21
			vmovapd		ZMMWORD PTR[r13 + 0 * ND1D*nVar * 8 + 2 * nVar * 8], zmm22
			vmovapd		ZMMWORD PTR[r13 + 0 * ND1D*nVar * 8 + 3 * nVar * 8], zmm23
			vmovapd		ZMMWORD PTR[r13 + 1 * ND1D*nVar * 8 + 0 * nVar * 8], zmm24
			vmovapd		ZMMWORD PTR[r13 + 1 * ND1D*nVar * 8 + 1 * nVar * 8], zmm25
			vmovapd		ZMMWORD PTR[r13 + 1 * ND1D*nVar * 8 + 2 * nVar * 8], zmm26
			vmovapd		ZMMWORD PTR[r13 + 1 * ND1D*nVar * 8 + 3 * nVar * 8], zmm27

			// second part 

			vmovapd		zmm4,  ZMMWORD PTR[r12 + 0 * ND3D*nVar * 8 + 2 * ND1D*nVar * 8 + 0 * nVar * 8]
			vmovapd		zmm8,  ZMMWORD PTR[r12 + 0 * ND3D*nVar * 8 + 2 * ND1D*nVar * 8 + 1 * nVar * 8]
			vmovapd		zmm12, ZMMWORD PTR[r12 + 0 * ND3D*nVar * 8 + 2 * ND1D*nVar * 8 + 2 * nVar * 8]
			vmovapd		zmm16, ZMMWORD PTR[r12 + 0 * ND3D*nVar * 8 + 2 * ND1D*nVar * 8 + 3 * nVar * 8]
			vmovapd		zmm5,  ZMMWORD PTR[r12 + 1 * ND3D*nVar * 8 + 2 * ND1D*nVar * 8 + 0 * nVar * 8]
			vmovapd		zmm9,  ZMMWORD PTR[r12 + 1 * ND3D*nVar * 8 + 2 * ND1D*nVar * 8 + 1 * nVar * 8]
			vmovapd		zmm13, ZMMWORD PTR[r12 + 1 * ND3D*nVar * 8 + 2 * ND1D*nVar * 8 + 2 * nVar * 8]
			vmovapd		zmm17, ZMMWORD PTR[r12 + 1 * ND3D*nVar * 8 + 2 * ND1D*nVar * 8 + 3 * nVar * 8]
			vmovapd		zmm6,  ZMMWORD PTR[r12 + 2 * ND3D*nVar * 8 + 2 * ND1D*nVar * 8 + 0 * nVar * 8]
			vmovapd		zmm10, ZMMWORD PTR[r12 + 2 * ND3D*nVar * 8 + 2 * ND1D*nVar * 8 + 1 * nVar * 8]
			vmovapd		zmm14, ZMMWORD PTR[r12 + 2 * ND3D*nVar * 8 + 2 * ND1D*nVar * 8 + 2 * nVar * 8]
			vmovapd		zmm18, ZMMWORD PTR[r12 + 2 * ND3D*nVar * 8 + 2 * ND1D*nVar * 8 + 3 * nVar * 8]
			vmovapd		zmm7,  ZMMWORD PTR[r12 + 3 * ND3D*nVar * 8 + 2 * ND1D*nVar * 8 + 0 * nVar * 8]
			vmovapd		zmm11, ZMMWORD PTR[r12 + 3 * ND3D*nVar * 8 + 2 * ND1D*nVar * 8 + 1 * nVar * 8]
			vmovapd		zmm15, ZMMWORD PTR[r12 + 3 * ND3D*nVar * 8 + 2 * ND1D*nVar * 8 + 2 * nVar * 8]
			vmovapd		zmm19, ZMMWORD PTR[r12 + 3 * ND3D*nVar * 8 + 2 * ND1D*nVar * 8 + 3 * nVar * 8]
			vmovapd		zmm28, ZMMWORD PTR[r12 + 0 * ND3D*nVar * 8 + 3 * ND1D*nVar * 8 + 3 * nVar * 8]
			vmovapd		zmm29, ZMMWORD PTR[r12 + 1 * ND3D*nVar * 8 + 3 * ND1D*nVar * 8 + 3 * nVar * 8]
			vmovapd		zmm30, ZMMWORD PTR[r12 + 2 * ND3D*nVar * 8 + 3 * ND1D*nVar * 8 + 3 * nVar * 8]
			vmovapd		zmm31, ZMMWORD PTR[r12 + 3 * ND3D*nVar * 8 + 3 * ND1D*nVar * 8 + 3 * nVar * 8]

			vmulpd		zmm20, zmm0, zmm4
			vmulpd		zmm21, zmm0, zmm8
			vmulpd		zmm22, zmm0, zmm12
			vmulpd		zmm23, zmm0, zmm16
			vmulpd      zmm24, zmm0, ZMMWORD PTR[r12 + 0 * ND3D*nVar * 8 + 3 * ND1D*nVar * 8 + 0 * nVar * 8]
			vmulpd      zmm25, zmm0, ZMMWORD PTR[r12 + 0 * ND3D*nVar * 8 + 3 * ND1D*nVar * 8 + 1 * nVar * 8]
			vmulpd      zmm26, zmm0, ZMMWORD PTR[r12 + 0 * ND3D*nVar * 8 + 3 * ND1D*nVar * 8 + 2 * nVar * 8]
			vmulpd      zmm27, zmm0, zmm28 

			prefetchw [r13 + 2 * ND1D*nVar * 8 + 0 * nVar * 8]
			prefetchw [r13 + 2 * ND1D*nVar * 8 + 1 * nVar * 8]
			prefetchw [r13 + 2 * ND1D*nVar * 8 + 2 * nVar * 8]
			prefetchw [r13 + 2 * ND1D*nVar * 8 + 3 * nVar * 8]
			prefetchw [r13 + 3 * ND1D*nVar * 8 + 0 * nVar * 8]
			prefetchw [r13 + 3 * ND1D*nVar * 8 + 1 * nVar * 8]
			prefetchw [r13 + 3 * ND1D*nVar * 8 + 2 * nVar * 8]
			prefetchw [r13 + 3 * ND1D*nVar * 8 + 3 * nVar * 8]

			vfmadd231pd zmm20, zmm1, zmm5
			vfmadd231pd zmm21, zmm1, zmm9
			vfmadd231pd zmm22, zmm1, zmm13
			vfmadd231pd zmm23, zmm1, zmm17
			vfmadd231pd zmm24, zmm1, ZMMWORD PTR[r12 + 1 * ND3D*nVar * 8 + 3 * ND1D*nVar * 8 + 0 * nVar * 8]
			vfmadd231pd zmm25, zmm1, ZMMWORD PTR[r12 + 1 * ND3D*nVar * 8 + 3 * ND1D*nVar * 8 + 1 * nVar * 8]
			vfmadd231pd zmm26, zmm1, ZMMWORD PTR[r12 + 1 * ND3D*nVar * 8 + 3 * ND1D*nVar * 8 + 2 * nVar * 8]
			vfmadd231pd zmm27, zmm1, zmm29 

			vfmadd231pd zmm20, zmm2, zmm6
			vfmadd231pd zmm21, zmm2, zmm10
			vfmadd231pd zmm22, zmm2, zmm14
			vfmadd231pd zmm23, zmm2, zmm18
			vfmadd231pd zmm24, zmm2, ZMMWORD PTR[r12 + 2 * ND3D*nVar * 8 + 3 * ND1D*nVar * 8 + 0 * nVar * 8]
			vfmadd231pd zmm25, zmm2, ZMMWORD PTR[r12 + 2 * ND3D*nVar * 8 + 3 * ND1D*nVar * 8 + 1 * nVar * 8]
			vfmadd231pd zmm26, zmm2, ZMMWORD PTR[r12 + 2 * ND3D*nVar * 8 + 3 * ND1D*nVar * 8 + 2 * nVar * 8]
			vfmadd231pd zmm27, zmm2, zmm30 

			vfmadd231pd zmm20, zmm3, zmm7
			vfmadd231pd zmm21, zmm3, zmm11
			vfmadd231pd zmm22, zmm3, zmm15
			vfmadd231pd zmm23, zmm3, zmm19
			vfmadd231pd zmm24, zmm3, ZMMWORD PTR[r12 + 3 * ND3D*nVar * 8 + 3 * ND1D*nVar * 8 + 0 * nVar * 8]
			vfmadd231pd zmm25, zmm3, ZMMWORD PTR[r12 + 3 * ND3D*nVar * 8 + 3 * ND1D*nVar * 8 + 1 * nVar * 8]
			vfmadd231pd zmm26, zmm3, ZMMWORD PTR[r12 + 3 * ND3D*nVar * 8 + 3 * ND1D*nVar * 8 + 2 * nVar * 8]
			vfmadd231pd zmm27, zmm3, zmm31

			vmovapd		ZMMWORD PTR[r13 + 2 * ND1D*nVar * 8 + 0 * nVar * 8], zmm20
			vmovapd		ZMMWORD PTR[r13 + 2 * ND1D*nVar * 8 + 1 * nVar * 8], zmm21
			vmovapd		ZMMWORD PTR[r13 + 2 * ND1D*nVar * 8 + 2 * nVar * 8], zmm22
			vmovapd		ZMMWORD PTR[r13 + 2 * ND1D*nVar * 8 + 3 * nVar * 8], zmm23
			vmovapd		ZMMWORD PTR[r13 + 3 * ND1D*nVar * 8 + 0 * nVar * 8], zmm24
			vmovapd		ZMMWORD PTR[r13 + 3 * ND1D*nVar * 8 + 1 * nVar * 8], zmm25
			vmovapd		ZMMWORD PTR[r13 + 3 * ND1D*nVar * 8 + 2 * nVar * 8], zmm26
			vmovapd		ZMMWORD PTR[r13 + 3 * ND1D*nVar * 8 + 3 * nVar * 8], zmm27


		add r15, VECTORLENGTH*8 
		cmp r15, nVar*8  
		jl loopvar

		add r10,   ND2D*nVar*8   
		cmp r10, 4*ND2D*nVar*8 
		jl loopz

	}

  /* 

	int i, j, k, ix, iy, iz;
	double weight[1];
	__m512d zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7;
	__m512d zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15;
	__m128d xmm;

	// load weights into the vector registers 
	xmm = _mm_load1_pd(&wGP[0]);
	zmm0 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&wGP[1]);
	zmm1 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&wGP[2]);
	zmm2 = _mm512_broadcastsd_pd(xmm);
	xmm = _mm_load1_pd(&wGP[3]);
	zmm3 = _mm512_broadcastsd_pd(xmm);

	for (iz = 0; iz < 4; iz++)
	{
		for (iy = 0; iy < 4; iy++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				// optimized AVX code for matmul on all tensor slices 

				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of matrix a on position i 
					// multiply the row of a with the weights wGP and add the result to c  
					// using FMA 
					zmm4 = _mm512_load_pd(&a[0][iz][iy][ix][i]);
					zmm9 = _mm512_mul_pd(zmm4, zmm0);
					zmm5 = _mm512_load_pd(&a[1][iz][iy][ix][i]);
					zmm9 = _mm512_fmadd_pd(zmm5, zmm1, zmm9);
					zmm6 = _mm512_load_pd(&a[2][iz][iy][ix][i]);
					zmm9 = _mm512_fmadd_pd(zmm6, zmm2, zmm9);
					zmm7 = _mm512_load_pd(&a[3][iz][iy][ix][i]);
					zmm9 = _mm512_fmadd_pd(zmm7, zmm3, zmm9);
					// store the result in c 
					_mm512_store_pd(&c[iz][iy][ix][i], zmm9);
				}
			}
		}
	}

	*/ 

#else
	__assume_aligned(c, 32);
	__assume_aligned(a, 32);
	int i, j, k, ix, iy, iz;
	double weight[1];
	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	__m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;

	// load weights into the vector registers 
	ymm0 = _mm256_broadcast_sd(&wGP[0]);
	ymm1 = _mm256_broadcast_sd(&wGP[1]);
	ymm2 = _mm256_broadcast_sd(&wGP[2]);
	ymm3 = _mm256_broadcast_sd(&wGP[3]); 
	
	for (iz = 0; iz < 4; iz++)
	{
		for (iy = 0; iy < 4; iy++)
		{
			for (ix = 0; ix < 4; ix++)
			{
				// optimized AVX code for matmul on all tensor slices 

				for (i = 0; i < nVar; i = i + VECTORLENGTH)
				{
					// load VECTORLENGTH rows of matrix a on position i 
					ymm4 = _mm256_load_pd(&a[0][iz][iy][ix][i]);
					ymm5 = _mm256_load_pd(&a[1][iz][iy][ix][i]);
					ymm6 = _mm256_load_pd(&a[2][iz][iy][ix][i]);
					ymm7 = _mm256_load_pd(&a[3][iz][iy][ix][i]);
					// set the result register to zero 
					ymm9 = _mm256_setzero_pd();
#ifdef AVX2 
					// multiply the row of a with the weights wGP and add the result to c  
					// using FMA 
					ymm9 = _mm256_fmadd_pd(ymm4, ymm0, ymm9);
					ymm9 = _mm256_fmadd_pd(ymm5, ymm1, ymm9);
					ymm9 = _mm256_fmadd_pd(ymm6, ymm2, ymm9);
					ymm9 = _mm256_fmadd_pd(ymm7, ymm3, ymm9);
#else
					// multiply the row of a with the weights wGP and add the result to c 
					// using vectorized multiply and add for CPU where FMA is not available 
					ymm8 = _mm256_mul_pd(ymm4, ymm0);
					ymm9 = _mm256_add_pd(ymm8, ymm9);
					ymm8 = _mm256_mul_pd(ymm5, ymm1);
					ymm9 = _mm256_add_pd(ymm8, ymm9);
					ymm8 = _mm256_mul_pd(ymm6, ymm2);
					ymm9 = _mm256_add_pd(ymm8, ymm9);
					ymm8 = _mm256_mul_pd(ymm7, ymm3);
					ymm9 = _mm256_add_pd(ymm8, ymm9);
#endif 
					// store the result in c 
					_mm256_store_pd(&c[iz][iy][ix][i], ymm9);
				}
			}
		}
	}
#endif 
}


