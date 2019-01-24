// Direct implementation of the PDE system with AVX intrinsics 

#define DONT_COUNT_FLOPS 
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

extern "C" void __cdecl pdefluxvectoravx( double F[3][nVar][VECTORLENGTH], double V[nVar][VECTORLENGTH], double Q[nVar][VECTORLENGTH] )
{
#ifdef AVX512

	__assume_aligned(F, 64);
	__assume_aligned(Q, 64);

	int i, j, k;
	__m512d zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7;
	__m512d zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15;
	__m128d xmm; 

	double gamma[1], gamma1[1];
	gamma[0] = 1.4;
	gamma1[0] = gamma[0] - 1.0;			// precompute gamma - 1 

	zmm0 = _mm512_load_pd(&Q[0][0]);		// load density
	zmm1 = _mm512_set1_pd(1.0);				// set unit 
	zmm2 = _mm512_div_pd(zmm1, zmm0);		// compute density^{-1} 

	zmm3 = _mm512_load_pd(&Q[1][0]);		// load x momentum 
	zmm4 = _mm512_mul_pd(zmm3, zmm2);		// compute x velocity 
	zmm5 = _mm512_load_pd(&Q[2][0]);		// load y momentum 
	zmm6 = _mm512_mul_pd(zmm5, zmm2);		// compute y velocity 
	zmm7 = _mm512_load_pd(&Q[3][0]);		// load z momentum 
	zmm8 = _mm512_mul_pd(zmm7, zmm2);		// compute z velocity 

	zmm9 = _mm512_load_pd(&Q[4][0]);		 // load total energy 
	xmm = _mm_load1_pd(&gamma1[0]);
	zmm1 = _mm512_broadcastsd_pd(xmm);

	// zmm1 = _mm512_broadcast_sd(&gamma1[0]); // broadcast gamma-1 
	zmm15 = _mm512_set1_pd(-0.5);			 // set -1/2 
											 // 
											 // compute the pressure via AVX instructions (AVX2 does not make sense here, since we also need the intermediate result of the products 
											 // p(:) = (EQN%gamma - 1)*(Q(:, 5) - 0.5*Q(:, 1)*(uu(:)*uu(:) + vv(:)*vv(:) + ww(:)*ww(:)))
											 // 
	zmm11 = _mm512_mul_pd(zmm3, zmm4);				// compute rho*u*u 
	zmm12 = _mm512_mul_pd(zmm5, zmm6);				// compute rho*v*v 
	zmm13 = _mm512_mul_pd(zmm7, zmm8);				// compute rho*w*w 
	zmm10 = _mm512_add_pd(zmm11, zmm12);				// compute the kinetic energy * 2 
	zmm10 = _mm512_add_pd(zmm10, zmm13);
	zmm10 = _mm512_mul_pd(zmm10, zmm15);			// divide the result by 2 and invert sign 
	zmm10 = _mm512_add_pd(zmm10, zmm9);				// subtract the kinetic energy density from the total energy density to get the internal energy density 
	zmm10 = _mm512_mul_pd(zmm10, zmm1);				// multiply the internal energy density with gamma-1 to get the pressure 
													// now compute the tensor of the nonlinear convective terms via AVX instructions 
	zmm0 = _mm512_mul_pd(zmm3, zmm6);				// compute rho*u*v  
	zmm1 = _mm512_mul_pd(zmm3, zmm8);				// compute rho*u*w  
	zmm2 = _mm512_mul_pd(zmm5, zmm8);				// compute rho*v*w 
													// add the pressure to the diagonal 
	zmm11 = _mm512_add_pd(zmm11, zmm10);			// compute rho*u*u + p 
	zmm12 = _mm512_add_pd(zmm12, zmm10);			// compute rho*v*v + p 
	zmm13 = _mm512_add_pd(zmm13, zmm10);			// compute rho*w*w + p 
	zmm14 = _mm512_add_pd(zmm9, zmm10);				// compute rhoE + p 

													// store the results in F 
	zmm15 = _mm512_mul_pd(zmm14, zmm4);				// compute u*(rhoE+p) 
	_mm512_store_pd(&F[0][0][0], zmm3);				// rho*u 
	_mm512_store_pd(&F[0][1][0], zmm11);			// rho*u*u + p   
	_mm512_store_pd(&F[0][2][0], zmm0);				// rho*u*v 
	_mm512_store_pd(&F[0][3][0], zmm1);				// rho*u*w 
	_mm512_store_pd(&F[0][4][0], zmm15);			// u*(rhoE+p)  

	zmm15 = _mm512_mul_pd(zmm14, zmm6);				// compute v*(rhoE+p) 
	_mm512_store_pd(&F[1][0][0], zmm5);				// rho*v 
	_mm512_store_pd(&F[1][1][0], zmm0);				// rho*v*u 
	_mm512_store_pd(&F[1][2][0], zmm12);			// rho*v*v + p  
	_mm512_store_pd(&F[1][3][0], zmm2);				// rho*v*w 
	_mm512_store_pd(&F[1][4][0], zmm15);			// v*(rhoE+p)  

	zmm15 = _mm512_mul_pd(zmm14, zmm8);				// compute w*(rhoE+p) 
	_mm512_store_pd(&F[2][0][0], zmm7);				// rho*w 
	_mm512_store_pd(&F[2][1][0], zmm1);				// rho*w*u 
	_mm512_store_pd(&F[2][2][0], zmm2);				// rho*w*v   
	_mm512_store_pd(&F[2][3][0], zmm13);			// rho*w*w + p 
	_mm512_store_pd(&F[2][4][0], zmm15);			// v*(rhoE+p)  

#else 
	__assume_aligned(F, 32);
	__assume_aligned(Q, 32);
	
	int i, j, k;
	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	__m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;

	double gamma[1], gamma1[1]; 
	gamma[0]   = 1.4; 
	gamma1[0]  = gamma[0] - 1.0;			// precompute gamma - 1 

	ymm0 = _mm256_load_pd(&Q[0][0]);		// load density
	ymm1 = _mm256_set1_pd(1.0);				// set unit 
	ymm2 = _mm256_div_pd(ymm1, ymm0);		// compute density^{-1} 

	ymm3 = _mm256_load_pd(&Q[1][0]);		// load x momentum 
	ymm4 = _mm256_mul_pd(ymm3, ymm2);		// compute x velocity 
	ymm5 = _mm256_load_pd(&Q[2][0]);		// load y momentum 
	ymm6 = _mm256_mul_pd(ymm5, ymm2);		// compute y velocity 
	ymm7 = _mm256_load_pd(&Q[3][0]);		// load z momentum 
	ymm8 = _mm256_mul_pd(ymm7, ymm2);		// compute z velocity 
	
	ymm9  = _mm256_load_pd(&Q[4][0]);		 // load total energy 
	ymm1  = _mm256_broadcast_sd(&gamma1[0]); // broadcast gamma-1 
	ymm15 = _mm256_set1_pd(-0.5);			 // set -1/2 
	// 
	// compute the pressure via AVX instructions (AVX2 does not make sense here, since we also need the intermediate result of the products 
	// p(:) = (EQN%gamma - 1)*(Q(:, 5) - 0.5*Q(:, 1)*(uu(:)*uu(:) + vv(:)*vv(:) + ww(:)*ww(:)))
	// 
	ymm11 = _mm256_mul_pd(ymm3, ymm4);				// compute rho*u*u 
	ymm12 = _mm256_mul_pd(ymm5, ymm6);				// compute rho*v*v 
	ymm13 = _mm256_mul_pd(ymm7, ymm8);				// compute rho*w*w 
	ymm10 = _mm256_add_pd(ymm11,ymm12);				// compute the kinetic energy * 2 
	ymm10 = _mm256_add_pd(ymm10,ymm13); 
	ymm10 = _mm256_mul_pd(ymm10, ymm15);			// divide the result by 2 and invert sign 
	ymm10 = _mm256_add_pd(ymm10, ymm9);				// subtract the kinetic energy density from the total energy density to get the internal energy density 
	ymm10 = _mm256_mul_pd(ymm10, ymm1);				// multiply the internal energy density with gamma-1 to get the pressure 
	// now compute the tensor of the nonlinear convective terms via AVX instructions 
	ymm0  = _mm256_mul_pd(ymm3, ymm6);				// compute rho*u*v  
	ymm1  = _mm256_mul_pd(ymm3, ymm8);				// compute rho*u*w  
	ymm2  = _mm256_mul_pd(ymm5, ymm8);				// compute rho*v*w 
	// add the pressure to the diagonal 
	ymm11 = _mm256_add_pd(ymm11, ymm10);			// compute rho*u*u + p 
	ymm12 = _mm256_add_pd(ymm12, ymm10);			// compute rho*v*v + p 
	ymm13 = _mm256_add_pd(ymm13, ymm10);			// compute rho*w*w + p 
	ymm14 = _mm256_add_pd(ymm9, ymm10);				// compute rhoE + p 

	// store the results in F 
	ymm15 = _mm256_mul_pd(ymm14, ymm4);				// compute u*(rhoE+p) 
	_mm256_store_pd(&F[0][0][0], ymm3);				// rho*u 
	_mm256_store_pd(&F[0][1][0], ymm11);			// rho*u*u + p   
	_mm256_store_pd(&F[0][2][0], ymm0);				// rho*u*v 
	_mm256_store_pd(&F[0][3][0], ymm1);				// rho*u*w 
	_mm256_store_pd(&F[0][4][0], ymm15);			// u*(rhoE+p)  

	ymm15 = _mm256_mul_pd(ymm14, ymm6);				// compute v*(rhoE+p) 
	_mm256_store_pd(&F[1][0][0], ymm5);				// rho*v 
	_mm256_store_pd(&F[1][1][0], ymm0);				// rho*v*u 
	_mm256_store_pd(&F[1][2][0], ymm12);			// rho*v*v + p  
	_mm256_store_pd(&F[1][3][0], ymm2);				// rho*v*w 
	_mm256_store_pd(&F[1][4][0], ymm15);			// v*(rhoE+p)  

	ymm15 = _mm256_mul_pd(ymm14, ymm8);				// compute w*(rhoE+p) 
	_mm256_store_pd(&F[2][0][0], ymm7);				// rho*w 
	_mm256_store_pd(&F[2][1][0], ymm1);				// rho*w*u 
	_mm256_store_pd(&F[2][2][0], ymm2);				// rho*w*v   
	_mm256_store_pd(&F[2][3][0], ymm13);			// rho*w*w + p 
	_mm256_store_pd(&F[2][4][0], ymm15);			// v*(rhoE+p)  

#endif 

}


extern "C" void __cdecl func_c2p_rmhd1_avx(double x[VECTORLENGTH], double f[VECTORLENGTH], 
											   double df[VECTORLENGTH], double gam[VECTORLENGTH],
	                                           double dd[VECTORLENGTH], double e[VECTORLENGTH],
											   double s2[VECTORLENGTH], double b2[VECTORLENGTH],
											   double sb2[VECTORLENGTH], double w[VECTORLENGTH] ) 
{
#ifdef AVX512

	__assume_aligned(x,   64);
	__assume_aligned(f,   64);
	__assume_aligned(df,  64);
	__assume_aligned(gam, 64);
	__assume_aligned(dd,  64);
	__assume_aligned(e,   64);
	__assume_aligned(s2,  64);
	__assume_aligned(b2,  64);
	__assume_aligned(sb2, 64);
	__assume_aligned(w,   64);

	__m512d zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7;
	__m512d zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15;
	int iNewton;

	zmm0 = _mm512_load_pd(&x[0]);       // v2 
										//
	zmm1 = _mm512_set1_pd(1.0);
	zmm2 = _mm512_sub_pd(zmm1, zmm0);   // 1 - v2 
	zmm3 = _mm512_add_pd(zmm1, zmm0);   // 1 + v2 
	zmm4 = _mm512_sqrt_pd(zmm2);
	zmm15 = _mm512_load_pd(&dd[0]);
	zmm5 = _mm512_mul_pd(zmm15, zmm4);   // rho(:) = d(:)*sqrt(1. - v2(:))
	zmm6 = _mm512_load_pd(&gam[0]);
	zmm7 = _mm512_set1_pd(-1.0);
	zmm8 = _mm512_mul_pd(zmm7, zmm6);   // -gam 
	zmm9 = _mm512_set1_pd(1.0);
	zmm9 = _mm512_fmadd_pd(zmm8, zmm2, zmm9); // c3(:) = 1.0 - gam * (1. - v2(:)) 
	zmm10 = _mm512_load_pd(&e[0]);
	zmm10 = _mm512_fmsub_pd(zmm6, zmm5, zmm10);
	zmm11 = _mm512_load_pd(&b2[0]);
	zmm12 = _mm512_set1_pd(0.5);
	zmm13 = _mm512_mul_pd(zmm11, zmm12);
	zmm10 = _mm512_fmadd_pd(zmm13, zmm3, zmm10); // c2(:) = gam(:)*rho(:) + 0.5*b2(:)*(1.0 + v2(:)) - e(:) 
	zmm12 = _mm512_set1_pd(-0.5);
	zmm11 = _mm512_load_pd(&sb2[0]);
	zmm13 = _mm512_mul_pd(zmm12, zmm11);   // c0(:) = -0.5*sb2(:) 
	zmm14 = _mm512_div_pd(zmm7, zmm9);     // -1/c3 
	zmm15 = _mm512_mul_pd(zmm10, zmm14);   // -c2/c3 
	zmm6 = _mm512_mul_pd(zmm13, zmm14);   // -c0/c3 
	zmm12 = _mm512_set1_pd(1.0 / 3.0);
	zmm8 = _mm512_pow_pd(zmm6, zmm12);     // (-c0/c3)**1./3. 
	zmm12 = _mm512_max_pd(zmm15, zmm8);     // w 	
	zmm11 = _mm512_set1_pd(2.0 / 3.0);

	zmm2 = zmm10;                             // save c2 
	zmm3 = zmm13;                             // save c0 
	zmm8 = _mm512_mul_pd(zmm11, zmm10);
	zmm8 = _mm512_fmadd_pd(zmm9, zmm12, zmm8);  // denominator 
	zmm10 = _mm512_fmadd_pd(zmm9, zmm12, zmm10);
	zmm15 = _mm512_mul_pd(zmm12, zmm12);      // w**2
	zmm13 = _mm512_fmadd_pd(zmm15, zmm10, zmm13); // numerator 
	zmm5 = _mm512_set1_pd(3.0);
	zmm6 = _mm512_mul_pd(zmm5, zmm12);        // 3*w 
	zmm4 = _mm512_mul_pd(zmm8, zmm6);
	zmm15 = _mm512_div_pd(zmm13, zmm4);       // dw 
	zmm8 = _mm512_set1_pd(-1.0);
	zmm12 = _mm512_fmadd_pd(zmm8, zmm15, zmm12); // w = w - dw 

	zmm10 = zmm2;                               // get c2 
	zmm13 = zmm3;                             // get c0 
	zmm8 = _mm512_mul_pd(zmm11, zmm10);
	zmm8 = _mm512_fmadd_pd(zmm9, zmm12, zmm8);  // denominator 
	zmm10 = _mm512_fmadd_pd(zmm9, zmm12, zmm10);
	zmm15 = _mm512_mul_pd(zmm12, zmm12);      // w**2
	zmm13 = _mm512_fmadd_pd(zmm15, zmm10, zmm13); // numerator 
	zmm5 = _mm512_set1_pd(3.0);
	zmm6 = _mm512_mul_pd(zmm5, zmm12);        // 3*w 
	zmm4 = _mm512_mul_pd(zmm8, zmm6);
	zmm15 = _mm512_div_pd(zmm13, zmm4);       // dw 
	zmm8 = _mm512_set1_pd(-1.0);
	zmm12 = _mm512_fmadd_pd(zmm8, zmm15, zmm12); // w = w - dw 

	zmm10 = zmm2;                               // get c2 
	zmm13 = zmm3;                             // get c0 
	zmm8 = _mm512_mul_pd(zmm11, zmm10);
	zmm8 = _mm512_fmadd_pd(zmm9, zmm12, zmm8);  // denominator 
	zmm10 = _mm512_fmadd_pd(zmm9, zmm12, zmm10);
	zmm15 = _mm512_mul_pd(zmm12, zmm12);      // w**2
	zmm13 = _mm512_fmadd_pd(zmm15, zmm10, zmm13); // numerator 
	zmm5 = _mm512_set1_pd(3.0);
	zmm6 = _mm512_mul_pd(zmm5, zmm12);        // 3*w 
	zmm4 = _mm512_mul_pd(zmm8, zmm6);
	zmm15 = _mm512_div_pd(zmm13, zmm4);       // dw 
	zmm8 = _mm512_set1_pd(-1.0);
	zmm12 = _mm512_fmadd_pd(zmm8, zmm15, zmm12); // w = w - dw 

	zmm10 = zmm2;                               // get c2 
	zmm13 = zmm3;                             // get c0 
	zmm8 = _mm512_mul_pd(zmm11, zmm10);
	zmm8 = _mm512_fmadd_pd(zmm9, zmm12, zmm8);  // denominator 
	zmm10 = _mm512_fmadd_pd(zmm9, zmm12, zmm10);
	zmm15 = _mm512_mul_pd(zmm12, zmm12);      // w**2
	zmm13 = _mm512_fmadd_pd(zmm15, zmm10, zmm13); // numerator 
	zmm5 = _mm512_set1_pd(3.0);
	zmm6 = _mm512_mul_pd(zmm5, zmm12);        // 3*w 
	zmm4 = _mm512_mul_pd(zmm8, zmm6);
	zmm15 = _mm512_div_pd(zmm13, zmm4);       // dw 
	zmm8 = _mm512_set1_pd(-1.0);
	zmm12 = _mm512_fmadd_pd(zmm8, zmm15, zmm12); // w = w - dw 

	zmm10 = zmm2;                               // get c2 
	zmm13 = zmm3;                             // get c0 
	zmm8 = _mm512_mul_pd(zmm11, zmm10);
	zmm8 = _mm512_fmadd_pd(zmm9, zmm12, zmm8);  // denominator 
	zmm10 = _mm512_fmadd_pd(zmm9, zmm12, zmm10);
	zmm15 = _mm512_mul_pd(zmm12, zmm12);      // w**2
	zmm13 = _mm512_fmadd_pd(zmm15, zmm10, zmm13); // numerator 
	zmm5 = _mm512_set1_pd(3.0);
	zmm6 = _mm512_mul_pd(zmm5, zmm12);        // 3*w 
	zmm4 = _mm512_mul_pd(zmm8, zmm6);
	zmm15 = _mm512_div_pd(zmm13, zmm4);       // dw 
	zmm8 = _mm512_set1_pd(-1.0);
	zmm12 = _mm512_fmadd_pd(zmm8, zmm15, zmm12); // w = w - dw 

	zmm10 = zmm2;                               // get c2 
	zmm13 = zmm3;                             // get c0 
	zmm8 = _mm512_mul_pd(zmm11, zmm10);
	zmm8 = _mm512_fmadd_pd(zmm9, zmm12, zmm8);  // denominator 
	zmm10 = _mm512_fmadd_pd(zmm9, zmm12, zmm10);
	zmm15 = _mm512_mul_pd(zmm12, zmm12);      // w**2
	zmm13 = _mm512_fmadd_pd(zmm15, zmm10, zmm13); // numerator 
	zmm5 = _mm512_set1_pd(3.0);
	zmm6 = _mm512_mul_pd(zmm5, zmm12);        // 3*w 
	zmm4 = _mm512_mul_pd(zmm8, zmm6);
	zmm15 = _mm512_div_pd(zmm13, zmm4);       // dw 
	zmm8 = _mm512_set1_pd(-1.0);
	zmm12 = _mm512_fmadd_pd(zmm8, zmm15, zmm12); // w = w - dw 

	_mm512_store_pd(&w[0], zmm12);

	zmm15 = zmm2;                       // c2 
	zmm1 = _mm512_set1_pd(1.0);
	zmm2 = _mm512_sub_pd(zmm1, zmm0);   // 1 - v2 
	zmm4 = _mm512_sqrt_pd(zmm2);
	zmm8 = _mm512_load_pd(&dd[0]);
	zmm5 = _mm512_mul_pd(zmm8, zmm4);   // rho(:) = d(:)*sqrt(1. - v2(:))
	zmm6 = _mm512_load_pd(&gam[0]);
	zmm8 = _mm512_set1_pd(-1.0);
	zmm11 = _mm512_mul_pd(zmm8, zmm6);  // -gam 
	zmm7 = _mm512_div_pd(zmm11, zmm2);  // -gam/(1-v2) 
	zmm8 = _mm512_load_pd(&b2[0]);

	zmm2 = zmm8;   // b2 
	zmm8 = _mm512_fmadd_pd(zmm7, zmm5, zmm8);
	zmm7 = _mm512_set1_pd(0.5);
	zmm10 = _mm512_mul_pd(zmm7, zmm8);   // dc2 
	zmm14 = _mm512_set1_pd(2. / 3.);
	zmm13 = _mm512_mul_pd(zmm14, zmm15); // 2./3*c2 
	zmm13 = _mm512_fmadd_pd(zmm9, zmm12, zmm13);    // denom  
	zmm7 = _mm512_set1_pd(-3.0);
	zmm14 = _mm512_mul_pd(zmm7, zmm13); // -3*denom 
	zmm10 = _mm512_fmadd_pd(zmm12, zmm6, zmm10); // nominator 
	zmm15 = _mm512_div_pd(zmm10, zmm14);   // dlogw 
	zmm3 = _mm512_add_pd(zmm12, zmm2);     // wb 
	zmm4 = _mm512_mul_pd(zmm3, zmm3);      // wb**2 
	zmm13 = _mm512_mul_pd(zmm12, zmm12);   // w**2
	zmm8 = _mm512_load_pd(&sb2[0]);        // sb2 
	zmm9 = _mm512_div_pd(zmm8, zmm13);     // vb2 
	zmm7 = _mm512_set1_pd(-2.0);
	zmm2 = _mm512_fmsub_pd(zmm7, zmm12, zmm2);  // -2*w - b2 
	zmm8 = _mm512_load_pd(&s2[0]);
	zmm8 = _mm512_fmsub_pd(zmm2, zmm9, zmm8);
	zmm8 = _mm512_fmadd_pd(zmm4, zmm0, zmm8);   // f 
	_mm512_store_pd(&f[0], zmm8);


	zmm7 = _mm512_set1_pd(2.0);
	zmm10 = _mm512_mul_pd(zmm7, zmm15);          // 2*dlogw 
	zmm1 = zmm3;                                // wb 
	zmm9 = _mm512_fmadd_pd(zmm12, zmm0, zmm9);
	zmm3 = _mm512_fmadd_pd(zmm10, zmm9, zmm3);
	zmm2 = _mm512_mul_pd(zmm3, zmm1);           // df 
	_mm512_store_pd(&df[0], zmm2);

	//zmm7 = _mm256_set1_pd(-1.0);
	//zmm9 = _mm256_div_pd(zmm7, zmm2);
	//zmm0 = _mm256_fmadd_pd(zmm9, zmm8, zmm0);  // v2(:) = v2(:) - F(:) / DF(:)

	// _mm256_store_pd(&x[0], zmm0); 



	// asm(nop); 

	// f(:) = wb(:)**2 * v2(:) + (-2.0 * w(:) - b2(:)) * vb2(:) - s2(:)
	// df(:) = wb(:) * (wb(:) + 2.0 * dlogw(:) * (w(:) * v2(:) + vb2(:)))


#else 
	__assume_aligned(x, 32);
	__assume_aligned(f, 32);
	__assume_aligned(df, 32);
	__assume_aligned(gam, 32);
	__assume_aligned(dd, 32);
	__assume_aligned(e, 32);
	__assume_aligned(s2, 32);
	__assume_aligned(b2, 32);
	__assume_aligned(sb2, 32);
	__assume_aligned(w, 32);

	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	__m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;
	int iNewton; 

	ymm0 = _mm256_load_pd(&x[0]);       // v2 
	//
		ymm1 = _mm256_set1_pd(1.0);
		ymm2 = _mm256_sub_pd(ymm1, ymm0);   // 1 - v2 
		ymm3 = _mm256_add_pd(ymm1, ymm0);   // 1 + v2 
		ymm4 = _mm256_sqrt_pd(ymm2);
		ymm15 = _mm256_load_pd(&dd[0]);
		ymm5 = _mm256_mul_pd(ymm15, ymm4);   // rho(:) = d(:)*sqrt(1. - v2(:))
		ymm6 = _mm256_load_pd(&gam[0]);
		ymm7 = _mm256_set1_pd(-1.0);
		ymm8 = _mm256_mul_pd(ymm7, ymm6);   // -gam 
		ymm9 = _mm256_set1_pd(1.0);
		ymm9 = _mm256_fmadd_pd(ymm8, ymm2, ymm9); // c3(:) = 1.0 - gam * (1. - v2(:)) 
		ymm10 = _mm256_load_pd(&e[0]);
		ymm10 = _mm256_fmsub_pd(ymm6, ymm5, ymm10);
		ymm11 = _mm256_load_pd(&b2[0]);
		ymm12 = _mm256_set1_pd(0.5);
		ymm13 = _mm256_mul_pd(ymm11, ymm12);
		ymm10 = _mm256_fmadd_pd(ymm13, ymm3, ymm10); // c2(:) = gam(:)*rho(:) + 0.5*b2(:)*(1.0 + v2(:)) - e(:) 
		ymm12 = _mm256_set1_pd(-0.5);
		ymm11 = _mm256_load_pd(&sb2[0]);
		ymm13 = _mm256_mul_pd(ymm12, ymm11);   // c0(:) = -0.5*sb2(:) 
		ymm14 = _mm256_div_pd(ymm7, ymm9);     // -1/c3 
		ymm15 = _mm256_mul_pd(ymm10, ymm14);   // -c2/c3 
		ymm6 = _mm256_mul_pd(ymm13, ymm14);   // -c0/c3 
		ymm12 = _mm256_set1_pd(1.0 / 3.0);
		ymm8 = _mm256_pow_pd(ymm6, ymm12);     // (-c0/c3)**1./3. 
		ymm12 = _mm256_max_pd(ymm15, ymm8);     // w 	
		ymm11 = _mm256_set1_pd(2.0 / 3.0);

		ymm2 = ymm10;                             // save c2 
		ymm3 = ymm13;                             // save c0 
		ymm8 = _mm256_mul_pd(ymm11, ymm10);
		ymm8 = _mm256_fmadd_pd(ymm9, ymm12, ymm8);  // denominator 
		ymm10 = _mm256_fmadd_pd(ymm9, ymm12, ymm10);
		ymm15 = _mm256_mul_pd(ymm12, ymm12);      // w**2
		ymm13 = _mm256_fmadd_pd(ymm15, ymm10, ymm13); // numerator 
		ymm5 = _mm256_set1_pd(3.0);
		ymm6 = _mm256_mul_pd(ymm5, ymm12);        // 3*w 
		ymm4 = _mm256_mul_pd(ymm8, ymm6);
		ymm15 = _mm256_div_pd(ymm13, ymm4);       // dw 
		ymm8 = _mm256_set1_pd(-1.0);
		ymm12 = _mm256_fmadd_pd(ymm8, ymm15, ymm12); // w = w - dw 

		ymm10 = ymm2;                               // get c2 
		ymm13 = ymm3;                             // get c0 
		ymm8 = _mm256_mul_pd(ymm11, ymm10);
		ymm8 = _mm256_fmadd_pd(ymm9, ymm12, ymm8);  // denominator 
		ymm10 = _mm256_fmadd_pd(ymm9, ymm12, ymm10);
		ymm15 = _mm256_mul_pd(ymm12, ymm12);      // w**2
		ymm13 = _mm256_fmadd_pd(ymm15, ymm10, ymm13); // numerator 
		ymm5 = _mm256_set1_pd(3.0);
		ymm6 = _mm256_mul_pd(ymm5, ymm12);        // 3*w 
		ymm4 = _mm256_mul_pd(ymm8, ymm6);
		ymm15 = _mm256_div_pd(ymm13, ymm4);       // dw 
		ymm8 = _mm256_set1_pd(-1.0);
		ymm12 = _mm256_fmadd_pd(ymm8, ymm15, ymm12); // w = w - dw 

		ymm10 = ymm2;                               // get c2 
		ymm13 = ymm3;                             // get c0 
		ymm8 = _mm256_mul_pd(ymm11, ymm10);
		ymm8 = _mm256_fmadd_pd(ymm9, ymm12, ymm8);  // denominator 
		ymm10 = _mm256_fmadd_pd(ymm9, ymm12, ymm10);
		ymm15 = _mm256_mul_pd(ymm12, ymm12);      // w**2
		ymm13 = _mm256_fmadd_pd(ymm15, ymm10, ymm13); // numerator 
		ymm5 = _mm256_set1_pd(3.0);
		ymm6 = _mm256_mul_pd(ymm5, ymm12);        // 3*w 
		ymm4 = _mm256_mul_pd(ymm8, ymm6);
		ymm15 = _mm256_div_pd(ymm13, ymm4);       // dw 
		ymm8 = _mm256_set1_pd(-1.0);
		ymm12 = _mm256_fmadd_pd(ymm8, ymm15, ymm12); // w = w - dw 

		ymm10 = ymm2;                               // get c2 
		ymm13 = ymm3;                             // get c0 
		ymm8 = _mm256_mul_pd(ymm11, ymm10);
		ymm8 = _mm256_fmadd_pd(ymm9, ymm12, ymm8);  // denominator 
		ymm10 = _mm256_fmadd_pd(ymm9, ymm12, ymm10);
		ymm15 = _mm256_mul_pd(ymm12, ymm12);      // w**2
		ymm13 = _mm256_fmadd_pd(ymm15, ymm10, ymm13); // numerator 
		ymm5 = _mm256_set1_pd(3.0);
		ymm6 = _mm256_mul_pd(ymm5, ymm12);        // 3*w 
		ymm4 = _mm256_mul_pd(ymm8, ymm6);
		ymm15 = _mm256_div_pd(ymm13, ymm4);       // dw 
		ymm8 = _mm256_set1_pd(-1.0);
		ymm12 = _mm256_fmadd_pd(ymm8, ymm15, ymm12); // w = w - dw 

		ymm10 = ymm2;                               // get c2 
		ymm13 = ymm3;                             // get c0 
		ymm8 = _mm256_mul_pd(ymm11, ymm10);
		ymm8 = _mm256_fmadd_pd(ymm9, ymm12, ymm8);  // denominator 
		ymm10 = _mm256_fmadd_pd(ymm9, ymm12, ymm10);
		ymm15 = _mm256_mul_pd(ymm12, ymm12);      // w**2
		ymm13 = _mm256_fmadd_pd(ymm15, ymm10, ymm13); // numerator 
		ymm5 = _mm256_set1_pd(3.0);
		ymm6 = _mm256_mul_pd(ymm5, ymm12);        // 3*w 
		ymm4 = _mm256_mul_pd(ymm8, ymm6);
		ymm15 = _mm256_div_pd(ymm13, ymm4);       // dw 
		ymm8 = _mm256_set1_pd(-1.0);
		ymm12 = _mm256_fmadd_pd(ymm8, ymm15, ymm12); // w = w - dw 

		ymm10 = ymm2;                               // get c2 
		ymm13 = ymm3;                             // get c0 
		ymm8 = _mm256_mul_pd(ymm11, ymm10);
		ymm8 = _mm256_fmadd_pd(ymm9, ymm12, ymm8);  // denominator 
		ymm10 = _mm256_fmadd_pd(ymm9, ymm12, ymm10);
		ymm15 = _mm256_mul_pd(ymm12, ymm12);      // w**2
		ymm13 = _mm256_fmadd_pd(ymm15, ymm10, ymm13); // numerator 
		ymm5 = _mm256_set1_pd(3.0);
		ymm6 = _mm256_mul_pd(ymm5, ymm12);        // 3*w 
		ymm4 = _mm256_mul_pd(ymm8, ymm6);
		ymm15 = _mm256_div_pd(ymm13, ymm4);       // dw 
		ymm8 = _mm256_set1_pd(-1.0);
		ymm12 = _mm256_fmadd_pd(ymm8, ymm15, ymm12); // w = w - dw 

		_mm256_store_pd(&w[0], ymm12);

		ymm15 = ymm2;                       // c2 
		ymm1 = _mm256_set1_pd(1.0);
		ymm2 = _mm256_sub_pd(ymm1, ymm0);   // 1 - v2 
		ymm4 = _mm256_sqrt_pd(ymm2);
		ymm8 = _mm256_load_pd(&dd[0]);
		ymm5 = _mm256_mul_pd(ymm8, ymm4);   // rho(:) = d(:)*sqrt(1. - v2(:))
		ymm6 = _mm256_load_pd(&gam[0]);
		ymm8 = _mm256_set1_pd(-1.0);
		ymm11 = _mm256_mul_pd(ymm8, ymm6);  // -gam 
		ymm7 = _mm256_div_pd(ymm11, ymm2);  // -gam/(1-v2) 
		ymm8 = _mm256_load_pd(&b2[0]);

		ymm2 = ymm8;   // b2 
		ymm8 = _mm256_fmadd_pd(ymm7, ymm5, ymm8);
		ymm7 = _mm256_set1_pd(0.5);
		ymm10 = _mm256_mul_pd(ymm7, ymm8);   // dc2 
		ymm14 = _mm256_set1_pd(2. / 3.);
		ymm13 = _mm256_mul_pd(ymm14, ymm15); // 2./3*c2 
		ymm13 = _mm256_fmadd_pd(ymm9, ymm12, ymm13);    // denom  
		ymm7 = _mm256_set1_pd(-3.0);
		ymm14 = _mm256_mul_pd(ymm7, ymm13); // -3*denom 
		ymm10 = _mm256_fmadd_pd(ymm12, ymm6, ymm10); // nominator 
		ymm15 = _mm256_div_pd(ymm10, ymm14);   // dlogw 
		ymm3 = _mm256_add_pd(ymm12, ymm2);     // wb 
		ymm4 = _mm256_mul_pd(ymm3, ymm3);      // wb**2 
		ymm13 = _mm256_mul_pd(ymm12, ymm12);   // w**2
		ymm8 = _mm256_load_pd(&sb2[0]);        // sb2 
		ymm9 = _mm256_div_pd(ymm8, ymm13);     // vb2 
		ymm7 = _mm256_set1_pd(-2.0);
		ymm2 = _mm256_fmsub_pd(ymm7, ymm12, ymm2);  // -2*w - b2 
		ymm8 = _mm256_load_pd(&s2[0]);
		ymm8 = _mm256_fmsub_pd(ymm2, ymm9, ymm8);
		ymm8 = _mm256_fmadd_pd(ymm4, ymm0, ymm8);   // f 
		_mm256_store_pd(&f[0], ymm8); 


		ymm7 = _mm256_set1_pd(2.0);
		ymm10 = _mm256_mul_pd(ymm7, ymm15);          // 2*dlogw 
		ymm1 = ymm3;                                // wb 
		ymm9 = _mm256_fmadd_pd(ymm12, ymm0, ymm9);
		ymm3 = _mm256_fmadd_pd(ymm10, ymm9, ymm3);
		ymm2 = _mm256_mul_pd(ymm3, ymm1);           // df 
		_mm256_store_pd(&df[0], ymm2);

		//ymm7 = _mm256_set1_pd(-1.0);
		//ymm9 = _mm256_div_pd(ymm7, ymm2);
		//ymm0 = _mm256_fmadd_pd(ymm9, ymm8, ymm0);  // v2(:) = v2(:) - F(:) / DF(:)

		// _mm256_store_pd(&x[0], ymm0); 
	    
     

	// asm(nop); 

	// f(:) = wb(:)**2 * v2(:) + (-2.0 * w(:) - b2(:)) * vb2(:) - s2(:)
	// df(:) = wb(:) * (wb(:) + 2.0 * dlogw(:) * (w(:) * v2(:) + vb2(:)))


#endif 

}


extern "C" void __cdecl ccz4testc(double dChristoffel_tildeNCP[3][3][3][3][VECTORLENGTH],    
								  double Z[3][VECTORLENGTH], double Zup[3][VECTORLENGTH],
								  double Gtilde[3][VECTORLENGTH], 
								  double Christoffel[3][3][3][VECTORLENGTH], 
								  double Christoffel_tilde[3][3][3][VECTORLENGTH], 
								  double dgup[3][3][3][VECTORLENGTH], 
								  double Q[nVar][VECTORLENGTH], 
								  double Qx[nVar][VECTORLENGTH],
								  double Qy[nVar][VECTORLENGTH],
								  double Qz[nVar][VECTORLENGTH], double ltime[1])
{
#ifdef AVX512 

	__assume_aligned(dChristoffel_tildeNCP, 64);
	__assume_aligned(Z, 64);
	__assume_aligned(Zup, 64);
	__assume_aligned(Gtilde, 64);
	__assume_aligned(Christoffel, 64); 
	__assume_aligned(Christoffel_tilde, 64);
	__assume_aligned(dgup, 64);
	__assume_aligned(Q, 64);
	__assume_aligned(Qx, 64);
	__assume_aligned(Qy, 64);
	__assume_aligned(Qz, 64);

	int i, j, k, v, nAdd, nMul;

#ifdef COUNT_FLOPS 
	nAdd = 0;
	nMul = 0; 
#endif 

	// We have 20 free registers. The other 12 are used for permanently storing the covariant and contravariant metric 
	__m512d zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7, zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15, zmm16, zmm17, zmm18, zmm19; 
	register __m512d gxx, gxy, gxz, gyy, gyz, gzz;
	register __m512d igxx, igxy, igxz, igyy, igyz, igzz;
	__m128d xmm;  

	_declspec(align(64)) double g_cov[3][3][VECTORLENGTH], g_contr[3][3][VECTORLENGTH], det[VECTORLENGTH], idet[VECTORLENGTH], DD[3][3][3][VECTORLENGTH];
	_declspec(align(64)) double dDDsymm[3][3][3][3][VECTORLENGTH];
	//_declspec(align(64)) double dChristoffel_tildeNCP[3][3][3][3][VECTORLENGTH];

	double s1[1], s2[1], s3[1]; 
	s1[0] = -2.0; 
	s2[0] = -1.0; 
	s3[0] =  0.5; 

	// shuffle load operations with computations 
	// first sub-determinant 
	gyz = _mm512_load_pd(&Q[4][0]); 
	zmm0 = _mm512_mul_pd(gyz, gyz);
	gyy = _mm512_load_pd(&Q[3][0]);
	gzz = _mm512_load_pd(&Q[5][0]);
	zmm0 = _mm512_fmsub_pd(gyy, gzz, zmm0);
	// second sub-determinant (order changed since the factor in front is negative) 
	gxy = _mm512_load_pd(&Q[1][0]);
	zmm1 = _mm512_mul_pd(gxy, gzz);
	gxz = _mm512_load_pd(&Q[2][0]);
	zmm1 = _mm512_fmsub_pd(gxz, gyz, zmm1);
	// third sub-determinant 
	zmm2 = _mm512_mul_pd(gxz, gyy);
	gxx = _mm512_load_pd(&Q[0][0]);
	zmm2 = _mm512_fmsub_pd(gxy, gyz, zmm2);
	// Now assemble the final determinant 
	zmm18 = _mm512_mul_pd(gxx, zmm0);
	zmm18 = _mm512_fmadd_pd(gxy, zmm1, zmm18);
	zmm18 = _mm512_fmadd_pd(gxz, zmm2, zmm18);
	// Compute 1/det 
	zmm10 = _mm512_set1_pd(1.0);	    		// set unit 
	zmm19 = _mm512_div_pd(zmm10, zmm18);		// compute 1/det 

#ifdef COUNT_FLOPS  
	nMul = nMul + VECTORLENGTH * (4 + 5);
	nAdd = nAdd + VECTORLENGTH * 5;
#endif 

	// now compute the independent components of the contravariant metric tensor 
	igxx = _mm512_mul_pd(zmm0, zmm19); 
	igxy = _mm512_mul_pd(zmm1, zmm19);
	igxz = _mm512_mul_pd(zmm2, zmm19);
	zmm3 = _mm512_mul_pd(gxz, gxz);
	zmm3 = _mm512_fmsub_pd(gxx, gzz, zmm3);
	igyy = _mm512_mul_pd(zmm3, zmm19);
	zmm4 = _mm512_mul_pd(gxx, gyz);
	zmm4 = _mm512_fmsub_pd(gxy, gxz, zmm4);
	igyz = _mm512_mul_pd(zmm4, zmm19);
	zmm5 = _mm512_mul_pd(gxy, gxy);
	zmm5 = _mm512_fmsub_pd(gxx, gyy, zmm5);
	igzz = _mm512_mul_pd(zmm5, zmm19);

#ifdef COUNT_FLOPS  
	nMul = nMul + VECTORLENGTH * 12;
	nAdd = nAdd + VECTORLENGTH * 3;
#endif 

	// be careful, here the order of the indices of D_kij is different than in Fortran. The left-most index corresponds to k, since it should be the slowest one. 

	// load scalar ( -1.0 ) and multiply with igxx, igxy and igxz   
	xmm = _mm_load1_pd(&s2[0]);
	zmm0  = _mm512_broadcastsd_pd(xmm);
	zmm17 = _mm512_mul_pd(igxx, zmm0);
	zmm18 = _mm512_mul_pd(igxy, zmm0);
	zmm19 = _mm512_mul_pd(igxz, zmm0);
	// load scalar ( -2.0 )  
	xmm = _mm_load1_pd(&s1[0]);
	zmm16 = _mm512_broadcastsd_pd(xmm);
#ifdef COUNT_FLOPS  
	nMul = nMul + VECTORLENGTH * 3 * 3;   // this is counted for all three slices. 
#endif 

	// load first slice of D_kij  (D_1ij) 
	zmm10 = _mm512_load_pd(&Q[35][0]);
	zmm11 = _mm512_load_pd(&Q[36][0]);
	zmm12 = _mm512_load_pd(&Q[37][0]);
	zmm13 = _mm512_load_pd(&Q[38][0]);
	zmm14 = _mm512_load_pd(&Q[39][0]);
	zmm15 = _mm512_load_pd(&Q[40][0]);

	// for Christoffel tilde, multiply with g_contr 
	zmm0 = _mm512_mul_pd(igxx, zmm10); 
	zmm0 = _mm512_fmadd_pd(igxy, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igxy, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igxz, zmm12, zmm0);
	zmm0 = _mm512_fmadd_pd(igxz, zmm12, zmm0);
	zmm1 = _mm512_mul_pd(igxy, zmm13); 
	zmm1 = _mm512_fmadd_pd(igxz, zmm14, zmm1);
	zmm2 = _mm512_mul_pd(igxy, zmm14);
	zmm2 = _mm512_fmadd_pd(igxz, zmm15, zmm2); 
	zmm3 = _mm512_mul_pd(zmm17, zmm13); 
	zmm4 = _mm512_mul_pd(zmm17, zmm14);
	zmm5 = _mm512_mul_pd(zmm17, zmm15);
	_mm512_store_pd(&Christoffel_tilde[0][0][0][0], zmm0); 
	_mm512_store_pd(&Christoffel_tilde[0][0][1][0], zmm1);
	_mm512_store_pd(&Christoffel_tilde[0][0][2][0], zmm2);
	_mm512_store_pd(&Christoffel_tilde[0][1][1][0], zmm3);
	_mm512_store_pd(&Christoffel_tilde[0][1][2][0], zmm4);
	_mm512_store_pd(&Christoffel_tilde[0][2][2][0], zmm5);
	// 
	zmm0 = _mm512_mul_pd(igxy, zmm10);
	zmm0 = _mm512_fmadd_pd(igyy, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igyy, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igyz, zmm12, zmm0);
	zmm0 = _mm512_fmadd_pd(igyz, zmm12, zmm0);
	zmm1 = _mm512_mul_pd(igyy, zmm13);
	zmm1 = _mm512_fmadd_pd(igyz, zmm14, zmm1);
	zmm2 = _mm512_mul_pd(igyy, zmm14);
	zmm2 = _mm512_fmadd_pd(igyz, zmm15, zmm2);
	zmm3 = _mm512_mul_pd(zmm18, zmm13);
	zmm4 = _mm512_mul_pd(zmm18, zmm14);
	zmm5 = _mm512_mul_pd(zmm18, zmm15);
	_mm512_store_pd(&Christoffel_tilde[1][0][0][0], zmm0);
	_mm512_store_pd(&Christoffel_tilde[1][0][1][0], zmm1);
	_mm512_store_pd(&Christoffel_tilde[1][0][2][0], zmm2);
	_mm512_store_pd(&Christoffel_tilde[1][1][1][0], zmm3);
	_mm512_store_pd(&Christoffel_tilde[1][1][2][0], zmm4);
	_mm512_store_pd(&Christoffel_tilde[1][2][2][0], zmm5);
	// 
	zmm0 = _mm512_mul_pd(igxz, zmm10);
	zmm0 = _mm512_fmadd_pd(igyz, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igyz, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igzz, zmm12, zmm0);
	zmm0 = _mm512_fmadd_pd(igzz, zmm12, zmm0);
	zmm1 = _mm512_mul_pd(igyz, zmm13);
	zmm1 = _mm512_fmadd_pd(igzz, zmm14, zmm1);
	zmm2 = _mm512_mul_pd(igyz, zmm14);
	zmm2 = _mm512_fmadd_pd(igzz, zmm15, zmm2);
	zmm3 = _mm512_mul_pd(zmm19, zmm13);
	zmm4 = _mm512_mul_pd(zmm19, zmm14);
	zmm5 = _mm512_mul_pd(zmm19, zmm15);
	_mm512_store_pd(&Christoffel_tilde[2][0][0][0], zmm0);
	_mm512_store_pd(&Christoffel_tilde[2][0][1][0], zmm1);
	_mm512_store_pd(&Christoffel_tilde[2][0][2][0], zmm2);
	_mm512_store_pd(&Christoffel_tilde[2][1][1][0], zmm3);
	_mm512_store_pd(&Christoffel_tilde[2][1][2][0], zmm4);
	_mm512_store_pd(&Christoffel_tilde[2][2][2][0], zmm5);

	// multiply with the scalar for dgup 
	zmm0 = _mm512_mul_pd(zmm10, zmm16);
	zmm1 = _mm512_mul_pd(zmm11, zmm16);
	zmm2 = _mm512_mul_pd(zmm12, zmm16);
	zmm3 = _mm512_mul_pd(zmm13, zmm16);
	zmm4 = _mm512_mul_pd(zmm14, zmm16);
	zmm5 = _mm512_mul_pd(zmm15, zmm16);

	// save into DD 
	//_mm512_store_pd(&DD[0][1][0][0], zmm11);
	//_mm512_store_pd(&DD[0][2][0][0], zmm12);
	//_mm512_store_pd(&DD[0][2][1][0], zmm14); 

	// raise first index of the first block 
	zmm6 = _mm512_mul_pd(igxx, zmm0); 
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	_mm512_store_pd(&DD[0][0][0][0], zmm10);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm7 = _mm512_mul_pd(igxx, zmm1);
	_mm512_store_pd(&DD[0][0][1][0], zmm11);
	zmm7 = _mm512_fmadd_pd(igxy, zmm3, zmm7);
	_mm512_store_pd(&DD[0][0][2][0], zmm12);
	zmm7 = _mm512_fmadd_pd(igxz, zmm4, zmm7);
	_mm512_store_pd(&DD[0][1][1][0], zmm13);
	zmm8 = _mm512_mul_pd(igxx, zmm2);
	_mm512_store_pd(&DD[0][1][2][0], zmm14);
	zmm8 = _mm512_fmadd_pd(igxy, zmm4, zmm8);
	_mm512_store_pd(&DD[0][2][2][0], zmm15);
	zmm8 = _mm512_fmadd_pd(igxz, zmm5, zmm8);

	zmm9 = _mm512_mul_pd(igxy, zmm0);
	zmm9 = _mm512_fmadd_pd(igyy, zmm1, zmm9);
	zmm9 = _mm512_fmadd_pd(igyz, zmm2, zmm9);
	zmm10 = _mm512_mul_pd(igxy, zmm1);
	zmm10 = _mm512_fmadd_pd(igyy, zmm3, zmm10);
	zmm10 = _mm512_fmadd_pd(igyz, zmm4, zmm10);
	zmm11 = _mm512_mul_pd(igxy, zmm2);
	zmm11 = _mm512_fmadd_pd(igyy, zmm4, zmm11);
	zmm11 = _mm512_fmadd_pd(igyz, zmm5, zmm11); 

	zmm12 = _mm512_mul_pd(igxz, zmm0);
	zmm12 = _mm512_fmadd_pd(igyz, zmm1, zmm12);
	zmm12 = _mm512_fmadd_pd(igzz, zmm2, zmm12);
	zmm13 = _mm512_mul_pd(igxz, zmm1);
	zmm13 = _mm512_fmadd_pd(igyz, zmm3, zmm13);
	zmm13 = _mm512_fmadd_pd(igzz, zmm4, zmm13);
	zmm14 = _mm512_mul_pd(igxz, zmm2);
	zmm14 = _mm512_fmadd_pd(igyz, zmm4, zmm14);
	zmm14 = _mm512_fmadd_pd(igzz, zmm5, zmm14); 

	// raise the second index. The resulting object (dgup) is symmetric in the ij indices. 

	zmm0 = _mm512_mul_pd(igxx, zmm6);
	zmm0 = _mm512_fmadd_pd(igxy, zmm7, zmm0);
	zmm0 = _mm512_fmadd_pd(igxz, zmm8, zmm0);
	_mm512_store_pd(&dgup[0][0][0][0], zmm0);
	zmm1 = _mm512_mul_pd(igxy, zmm6);
	zmm1 = _mm512_fmadd_pd(igyy, zmm7, zmm1);
	zmm1 = _mm512_fmadd_pd(igyz, zmm8, zmm1);
	_mm512_store_pd(&dgup[0][0][1][0], zmm1);
	zmm2 = _mm512_mul_pd(igxz, zmm6);
	zmm2 = _mm512_fmadd_pd(igyz, zmm7, zmm2);
	zmm2 = _mm512_fmadd_pd(igzz, zmm8, zmm2);
	_mm512_store_pd(&dgup[0][0][2][0], zmm2);
	zmm3 = _mm512_mul_pd(igxy, zmm9);
	zmm3 = _mm512_fmadd_pd(igyy, zmm10, zmm3);
	zmm3 = _mm512_fmadd_pd(igyz, zmm11, zmm3);
	_mm512_store_pd(&dgup[0][1][1][0], zmm3);
	zmm4 = _mm512_mul_pd(igxz, zmm9);
	zmm4 = _mm512_fmadd_pd(igyz, zmm10, zmm4);
	zmm4 = _mm512_fmadd_pd(igzz, zmm11, zmm4);
	_mm512_store_pd(&dgup[0][1][2][0], zmm4);
	zmm5 = _mm512_mul_pd(igxz, zmm12);
	zmm5 = _mm512_fmadd_pd(igyz, zmm13, zmm5);
	zmm5 = _mm512_fmadd_pd(igzz, zmm14, zmm5); 

	// save the result into dgup 
	//_mm512_store_pd(&dgup[0][1][0][0], zmm1);
	//_mm512_store_pd(&dgup[0][2][0][0], zmm2);
	//_mm512_store_pd(&dgup[0][2][1][0], zmm4);
	_mm512_store_pd(&dgup[0][2][2][0], zmm5);


	// load scalar ( -1.0 ) and multiply with igxx, igxy and igxz   
	xmm = _mm_load1_pd(&s2[0]);
	zmm0 = _mm512_broadcastsd_pd(xmm);
	zmm17 = _mm512_mul_pd(igxy, zmm0);
	zmm18 = _mm512_mul_pd(igyy, zmm0);
	zmm19 = _mm512_mul_pd(igyz, zmm0);
	// load scalar ( -2.0 )  
	xmm = _mm_load1_pd(&s1[0]);
	zmm16 = _mm512_broadcastsd_pd(xmm);

	// load second slice of D_kij 
	zmm10 = _mm512_load_pd(&Q[41][0]);
	zmm11 = _mm512_load_pd(&Q[42][0]);
	zmm12 = _mm512_load_pd(&Q[43][0]);
	zmm13 = _mm512_load_pd(&Q[44][0]);
	zmm14 = _mm512_load_pd(&Q[45][0]);
	zmm15 = _mm512_load_pd(&Q[46][0]);

	// for Christoffel tilde, load the result of the previous slice, multiply DD properly with g_contr and add contribution 
	zmm0 = _mm512_load_pd(&Christoffel_tilde[0][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm17, zmm10, zmm0); 
	zmm1 = _mm512_load_pd(&Christoffel_tilde[0][0][1][0]);
	zmm1 = _mm512_fmadd_pd(igxx, zmm10, zmm1); 
	zmm1 = _mm512_fmadd_pd(igxz, zmm12, zmm1);
	zmm2 = _mm512_load_pd(&Christoffel_tilde[0][0][2][0]); 
	zmm2 = _mm512_fmadd_pd(zmm17, zmm12, zmm2); 
	zmm3 = _mm512_load_pd(&Christoffel_tilde[0][1][1][0]);
	zmm3 = _mm512_fmadd_pd(igxx, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igxx, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igxy, zmm13, zmm3);
	zmm3 = _mm512_fmadd_pd(igxz, zmm14, zmm3);
	zmm3 = _mm512_fmadd_pd(igxz, zmm14, zmm3);
	zmm4 = _mm512_load_pd(&Christoffel_tilde[0][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxx, zmm12, zmm4);
	zmm4 = _mm512_fmadd_pd(igxz, zmm15, zmm4);
	zmm5 = _mm512_load_pd(&Christoffel_tilde[0][2][2][0]);
	zmm5 = _mm512_fmadd_pd(zmm17, zmm15, zmm5); 
	_mm512_store_pd(&Christoffel_tilde[0][0][0][0], zmm0);
	_mm512_store_pd(&Christoffel_tilde[0][0][1][0], zmm1);
	_mm512_store_pd(&Christoffel_tilde[0][0][2][0], zmm2);
	_mm512_store_pd(&Christoffel_tilde[0][1][1][0], zmm3);
	_mm512_store_pd(&Christoffel_tilde[0][1][2][0], zmm4);
	_mm512_store_pd(&Christoffel_tilde[0][2][2][0], zmm5);
	// 
	zmm0 = _mm512_load_pd(&Christoffel_tilde[1][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm18, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&Christoffel_tilde[1][0][1][0]);
	zmm1 = _mm512_fmadd_pd(igxy, zmm10, zmm1);
	zmm1 = _mm512_fmadd_pd(igyz, zmm12, zmm1);
	zmm2 = _mm512_load_pd(&Christoffel_tilde[1][0][2][0]);
	zmm2 = _mm512_fmadd_pd(zmm18, zmm12, zmm2);
	zmm3 = _mm512_load_pd(&Christoffel_tilde[1][1][1][0]);
	zmm3 = _mm512_fmadd_pd(igxy, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igxy, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igyy, zmm13, zmm3);
	zmm3 = _mm512_fmadd_pd(igyz, zmm14, zmm3);
	zmm3 = _mm512_fmadd_pd(igyz, zmm14, zmm3);
	zmm4 = _mm512_load_pd(&Christoffel_tilde[1][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxy, zmm12, zmm4);
	zmm4 = _mm512_fmadd_pd(igyz, zmm15, zmm4);
	zmm5 = _mm512_load_pd(&Christoffel_tilde[1][2][2][0]);
	zmm5 = _mm512_fmadd_pd(zmm18, zmm15, zmm5);
	_mm512_store_pd(&Christoffel_tilde[1][0][0][0], zmm0);
	_mm512_store_pd(&Christoffel_tilde[1][0][1][0], zmm1);
	_mm512_store_pd(&Christoffel_tilde[1][0][2][0], zmm2);
	_mm512_store_pd(&Christoffel_tilde[1][1][1][0], zmm3);
	_mm512_store_pd(&Christoffel_tilde[1][1][2][0], zmm4);
	_mm512_store_pd(&Christoffel_tilde[1][2][2][0], zmm5);
	// 
	zmm0 = _mm512_load_pd(&Christoffel_tilde[2][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm19, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&Christoffel_tilde[2][0][1][0]);
	zmm1 = _mm512_fmadd_pd(igxz, zmm10, zmm1);
	zmm1 = _mm512_fmadd_pd(igzz, zmm12, zmm1);
	zmm2 = _mm512_load_pd(&Christoffel_tilde[2][0][2][0]);
	zmm2 = _mm512_fmadd_pd(zmm19, zmm12, zmm2);
	zmm3 = _mm512_load_pd(&Christoffel_tilde[2][1][1][0]);
	zmm3 = _mm512_fmadd_pd(igxz, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igxz, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igyz, zmm13, zmm3);
	zmm3 = _mm512_fmadd_pd(igzz, zmm14, zmm3);
	zmm3 = _mm512_fmadd_pd(igzz, zmm14, zmm3);
	zmm4 = _mm512_load_pd(&Christoffel_tilde[2][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxz, zmm12, zmm4);
	zmm4 = _mm512_fmadd_pd(igzz, zmm15, zmm4);
	zmm5 = _mm512_load_pd(&Christoffel_tilde[2][2][2][0]);
	zmm5 = _mm512_fmadd_pd(zmm19, zmm15, zmm5);
	_mm512_store_pd(&Christoffel_tilde[2][0][0][0], zmm0);
	_mm512_store_pd(&Christoffel_tilde[2][0][1][0], zmm1);
	_mm512_store_pd(&Christoffel_tilde[2][0][2][0], zmm2);
	_mm512_store_pd(&Christoffel_tilde[2][1][1][0], zmm3);
	_mm512_store_pd(&Christoffel_tilde[2][1][2][0], zmm4);
	_mm512_store_pd(&Christoffel_tilde[2][2][2][0], zmm5);

	// for dgup, multiply with the scalar (-2) 
	zmm0 = _mm512_mul_pd(zmm10, zmm16);
	zmm1 = _mm512_mul_pd(zmm11, zmm16);
	zmm2 = _mm512_mul_pd(zmm12, zmm16);
	zmm3 = _mm512_mul_pd(zmm13, zmm16);
	zmm4 = _mm512_mul_pd(zmm14, zmm16);
	zmm5 = _mm512_mul_pd(zmm15, zmm16);
	// save into DD 
	//_mm512_store_pd(&DD[1][1][0][0], zmm11);
	//_mm512_store_pd(&DD[1][2][0][0], zmm12);
	//_mm512_store_pd(&DD[1][2][1][0], zmm14);

	// raise first index of the first block 
	zmm6 = _mm512_mul_pd(igxx, zmm0);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	_mm512_store_pd(&DD[1][0][0][0], zmm10);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	_mm512_store_pd(&DD[1][0][1][0], zmm11);
	zmm7 = _mm512_mul_pd(igxx, zmm1);
	_mm512_store_pd(&DD[1][0][2][0], zmm12);
	zmm7 = _mm512_fmadd_pd(igxy, zmm3, zmm7);
	_mm512_store_pd(&DD[1][1][1][0], zmm13);
	zmm7 = _mm512_fmadd_pd(igxz, zmm4, zmm7);
	_mm512_store_pd(&DD[1][1][2][0], zmm14);
	zmm8 = _mm512_mul_pd(igxx, zmm2);
	_mm512_store_pd(&DD[1][2][2][0], zmm15);
	zmm8 = _mm512_fmadd_pd(igxy, zmm4, zmm8);
	zmm8 = _mm512_fmadd_pd(igxz, zmm5, zmm8);

	zmm9 = _mm512_mul_pd(igxy, zmm0);
	zmm9 = _mm512_fmadd_pd(igyy, zmm1, zmm9);
	zmm9 = _mm512_fmadd_pd(igyz, zmm2, zmm9);
	zmm10 = _mm512_mul_pd(igxy, zmm1);
	zmm10 = _mm512_fmadd_pd(igyy, zmm3, zmm10);
	zmm10 = _mm512_fmadd_pd(igyz, zmm4, zmm10);
	zmm11 = _mm512_mul_pd(igxy, zmm2);
	zmm11 = _mm512_fmadd_pd(igyy, zmm4, zmm11);
	zmm11 = _mm512_fmadd_pd(igyz, zmm5, zmm11);

	zmm12 = _mm512_mul_pd(igxz, zmm0);
	zmm12 = _mm512_fmadd_pd(igyz, zmm1, zmm12);
	zmm12 = _mm512_fmadd_pd(igzz, zmm2, zmm12);
	zmm13 = _mm512_mul_pd(igxz, zmm1);
	zmm13 = _mm512_fmadd_pd(igyz, zmm3, zmm13);
	zmm13 = _mm512_fmadd_pd(igzz, zmm4, zmm13);
	zmm14 = _mm512_mul_pd(igxz, zmm2);
	zmm14 = _mm512_fmadd_pd(igyz, zmm4, zmm14);
	zmm14 = _mm512_fmadd_pd(igzz, zmm5, zmm14);

	// raise the second index. The resulting object (dgup) is symmetric in the ij indices. 

	zmm0 = _mm512_mul_pd(igxx, zmm6);
	zmm0 = _mm512_fmadd_pd(igxy, zmm7, zmm0);
	zmm0 = _mm512_fmadd_pd(igxz, zmm8, zmm0);
	_mm512_store_pd(&dgup[1][0][0][0], zmm0);
	zmm1 = _mm512_mul_pd(igxy, zmm6);
	zmm1 = _mm512_fmadd_pd(igyy, zmm7, zmm1);
	zmm1 = _mm512_fmadd_pd(igyz, zmm8, zmm1);
	_mm512_store_pd(&dgup[1][0][1][0], zmm1);
	zmm2 = _mm512_mul_pd(igxz, zmm6);
	zmm2 = _mm512_fmadd_pd(igyz, zmm7, zmm2);
	zmm2 = _mm512_fmadd_pd(igzz, zmm8, zmm2);
	_mm512_store_pd(&dgup[1][0][2][0], zmm2);
	zmm3 = _mm512_mul_pd(igxy, zmm9);
	zmm3 = _mm512_fmadd_pd(igyy, zmm10, zmm3);
	zmm3 = _mm512_fmadd_pd(igyz, zmm11, zmm3);
	_mm512_store_pd(&dgup[1][1][1][0], zmm3);
	zmm4 = _mm512_mul_pd(igxz, zmm9);
	zmm4 = _mm512_fmadd_pd(igyz, zmm10, zmm4);
	zmm4 = _mm512_fmadd_pd(igzz, zmm11, zmm4);
	_mm512_store_pd(&dgup[1][1][2][0], zmm4);
	zmm5 = _mm512_mul_pd(igxz, zmm12);
	zmm5 = _mm512_fmadd_pd(igyz, zmm13, zmm5);
	zmm5 = _mm512_fmadd_pd(igzz, zmm14, zmm5);

	// save the result into dgup 
	_mm512_store_pd(&dgup[1][2][2][0], zmm5);
	//_mm512_store_pd(&dgup[1][1][0][0], zmm1);
	//_mm512_store_pd(&dgup[1][2][0][0], zmm2);
	//_mm512_store_pd(&dgup[1][2][1][0], zmm4);


	// load scalar ( -1.0 ) and multiply with igxx, igxy and igxz   
	xmm = _mm_load1_pd(&s2[0]);
	zmm0 = _mm512_broadcastsd_pd(xmm);
	zmm17 = _mm512_mul_pd(igxz, zmm0);
	zmm18 = _mm512_mul_pd(igyz, zmm0);
	zmm19 = _mm512_mul_pd(igzz, zmm0);
	// load scalar ( -2.0 )  
	xmm = _mm_load1_pd(&s1[0]);
	zmm16 = _mm512_broadcastsd_pd(xmm);


	// load third slice of D_kij  
	
	zmm10 = _mm512_load_pd(&Q[47][0]);
	zmm11 = _mm512_load_pd(&Q[48][0]);
	zmm12 = _mm512_load_pd(&Q[49][0]);
	zmm13 = _mm512_load_pd(&Q[50][0]);
	zmm14 = _mm512_load_pd(&Q[51][0]);
	zmm15 = _mm512_load_pd(&Q[52][0]);

	// for Christoffel tilde, load the result of the previous slice, multiply DD properly with g_contr and add contribution 
	zmm0 = _mm512_load_pd(&Christoffel_tilde[0][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm17, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&Christoffel_tilde[0][0][1][0]);
	zmm1 = _mm512_fmadd_pd(zmm17, zmm11, zmm1);
	zmm2 = _mm512_load_pd(&Christoffel_tilde[0][0][2][0]);
	zmm2 = _mm512_fmadd_pd(igxx, zmm10, zmm2);
	zmm2 = _mm512_fmadd_pd(igxy, zmm11, zmm2);
	zmm3 = _mm512_load_pd(&Christoffel_tilde[0][1][1][0]);
	zmm3 = _mm512_fmadd_pd(zmm17, zmm13, zmm3);
	zmm4 = _mm512_load_pd(&Christoffel_tilde[0][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxx, zmm11, zmm4);
	zmm4 = _mm512_fmadd_pd(igxy, zmm13, zmm4);
	zmm5 = _mm512_load_pd(&Christoffel_tilde[0][2][2][0]);
	zmm5 = _mm512_fmadd_pd(igxx, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igxx, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igxz, zmm15, zmm5);
	zmm5 = _mm512_fmadd_pd(igxy, zmm14, zmm5);
	zmm5 = _mm512_fmadd_pd(igxy, zmm14, zmm5);
	_mm512_store_pd(&Christoffel_tilde[0][0][0][0], zmm0);
	_mm512_store_pd(&Christoffel_tilde[0][0][1][0], zmm1);
	_mm512_store_pd(&Christoffel_tilde[0][0][2][0], zmm2);
	_mm512_store_pd(&Christoffel_tilde[0][1][1][0], zmm3);
	_mm512_store_pd(&Christoffel_tilde[0][1][2][0], zmm4);
	_mm512_store_pd(&Christoffel_tilde[0][2][2][0], zmm5);
	// we are now on the last slice and Christoffel_tilde is complete. Immediately use this result to compute Gtilde^1  
	zmm6 = _mm512_mul_pd(igxx, zmm0); 
	zmm6 = _mm512_fmadd_pd(igxy, zmm1,zmm6); 
	zmm6 = _mm512_fmadd_pd(igxy, zmm1,zmm6); 
	zmm6 = _mm512_fmadd_pd(igxz, zmm2,zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2,zmm6);
	zmm6 = _mm512_fmadd_pd(igyy, zmm3,zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4,zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4,zmm6);
	zmm6 = _mm512_fmadd_pd(igzz, zmm5,zmm6);
	_mm512_store_pd(&Gtilde[0][0], zmm6); 
	// 
	zmm0 = _mm512_load_pd(&Christoffel_tilde[1][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm18, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&Christoffel_tilde[1][0][1][0]);
	zmm1 = _mm512_fmadd_pd(zmm18, zmm11, zmm1);
	zmm2 = _mm512_load_pd(&Christoffel_tilde[1][0][2][0]);
	zmm2 = _mm512_fmadd_pd(igxy, zmm10, zmm2);
	zmm2 = _mm512_fmadd_pd(igyy, zmm11, zmm2);
	zmm3 = _mm512_load_pd(&Christoffel_tilde[1][1][1][0]);
	zmm3 = _mm512_fmadd_pd(zmm18, zmm13, zmm3);
	zmm4 = _mm512_load_pd(&Christoffel_tilde[1][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxy, zmm11, zmm4);
	zmm4 = _mm512_fmadd_pd(igyy, zmm13, zmm4);
	zmm5 = _mm512_load_pd(&Christoffel_tilde[1][2][2][0]);
	zmm5 = _mm512_fmadd_pd(igxy, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igxy, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igyz, zmm15, zmm5);
	zmm5 = _mm512_fmadd_pd(igyy, zmm14, zmm5);
	zmm5 = _mm512_fmadd_pd(igyy, zmm14, zmm5);
	_mm512_store_pd(&Christoffel_tilde[1][0][0][0], zmm0);
	_mm512_store_pd(&Christoffel_tilde[1][0][1][0], zmm1);
	_mm512_store_pd(&Christoffel_tilde[1][0][2][0], zmm2);
	_mm512_store_pd(&Christoffel_tilde[1][1][1][0], zmm3);
	_mm512_store_pd(&Christoffel_tilde[1][1][2][0], zmm4);
	_mm512_store_pd(&Christoffel_tilde[1][2][2][0], zmm5);
	// we are now on the last slice and Christoffel_tilde is complete. Immediately use this result to compute Gtilde^2  
	zmm6 = _mm512_mul_pd(igxx, zmm0);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igyy, zmm3, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igzz, zmm5, zmm6);
	_mm512_store_pd(&Gtilde[1][0], zmm6);

	// 
	zmm0 = _mm512_load_pd(&Christoffel_tilde[2][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm19, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&Christoffel_tilde[2][0][1][0]);
	zmm1 = _mm512_fmadd_pd(zmm19, zmm11, zmm1);
	zmm2 = _mm512_load_pd(&Christoffel_tilde[2][0][2][0]);
	zmm2 = _mm512_fmadd_pd(igxz, zmm10, zmm2);
	zmm2 = _mm512_fmadd_pd(igyz, zmm11, zmm2);
	zmm3 = _mm512_load_pd(&Christoffel_tilde[2][1][1][0]);
	zmm3 = _mm512_fmadd_pd(zmm19, zmm13, zmm3);
	zmm4 = _mm512_load_pd(&Christoffel_tilde[2][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxz, zmm11, zmm4);
	zmm4 = _mm512_fmadd_pd(igyz, zmm13, zmm4);
	zmm5 = _mm512_load_pd(&Christoffel_tilde[2][2][2][0]);
	zmm5 = _mm512_fmadd_pd(igxz, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igxz, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igzz, zmm15, zmm5);
	zmm5 = _mm512_fmadd_pd(igyz, zmm14, zmm5);
	zmm5 = _mm512_fmadd_pd(igyz, zmm14, zmm5);
	_mm512_store_pd(&Christoffel_tilde[2][0][0][0], zmm0);
	_mm512_store_pd(&Christoffel_tilde[2][0][1][0], zmm1);
	_mm512_store_pd(&Christoffel_tilde[2][0][2][0], zmm2);
	_mm512_store_pd(&Christoffel_tilde[2][1][1][0], zmm3);
	_mm512_store_pd(&Christoffel_tilde[2][1][2][0], zmm4);
	_mm512_store_pd(&Christoffel_tilde[2][2][2][0], zmm5);
	// we are now on the last slice and Christoffel_tilde is complete. Immediately use this result to compute Gtilde^2  
	zmm6 = _mm512_mul_pd(igxx, zmm0);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igyy, zmm3, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igzz, zmm5, zmm6);
	_mm512_store_pd(&Gtilde[2][0], zmm6);


	// for dgup multiply with the scalar (-2.0) 
	zmm0 = _mm512_mul_pd(zmm10, zmm16);
	zmm1 = _mm512_mul_pd(zmm11, zmm16);
	zmm2 = _mm512_mul_pd(zmm12, zmm16);
	zmm3 = _mm512_mul_pd(zmm13, zmm16);
	zmm4 = _mm512_mul_pd(zmm14, zmm16);
	zmm5 = _mm512_mul_pd(zmm15, zmm16);
	// save into DD 
	//_mm512_store_pd(&DD[2][1][0][0], zmm11);
	//_mm512_store_pd(&DD[2][2][0][0], zmm12);
	//_mm512_store_pd(&DD[2][2][1][0], zmm14);

	// raise first index of the first block 
	zmm6 = _mm512_mul_pd(igxx, zmm0);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	_mm512_store_pd(&DD[2][0][0][0], zmm10);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	_mm512_store_pd(&DD[2][0][1][0], zmm11);
	zmm7 = _mm512_mul_pd(igxx, zmm1);
	_mm512_store_pd(&DD[2][0][2][0], zmm12);
	zmm7 = _mm512_fmadd_pd(igxy, zmm3, zmm7);
	_mm512_store_pd(&DD[2][1][1][0], zmm13);
	zmm7 = _mm512_fmadd_pd(igxz, zmm4, zmm7);
	zmm8 = _mm512_mul_pd(igxx, zmm2);
	_mm512_store_pd(&DD[2][1][2][0], zmm14);
	zmm8 = _mm512_fmadd_pd(igxy, zmm4, zmm8);
	zmm8 = _mm512_fmadd_pd(igxz, zmm5, zmm8);
	_mm512_store_pd(&DD[2][2][2][0], zmm15);

	zmm9 = _mm512_mul_pd(igxy, zmm0);
	zmm9 = _mm512_fmadd_pd(igyy, zmm1, zmm9);
	zmm9 = _mm512_fmadd_pd(igyz, zmm2, zmm9);
	zmm10 = _mm512_mul_pd(igxy, zmm1);
	zmm10 = _mm512_fmadd_pd(igyy, zmm3, zmm10);
	zmm10 = _mm512_fmadd_pd(igyz, zmm4, zmm10);
	zmm11 = _mm512_mul_pd(igxy, zmm2);
	zmm11 = _mm512_fmadd_pd(igyy, zmm4, zmm11);
	zmm11 = _mm512_fmadd_pd(igyz, zmm5, zmm11);

	zmm12 = _mm512_mul_pd(igxz, zmm0);
	zmm12 = _mm512_fmadd_pd(igyz, zmm1, zmm12);
	zmm12 = _mm512_fmadd_pd(igzz, zmm2, zmm12);
	zmm13 = _mm512_mul_pd(igxz, zmm1);
	zmm13 = _mm512_fmadd_pd(igyz, zmm3, zmm13);
	zmm13 = _mm512_fmadd_pd(igzz, zmm4, zmm13);
	zmm14 = _mm512_mul_pd(igxz, zmm2);
	zmm14 = _mm512_fmadd_pd(igyz, zmm4, zmm14);
	zmm14 = _mm512_fmadd_pd(igzz, zmm5, zmm14);

	// raise the second index. The resulting object (dgup) is symmetric in the ij indices. 

	zmm0 = _mm512_mul_pd(igxx, zmm6);
	zmm0 = _mm512_fmadd_pd(igxy, zmm7, zmm0);
	zmm0 = _mm512_fmadd_pd(igxz, zmm8, zmm0);
	_mm512_store_pd(&dgup[2][0][0][0], zmm0);
	zmm1 = _mm512_mul_pd(igxy, zmm6);
	zmm1 = _mm512_fmadd_pd(igyy, zmm7, zmm1);
	zmm1 = _mm512_fmadd_pd(igyz, zmm8, zmm1);
	_mm512_store_pd(&dgup[2][0][1][0], zmm1);
	zmm2 = _mm512_mul_pd(igxz, zmm6);
	zmm2 = _mm512_fmadd_pd(igyz, zmm7, zmm2);
	zmm2 = _mm512_fmadd_pd(igzz, zmm8, zmm2);
	_mm512_store_pd(&dgup[2][0][2][0], zmm2);
	zmm3 = _mm512_mul_pd(igxy, zmm9);
	zmm3 = _mm512_fmadd_pd(igyy, zmm10, zmm3);
	zmm3 = _mm512_fmadd_pd(igyz, zmm11, zmm3);
	_mm512_store_pd(&dgup[2][1][1][0], zmm3);
	zmm4 = _mm512_mul_pd(igxz, zmm9);
	zmm4 = _mm512_fmadd_pd(igyz, zmm10, zmm4);
	zmm4 = _mm512_fmadd_pd(igzz, zmm11, zmm4);
	_mm512_store_pd(&dgup[2][1][2][0], zmm4);
	zmm5 = _mm512_mul_pd(igxz, zmm12);
	zmm5 = _mm512_fmadd_pd(igyz, zmm13, zmm5);
	zmm5 = _mm512_fmadd_pd(igzz, zmm14, zmm5);
	_mm512_store_pd(&dgup[2][2][2][0], zmm5);

#ifdef COUNT_FLOPS  
	// cost for the computation of dgup 
	nMul = nMul + VECTORLENGTH * (21 + 30 ) * 3;
	nAdd = nAdd + VECTORLENGTH * 30 * 3;
	// cost for the computation of Christoffel_tilde 
	nMul = nMul + VECTORLENGTH * 12 * 3 * 3;
	nAdd = nAdd + VECTORLENGTH * (12 * 3 * 3 - 6 * 3);
	// cost for the computation of Gtilde 
	nMul = nMul + VECTORLENGTH * 9 * 3;
	nAdd = nAdd + VECTORLENGTH * 8 * 3;
#endif 

	// save the result into dgup 
	//_mm512_store_pd(&dgup[2][1][0][0], zmm1);
	//_mm512_store_pd(&dgup[2][2][0][0], zmm2);
	//_mm512_store_pd(&dgup[2][2][1][0], zmm4);

	// compute the true Christoffel symbol

	// load scalar ( -1.0 )  
	xmm = _mm_load1_pd(&s2[0]); 
	zmm19 = _mm512_broadcastsd_pd(xmm); 

	// we now need to load P_l  
	zmm0 = _mm512_load_pd(&Q[55][0]);
	zmm1 = _mm512_load_pd(&Q[56][0]);
	zmm2 = _mm512_load_pd(&Q[57][0]); 
	// compute Pup_k = g_contr(k,l)*P(l) 
	zmm3 = _mm512_mul_pd(igxx, zmm0); 
	zmm3 = _mm512_fmadd_pd(igxy, zmm1, zmm3);
	zmm3 = _mm512_fmadd_pd(igxz, zmm2, zmm3);
	zmm4 = _mm512_mul_pd(igxy, zmm0);
	zmm4 = _mm512_fmadd_pd(igyy, zmm1, zmm4);
	zmm4 = _mm512_fmadd_pd(igyz, zmm2, zmm4);
	zmm5 = _mm512_mul_pd(igxz, zmm0);
	zmm5 = _mm512_fmadd_pd(igyz, zmm1, zmm5);
	zmm5 = _mm512_fmadd_pd(igzz, zmm2, zmm5);

	// we then load the Christoffel tilde, which should still be in the L1 cache 
	zmm7  = _mm512_load_pd(&Christoffel_tilde[0][0][0][0]);
	zmm7  = _mm512_fmadd_pd(zmm16, zmm0, zmm7);   // multiply with scalar (-2.0) 
	zmm7  = _mm512_fmadd_pd(gxx, zmm3, zmm7);     // add g_cov * Pup 
	_mm512_store_pd(&Christoffel[0][0][0][0], zmm7); 
	zmm7  = _mm512_load_pd(&Christoffel_tilde[0][0][1][0]);
	zmm7 = _mm512_fmadd_pd(zmm19, zmm1, zmm7);   // multiply with scalar (-1.0) 
	zmm7 = _mm512_fmadd_pd(gxy, zmm3, zmm7);     // add g_cov * Pup  
	_mm512_store_pd(&Christoffel[0][0][1][0], zmm7);
	zmm7 = _mm512_load_pd(&Christoffel_tilde[0][0][2][0]);
	zmm7 = _mm512_fmadd_pd(zmm19, zmm2, zmm7);   // multiply with scalar (-1.0) 
	zmm7 = _mm512_fmadd_pd(gxz, zmm3, zmm7);     // add g_cov * Pup  
	_mm512_store_pd(&Christoffel[0][0][2][0], zmm7); 
	zmm7 = _mm512_load_pd(&Christoffel_tilde[0][1][1][0]); 
	zmm7 = _mm512_fmadd_pd(gyy, zmm3, zmm7);     // add g_cov * Pup  
	_mm512_store_pd(&Christoffel[0][1][1][0], zmm7);
	zmm7 = _mm512_load_pd(&Christoffel_tilde[0][1][2][0]);
	zmm7 = _mm512_fmadd_pd(gyz, zmm3, zmm7);     // add g_cov * Pup  
	_mm512_store_pd(&Christoffel[0][1][2][0], zmm7);
	zmm7 = _mm512_load_pd(&Christoffel_tilde[0][2][2][0]);
	zmm7 = _mm512_fmadd_pd(gzz, zmm3, zmm7);     // add g_cov * Pup  
	_mm512_store_pd(&Christoffel[0][2][2][0], zmm7);
	// 
	zmm7 = _mm512_load_pd(&Christoffel_tilde[1][0][0][0]);
	zmm7 = _mm512_fmadd_pd(gxx, zmm4, zmm7);     // add g_cov * Pup  
	_mm512_store_pd(&Christoffel[1][0][0][0], zmm7);
	zmm7 = _mm512_load_pd(&Christoffel_tilde[1][0][1][0]);
	zmm7 = _mm512_fmadd_pd(zmm19, zmm0, zmm7);   // multiply with scalar (-1.0) 
	zmm7 = _mm512_fmadd_pd(gxy, zmm4, zmm7);     // add g_cov * Pup  
	_mm512_store_pd(&Christoffel[1][0][1][0], zmm7);
	zmm7 = _mm512_load_pd(&Christoffel_tilde[1][0][2][0]);
	zmm7 = _mm512_fmadd_pd(gxz, zmm4, zmm7);     // add g_cov * Pup  
	_mm512_store_pd(&Christoffel[1][0][2][0], zmm7);
	zmm7 = _mm512_load_pd(&Christoffel_tilde[1][1][1][0]);
	zmm7 = _mm512_fmadd_pd(zmm16, zmm1, zmm7);   // multiply with scalar (-2.0) 
	zmm7 = _mm512_fmadd_pd(gyy, zmm4, zmm7);     // add g_cov * Pup 
	_mm512_store_pd(&Christoffel[1][1][1][0], zmm7);
	zmm7 = _mm512_load_pd(&Christoffel_tilde[1][1][2][0]);
	zmm7 = _mm512_fmadd_pd(zmm19, zmm2, zmm7);   // multiply with scalar (-1.0) 
	zmm7 = _mm512_fmadd_pd(gyz, zmm4, zmm7);     // add g_cov * Pup  
	_mm512_store_pd(&Christoffel[1][1][2][0], zmm7);
	zmm7 = _mm512_load_pd(&Christoffel_tilde[1][2][2][0]);
	zmm7 = _mm512_fmadd_pd(gzz, zmm4, zmm7);     // add g_cov * Pup  
	_mm512_store_pd(&Christoffel[1][2][2][0], zmm7);
	// 
	zmm7 = _mm512_load_pd(&Christoffel_tilde[2][0][0][0]);
	zmm7 = _mm512_fmadd_pd(gxx, zmm5, zmm7);     // add g_cov * Pup  
	_mm512_store_pd(&Christoffel[2][0][0][0], zmm7);
	zmm7 = _mm512_load_pd(&Christoffel_tilde[2][0][1][0]);
	zmm7 = _mm512_fmadd_pd(gxy, zmm5, zmm7);     // add g_cov * Pup  
	_mm512_store_pd(&Christoffel[2][0][1][0], zmm7);
	zmm7 = _mm512_load_pd(&Christoffel_tilde[2][0][2][0]);
	zmm7 = _mm512_fmadd_pd(zmm19, zmm0, zmm7);   // multiply with scalar (-1.0) 
	zmm7 = _mm512_fmadd_pd(gxz, zmm5, zmm7);     // add g_cov * Pup  
	_mm512_store_pd(&Christoffel[2][0][2][0], zmm7);
	zmm7 = _mm512_load_pd(&Christoffel_tilde[2][1][1][0]);
	zmm7 = _mm512_fmadd_pd(gyy, zmm5, zmm7);     // add g_cov * Pup  
	_mm512_store_pd(&Christoffel[2][1][1][0], zmm7);
	zmm7 = _mm512_load_pd(&Christoffel_tilde[2][1][2][0]);
	zmm7 = _mm512_fmadd_pd(zmm19, zmm1, zmm7);   // multiply with scalar (-1.0) 
	zmm7 = _mm512_fmadd_pd(gyz, zmm5, zmm7);     // add g_cov * Pup  
	_mm512_store_pd(&Christoffel[2][1][2][0], zmm7);
	zmm7 = _mm512_load_pd(&Christoffel_tilde[2][2][2][0]);
	zmm7 = _mm512_fmadd_pd(zmm16, zmm2, zmm7);   // multiply with scalar (-2.0) 
	zmm7 = _mm512_fmadd_pd(gzz, zmm5, zmm7);     // add g_cov * Pup 
	_mm512_store_pd(&Christoffel[2][2][2][0], zmm7);

#ifdef COUNT_FLOPS
	nMul = nMul + VECTORLENGTH * (3 + 33);
	nAdd = nAdd + VECTORLENGTH * 33;
#endif


	// compute Z and Zup 
	zmm14 = _mm512_load_pd(&Q[54][0]);    // load log phi 
	zmm15 = _mm512_set1_pd(1.0); // _mm512_exp_pd(zmm14);         // compute phi via exponential 
	// zmm15 = _mm512_exp_pd(zmm14);         // compute phi via exponential  
	zmm16 = _mm512_mul_pd(zmm15, zmm15);  // phi*phi  
	zmm0 = _mm512_set1_pd(1.0);			  // set unit 
	zmm15 = _mm512_div_pd(zmm0,zmm16);	  //  1 / phi**2 

	xmm = _mm_load1_pd(&s3[0]);		      // scalar 0.5 
	zmm13 = _mm512_broadcastsd_pd(xmm);   
	zmm14 = _mm512_mul_pd(zmm13, zmm16);  // 0.5*phi*phi 

	zmm0 = _mm512_load_pd(&Gtilde[0][0]); // load Gtilde and multiply with 0.5*phi**2  
	zmm17 = _mm512_mul_pd(zmm0, zmm14); 
	zmm1 = _mm512_load_pd(&Gtilde[1][0]);
	zmm18 = _mm512_mul_pd(zmm1, zmm14);
	zmm2 = _mm512_load_pd(&Gtilde[2][0]);
	zmm19 = _mm512_mul_pd(zmm2, zmm14);


	zmm0 = _mm512_load_pd(&Q[13][0]);     // load Ghat 
	zmm17 = _mm512_fmsub_pd(zmm0, zmm14, zmm17); // multiply Ghat with 0.5 and subtract 0.5*phi**2*Gtilde from the result 
	_mm512_store_pd(&Zup[0][0], zmm17);
	zmm1 = _mm512_load_pd(&Q[14][0]);
	zmm18 = _mm512_fmsub_pd(zmm1, zmm14, zmm18);
	_mm512_store_pd(&Zup[1][0], zmm18);
	zmm2 = _mm512_load_pd(&Q[15][0]);
	zmm19 = _mm512_fmsub_pd(zmm2, zmm14, zmm19);
	_mm512_store_pd(&Zup[2][0], zmm19);


	zmm7 = _mm512_mul_pd(zmm17, zmm15);  // scale Zup with 1/phi**2 
	zmm8 = _mm512_mul_pd(zmm18, zmm15);
	zmm9 = _mm512_mul_pd(zmm19, zmm15);

	// multiply with the covariant metric and save result in Z 
	zmm0 = _mm512_mul_pd(gxx, zmm7); 
	zmm0 = _mm512_fmadd_pd(gxy, zmm8, zmm0); 
	zmm0 = _mm512_fmadd_pd(gxz, zmm9, zmm0); 
	_mm512_store_pd(&Z[0][0],zmm0); 

	zmm1 = _mm512_mul_pd(gxy, zmm7);
	zmm1 = _mm512_fmadd_pd(gyy, zmm8, zmm1);
	zmm1 = _mm512_fmadd_pd(gyz, zmm9, zmm1);
	_mm512_store_pd(&Z[1][0],zmm1);

	zmm2 = _mm512_mul_pd(gxz, zmm7);
	zmm2 = _mm512_fmadd_pd(gyz, zmm8, zmm2);
	zmm2 = _mm512_fmadd_pd(gzz, zmm9, zmm2);
	_mm512_store_pd(&Z[2][0],zmm2);


#ifdef COUNT_FLOPS  
	nMul = nMul + VECTORLENGTH * (11 + 9);
	nAdd = nAdd + VECTORLENGTH * 9;
#endif 

	// compute DDsymm 
	// diagonal elements. This is only memory shuffling, no single flop is done here, which is very bad. 
	// xx 
	zmm0 = _mm512_load_pd(&Qx[35][0]);
	_mm512_store_pd(&dDDsymm[0][0][0][0][0], zmm0);
	zmm1 = _mm512_load_pd(&Qx[36][0]);
	_mm512_store_pd(&dDDsymm[0][0][0][1][0], zmm1);
	zmm2 = _mm512_load_pd(&Qx[37][0]);
	_mm512_store_pd(&dDDsymm[0][0][0][2][0], zmm2);
	zmm3 = _mm512_load_pd(&Qx[38][0]);
	_mm512_store_pd(&dDDsymm[0][0][1][1][0], zmm3);
	zmm4 = _mm512_load_pd(&Qx[39][0]);
	_mm512_store_pd(&dDDsymm[0][0][1][2][0], zmm4);
	zmm5 = _mm512_load_pd(&Qx[40][0]);
	_mm512_store_pd(&dDDsymm[0][0][2][2][0], zmm5);
	// yy 
	zmm0 = _mm512_load_pd(&Qy[41][0]);
	_mm512_store_pd(&dDDsymm[1][1][0][0][0], zmm0);
	zmm1 = _mm512_load_pd(&Qy[42][0]);
	_mm512_store_pd(&dDDsymm[1][1][0][1][0], zmm1);
	zmm2 = _mm512_load_pd(&Qy[43][0]);
	_mm512_store_pd(&dDDsymm[1][1][0][2][0], zmm2);
	zmm3 = _mm512_load_pd(&Qy[44][0]);
	_mm512_store_pd(&dDDsymm[1][1][1][1][0], zmm3);
	zmm4 = _mm512_load_pd(&Qy[45][0]);
	_mm512_store_pd(&dDDsymm[1][1][1][2][0], zmm4);
	zmm5 = _mm512_load_pd(&Qy[46][0]);
	_mm512_store_pd(&dDDsymm[1][1][2][2][0], zmm5);
	// zz 
	zmm0 = _mm512_load_pd(&Qz[47][0]);
	_mm512_store_pd(&dDDsymm[2][2][0][0][0], zmm0);
	zmm1 = _mm512_load_pd(&Qz[48][0]);
	_mm512_store_pd(&dDDsymm[2][2][0][1][0], zmm1);
	zmm2 = _mm512_load_pd(&Qz[49][0]);
	_mm512_store_pd(&dDDsymm[2][2][0][2][0], zmm2);
	zmm3 = _mm512_load_pd(&Qz[50][0]);
	_mm512_store_pd(&dDDsymm[2][2][1][1][0], zmm3);
	zmm4 = _mm512_load_pd(&Qz[51][0]);
	_mm512_store_pd(&dDDsymm[2][2][1][2][0], zmm4);
	zmm5 = _mm512_load_pd(&Qz[52][0]);
	_mm512_store_pd(&dDDsymm[2][2][2][2][0], zmm5);
	// mixed elements 
	xmm = _mm_load1_pd(&s3[0]);		      // scalar 0.5 
	zmm19 = _mm512_broadcastsd_pd(xmm);
	// xy 
	zmm0 = _mm512_load_pd(&Qy[35][0]);
	zmm12 = _mm512_mul_pd(zmm0, zmm19); 
	zmm1 = _mm512_load_pd(&Qy[36][0]);
	zmm13 = _mm512_mul_pd(zmm1, zmm19); 
	zmm2 = _mm512_load_pd(&Qy[37][0]);
	zmm14 = _mm512_mul_pd(zmm2, zmm19);
	zmm3 = _mm512_load_pd(&Qy[38][0]);
	zmm15 = _mm512_mul_pd(zmm3, zmm19);
	zmm4 = _mm512_load_pd(&Qy[39][0]);
	zmm16 = _mm512_mul_pd(zmm4, zmm19);
	zmm5 = _mm512_load_pd(&Qy[40][0]);
	zmm17 = _mm512_mul_pd(zmm5, zmm19); 
	zmm6 = _mm512_load_pd(&Qx[41][0]);
	zmm12 = _mm512_fmadd_pd(zmm6, zmm19, zmm12);
	_mm512_store_pd(&dDDsymm[0][1][0][0][0], zmm12);
	zmm7 = _mm512_load_pd(&Qx[42][0]);
	zmm13 = _mm512_fmadd_pd(zmm7, zmm19, zmm13);
	_mm512_store_pd(&dDDsymm[0][1][0][1][0], zmm13);
	zmm8 = _mm512_load_pd(&Qx[43][0]);
	zmm14 = _mm512_fmadd_pd(zmm8, zmm19, zmm14);
	_mm512_store_pd(&dDDsymm[0][1][0][2][0], zmm14);
	zmm9 = _mm512_load_pd(&Qx[44][0]);
	zmm15 = _mm512_fmadd_pd(zmm9, zmm19, zmm15);
	_mm512_store_pd(&dDDsymm[0][1][1][1][0], zmm15);
	zmm10 = _mm512_load_pd(&Qx[45][0]);
	zmm16 = _mm512_fmadd_pd(zmm10, zmm19, zmm16);
	_mm512_store_pd(&dDDsymm[0][1][1][2][0], zmm16);
	zmm11 = _mm512_load_pd(&Qx[46][0]);
	zmm17 = _mm512_fmadd_pd(zmm11, zmm19, zmm17);
	_mm512_store_pd(&dDDsymm[0][1][2][2][0], zmm17);
	// xz 
	zmm0 = _mm512_load_pd(&Qz[35][0]);
	zmm12 = _mm512_mul_pd(zmm0, zmm19);
	zmm1 = _mm512_load_pd(&Qz[36][0]);
	zmm13 = _mm512_mul_pd(zmm1, zmm19);
	zmm2 = _mm512_load_pd(&Qz[37][0]);
	zmm14 = _mm512_mul_pd(zmm2, zmm19);
	zmm3 = _mm512_load_pd(&Qz[38][0]);
	zmm15 = _mm512_mul_pd(zmm3, zmm19);
	zmm4 = _mm512_load_pd(&Qz[39][0]);
	zmm16 = _mm512_mul_pd(zmm4, zmm19);
	zmm5 = _mm512_load_pd(&Qz[40][0]);
	zmm17 = _mm512_mul_pd(zmm5, zmm19);
	zmm6 = _mm512_load_pd(&Qx[47][0]);
	zmm12 = _mm512_fmadd_pd(zmm6, zmm19, zmm12);
	_mm512_store_pd(&dDDsymm[0][2][0][0][0], zmm12);
	zmm7 = _mm512_load_pd(&Qx[48][0]);
	zmm13 = _mm512_fmadd_pd(zmm7, zmm19, zmm13);
	_mm512_store_pd(&dDDsymm[0][2][0][1][0], zmm13);
	zmm8 = _mm512_load_pd(&Qx[49][0]);
	zmm14 = _mm512_fmadd_pd(zmm8, zmm19, zmm14);
	_mm512_store_pd(&dDDsymm[0][2][0][2][0], zmm14);
	zmm9 = _mm512_load_pd(&Qx[50][0]);
	zmm15 = _mm512_fmadd_pd(zmm9, zmm19, zmm15);
	_mm512_store_pd(&dDDsymm[0][2][1][1][0], zmm15);
	zmm10 = _mm512_load_pd(&Qx[51][0]);
	zmm16 = _mm512_fmadd_pd(zmm10, zmm19, zmm16);
	_mm512_store_pd(&dDDsymm[0][2][1][2][0], zmm16);
	zmm11 = _mm512_load_pd(&Qx[52][0]);
	zmm17 = _mm512_fmadd_pd(zmm11, zmm19, zmm17);
	_mm512_store_pd(&dDDsymm[0][2][2][2][0], zmm17);
	// yz  
	zmm0 = _mm512_load_pd(&Qz[41][0]);
	zmm12 = _mm512_mul_pd(zmm0, zmm19);
	zmm1 = _mm512_load_pd(&Qz[42][0]);
	zmm13 = _mm512_mul_pd(zmm1, zmm19);
	zmm2 = _mm512_load_pd(&Qz[43][0]);
	zmm14 = _mm512_mul_pd(zmm2, zmm19);
	zmm3 = _mm512_load_pd(&Qz[44][0]);
	zmm15 = _mm512_mul_pd(zmm3, zmm19);
	zmm4 = _mm512_load_pd(&Qz[45][0]);
	zmm16 = _mm512_mul_pd(zmm4, zmm19);
	zmm5 = _mm512_load_pd(&Qz[46][0]);
	zmm17 = _mm512_mul_pd(zmm5, zmm19);
	zmm6 = _mm512_load_pd(&Qy[47][0]);
	zmm12 = _mm512_fmadd_pd(zmm6, zmm19, zmm12);
	_mm512_store_pd(&dDDsymm[1][2][0][0][0], zmm12);
	zmm7 = _mm512_load_pd(&Qy[48][0]);
	zmm13 = _mm512_fmadd_pd(zmm7, zmm19, zmm13);
	_mm512_store_pd(&dDDsymm[1][2][0][1][0], zmm13);
	zmm8 = _mm512_load_pd(&Qy[49][0]);
	zmm14 = _mm512_fmadd_pd(zmm8, zmm19, zmm14);
	_mm512_store_pd(&dDDsymm[1][2][0][2][0], zmm14);
	zmm9 = _mm512_load_pd(&Qy[50][0]);
	zmm15 = _mm512_fmadd_pd(zmm9, zmm19, zmm15);
	_mm512_store_pd(&dDDsymm[1][2][1][1][0], zmm15);
	zmm10 = _mm512_load_pd(&Qy[51][0]);
	zmm16 = _mm512_fmadd_pd(zmm10, zmm19, zmm16);
	_mm512_store_pd(&dDDsymm[1][2][1][2][0], zmm16);
	zmm11 = _mm512_load_pd(&Qy[52][0]);
	zmm17 = _mm512_fmadd_pd(zmm11, zmm19, zmm17);
	_mm512_store_pd(&dDDsymm[1][2][2][2][0], zmm17);

#ifdef COUNT_FLOPS
	nMul = nMul + VECTORLENGTH * 3 * 12;
	nAdd = nAdd + VECTORLENGTH * 3 * 6;
#endif

	// Compute the spatial derivative of the Christoffel symbol

	// **** 1/3 **** 
	// load scalar ( -1.0 ) and multiply with igxx, igxy and igxz   
	xmm = _mm_load1_pd(&s2[0]);
	zmm0 = _mm512_broadcastsd_pd(xmm);
	zmm17 = _mm512_mul_pd(igxx, zmm0);
	zmm18 = _mm512_mul_pd(igxy, zmm0);
	zmm19 = _mm512_mul_pd(igxz, zmm0);
#ifdef COUNT_FLOPS  
	nMul = nMul + VECTORLENGTH * 3 * 3;   // this is counted for all three slices. 
#endif 

	// load first slice of dDDsymm_klij  (D_1ij) 
	zmm10 = _mm512_load_pd(&dDDsymm[0][0][0][0][0]);
	zmm11 = _mm512_load_pd(&dDDsymm[0][0][0][1][0]);
	zmm12 = _mm512_load_pd(&dDDsymm[0][0][0][2][0]);
	zmm13 = _mm512_load_pd(&dDDsymm[0][0][1][1][0]);
	zmm14 = _mm512_load_pd(&dDDsymm[0][0][1][2][0]);
	zmm15 = _mm512_load_pd(&dDDsymm[0][0][2][2][0]);

	// for dChristoffel tildeNCP, multiply with g_contr 
	zmm0 = _mm512_mul_pd(igxx, zmm10);
	zmm0 = _mm512_fmadd_pd(igxy, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igxy, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igxz, zmm12, zmm0);
	zmm0 = _mm512_fmadd_pd(igxz, zmm12, zmm0);
	zmm1 = _mm512_mul_pd(igxy, zmm13);
	zmm1 = _mm512_fmadd_pd(igxz, zmm14, zmm1);
	zmm2 = _mm512_mul_pd(igxy, zmm14);
	zmm2 = _mm512_fmadd_pd(igxz, zmm15, zmm2);
	zmm3 = _mm512_mul_pd(zmm17, zmm13);
	zmm4 = _mm512_mul_pd(zmm17, zmm14);
	zmm5 = _mm512_mul_pd(zmm17, zmm15);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][0][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][0][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][0][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][0][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][0][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][0][2][2][0], zmm5);
	// 
	zmm0 = _mm512_mul_pd(igxy, zmm10);
	zmm0 = _mm512_fmadd_pd(igyy, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igyy, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igyz, zmm12, zmm0);
	zmm0 = _mm512_fmadd_pd(igyz, zmm12, zmm0);
	zmm1 = _mm512_mul_pd(igyy, zmm13);
	zmm1 = _mm512_fmadd_pd(igyz, zmm14, zmm1);
	zmm2 = _mm512_mul_pd(igyy, zmm14);
	zmm2 = _mm512_fmadd_pd(igyz, zmm15, zmm2);
	zmm3 = _mm512_mul_pd(zmm18, zmm13);
	zmm4 = _mm512_mul_pd(zmm18, zmm14);
	zmm5 = _mm512_mul_pd(zmm18, zmm15);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][1][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][1][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][1][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][1][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][1][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][1][2][2][0], zmm5);
	// 
	zmm0 = _mm512_mul_pd(igxz, zmm10);
	zmm0 = _mm512_fmadd_pd(igyz, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igyz, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igzz, zmm12, zmm0);
	zmm0 = _mm512_fmadd_pd(igzz, zmm12, zmm0);
	zmm1 = _mm512_mul_pd(igyz, zmm13);
	zmm1 = _mm512_fmadd_pd(igzz, zmm14, zmm1);
	zmm2 = _mm512_mul_pd(igyz, zmm14);
	zmm2 = _mm512_fmadd_pd(igzz, zmm15, zmm2);
	zmm3 = _mm512_mul_pd(zmm19, zmm13);
	zmm4 = _mm512_mul_pd(zmm19, zmm14);
	zmm5 = _mm512_mul_pd(zmm19, zmm15);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][2][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][2][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][2][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][2][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][2][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][2][2][2][0], zmm5);

	// load scalar ( -1.0 ) and multiply with igxx, igxy and igxz   
	xmm = _mm_load1_pd(&s2[0]);
	zmm0 = _mm512_broadcastsd_pd(xmm);
	zmm17 = _mm512_mul_pd(igxy, zmm0);
	zmm18 = _mm512_mul_pd(igyy, zmm0);
	zmm19 = _mm512_mul_pd(igyz, zmm0);

	// load second slice of dDDsymm_klij 
	zmm10 = _mm512_load_pd(&dDDsymm[0][1][0][0][0]); 
	zmm11 = _mm512_load_pd(&dDDsymm[0][1][0][1][0]);
	zmm12 = _mm512_load_pd(&dDDsymm[0][1][0][2][0]);
	zmm13 = _mm512_load_pd(&dDDsymm[0][1][1][1][0]);
	zmm14 = _mm512_load_pd(&dDDsymm[0][1][1][2][0]);
	zmm15 = _mm512_load_pd(&dDDsymm[0][1][2][2][0]);

	// for Christoffel tilde, load the result of the previous slice, multiply DD properly with g_contr and add contribution 
	zmm0 = _mm512_load_pd(&dChristoffel_tildeNCP[0][0][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm17, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&dChristoffel_tildeNCP[0][0][0][1][0]);
	zmm1 = _mm512_fmadd_pd(igxx, zmm10, zmm1);
	zmm1 = _mm512_fmadd_pd(igxz, zmm12, zmm1);
	zmm2 = _mm512_load_pd(&dChristoffel_tildeNCP[0][0][0][2][0]);
	zmm2 = _mm512_fmadd_pd(zmm17, zmm12, zmm2);
	zmm3 = _mm512_load_pd(&dChristoffel_tildeNCP[0][0][1][1][0]);
	zmm3 = _mm512_fmadd_pd(igxx, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igxx, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igxy, zmm13, zmm3);
	zmm3 = _mm512_fmadd_pd(igxz, zmm14, zmm3);
	zmm3 = _mm512_fmadd_pd(igxz, zmm14, zmm3);
	zmm4 = _mm512_load_pd(&dChristoffel_tildeNCP[0][0][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxx, zmm12, zmm4);
	zmm4 = _mm512_fmadd_pd(igxz, zmm15, zmm4);
	zmm5 = _mm512_load_pd(&dChristoffel_tildeNCP[0][0][2][2][0]);
	zmm5 = _mm512_fmadd_pd(zmm17, zmm15, zmm5);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][0][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][0][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][0][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][0][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][0][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][0][2][2][0], zmm5);
	// 
	zmm0 = _mm512_load_pd(&dChristoffel_tildeNCP[0][1][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm18, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&dChristoffel_tildeNCP[0][1][0][1][0]);
	zmm1 = _mm512_fmadd_pd(igxy, zmm10, zmm1);
	zmm1 = _mm512_fmadd_pd(igyz, zmm12, zmm1);
	zmm2 = _mm512_load_pd(&dChristoffel_tildeNCP[0][1][0][2][0]);
	zmm2 = _mm512_fmadd_pd(zmm18, zmm12, zmm2);
	zmm3 = _mm512_load_pd(&dChristoffel_tildeNCP[0][1][1][1][0]);
	zmm3 = _mm512_fmadd_pd(igxy, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igxy, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igyy, zmm13, zmm3);
	zmm3 = _mm512_fmadd_pd(igyz, zmm14, zmm3);
	zmm3 = _mm512_fmadd_pd(igyz, zmm14, zmm3);
	zmm4 = _mm512_load_pd(&dChristoffel_tildeNCP[0][1][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxy, zmm12, zmm4);
	zmm4 = _mm512_fmadd_pd(igyz, zmm15, zmm4);
	zmm5 = _mm512_load_pd(&dChristoffel_tildeNCP[0][1][2][2][0]);
	zmm5 = _mm512_fmadd_pd(zmm18, zmm15, zmm5);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][1][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][1][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][1][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][1][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][1][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][1][2][2][0], zmm5);
	// 
	zmm0 = _mm512_load_pd(&dChristoffel_tildeNCP[0][2][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm19, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&dChristoffel_tildeNCP[0][2][0][1][0]);
	zmm1 = _mm512_fmadd_pd(igxz, zmm10, zmm1);
	zmm1 = _mm512_fmadd_pd(igzz, zmm12, zmm1);
	zmm2 = _mm512_load_pd(&dChristoffel_tildeNCP[0][2][0][2][0]);
	zmm2 = _mm512_fmadd_pd(zmm19, zmm12, zmm2);
	zmm3 = _mm512_load_pd(&dChristoffel_tildeNCP[0][2][1][1][0]);
	zmm3 = _mm512_fmadd_pd(igxz, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igxz, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igyz, zmm13, zmm3);
	zmm3 = _mm512_fmadd_pd(igzz, zmm14, zmm3);
	zmm3 = _mm512_fmadd_pd(igzz, zmm14, zmm3);
	zmm4 = _mm512_load_pd(&dChristoffel_tildeNCP[0][2][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxz, zmm12, zmm4);
	zmm4 = _mm512_fmadd_pd(igzz, zmm15, zmm4);
	zmm5 = _mm512_load_pd(&dChristoffel_tildeNCP[0][2][2][2][0]);
	zmm5 = _mm512_fmadd_pd(zmm19, zmm15, zmm5);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][2][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][2][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][2][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][2][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][2][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][2][2][2][0], zmm5);


	// load scalar ( -1.0 ) and multiply with igxx, igxy and igxz   
	xmm = _mm_load1_pd(&s2[0]);
	zmm0 = _mm512_broadcastsd_pd(xmm);
	zmm17 = _mm512_mul_pd(igxz, zmm0);
	zmm18 = _mm512_mul_pd(igyz, zmm0);
	zmm19 = _mm512_mul_pd(igzz, zmm0);

	// load third slice of D_kij  
	zmm10 = _mm512_load_pd(&dDDsymm[0][2][0][0][0]);
	zmm11 = _mm512_load_pd(&dDDsymm[0][2][0][1][0]);
	zmm12 = _mm512_load_pd(&dDDsymm[0][2][0][2][0]);
	zmm13 = _mm512_load_pd(&dDDsymm[0][2][1][1][0]);
	zmm14 = _mm512_load_pd(&dDDsymm[0][2][1][2][0]);
	zmm15 = _mm512_load_pd(&dDDsymm[0][2][2][2][0]);

	// for Christoffel tilde, load the result of the previous slice, multiply DD properly with g_contr and add contribution 
	zmm0 = _mm512_load_pd(&dChristoffel_tildeNCP[0][0][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm17, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&dChristoffel_tildeNCP[0][0][0][1][0]);
	zmm1 = _mm512_fmadd_pd(zmm17, zmm11, zmm1);
	zmm2 = _mm512_load_pd(&dChristoffel_tildeNCP[0][0][0][2][0]);
	zmm2 = _mm512_fmadd_pd(igxx, zmm10, zmm2);
	zmm2 = _mm512_fmadd_pd(igxy, zmm11, zmm2);
	zmm3 = _mm512_load_pd(&dChristoffel_tildeNCP[0][0][1][1][0]);
	zmm3 = _mm512_fmadd_pd(zmm17, zmm13, zmm3);
	zmm4 = _mm512_load_pd(&dChristoffel_tildeNCP[0][0][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxx, zmm11, zmm4);
	zmm4 = _mm512_fmadd_pd(igxy, zmm13, zmm4);
	zmm5 = _mm512_load_pd(&dChristoffel_tildeNCP[0][0][2][2][0]);
	zmm5 = _mm512_fmadd_pd(igxx, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igxx, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igxz, zmm15, zmm5);
	zmm5 = _mm512_fmadd_pd(igxy, zmm14, zmm5);
	zmm5 = _mm512_fmadd_pd(igxy, zmm14, zmm5);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][0][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][0][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][0][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][0][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][0][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][0][2][2][0], zmm5);
	/* we are now on the last slice and Christoffel_tilde is complete. Immediately use this result to compute Gtilde^1  
	zmm6 = _mm512_mul_pd(igxx, zmm0);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igyy, zmm3, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igzz, zmm5, zmm6);
	_mm512_store_pd(&Gtilde[0][0], zmm6); */ 
	// 
	zmm0 = _mm512_load_pd(&dChristoffel_tildeNCP[0][1][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm18, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&dChristoffel_tildeNCP[0][1][0][1][0]);
	zmm1 = _mm512_fmadd_pd(zmm18, zmm11, zmm1);
	zmm2 = _mm512_load_pd(&dChristoffel_tildeNCP[0][1][0][2][0]);
	zmm2 = _mm512_fmadd_pd(igxy, zmm10, zmm2);
	zmm2 = _mm512_fmadd_pd(igyy, zmm11, zmm2);
	zmm3 = _mm512_load_pd(&dChristoffel_tildeNCP[0][1][1][1][0]);
	zmm3 = _mm512_fmadd_pd(zmm18, zmm13, zmm3);
	zmm4 = _mm512_load_pd(&dChristoffel_tildeNCP[0][1][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxy, zmm11, zmm4);
	zmm4 = _mm512_fmadd_pd(igyy, zmm13, zmm4);
	zmm5 = _mm512_load_pd(&dChristoffel_tildeNCP[0][1][2][2][0]);
	zmm5 = _mm512_fmadd_pd(igxy, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igxy, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igyz, zmm15, zmm5);
	zmm5 = _mm512_fmadd_pd(igyy, zmm14, zmm5);
	zmm5 = _mm512_fmadd_pd(igyy, zmm14, zmm5);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][1][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][1][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][1][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][1][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][1][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][1][2][2][0], zmm5);
	/* we are now on the last slice and Christoffel_tilde is complete. Immediately use this result to compute Gtilde^2  
	zmm6 = _mm512_mul_pd(igxx, zmm0);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igyy, zmm3, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igzz, zmm5, zmm6);
	_mm512_store_pd(&Gtilde[1][0], zmm6); */ 

	// 
	zmm0 = _mm512_load_pd(&dChristoffel_tildeNCP[0][2][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm19, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&dChristoffel_tildeNCP[0][2][0][1][0]);
	zmm1 = _mm512_fmadd_pd(zmm19, zmm11, zmm1);
	zmm2 = _mm512_load_pd(&dChristoffel_tildeNCP[0][2][0][2][0]);
	zmm2 = _mm512_fmadd_pd(igxz, zmm10, zmm2);
	zmm2 = _mm512_fmadd_pd(igyz, zmm11, zmm2);
	zmm3 = _mm512_load_pd(&dChristoffel_tildeNCP[0][2][1][1][0]);
	zmm3 = _mm512_fmadd_pd(zmm19, zmm13, zmm3);
	zmm4 = _mm512_load_pd(&dChristoffel_tildeNCP[0][2][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxz, zmm11, zmm4);
	zmm4 = _mm512_fmadd_pd(igyz, zmm13, zmm4);
	zmm5 = _mm512_load_pd(&dChristoffel_tildeNCP[0][2][2][2][0]);
	zmm5 = _mm512_fmadd_pd(igxz, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igxz, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igzz, zmm15, zmm5);
	zmm5 = _mm512_fmadd_pd(igyz, zmm14, zmm5);
	zmm5 = _mm512_fmadd_pd(igyz, zmm14, zmm5);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][2][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][2][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][2][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][2][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][2][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[0][2][2][2][0], zmm5);
	/* we are now on the last slice and Christoffel_tilde is complete. Immediately use this result to compute Gtilde^2  
	zmm6 = _mm512_mul_pd(igxx, zmm0);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igyy, zmm3, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igzz, zmm5, zmm6);
	_mm512_store_pd(&Gtilde[2][0], zmm6); */ 

	// **** 2/3 **** 
	// load scalar ( -1.0 ) and multiply with igxx, igxy and igxz   
	xmm = _mm_load1_pd(&s2[0]);
	zmm0 = _mm512_broadcastsd_pd(xmm);
	zmm17 = _mm512_mul_pd(igxx, zmm0);
	zmm18 = _mm512_mul_pd(igxy, zmm0);
	zmm19 = _mm512_mul_pd(igxz, zmm0);
#ifdef COUNT_FLOPS  
	nMul = nMul + VECTORLENGTH * 3 * 3;   // this is counted for all three slices. 
#endif 

										  // load first slice of dDDsymm_klij  (D_1ij) 
	zmm10 = _mm512_load_pd(&dDDsymm[0][1][0][0][0]);
	zmm11 = _mm512_load_pd(&dDDsymm[0][1][0][1][0]);
	zmm12 = _mm512_load_pd(&dDDsymm[0][1][0][2][0]);
	zmm13 = _mm512_load_pd(&dDDsymm[0][1][1][1][0]);
	zmm14 = _mm512_load_pd(&dDDsymm[0][1][1][2][0]);
	zmm15 = _mm512_load_pd(&dDDsymm[0][1][2][2][0]);

	// for dChristoffel tildeNCP, multiply with g_contr 
	zmm0 = _mm512_mul_pd(igxx, zmm10);
	zmm0 = _mm512_fmadd_pd(igxy, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igxy, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igxz, zmm12, zmm0);
	zmm0 = _mm512_fmadd_pd(igxz, zmm12, zmm0);
	zmm1 = _mm512_mul_pd(igxy, zmm13);
	zmm1 = _mm512_fmadd_pd(igxz, zmm14, zmm1);
	zmm2 = _mm512_mul_pd(igxy, zmm14);
	zmm2 = _mm512_fmadd_pd(igxz, zmm15, zmm2);
	zmm3 = _mm512_mul_pd(zmm17, zmm13);
	zmm4 = _mm512_mul_pd(zmm17, zmm14);
	zmm5 = _mm512_mul_pd(zmm17, zmm15);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][0][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][0][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][0][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][0][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][0][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][0][2][2][0], zmm5);
	// 
	zmm0 = _mm512_mul_pd(igxy, zmm10);
	zmm0 = _mm512_fmadd_pd(igyy, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igyy, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igyz, zmm12, zmm0);
	zmm0 = _mm512_fmadd_pd(igyz, zmm12, zmm0);
	zmm1 = _mm512_mul_pd(igyy, zmm13);
	zmm1 = _mm512_fmadd_pd(igyz, zmm14, zmm1);
	zmm2 = _mm512_mul_pd(igyy, zmm14);
	zmm2 = _mm512_fmadd_pd(igyz, zmm15, zmm2);
	zmm3 = _mm512_mul_pd(zmm18, zmm13);
	zmm4 = _mm512_mul_pd(zmm18, zmm14);
	zmm5 = _mm512_mul_pd(zmm18, zmm15);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][1][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][1][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][1][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][1][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][1][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][1][2][2][0], zmm5);
	// 
	zmm0 = _mm512_mul_pd(igxz, zmm10);
	zmm0 = _mm512_fmadd_pd(igyz, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igyz, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igzz, zmm12, zmm0);
	zmm0 = _mm512_fmadd_pd(igzz, zmm12, zmm0);
	zmm1 = _mm512_mul_pd(igyz, zmm13);
	zmm1 = _mm512_fmadd_pd(igzz, zmm14, zmm1);
	zmm2 = _mm512_mul_pd(igyz, zmm14);
	zmm2 = _mm512_fmadd_pd(igzz, zmm15, zmm2);
	zmm3 = _mm512_mul_pd(zmm19, zmm13);
	zmm4 = _mm512_mul_pd(zmm19, zmm14);
	zmm5 = _mm512_mul_pd(zmm19, zmm15);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][2][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][2][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][2][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][2][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][2][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][2][2][2][0], zmm5);

	// load scalar ( -1.0 ) and multiply with igxx, igxy and igxz   
	xmm = _mm_load1_pd(&s2[0]);
	zmm0 = _mm512_broadcastsd_pd(xmm);
	zmm17 = _mm512_mul_pd(igxy, zmm0);
	zmm18 = _mm512_mul_pd(igyy, zmm0);
	zmm19 = _mm512_mul_pd(igyz, zmm0);

	// load second slice of dDDsymm_klij 
	zmm10 = _mm512_load_pd(&dDDsymm[1][1][0][0][0]);
	zmm11 = _mm512_load_pd(&dDDsymm[1][1][0][1][0]);
	zmm12 = _mm512_load_pd(&dDDsymm[1][1][0][2][0]);
	zmm13 = _mm512_load_pd(&dDDsymm[1][1][1][1][0]);
	zmm14 = _mm512_load_pd(&dDDsymm[1][1][1][2][0]);
	zmm15 = _mm512_load_pd(&dDDsymm[1][1][2][2][0]);

	// for Christoffel tilde, load the result of the previous slice, multiply DD properly with g_contr and add contribution 
	zmm0 = _mm512_load_pd(&dChristoffel_tildeNCP[1][0][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm17, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&dChristoffel_tildeNCP[1][0][0][1][0]);
	zmm1 = _mm512_fmadd_pd(igxx, zmm10, zmm1);
	zmm1 = _mm512_fmadd_pd(igxz, zmm12, zmm1);
	zmm2 = _mm512_load_pd(&dChristoffel_tildeNCP[1][0][0][2][0]);
	zmm2 = _mm512_fmadd_pd(zmm17, zmm12, zmm2);
	zmm3 = _mm512_load_pd(&dChristoffel_tildeNCP[1][0][1][1][0]);
	zmm3 = _mm512_fmadd_pd(igxx, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igxx, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igxy, zmm13, zmm3);
	zmm3 = _mm512_fmadd_pd(igxz, zmm14, zmm3);
	zmm3 = _mm512_fmadd_pd(igxz, zmm14, zmm3);
	zmm4 = _mm512_load_pd(&dChristoffel_tildeNCP[1][0][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxx, zmm12, zmm4);
	zmm4 = _mm512_fmadd_pd(igxz, zmm15, zmm4);
	zmm5 = _mm512_load_pd(&dChristoffel_tildeNCP[1][0][2][2][0]);
	zmm5 = _mm512_fmadd_pd(zmm17, zmm15, zmm5);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][0][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][0][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][0][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][0][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][0][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][0][2][2][0], zmm5);
	// 
	zmm0 = _mm512_load_pd(&dChristoffel_tildeNCP[1][1][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm18, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&dChristoffel_tildeNCP[1][1][0][1][0]);
	zmm1 = _mm512_fmadd_pd(igxy, zmm10, zmm1);
	zmm1 = _mm512_fmadd_pd(igyz, zmm12, zmm1);
	zmm2 = _mm512_load_pd(&dChristoffel_tildeNCP[1][1][0][2][0]);
	zmm2 = _mm512_fmadd_pd(zmm18, zmm12, zmm2);
	zmm3 = _mm512_load_pd(&dChristoffel_tildeNCP[1][1][1][1][0]);
	zmm3 = _mm512_fmadd_pd(igxy, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igxy, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igyy, zmm13, zmm3);
	zmm3 = _mm512_fmadd_pd(igyz, zmm14, zmm3);
	zmm3 = _mm512_fmadd_pd(igyz, zmm14, zmm3);
	zmm4 = _mm512_load_pd(&dChristoffel_tildeNCP[1][1][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxy, zmm12, zmm4);
	zmm4 = _mm512_fmadd_pd(igyz, zmm15, zmm4);
	zmm5 = _mm512_load_pd(&dChristoffel_tildeNCP[1][1][2][2][0]);
	zmm5 = _mm512_fmadd_pd(zmm18, zmm15, zmm5);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][1][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][1][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][1][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][1][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][1][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][1][2][2][0], zmm5);
	// 
	zmm0 = _mm512_load_pd(&dChristoffel_tildeNCP[1][2][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm19, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&dChristoffel_tildeNCP[1][2][0][1][0]);
	zmm1 = _mm512_fmadd_pd(igxz, zmm10, zmm1);
	zmm1 = _mm512_fmadd_pd(igzz, zmm12, zmm1);
	zmm2 = _mm512_load_pd(&dChristoffel_tildeNCP[1][2][0][2][0]);
	zmm2 = _mm512_fmadd_pd(zmm19, zmm12, zmm2);
	zmm3 = _mm512_load_pd(&dChristoffel_tildeNCP[1][2][1][1][0]);
	zmm3 = _mm512_fmadd_pd(igxz, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igxz, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igyz, zmm13, zmm3);
	zmm3 = _mm512_fmadd_pd(igzz, zmm14, zmm3);
	zmm3 = _mm512_fmadd_pd(igzz, zmm14, zmm3);
	zmm4 = _mm512_load_pd(&dChristoffel_tildeNCP[1][2][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxz, zmm12, zmm4);
	zmm4 = _mm512_fmadd_pd(igzz, zmm15, zmm4);
	zmm5 = _mm512_load_pd(&dChristoffel_tildeNCP[1][2][2][2][0]);
	zmm5 = _mm512_fmadd_pd(zmm19, zmm15, zmm5);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][2][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][2][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][2][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][2][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][2][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][2][2][2][0], zmm5);


	// load scalar ( -1.0 ) and multiply with igxx, igxy and igxz   
	xmm = _mm_load1_pd(&s2[0]);
	zmm0 = _mm512_broadcastsd_pd(xmm);
	zmm17 = _mm512_mul_pd(igxz, zmm0);
	zmm18 = _mm512_mul_pd(igyz, zmm0);
	zmm19 = _mm512_mul_pd(igzz, zmm0);

	// load third slice of D_kij  
	zmm10 = _mm512_load_pd(&dDDsymm[1][2][0][0][0]);
	zmm11 = _mm512_load_pd(&dDDsymm[1][2][0][1][0]);
	zmm12 = _mm512_load_pd(&dDDsymm[1][2][0][2][0]);
	zmm13 = _mm512_load_pd(&dDDsymm[1][2][1][1][0]);
	zmm14 = _mm512_load_pd(&dDDsymm[1][2][1][2][0]);
	zmm15 = _mm512_load_pd(&dDDsymm[1][2][2][2][0]);

	// for Christoffel tilde, load the result of the previous slice, multiply DD properly with g_contr and add contribution 
	zmm0 = _mm512_load_pd(&dChristoffel_tildeNCP[1][0][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm17, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&dChristoffel_tildeNCP[1][0][0][1][0]);
	zmm1 = _mm512_fmadd_pd(zmm17, zmm11, zmm1);
	zmm2 = _mm512_load_pd(&dChristoffel_tildeNCP[1][0][0][2][0]);
	zmm2 = _mm512_fmadd_pd(igxx, zmm10, zmm2);
	zmm2 = _mm512_fmadd_pd(igxy, zmm11, zmm2);
	zmm3 = _mm512_load_pd(&dChristoffel_tildeNCP[1][0][1][1][0]);
	zmm3 = _mm512_fmadd_pd(zmm17, zmm13, zmm3);
	zmm4 = _mm512_load_pd(&dChristoffel_tildeNCP[1][0][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxx, zmm11, zmm4);
	zmm4 = _mm512_fmadd_pd(igxy, zmm13, zmm4);
	zmm5 = _mm512_load_pd(&dChristoffel_tildeNCP[1][0][2][2][0]);
	zmm5 = _mm512_fmadd_pd(igxx, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igxx, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igxz, zmm15, zmm5);
	zmm5 = _mm512_fmadd_pd(igxy, zmm14, zmm5);
	zmm5 = _mm512_fmadd_pd(igxy, zmm14, zmm5);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][0][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][0][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][0][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][0][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][0][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][0][2][2][0], zmm5);
	/* we are now on the last slice and Christoffel_tilde is complete. Immediately use this result to compute Gtilde^1
	zmm6 = _mm512_mul_pd(igxx, zmm0);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igyy, zmm3, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igzz, zmm5, zmm6);
	_mm512_store_pd(&Gtilde[0][0], zmm6); */
	// 
	zmm0 = _mm512_load_pd(&dChristoffel_tildeNCP[1][1][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm18, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&dChristoffel_tildeNCP[1][1][0][1][0]);
	zmm1 = _mm512_fmadd_pd(zmm18, zmm11, zmm1);
	zmm2 = _mm512_load_pd(&dChristoffel_tildeNCP[1][1][0][2][0]);
	zmm2 = _mm512_fmadd_pd(igxy, zmm10, zmm2);
	zmm2 = _mm512_fmadd_pd(igyy, zmm11, zmm2);
	zmm3 = _mm512_load_pd(&dChristoffel_tildeNCP[1][1][1][1][0]);
	zmm3 = _mm512_fmadd_pd(zmm18, zmm13, zmm3);
	zmm4 = _mm512_load_pd(&dChristoffel_tildeNCP[1][1][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxy, zmm11, zmm4);
	zmm4 = _mm512_fmadd_pd(igyy, zmm13, zmm4);
	zmm5 = _mm512_load_pd(&dChristoffel_tildeNCP[1][1][2][2][0]);
	zmm5 = _mm512_fmadd_pd(igxy, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igxy, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igyz, zmm15, zmm5);
	zmm5 = _mm512_fmadd_pd(igyy, zmm14, zmm5);
	zmm5 = _mm512_fmadd_pd(igyy, zmm14, zmm5);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][1][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][1][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][1][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][1][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][1][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][1][2][2][0], zmm5);
	/* we are now on the last slice and Christoffel_tilde is complete. Immediately use this result to compute Gtilde^2
	zmm6 = _mm512_mul_pd(igxx, zmm0);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igyy, zmm3, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igzz, zmm5, zmm6);
	_mm512_store_pd(&Gtilde[1][0], zmm6); */

	// 
	zmm0 = _mm512_load_pd(&dChristoffel_tildeNCP[1][2][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm19, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&dChristoffel_tildeNCP[1][2][0][1][0]);
	zmm1 = _mm512_fmadd_pd(zmm19, zmm11, zmm1);
	zmm2 = _mm512_load_pd(&dChristoffel_tildeNCP[1][2][0][2][0]);
	zmm2 = _mm512_fmadd_pd(igxz, zmm10, zmm2);
	zmm2 = _mm512_fmadd_pd(igyz, zmm11, zmm2);
	zmm3 = _mm512_load_pd(&dChristoffel_tildeNCP[1][2][1][1][0]);
	zmm3 = _mm512_fmadd_pd(zmm19, zmm13, zmm3);
	zmm4 = _mm512_load_pd(&dChristoffel_tildeNCP[1][2][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxz, zmm11, zmm4);
	zmm4 = _mm512_fmadd_pd(igyz, zmm13, zmm4);
	zmm5 = _mm512_load_pd(&dChristoffel_tildeNCP[1][2][2][2][0]);
	zmm5 = _mm512_fmadd_pd(igxz, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igxz, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igzz, zmm15, zmm5);
	zmm5 = _mm512_fmadd_pd(igyz, zmm14, zmm5);
	zmm5 = _mm512_fmadd_pd(igyz, zmm14, zmm5);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][2][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][2][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][2][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][2][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][2][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[1][2][2][2][0], zmm5);
	/* we are now on the last slice and Christoffel_tilde is complete. Immediately use this result to compute Gtilde^2
	zmm6 = _mm512_mul_pd(igxx, zmm0);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igyy, zmm3, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igzz, zmm5, zmm6);
	_mm512_store_pd(&Gtilde[2][0], zmm6); */

	// **** 3/3 **** 
	// load scalar ( -1.0 ) and multiply with igxx, igxy and igxz   
	xmm = _mm_load1_pd(&s2[0]);
	zmm0 = _mm512_broadcastsd_pd(xmm);
	zmm17 = _mm512_mul_pd(igxx, zmm0);
	zmm18 = _mm512_mul_pd(igxy, zmm0);
	zmm19 = _mm512_mul_pd(igxz, zmm0);
#ifdef COUNT_FLOPS  
	nMul = nMul + VECTORLENGTH * 3 * 3;   // this is counted for all three slices. 
#endif 

										  // load first slice of dDDsymm_klij  (D_1ij) 
	zmm10 = _mm512_load_pd(&dDDsymm[0][2][0][0][0]);
	zmm11 = _mm512_load_pd(&dDDsymm[0][2][0][1][0]);
	zmm12 = _mm512_load_pd(&dDDsymm[0][2][0][2][0]);
	zmm13 = _mm512_load_pd(&dDDsymm[0][2][1][1][0]);
	zmm14 = _mm512_load_pd(&dDDsymm[0][2][1][2][0]);
	zmm15 = _mm512_load_pd(&dDDsymm[0][2][2][2][0]);

	// for dChristoffel tildeNCP, multiply with g_contr 
	zmm0 = _mm512_mul_pd(igxx, zmm10);
	zmm0 = _mm512_fmadd_pd(igxy, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igxy, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igxz, zmm12, zmm0);
	zmm0 = _mm512_fmadd_pd(igxz, zmm12, zmm0);
	zmm1 = _mm512_mul_pd(igxy, zmm13);
	zmm1 = _mm512_fmadd_pd(igxz, zmm14, zmm1);
	zmm2 = _mm512_mul_pd(igxy, zmm14);
	zmm2 = _mm512_fmadd_pd(igxz, zmm15, zmm2);
	zmm3 = _mm512_mul_pd(zmm17, zmm13);
	zmm4 = _mm512_mul_pd(zmm17, zmm14);
	zmm5 = _mm512_mul_pd(zmm17, zmm15);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][0][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][0][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][0][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][0][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][0][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][0][2][2][0], zmm5);
	// 
	zmm0 = _mm512_mul_pd(igxy, zmm10);
	zmm0 = _mm512_fmadd_pd(igyy, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igyy, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igyz, zmm12, zmm0);
	zmm0 = _mm512_fmadd_pd(igyz, zmm12, zmm0);
	zmm1 = _mm512_mul_pd(igyy, zmm13);
	zmm1 = _mm512_fmadd_pd(igyz, zmm14, zmm1);
	zmm2 = _mm512_mul_pd(igyy, zmm14);
	zmm2 = _mm512_fmadd_pd(igyz, zmm15, zmm2);
	zmm3 = _mm512_mul_pd(zmm18, zmm13);
	zmm4 = _mm512_mul_pd(zmm18, zmm14);
	zmm5 = _mm512_mul_pd(zmm18, zmm15);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][1][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][1][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][1][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][1][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][1][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][1][2][2][0], zmm5);
	// 
	zmm0 = _mm512_mul_pd(igxz, zmm10);
	zmm0 = _mm512_fmadd_pd(igyz, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igyz, zmm11, zmm0);
	zmm0 = _mm512_fmadd_pd(igzz, zmm12, zmm0);
	zmm0 = _mm512_fmadd_pd(igzz, zmm12, zmm0);
	zmm1 = _mm512_mul_pd(igyz, zmm13);
	zmm1 = _mm512_fmadd_pd(igzz, zmm14, zmm1);
	zmm2 = _mm512_mul_pd(igyz, zmm14);
	zmm2 = _mm512_fmadd_pd(igzz, zmm15, zmm2);
	zmm3 = _mm512_mul_pd(zmm19, zmm13);
	zmm4 = _mm512_mul_pd(zmm19, zmm14);
	zmm5 = _mm512_mul_pd(zmm19, zmm15);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][2][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][2][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][2][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][2][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][2][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][2][2][2][0], zmm5);

	// load scalar ( -1.0 ) and multiply with igxx, igxy and igxz   
	xmm = _mm_load1_pd(&s2[0]);
	zmm0 = _mm512_broadcastsd_pd(xmm);
	zmm17 = _mm512_mul_pd(igxy, zmm0);
	zmm18 = _mm512_mul_pd(igyy, zmm0);
	zmm19 = _mm512_mul_pd(igyz, zmm0);

	// load second slice of dDDsymm_klij 
	zmm10 = _mm512_load_pd(&dDDsymm[1][2][0][0][0]);
	zmm11 = _mm512_load_pd(&dDDsymm[1][2][0][1][0]);
	zmm12 = _mm512_load_pd(&dDDsymm[1][2][0][2][0]);
	zmm13 = _mm512_load_pd(&dDDsymm[1][2][1][1][0]);
	zmm14 = _mm512_load_pd(&dDDsymm[1][2][1][2][0]);
	zmm15 = _mm512_load_pd(&dDDsymm[1][2][2][2][0]);

	// for Christoffel tilde, load the result of the previous slice, multiply DD properly with g_contr and add contribution 
	zmm0 = _mm512_load_pd(&dChristoffel_tildeNCP[2][0][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm17, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&dChristoffel_tildeNCP[2][0][0][1][0]);
	zmm1 = _mm512_fmadd_pd(igxx, zmm10, zmm1);
	zmm1 = _mm512_fmadd_pd(igxz, zmm12, zmm1);
	zmm2 = _mm512_load_pd(&dChristoffel_tildeNCP[2][0][0][2][0]);
	zmm2 = _mm512_fmadd_pd(zmm17, zmm12, zmm2);
	zmm3 = _mm512_load_pd(&dChristoffel_tildeNCP[2][0][1][1][0]);
	zmm3 = _mm512_fmadd_pd(igxx, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igxx, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igxy, zmm13, zmm3);
	zmm3 = _mm512_fmadd_pd(igxz, zmm14, zmm3);
	zmm3 = _mm512_fmadd_pd(igxz, zmm14, zmm3);
	zmm4 = _mm512_load_pd(&dChristoffel_tildeNCP[2][0][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxx, zmm12, zmm4);
	zmm4 = _mm512_fmadd_pd(igxz, zmm15, zmm4);
	zmm5 = _mm512_load_pd(&dChristoffel_tildeNCP[2][0][2][2][0]);
	zmm5 = _mm512_fmadd_pd(zmm17, zmm15, zmm5);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][0][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][0][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][0][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][0][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][0][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][0][2][2][0], zmm5);
	// 
	zmm0 = _mm512_load_pd(&dChristoffel_tildeNCP[2][1][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm18, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&dChristoffel_tildeNCP[2][1][0][1][0]);
	zmm1 = _mm512_fmadd_pd(igxy, zmm10, zmm1);
	zmm1 = _mm512_fmadd_pd(igyz, zmm12, zmm1);
	zmm2 = _mm512_load_pd(&dChristoffel_tildeNCP[2][1][0][2][0]);
	zmm2 = _mm512_fmadd_pd(zmm18, zmm12, zmm2);
	zmm3 = _mm512_load_pd(&dChristoffel_tildeNCP[2][1][1][1][0]);
	zmm3 = _mm512_fmadd_pd(igxy, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igxy, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igyy, zmm13, zmm3);
	zmm3 = _mm512_fmadd_pd(igyz, zmm14, zmm3);
	zmm3 = _mm512_fmadd_pd(igyz, zmm14, zmm3);
	zmm4 = _mm512_load_pd(&dChristoffel_tildeNCP[2][1][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxy, zmm12, zmm4);
	zmm4 = _mm512_fmadd_pd(igyz, zmm15, zmm4);
	zmm5 = _mm512_load_pd(&dChristoffel_tildeNCP[2][1][2][2][0]);
	zmm5 = _mm512_fmadd_pd(zmm18, zmm15, zmm5);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][1][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][1][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][1][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][1][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][1][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][1][2][2][0], zmm5);
	// 
	zmm0 = _mm512_load_pd(&dChristoffel_tildeNCP[2][2][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm19, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&dChristoffel_tildeNCP[2][2][0][1][0]);
	zmm1 = _mm512_fmadd_pd(igxz, zmm10, zmm1);
	zmm1 = _mm512_fmadd_pd(igzz, zmm12, zmm1);
	zmm2 = _mm512_load_pd(&dChristoffel_tildeNCP[2][2][0][2][0]);
	zmm2 = _mm512_fmadd_pd(zmm19, zmm12, zmm2);
	zmm3 = _mm512_load_pd(&dChristoffel_tildeNCP[2][2][1][1][0]);
	zmm3 = _mm512_fmadd_pd(igxz, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igxz, zmm11, zmm3);
	zmm3 = _mm512_fmadd_pd(igyz, zmm13, zmm3);
	zmm3 = _mm512_fmadd_pd(igzz, zmm14, zmm3);
	zmm3 = _mm512_fmadd_pd(igzz, zmm14, zmm3);
	zmm4 = _mm512_load_pd(&dChristoffel_tildeNCP[2][2][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxz, zmm12, zmm4);
	zmm4 = _mm512_fmadd_pd(igzz, zmm15, zmm4);
	zmm5 = _mm512_load_pd(&dChristoffel_tildeNCP[2][2][2][2][0]);
	zmm5 = _mm512_fmadd_pd(zmm19, zmm15, zmm5);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][2][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][2][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][2][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][2][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][2][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][2][2][2][0], zmm5);


	// load scalar ( -1.0 ) and multiply with igxx, igxy and igxz   
	xmm = _mm_load1_pd(&s2[0]);
	zmm0 = _mm512_broadcastsd_pd(xmm);
	zmm17 = _mm512_mul_pd(igxz, zmm0);
	zmm18 = _mm512_mul_pd(igyz, zmm0);
	zmm19 = _mm512_mul_pd(igzz, zmm0);

	// load third slice of D_kij  
	zmm10 = _mm512_load_pd(&dDDsymm[2][2][0][0][0]);
	zmm11 = _mm512_load_pd(&dDDsymm[2][2][0][1][0]);
	zmm12 = _mm512_load_pd(&dDDsymm[2][2][0][2][0]);
	zmm13 = _mm512_load_pd(&dDDsymm[2][2][1][1][0]);
	zmm14 = _mm512_load_pd(&dDDsymm[2][2][1][2][0]);
	zmm15 = _mm512_load_pd(&dDDsymm[2][2][2][2][0]);

	// for Christoffel tilde, load the result of the previous slice, multiply DD properly with g_contr and add contribution 
	zmm0 = _mm512_load_pd(&dChristoffel_tildeNCP[2][0][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm17, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&dChristoffel_tildeNCP[2][0][0][1][0]);
	zmm1 = _mm512_fmadd_pd(zmm17, zmm11, zmm1);
	zmm2 = _mm512_load_pd(&dChristoffel_tildeNCP[2][0][0][2][0]);
	zmm2 = _mm512_fmadd_pd(igxx, zmm10, zmm2);
	zmm2 = _mm512_fmadd_pd(igxy, zmm11, zmm2);
	zmm3 = _mm512_load_pd(&dChristoffel_tildeNCP[2][0][1][1][0]);
	zmm3 = _mm512_fmadd_pd(zmm17, zmm13, zmm3);
	zmm4 = _mm512_load_pd(&dChristoffel_tildeNCP[2][0][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxx, zmm11, zmm4);
	zmm4 = _mm512_fmadd_pd(igxy, zmm13, zmm4);
	zmm5 = _mm512_load_pd(&dChristoffel_tildeNCP[2][0][2][2][0]);
	zmm5 = _mm512_fmadd_pd(igxx, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igxx, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igxz, zmm15, zmm5);
	zmm5 = _mm512_fmadd_pd(igxy, zmm14, zmm5);
	zmm5 = _mm512_fmadd_pd(igxy, zmm14, zmm5);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][0][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][0][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][0][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][0][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][0][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][0][2][2][0], zmm5);
	/* we are now on the last slice and Christoffel_tilde is complete. Immediately use this result to compute Gtilde^1
	zmm6 = _mm512_mul_pd(igxx, zmm0);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igyy, zmm3, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igzz, zmm5, zmm6);
	_mm512_store_pd(&Gtilde[0][0], zmm6); */
	// 
	zmm0 = _mm512_load_pd(&dChristoffel_tildeNCP[2][1][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm18, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&dChristoffel_tildeNCP[2][1][0][1][0]);
	zmm1 = _mm512_fmadd_pd(zmm18, zmm11, zmm1);
	zmm2 = _mm512_load_pd(&dChristoffel_tildeNCP[2][1][0][2][0]);
	zmm2 = _mm512_fmadd_pd(igxy, zmm10, zmm2);
	zmm2 = _mm512_fmadd_pd(igyy, zmm11, zmm2);
	zmm3 = _mm512_load_pd(&dChristoffel_tildeNCP[2][1][1][1][0]);
	zmm3 = _mm512_fmadd_pd(zmm18, zmm13, zmm3);
	zmm4 = _mm512_load_pd(&dChristoffel_tildeNCP[2][1][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxy, zmm11, zmm4);
	zmm4 = _mm512_fmadd_pd(igyy, zmm13, zmm4);
	zmm5 = _mm512_load_pd(&dChristoffel_tildeNCP[2][1][2][2][0]);
	zmm5 = _mm512_fmadd_pd(igxy, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igxy, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igyz, zmm15, zmm5);
	zmm5 = _mm512_fmadd_pd(igyy, zmm14, zmm5);
	zmm5 = _mm512_fmadd_pd(igyy, zmm14, zmm5);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][1][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][1][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][1][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][1][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][1][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][1][2][2][0], zmm5);
	/* we are now on the last slice and Christoffel_tilde is complete. Immediately use this result to compute Gtilde^2
	zmm6 = _mm512_mul_pd(igxx, zmm0);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igyy, zmm3, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igzz, zmm5, zmm6);
	_mm512_store_pd(&Gtilde[1][0], zmm6); */

	// 
	zmm0 = _mm512_load_pd(&dChristoffel_tildeNCP[2][2][0][0][0]);
	zmm0 = _mm512_fmadd_pd(zmm19, zmm10, zmm0);
	zmm1 = _mm512_load_pd(&dChristoffel_tildeNCP[2][2][0][1][0]);
	zmm1 = _mm512_fmadd_pd(zmm19, zmm11, zmm1);
	zmm2 = _mm512_load_pd(&dChristoffel_tildeNCP[2][2][0][2][0]);
	zmm2 = _mm512_fmadd_pd(igxz, zmm10, zmm2);
	zmm2 = _mm512_fmadd_pd(igyz, zmm11, zmm2);
	zmm3 = _mm512_load_pd(&dChristoffel_tildeNCP[2][2][1][1][0]);
	zmm3 = _mm512_fmadd_pd(zmm19, zmm13, zmm3);
	zmm4 = _mm512_load_pd(&dChristoffel_tildeNCP[2][2][1][2][0]);
	zmm4 = _mm512_fmadd_pd(igxz, zmm11, zmm4);
	zmm4 = _mm512_fmadd_pd(igyz, zmm13, zmm4);
	zmm5 = _mm512_load_pd(&dChristoffel_tildeNCP[2][2][2][2][0]);
	zmm5 = _mm512_fmadd_pd(igxz, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igxz, zmm12, zmm5);
	zmm5 = _mm512_fmadd_pd(igzz, zmm15, zmm5);
	zmm5 = _mm512_fmadd_pd(igyz, zmm14, zmm5);
	zmm5 = _mm512_fmadd_pd(igyz, zmm14, zmm5);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][2][0][0][0], zmm0);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][2][0][1][0], zmm1);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][2][0][2][0], zmm2);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][2][1][1][0], zmm3);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][2][1][2][0], zmm4);
	_mm512_store_pd(&dChristoffel_tildeNCP[2][2][2][2][0], zmm5);
	/* we are now on the last slice and Christoffel_tilde is complete. Immediately use this result to compute Gtilde^2
	zmm6 = _mm512_mul_pd(igxx, zmm0);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxy, zmm1, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igxz, zmm2, zmm6);
	zmm6 = _mm512_fmadd_pd(igyy, zmm3, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igyz, zmm4, zmm6);
	zmm6 = _mm512_fmadd_pd(igzz, zmm5, zmm6);
	_mm512_store_pd(&Gtilde[2][0], zmm6); */

#ifdef COUNT_FLOPS  
	// cost for the computation of dChristoffel_tilde 
	nMul = nMul + VECTORLENGTH * 12 * 3 * 3 * 3;
	nAdd = nAdd + VECTORLENGTH * (12 * 3 * 3 - 6 * 3)*3;
#endif 

	// dChristoffel_tildeNCP(:, k, i, j, m) = dChristoffel_tildeNCP(:, k, i, j, m) + g_contr(:, m, l)*(dDDSymm(:, k, i, j, l) + dDDSymm(:, k, j, i, l) - dDDSymm(:, k, l, i, j))
	
	asm(nop); 


#else

#endif 

}

