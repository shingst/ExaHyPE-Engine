
// Independently:


#if defined(Dim2)
#define DIMENSIONS 2
#elif defined(Dim3)
#define DIMENSIONS 3
#else
#error You must define Dim2 or Dim3
#endif

#include <stdio.h>
#include "InitialData/InitialData.h"
#include "InitialData/AnalyticalID.h"
#include <cmath>
#include <algorithm>
using namespace std;
using namespace SVEC;

#include "PDE/PDE-GRMHD-ExaHyPE.h"
using namespace GRMHD;
using namespace ExaHyPEAdapter;

#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>

constexpr double SNAN = std::numeric_limits<double>::signaling_NaN();

// Fortran directly here to avoid double inclusion from similar exahype projects
extern "C" {
void alfenwave_(const double* x,const double* const t,  double* Q);
//void pdeprim2cons_(double* /* OUT */ Q, double* /* IN */ V);
//void pdecons2prim_(double* /* OUT */ V, double* /* IN */ Q, int* iErr);
void minimumtreedepth_(int* depth);
void hastoadjustsolution_(double* t, bool* refine);
void adjustedsolutionvalues_(const double* const x,const double* w,const double* t,const double* dt,double* Q);
void pdeflux_(double* Fx, double* Fy, double* Fz, const double* const Q);
void pdesource_(double* S, const double* const Q);
void pdeeigenvalues_(double* lambda, const double* const Q, double* nv);
void registerinitialdata_(const char* const id_name, int* id_name_len);

void pdencp_(double* BgradQ, const double* const Q, const double* const gradQ);
void pdematrixb_(double* Bn, const double* const Q, double* nv);

void pdecons2prim_(double* V, const double* const Q, int* iErr);
void pdeprim2cons_(double* Q, const double* const V);
}

int vacuumtest() {
	double x[3] = {0};
	double t = 0;
	constexpr int nCheck = 19;
	constexpr int nSize = 30;
	constexpr double eps = 1e-7;
	double QC[nSize];

	VacuumInitialData Q(QC);
	for(int j=0; j<nCheck; j++) {
		printf("QC[%d] = %e\n", j, QC[j]);
	}
}

constexpr double eps = 5e-9; // todo: raise
constexpr int nSize = 300; // array size
constexpr int nCheck = 19; // hydro vars
#define CHECK(A,B) check(#A, A, #B, B)
#define P(A) ((std::abs(A) < 1e-20) ? 0.0 : A)
bool check(const char* const labelA, const double* const A, const char* const labelB, const double* const B) {
	bool globalDiffer = false;
	for(int j=0; j<nCheck; j++) {
		double difference = std::abs(A[j] - B[j]);
		bool differ = ( difference > eps );
		if(differ) globalDiffer = true;
		printf("[%2d]: %s=%15e %s=%15e %s\t(diff=%e)\n", j, labelA, P(A[j]), labelB, P(B[j]), differ ? "differ" : "SAME", difference);
	}
	return globalDiffer;
}

void readArray(double* target, const std::string& text) {
	stringstream input(text);

	int counter=0;
	double current_number = 0;
	while (input >> current_number) {
		target[counter++] = current_number;
	}

	cout << "Read " << counter << " Numbers: ";
        for (int count = 0; count < counter; count++){
            cout << target[count] << " ";
        }
        cout << "\n";
}



// Does derivatives in three directions
template <typename IDFunc>
void interp_deriv(const double* const xc, const double t, double **gradQ) {
	constexpr double epsilon = 1e-7, eps4 = 1e-4;
	constexpr int stencil=4; // stencil size	
	double Q[stencil][nSize];
	// for debugging purposes, find errors in algorithm:
	for(int j=0;j<stencil;j++) for(int i=0;i<nSize;i++) Q[j][i] = NAN;
	vec::stored<3> x[stencil];
	constexpr int p=0, m=1, pp=2, mm=3;
	// Stencil indices: {p,m,pp,mm}]
	
	
	// Metric derivative computed with a fourth order central finite difference 
	// as done by Michael D
	double *Qd;
	for(int d=0; d<3; d++) {
		DFOR(i) for(int j=0;j<stencil;j++) x[j](i) = xc[i];

		x[p](d) += eps4;
		x[m](d) -= eps4;
		x[pp](d) += 2*eps4;
		x[mm](d) -= 2*eps4;
		
		for(int j=0; j<stencil;j++) IDFunc(x[j].data, t, Q[j]);
		
		for(int i=0; i<nSize; i++) {
			gradQ[d][i] = ( 8.0*Q[p][i] - 8.0*Q[m][i]  + Q[mm][i]   - Q[pp][i]  )/(12.0*eps4);
		}
	}
}

// enable nan tracker
#include <fenv.h>
#include <../../home/sven/numrel/exahype/Engine-ExaHyPE/ApplicationExamples/GRMHD/GRMHD_cpp/InitialData/AnalyticalID.h>


extern "C" {
void cmain_() {
	
	printf("Inditest, compiled for DIMENSIONS=%d\n", DIMENSIONS);
	
	// try NaN catcher
	feenableexcept(FE_INVALID | FE_OVERFLOW);  // Enable all floating point exceptions but FE_INEXACT
	
	double x[3] = {0};
	double t = 0;
	double dx = 0.3;
	int N = 1;
	
	double QC[nSize], QQC[nSize], QF[nSize], VC[nSize], VF[nSize];
	double VVC[nSize], VVF[nSize]; // Cons2Prim test
	fill_n(QC, nSize, 137);
	fill_n(QQC, nSize, 137);
	fill_n(QF, nSize, 137);
	fill_n(VC, nSize, 137);
	fill_n(VF, nSize, 137);
	InitialState q(QC); // NaN default
	
	for(int i=0; i<N; i++) {
		x[1] += dx;
		constexpr bool doC=true;
		constexpr bool doF=true;
		
		printf("Generating ID at x=[%f %f %f].\n", x[0], x[1], x[2]);
		
		if(doC) AlfenWave id(x,t,VC); // gives prims
		if(doF) alfenwave_(x,&t,VF);  // probably gives cons
		
		if(doC&&doF) printf("Initial Data (Primitives) from Fortran and C:\n");
		if(doC&&doF) CHECK(VC, VF);
		
		// understand whether raising works correctly
		//Hydro::Conserved::ConstShadowExtendable
		//public Magneto::ConstShadowExtendable, 
		/*
		ADMBase::Full adm(VC);
		Hydro::Primitives::ShadowExtendable prim(VC);
		Lo<vec::stored3> myvel;
		myvel.lo = 12345;
		*/
		
		//myvel.lo=0; DFOR(i) CONTRACT(j) myvel.lo(i) += prim.vel.up(j) * adm.gam.lo(i,j);
		//adm.gam.lower_vec(prim.vel, myvel);
		/*
		for(int i=0;i<6;i++) printf("gmunu.lo(%d)=%e\n",i,adm.gam.lo.data[i]);
		DFOR(i) DFOR(j) printf("g.lo(%d,%d)=%e\n", i,j, adm.gam.lo(i,j));
		for(int i=0;i<3; i++) printf(" prim.vel.up(%d)=%e\n",i,prim.vel.up(i));
		for(int i=0;i<3;i++) printf("myvel.vel.lo(%d)=%e\n",i,myvel.lo(i));
		
		return -1;
		*/
		
		// get the fortran cons from the gfortran executable
		std::string FortCons = " 1.1241719689735967       0.45685025174785676        2.1089451544283881      -0.58620660444159478        2.7714719547653637        1.0000000000000000      -0.96347212472303134       0.26780863481539874        0.0000000000000000        1.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        1.0000000000000000        0.0000000000000000        0.0000000000000000        1.0000000000000000        0.0000000000000000        1.0000000000000000";
		readArray(QF, FortCons);
		//printf("C vs F cons\n");
		//CHECK(QC,QF);
		
		// Do the c2p
		//if(doF) pdeprim2cons_(QF, VF);
		if(doC) GRMHD::Prim2Cons(QC, VC).copyFullStateVector();
		if(doF) { printf("Compare Fortran PRIM <-> Fortran COONS:\n"); CHECK(QF, VF); }
		
		printf("Orig prims and cons:\n");
		AlfenWaveCons id(x,t,QQC);
		CHECK(QQC,QC);

		
		// do again the comparison
		if(doC&&doF) printf("Prim2Cons:\n");
		if(doC&&doF) CHECK(QC, QF);
		
		// do the conversion back
		if(doF) { int iErr; pdecons2prim_(VVF, QF, &iErr); }
		if(doC) GRMHD::Cons2Prim(VVC, QC).copyFullStateVector();
		
		// Compare again:
		if(doC&&doF) printf("Cons2Prim:\n");
		if(doC&&doF) CHECK(VVC, VVF);
		
		if(doC) printf("Cons2Prim C outcome with original C prims (VERY GOOD):\n");
		if(doC) CHECK(VC, VVC);
		if(doF) printf("Cons2Prim F outcome with original F prims (BAD due to compiler bug):\n");
		if(doF) CHECK(VF, VVF);
		
		// Fluxes in one test direction
		constexpr int dir = 0;
		double FC[3][nSize], FF[3][nSize];
		//double *FCp[3] = {FC[0],FC[1],FC[2]};
		// to detect non set values
		std::fill_n(FC[0], 3*nSize, 123456789);
		std::fill_n(FF[0], 3*nSize, 123456789);

		//PDE(QC).flux(FCp);
		flux(QC, FC[0], FC[1], FC[2]).zeroMaterialFluxes();
		// pdeflux_(FF[0], FF[1], FF[2], QF); // will not work due to broken F C2P
		
		readArray(FF[0], std::string("1.8880309866282984E-012   1.1043560762900662       0.96347212472583399      -0.26780863481617778       0.45685025174596872        0.0000000000000000      -0.44016248272562447       0.12234844223395951        1.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000"));
		readArray(FF[1], std::string("0.49481832487215294       0.96347212472583399        1.8956439237513409        6.2438387793406491E-012   1.6141268295562352       0.44016248272562447        0.0000000000000000       -1.3877787807814457E-017 -0.96347212472303134        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000"));
		readArray(FF[2], std::string("-0.13754068920649698      -0.26780863481617778        6.2438942904918804E-012   1.8956439237720684      -0.44866591523509780      -0.12234844223395951        1.3877787807814457E-017   0.0000000000000000       0.26780863481539874        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000 "));

		DFOR(dir) {
		printf("Fluxes in direction dir=%i\n", dir);
		CHECK(FC[dir], FF[dir]);
		}
		
		// Source terms check
		double SC[nSize], SF[nSize];
		std::fill_n(SC, nSize, 123456789);
		double gradQC[3*nSize], gradQF[3*nSize];
		double *QkC[]={gradQC+0, gradQC+nSize, gradQC+2*nSize};
		double *QkF[]={gradQF+0, gradQF+nSize, gradQF+2*nSize};
		interp_deriv<AlfenWaveCons>(x, t, QkC);
		
		readArray(QkF[0], std::string("0.0000000000000000       -9.2518585385429703E-014   3.6832447239990942        13.250893207950771        0.0000000000000000        0.0000000000000000       -1.6826912794078626       -6.0536738979360027        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000"));
		readArray(QkF[1], std::string("0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000"));
		readArray(QkF[2], std::string("0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000"));
		
		DFOR(dir){
		printf("Derivative in direction dir=%i:\n",dir);
		CHECK(QkC[dir],QkF[dir]);
		}

		fusedSource(QC, QkC[0],QkC[1],QkC[2], SC).zero_adm(); // TODO: Zero Source terms for material parameters!
		//GRMHD::Shadow(SC).zero_adm();
		
		readArray(SF, std::string("-0.0000000000000000       -53.205687704723701       -0.0000000000000000       -0.0000000000000000       -6.0536738979362745       -0.0000000000000000       -0.0000000000000000       -0.0000000000000000       -0.0000000000000000       -0.0000000000000000       -0.0000000000000000       -0.0000000000000000       -0.0000000000000000       -0.0000000000000000       -0.0000000000000000       -0.0000000000000000       -0.0000000000000000       -0.0000000000000000       -0.0000000000000000"));
		
		printf("FusedSource terms:\n");
		CHECK(SC,SF);
		
		double QncpC[nSize];
		std::fill_n(QncpC, nSize, NAN);
		double *QkCpointer = QkC[0];
		// const Gradients g(QkCpointer); // this does not work as nSize >> nVar
		nonConservativeProduct(QC, QkC[0], QkC[1], QkC[2], QncpC).zero_adm();
		printf("NCP terms:\n");
		CHECK(QncpC,SC);
		
	}
}
} // ext C

int main() { cmain_(); }
