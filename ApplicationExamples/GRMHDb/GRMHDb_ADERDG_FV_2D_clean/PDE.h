// PDE.h

#ifndef __EXAHYPE_USER_PDE__
#define __EXAHYPE_USER_PDE__

#ifndef __C2P_MHD_CPP__
#define __C2P_MHD_CPP__



// Fortran functions:
extern "C" {
//void minimumtreedepth_(int* depth);
//void hastoadjustsolution_(double* t, bool* refine);
//void adjustedsolutionvalues_(const double* const x,const double* w,const double* t,const double* dt,double* Q);
void pdeprim2cons_(double* /* OUT */ Q, double* /* IN */ V);
void pdecons2prim_(double* /* OUT */ V, double* /* IN */ Q, int* iErr);
void pdeflux_(double* Fx, double* Fy, double* Fz, const double* const Q);
//void pdeflux_(double* const F, const double* const Q);
void pdesource_(double* S, const double* const Q);
void pdencp_(double* BgradQ, const double* const Q, const double* const gradQ);
void pdeeigenvectors_(double* R, double* lambda, double* iR, const double* const Q, const double* const nv);
void pdeeigenvalues_(double* lambda, const double* const Q, double* nv);
void pdevarname_(char* MyNameOUT, int* ind);
void pdeauxname_(char* AuxName, int* ind);
void pdeauxvar_(double* aux, double* V, double* x);
void hllemfluxfv_(double* FL, double* FR, const double* const  QL, const double* const  QR, const double* const  QavL, const double* const  QavR, const int* normalNonZeroIndex);
//void pdeeigenvectors_(double* L, double* R, double* iR, const double* const Q, const double* const nv);
//void hllemriemannsolverS_(const int* basisSize, const int* normalNonZeroIndex, double* FL, double* FR, const double* const  QL, const double* const  QR, const double* const  QavL, const double* const  QavR);

//void registerinitialdata_(const char* const id_name, int* id_name_len);
//void getnumericalsolution_(double* V,double* Q);
//void getexactsolution_(double* V,double* pos,double* timeStamp);
//void inittecplot_(const int* N_in,const int* M_in);
}/* extern "C" */

#endif /* __C2P_MHD_CPP__ */

#endif /* __EXAHYPE_USER_PDE__ */
