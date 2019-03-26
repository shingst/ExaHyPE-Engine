
#ifndef _EXAHYPE_SWE_MYSWESOLVER_ADERDG_KERNELS_ADERDG_QUADRATURE_H_
#define _EXAHYPE_SWE_MYSWESOLVER_ADERDG_KERNELS_ADERDG_QUADRATURE_H_

// Use Gauss-Legendre quadrature

namespace SWE {
namespace MySWESolver_ADERDG_kernels {
namespace aderdg {

void initQuadratureNodesAndWeights();
void freeQuadratureNodesAndWeights();

extern double *nodes;
extern double *weights1;
extern double *weights2;
extern double *weights3;
extern double *iweights3;

extern double* uh2lob;
extern double* dg2fv;
extern double* fv2dg;

}
}
}
#endif /* _EXAHYPE_SWE_MYSWESOLVER_ADERDG_KERNELS_ADERDG_QUADRATURE_H_ */