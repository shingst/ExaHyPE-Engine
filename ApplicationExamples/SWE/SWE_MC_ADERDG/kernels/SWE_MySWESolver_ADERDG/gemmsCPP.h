
#ifndef EXAHYPE_SWE_MYSWESOLVER_ADERDG_KERNELS_ADERDG_GEMMSCPP_H_
#define EXAHYPE_SWE_MYSWESOLVER_ADERDG_KERNELS_ADERDG_GEMMSCPP_H_


namespace SWE {
namespace MySWESolver_ADERDG_kernels {
namespace aderdg {

void gemm_4_2_2_face_Q_x(const double* A, const double* B, double* C); 
void gemm_4_2_2_face_F_x(const double* A, const double* B, double* C); 
void gemm_4_2_2_volume_x(const double* A, const double* B, double* C); 
void gemm_4_2_2_volume_y(const double* A, const double* B, double* C); 
void gemm_4_2_2_volume_y_add(const double* A, const double* B, double* C); 
void gemm_4_2_2_rhs_x(const double* A, const double* B, double* C); 
void gemm_4_2_2_rhs_y(const double* A, const double* B, double* C); 
void gemm_4_2_2_lduh_x(const double* A, const double* B, double* C); 
void gemm_4_2_2_lduh_y(const double* A, const double* B, double* C); 
void gemm_4_2_2_gradQ_x(const double* A, const double* B, double* C); 
void gemm_4_2_2_gradQ_y(const double* A, const double* B, double* C); 
void gemm_4_2_2_lqi(const double* A, const double* B, double* C); 
void gemm_4_3_2_dg2fv_x(const double* A, const double* B, double* C); 
void gemm_4_3_2_dg2fv_y(const double* A, const double* B, double* C); 
void gemm_4_2_3_fv2dg_x(const double* A, const double* B, double* C); 
void gemm_4_2_3_fv2dg_y(const double* A, const double* B, double* C); 
void gemm_4_2_2_uh2lob_x(const double* A, const double* B, double* C); 
void gemm_4_2_2_uh2lob_y_slice(const double* A, const double* B, double* C); 
void gemm_4_3_2_dg2fv_y_slice(const double* A, const double* B, double* C); 

}
}
}


#endif //EXAHYPE_SWE_MYSWESOLVER_ADERDG_KERNELS_ADERDG_GEMMSCPP_H_