
#ifndef EXAHYPE_SWE_MYSWESOLVER_ADERDG_KERNELS_ADERDG_CONVERTER_H_
#define EXAHYPE_SWE_MYSWESOLVER_ADERDG_KERNELS_ADERDG_CONVERTER_H_
#include <cstring>

#include "kernels/KernelUtils.h"

namespace SWE {
namespace MySWESolver_ADERDG_kernels {
namespace aderdg {
namespace converter {

  constexpr int getLduhGenArraySize()    { return 16;}
  constexpr int getLduhOptArraySize()    { return 16;}

  constexpr int getQhbndGenArraySize()   { return 32;}
  constexpr int getQhbndOptArraySize()   { return 32;}
  
  constexpr int getFhbndGenArraySize()   { return 32;}
  constexpr int getFhbndOptArraySize()   { return 32;}

  constexpr int getQFaceGenArraySize()   { return 8;}
  constexpr int getQFaceOptArraySize()   { return 8;}
  
  constexpr int getFFaceGenArraySize()   { return 8;}
  constexpr int getFFaceOptArraySize()   { return 8;}

  constexpr int getQhiGenArraySize()     { return 16;}
  constexpr int getQhiOptArraySize()     { return 16;}
  
  constexpr int getFhiGenArraySize()     { return 48;}
  constexpr int getFhiOptArraySize()     { return 48;}
  
  constexpr int getPSiGenArraySize()     { return 48;}
  constexpr int getPSiOptArraySize()     { return 48;}


  void Qhbnd_optimised2generic(const double* const opt, double* gen);
  void Qhbnd_generic2optimised(const double* const gen, double* opt);
  
  void Fhbnd_optimised2generic(const double* const opt, double* gen);
  void Fhbnd_generic2optimised(const double* const gen, double* opt);
  
  void QFace_optimised2generic(const double* const opt, double* gen);
  void QFace_generic2optimised(const double* const gen, double* opt);
  
  void FFace_optimised2generic(const double* const opt, double* gen);
  void FFace_generic2optimised(const double* const gen, double* opt);

  void Qhi_optimised2generic(const double* const opt, double* gen);
  void Qhi_generic2optimised(const double* const gen, double* opt);
  
  void Fhi_optimised2generic(const double* const opt, double* gen);
  void Fhi_generic2optimised(const double* const gen, double* opt);
    
  void PSi_optimised2generic(const double* const opt, double* gen);
  void PSi_generic2optimised(const double* const gen, double* opt);
  
  void Lduh_optimised2generic(const double* const opt, double* gen);
  void Lduh_generic2optimised(const double* const gen, double* opt);
}
}
}
}

#endif //EXAHYPE_SWE_MYSWESOLVER_ADERDG_KERNELS_ADERDG_CONVERTER_H_
