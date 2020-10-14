#ifndef __InitialData_CLASS_HEADER__
#define __InitialData_CLASS_HEADER__

#if defined(USE_ASAGI)
namespace easi {
  class YAMLParser;
  class Component;
  class ArraysAdapter;
};
class AsagiReader;
#endif

class InitialData {
 public:
  InitialData();
  InitialData(int scenario);  
  void getInitialData(const double* const x,double* Q);
  
 private:
  int scenario;
#if defined(USE_ASAGI)
  easi::YAMLParser* parser;
  easi::Component* model;
  AsagiReader* asagiReader;
#endif

  void ShockShockProblem(const double* const x,double* Q);
  void RareRareProblem(const double* const x,double* Q);
  void GaussFunctionProblem(const double* const x,double* Q);
  void ExpBreakProblem(const double* const x,double* Q);
  void DamBreakProblem(const double* const x,double* Q);
  void SeaAtRestProblem(const double* const x,double* Q);
  void SteadyRunUpLinear(const double* const x,double* Q);
  void RunUpLinear(const double* const x,double* Q);
  void SteadyRunUpShelf(const double* const x,double* Q);
  void RunUpShelf(const double* const x,double* Q);
  void WettingDryingProblem(const double* const x, double* Q);
  void RunUpTest(const double* const x, double* Q);
  void OscillatingLake(const double* const x, double* Q);
  void SolitaryWaveOnSimpleBeach(const double*const x, double* Q);
  void readAsagiData(const double* const x,double* Q);
};

#endif // __InitialData_CLASS_HEADER__
