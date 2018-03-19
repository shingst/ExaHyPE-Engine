#ifndef __InitialData_CLASS_HEADER__
#define __InitialData_CLASS_HEADER__

namespace SWE {
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
void initialData(const double* const x,double* Q);

}

#endif // __InitialData_CLASS_HEADER__
