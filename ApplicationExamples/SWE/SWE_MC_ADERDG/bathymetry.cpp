#include "bathymetry.h"
#include "MySWESolver_ADERDG_Variables.h"
#include "MySWESolver_FV_Variables.h"
#include "MySWESolver.h"
#include "readCsv.h"

#include <cmath>

using namespace std;
///// 2D /////

#ifdef Dim2


double linearInterpolation_1D(double x, double dx, int idx, int nelem){
    if(idx == 0){
        if(x < dx)
            return (1.0 - x/dx);
    }
    if(idx == nelem){
        if(x > 1.0 - dx)
            return (x/dx - (nelem-1.0));
    }
    for(int i = 1; i < nelem+1; i++){
        if(idx == i){
            //left
            if(x > (i-1)*dx && x <= i*dx)
                return (x/dx - (i-1.0));
            //right
            if(x > i*dx && x < (i+1)*dx)
                return (- x/dx + i+1.0);
        }
    }
    return 0.0;
}

double SWE::linearInterpolation(double x, double y, std::vector<double>& a){
    int nelem = (int)std::pow(a.size(),0.5)-1;
    assert((nelem+1)*(nelem+1) == a.size());
    double res = 0.0;
    //loop over "basis functions"
    for(int i=0; i < nelem + 1; i ++){
        for(int j=0; j < nelem + 1; j ++){
             res += a[(nelem+1)*i+j]*linearInterpolation_1D(x,1.0/nelem,i,nelem)*linearInterpolation_1D(y,1.0/nelem,j,nelem);
        }
    }
    return res;
}




#endif
