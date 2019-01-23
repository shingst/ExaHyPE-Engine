#include "bathymetry.h"
#include "MySWESolver_Variables.h"
#include "MySWESolver.h"

#include <cmath>

using namespace std;
///// 2D /////

extern const int nelem;

#ifdef Dim2


double linearInterpolation_1D(double x, double dx, int idx){
    if(idx == 0){
        if(x < dx)
            return (1.0 - x/dx);
    }
    if(idx == nelem){
        if(x > 1.0 - dx)
            return (x/dx - (nelem-1.0));
    }
    for(int i = 1; i < nelem; i++){
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

double SWE::linearInterpolation(double x, double y, double* a){
    double res = 0.0;
    //loop over "basis functions"
    for(int i=0; i < nelem + 1; i ++){
        for(int j=0; j < nelem + 1; j ++){
             res += a[(nelem+1)*i+j]*linearInterpolation_1D(x,1.0/nelem,i)*linearInterpolation_1D(y,1.0/nelem,j);
        }
    }
    return res;
}

double linearInterpolationCoarse(double x, double y, double* a){
    double res = 0.0;
    //loop over "basis functions"
    for(int i=1; i < 3 + 1; i ++){
        for(int j=1; j < 3 + 1; j ++){
             res += 0.02*a[(3+1)*(i-1)+(j-1)]*linearInterpolation_1D(x,1.0/4,i)*linearInterpolation_1D(y,1.0/4,j);
        }
    }
    return res;
}

double SWE::bathymetry(double x, double y) {
    double a[16] = {3.0,6.0,4.0,-1.0,-3.0,-6.0,18.0,1.0,-4.0,15.0,9.0,12.0,17.0,-18.0,24.0,6.0};
    return linearInterpolationCoarse(x,y,a);
}



#endif
