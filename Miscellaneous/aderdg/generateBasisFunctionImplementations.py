#!/usr/bin/env python3
"""
.. module:: main
  :platform: Unix, Windows, Mac
  :synopsis: Generates reference basis functions and their first and second derivatives up to a specific N.

:synopsis: Generates reference basis functions and their first and second derivatives up to a specific N.
"""
from lagrangeInterpolation import *
import mpmath as mp
import sympy
import re

def generateBasisFunctionDefinitions(outfile,rule,order,x,printPrec):
    mp2 = mp.mp.clone()
    mp2.dps=64
    x2 = x.copy();

    for i,xi in enumerate(x2): 
        x2[i] = mp2.mpf(xi)
    
    for m in range(0,order+1):
        s=sympy.symbols('s')
        outfile.write("double kernels::%s::basisFunction_%d_%d(const double& s) {" % (rule,order,m))
        refphi = LagrangBasisPoly(s,order,m,x2)


        ret = "  return %s;" % sympy.simplify(refphi)
        ret = re.sub(r"s\*\*([0-9]+)", r"pow\1(s)", ret)
        outfile.write(ret)
        outfile.write("}\n")

    for m in range(0,order+1):
        s=sympy.symbols('s')
        refphi     = LagrangBasisPoly(s,order,m,x2)
        dds_refphi = sympy.diff(refphi, s)
        outfile.write("double kernels::%s::basisFunctionFirstDerivative_%d_%d(const double& s) {" % (rule,order,m))
        ret = "  return %s;" % sympy.simplify(dds_refphi)
        ret = re.sub(r"s\*\*([0-9]+)", r"pow\1(s)", ret)
        outfile.write(ret)
        outfile.write("}\n")

    for m in range(0,order+1):
        s=sympy.symbols('s')
        refphi       = LagrangBasisPoly(s,order,m,x2)
        d2ds2_refphi = sympy.simplify(sympy.diff(refphi, s, s))
        outfile.write("double kernels::%s::basisFunctionSecondDerivative_%d_%d(const double& s) {" % (rule,order,m))
        ret = "  return %s;" % sympy.simplify(d2ds2_refphi)
        ret = re.sub(r"s\*\*([0-9]+)", r"pow\1(s)", ret)
        outfile.write(ret)
        outfile.write("}\n")
