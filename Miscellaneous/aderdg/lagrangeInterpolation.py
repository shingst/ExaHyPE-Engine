"""
.. module:: lagrangeinterp
  :platform: Unix, Windows, Mac
  :synopsis: Provides a routine to obtain symbolic expressions 
   for interpolating polynomials for Lagrange type.

  :synopsis: Provides a routine to obtain symbolic expressions 
   for interpolating polynomials for Lagrange type.
"""
import sympy

def LagrangBasisPoly(x,order,i,xi=None):
    """
    Returns the i-th Lagrange basis polynomial as symbolic expression.
    Algorithm was obtained from https://wanglongqi.github.io/python/2014/03/24/implement-of-lagrange-polynomial/.

    Parameters
    ----------
    x: symbolic variable
       Symbolic variable.
    order : number
       Order of the Lagrange basis polynomial.
    xi : number
       The node the basis polynomial should pass through
    i : number
       The node we are interested in.
    Returns
    -------
    Li: symbolic expression
       The i-th Lagrange basis polynomial as symbolic expression.
    """
    if xi==None:
        print("Please specify some quadrature nodes")
        raise 
    index = list(range(order+1))
    index.pop(i)
    symbolic_function=sympy.prod([(x-xi[j])/(xi[i]-xi[j]) for j in index])
    return symbolic_function


def LagrangePoly(x,Lx, Ly):
    """
    Returns the Lagrange interpolation polynomial corresponding to the support points Lx 
    and nodal values Ly as symbolic expression.
    Algorithm was obtained from http://stackoverflow.com/questions/27744475/lagrange-interpolation-in-python-as-a-result-matematical-formula

    Parameters
    ----------
    x: symbolic variable
       Symbolic variable.
    Lx : list
       The support points.
    Ly : list
       Nodal values associated with the support points.
    Returns
    -------
    L: symbolic expression
       The Lagrange interpolation polynomial as symbolic expression.
    """
    if  len(Lx)!= len(Ly):
        print("ERROR")
        return 1
    symbolic_function=0
    for k in range ( len(Lx) ):
        t=1
        for j in range ( len(Lx) ):
            if j != k:
                t=t* ( (x-Lx[j]) /(Lx[k]-Lx[j]) )
        symbolic_function+= t*Ly[k]
    return symbolic_function
