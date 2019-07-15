#!/usr/bin/python
import mpmath as mp

def gaulob(xL,xR,degree):
    mp.dps = 300

    eps = mp.power(10,-200)
    stop = False
    iter = 0
    diff = 0

    nodes = [mp.mpf(0.)] * degree
    weight = [mp.mpf(0.)] * degree
    nodes_old = [mp.mpf(0.)] * degree
    
    P = [mp.mpf(0.)] * (degree * degree)

    for i in range(degree):
        nodes[i] = mp.cos(mp.pi*mp.mpf(i)/mp.mpf((degree-1)))
    
    xL = mp.mpf(xL)
    xR = mp.mpf(xR)
    while(not stop):
        iter = iter + 1
        nodes_old = nodes[:]
        for i in range(degree):
            P[i] = mp.mpf(1)
            P[i+degree] = nodes[i]
        for k in range(2,degree):
            for i in range(degree):
                P[i+k*degree] = ((2*k-1)*nodes[i]*P[i+(k-1)*degree] - (k-1)*P[i+(k-2)*degree]) / k
        for i in range(degree):
            nodes[i] = nodes_old[i] - (nodes[i] * P[i+(degree-1)*degree] - P[i+(degree-2)*degree]) / (degree*P[i+(degree-1)*degree])
        if(max([abs(x-y) for x, y in zip(nodes,nodes_old)]) < eps):
            break
        
    for i in range(degree):
        weight[i] = 2. / (degree*(degree-1)*P[i+(degree-1)*degree]*P[i+(degree-1)*degree])
        weight[i] = 0.5*weight[i]*(xR-xL)
        nodes[i]  = xL+0.5*(nodes[i]+1)*(xR-xL)

    #print("p=%d: " % p);
    #print("sumWeights=%s " % sum(weight));
    #print("sumNodes=%s"    % sum(nodes));
    
    return nodes, weight

if __name__ == '__main__':
    from mpmath import mp
    import numpy,itertools

    printPrec=64
    mp.dps=256
    maxOrder = 15
    generateCPPArrays = False

    for p in range(1,maxOrder+1):
        x, w = gaulob(0,1,p+1)

        if generateCPPArrays:
            print("const double kernels::gaussLobattoWeights%d[%d]={\n  %s\n};"% (p,p+1,",\n  ".join([mp.nstr(i,printPrec) for i in w])))
            print("const double kernels::gaussLobattoNodes%d[%d]={\n  %s\n};"% (p,p+1,",\n  ".join([mp.nstr(i,printPrec) for i in x])))
        else:
            print("if nDof == %s:" % (p+1))
            print("    return [%s], [%s]" % (",".join([mp.nstr(i,printPrec) for i in w]),",".join([mp.nstr(i,printPrec) for i in x])))

        print("")
