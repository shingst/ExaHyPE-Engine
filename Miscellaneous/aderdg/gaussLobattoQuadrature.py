#!/usr/bin/python
import mpmath as mp

def gaulob(xL,xR,degree):
	mp.dps = 300

    eps = math.pow(10,-200)
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

    print("p=%d: " % p);
    print("sumWeights=%s " % sum(weight));
    print("sumNodes=%s"    % sum(nodes));
    
    return nodes, weight

if __name__ == '__main__':
    from mpmath import mp
    import numpy,itertools

    mp.dps=256
    maxOrder = 20
    generateCPPArrays = False

    for p in range(1,maxOrder+1):
        x, w = gaulob(0,1,p+1)

        if generateCPPArrays:
            for i,wi in enumerate(w):
                print("gaussLegendreWeights[%d][%d]=%s;"% (p,i,wi))
            for i,xi in enumerate(x):
                print("gaussLegendreNodes[%d][%d]=%s;" % (p,i,xi))
        else:
            print("if nDOF == %s:" % (p+1))
            print("    return [%s], [%s]" % (",".join([str(i) for i in x]),",".join([str(i) for i in w])))

        print("")