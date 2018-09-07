import numpy as np
import scipy.spatial.distance as scp
from math import hypot,sqrt
from PyGMO.util import hypervolume

def IGDmetric(Z,Zmin,A,function):

    Z = getParetoIntercepts(Z,Zmin,function)

    sumDist = 0.0
    
    for i in np.arange(Z.shape[0]):
        
        minDist = np.inf
        
        for j in np.arange(A.shape[0]):
            dist = scp.euclidean(Z[i],A[j])
            if(dist < minDist):
                minDist = dist

        sumDist = sumDist + minDist

    igd = sumDist/Z.shape[0]

    return igd,Z


def getParetoIntercepts(Z,Zmin,function):

    ZPar = np.copy(Z)

    if(function == 'DTLZ1'):
        t = (0.5 - sum(Zmin))/((Z-Zmin).sum(axis=1))
        for i in np.arange(Z.shape[0]):
            ZPar[i] = Zmin*(1-t[i]) + t[i]*Z[i]

    elif(function == 'DTLZ2'):
        a = ((Z-Zmin)**2).sum(axis=1)
        b = 2*(Zmin*(Z-Zmin)).sum(axis=1)
        c = (Zmin**2).sum() - 1
        t = (-b + np.sqrt((b**2)-4*a*c))/(2*a)
        for i in np.arange(Z.shape[0]):
            ZPar[i] = Zmin*(1-t[i]) + t[i]*Z[i]
            
    elif((function == 'ZDT1') or (function == 'ZDT4')):

        dy = Z[:,1] - Zmin[1]
        dx = Z[:,0] - Zmin[0]
        a = dy**2
        b = 2*Zmin[1]*dy - 2*dy - dx
        c = 1 - 2*Zmin[1] + Zmin[1]**2 - Zmin[0]
        delta = b**2-4*a*c
        ind = np.where(delta >= 0)[0]
        t = np.ones(Z.shape[0])
        t[ind] = (-b[ind]-np.sqrt(delta[ind]))/(2*a[ind])
        t[a == 0] = -c/b[a==0]
        for i in ind:
            ZPar[i] = Zmin*(1-t[i]) + t[i]*Z[i]
        ZPar = ZPar[ind]

    elif((function == 'ZDT2') or (function == 'ZDT6')):
        dy = Z[:,1] - Zmin[1]
        dx = Z[:,0] - Zmin[0]
        a = dx**2
        b = dy + 2*Zmin[0]*dx
        c = Zmin[1] - 1 + Zmin[0]**2
        delta = b**2-4*a*c
        ind = np.where(delta >= 0)[0]
        t = np.ones(Z.shape[0])
        t[ind] = (-b[ind]+np.sqrt(delta[ind]))/(2*a[ind])
        t[a == 0] = -c/b[a==0]
        for i in ind:
            ZPar[i] = Zmin*(1-t[i]) + t[i]*Z[i]
        ZPar = ZPar[ind]
              
    return ZPar

def normIGDmetric(Z,Zmin,a,A,function):

    Zint = getParetoIntercepts(Z,Zmin,function)
    Z = (Zint-Zmin)/a
    A = (A - Zmin)/a

    sumDist = 0.0
    
    for i in np.arange(Z.shape[0]):
        
        minDist = np.inf
        
        for j in np.arange(A.shape[0]):
            dist = scp.euclidean(Z[i],A[j])
            if(dist < minDist):
                minDist = dist

        sumDist = sumDist + minDist

    igd = sumDist/Z.shape[0]

    return igd,Zint


def distPareto(Pop, Z, Zmin, a):

    ZPar = getParetoIntercepts(Z,Zmin,Pop.function)
    Z = (ZPar - Zmin)/a
    points = Pop.pi[np.arange(Pop.NPop)]
    dist = (Z[points] - Pop.obj)**2
    dist = np.sqrt(dist.sum(axis=1))
    
    Pop.distPareto = dist

def hvContribution(Pop,Zmin,a):

    hvCont = np.zeros(Pop.NPop)
    indND = np.where(Pop.rank == 1)[0]
    NDobj = Pop.obj[indND]*a + Zmin
    ref = NDobj.max(axis=0)*1.1
    hv = hypervolume(NDobj.tolist())
    for i in np.arange(len(indND)):
        hvCont[indND[i]] = hv.exclusive(i,ref.tolist())

    Pop.hvCont = hvCont
        
def diversity(first_front, first, last):
    """Given a Pareto front `first_front` and the two extreme points of the 
    optimal Pareto front, this function returns a metric of the diversity 
    of the front as explained in the original NSGA-II article by K. Deb.
    The smaller the value is, the better the front is.
    """
    df = hypot(first_front[0,0] - first[0],
               first_front[0,1] - first[1])
    dl = hypot(first_front[-1,0] - last[0],
               first_front[-1,1] - last[1])
    dt = [hypot(first[0] - second[0],
                first[1] - second[1])
          for first, second in zip(first_front[:-1], first_front[1:])]

    if (first_front.shape[0]) == 1:
        return df + dl

    dm = sum(dt)/len(dt)
    di = sum(abs(d_i - dm) for d_i in dt)
    delta = (df + dl + di)/(df + dl + len(dt) * dm )
    return delta

def convergence(first_front, optimal_front):
    """Given a Pareto front `first_front` and the optimal Pareto front, 
    this function returns a metric of convergence
    of the front as explained in the original NSGA-II article by K. Deb.
    The smaller the value is, the closer the front is to the optimal one.
    """
    distances = []
    
    for ind in first_front:
        distances.append(float("inf"))
        for opt_ind in optimal_front:
            dist = 0.
            for i in xrange(len(opt_ind)):
                dist += (ind[i] - opt_ind[i])**2
            if dist < distances[-1]:
                distances[-1] = dist
        distances[-1] = sqrt(distances[-1])
        
    return sum(distances) / len(distances)
    
