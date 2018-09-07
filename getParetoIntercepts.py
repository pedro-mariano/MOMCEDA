from numpy import *
import matplotlib.pyplot as plt

def getParetoIntercepts(Z,objRec,function):

    a = (Z[:,1] - objRec.objIdeal[1])/(Z[:,0] - objRec.objIdeal[0])
    b = Z[:,1] - a*Z[:,0]
    
    f1 = ((sqrt(1-4*a*(b-1))-1)/(2*a))**2
    f1[a == inf] = objRec.objIdeal[0]
    f2 = 1 - sqrt(f1)

    parIn = c_[f1,f2]

    return parIn
