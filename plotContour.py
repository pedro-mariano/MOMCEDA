import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plotContour(Z,ZPar,NObj,p,ax):

    for i in np.arange(NObj):
        for j in np.arange(i+1,NObj):
            ext1 = np.where(Z[:,i] == 0)[0]
            ext2 = np.where(Z[:,j] == 0)[0]
            if(any(Z[ext2,i] != Z[ext1,j])):
                ext2 = ext2[::-1]
            for k in np.arange(p+1):
                coord = np.vstack((ZPar[ext1[k]],ZPar[ext2[k]]))
                ax.plot(coord[:,0],coord[:,1],coord[:,2],color='k',linewidth=0.5)
    
def plotContour2(ZPar,p,ax):

    begin = 0
    for line in np.arange(p):
        end = begin+p-line
        points = np.arange(begin,end+1)
        for i in points:
            if(i != end):
                coord = [i,i+1]
                ax.plot(ZPar[coord,0],ZPar[coord,1],ZPar[coord,2],color='k',linewidth=0.5)
                coord = [i,i+p+1-line]
                ax.plot(ZPar[coord,0],ZPar[coord,1],ZPar[coord,2],color='k',linewidth=0.5)
            if(i != begin):
                coord = [i,i+p-line]
                ax.plot(ZPar[coord,0],ZPar[coord,1],ZPar[coord,2],color='k',linewidth=0.5)
        begin = end + 1
