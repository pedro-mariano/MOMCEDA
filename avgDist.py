import numpy as np

def avgDist(M):

    NPop = M.shape[0]
    Csize = M.shape[1]
    ncomb = NPop*(NPop-1)/2
    dist = np.zeros(Csize)
    
    for i in np.arange(NPop):
        for j in np.arange(i+1,NPop):
            dist = dist + abs(M[i,:] - M[j,:])
    dist = dist/ncomb

##    for i in np.arange(NPop):
##        for j in np.arange(i+1,NPop):
##            nDist = abs(M[i,:] - M[j,:])
##            dist = np.vstack((dist,nDist)).min(axis=0)


    return dist

def nicheDist(Pop):

    nDir = len(Pop.rho)

    dist = np.zeros(Pop.Csize)

    for i in np.arange(nDir):
        members = Pop.members[Pop.pi == i]
        if(members.shape[0] > 1):
            dist = dist + avgDist(members)

    dist = dist/nDir

    return dist

def maxDist(M):

    NPop = M.shape[0]
    Csize = M.shape[1]
    ncomb = NPop*(NPop-1)/2
    dist = np.zeros(Csize)

    for i in np.arange(NPop):
        for j in np.arange(i+1,NPop):
            nDist = abs(M[i,:] - M[j,:])
            dist = np.vstack((dist,nDist)).max(axis=0)

    return dist

def minDist(M):

    NPop = M.shape[0]
    Csize = M.shape[1]
    ncomb = NPop*(NPop-1)/2
    dist = np.ones(Csize)*np.inf

    for i in np.arange(NPop):
        for j in np.arange(i+1,NPop):
            nDist = abs(M[i,:] - M[j,:])
            dist = np.vstack((dist,nDist)).min(axis=0)

    return dist

def rouletDist(Pop,coefGau):
    
    ind = np.argsort(Pop.globalRank)
    members = np.zeros((2,Pop.Csize))
    dist = np.zeros(Pop.Csize)
    
    for nreps in np.arange(Pop.NPop):
        for i in np.arange(2):

            pos = 0
            roulet = np.random.rand()

            sumCoef = coefGau[0]
            # Find member position on the roulet      
            while(roulet >= sumCoef):
                pos = pos + 1
                sumCoef = sumCoef + coefGau[pos]

            members[i] = Pop.members[ind[pos],:]

        dist = dist + abs(members[0,:] - members[1,:])
        
    dist = dist/Pop.NPop
    return dist

def randomDist(M):

    NPop = M.shape[0]
    Csize = M.shape[1]
    dist = np.zeros(Csize)
    
    for reps in np.arange(NPop):
        ind = np.random.randint(NPop,size=2)
        dist = dist + abs(M[ind[0]] - M[ind[1]])
    dist = dist/NPop

    return dist
        
