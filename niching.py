from numpy import *
from topsisEval import *

def niching(K,Pop,sizeEvol,weight,multiple=True):
#   Niching procedure

#    Input:
#    K: number of members to be added to the population
#    Pop: population
#    sizeEvol: evolution of population's size
#    weight: weight vector for TOPSIS ranking
#    multiple: use multiple criteria to select solutions
#              use only niche count if False

#    Output:
#    indRemove: indices of members to be removed from Pop

    F = arange(sizeEvol[-2],sizeEvol[-1]) # indices of the last front

    k = 0

    if(multiple):

        while (k < K):

            j = Pop.pi[F] # directions for members of the last front
            matrixEval = c_[Pop.domCount[F],Pop.d[F,j],Pop.rho[j],Pop.rank[F],-Pop.hvCont[F]]
            relSim = topsisEval(matrixEval,Pop.minim,weight) # TOPSIS ranking
            ind = F[argmax(relSim)] # Selected member to stay on population
            direct = Pop.pi[ind]
            Pop.rho[direct] = Pop.rho[direct] + 1 # Update niche count
            F = delete(F,where(F == ind)) # delete chosen individual from the list
            k = k + 1

        indRemove = F
        
    else:
    
        indList = zeros(K, dtype = int) # Members of the last front that will stay on Pop

        rho = Pop.rho # niche count for each reference direction
        pi = Pop.pi # index of associated reference direction for each member
        d = Pop.d # distance to the reference directions for each member

        validDir = ones(len(rho),dtype=bool) # All directions are valid in the beginning

        while (k < K):

            # Choose a random valid reference direction among the least populated
            minRho = min(rho[validDir])
            Jmin = where(rho == minRho)[0]
            j = random.choice(Jmin)
            Ij = F[pi[F] == j] # Members of the last front close to the selected direction

            if(size(Ij)>0):
                if(rho[j] == 0): # If the direction is not populated, add the closest individual
                    minD = d[Ij,j]
                    ind = argmin(minD)
                    indList[k] = Ij[ind]
                else: # Otherwise, choose a random individual to be added
                    indList[k] = random.choice(Ij)
                
                rho[j] = rho[j] + 1 # Update niche count
                F = delete(F,where(F == indList[k])) # delete individual from the last front
                k = k + 1
            else: # If there are no members in the last front close to the direction
                validDir[j] = False

        Pop.rho = rho

        indRemove = setdiff1d(arange(sizeEvol[-2],sizeEvol[-1]),indList,True)
    
    return indRemove
