from numpy import *

def associate( Pop, Zr, nextNPop ):
#associate 
#   Association of population members with reference directions

#    Input:
#    Pop: population
#    Zr : reference directions
#    nextNPop : number of members on the first fronts

    nDir = Zr.shape[0] # number of directions

    pi = zeros(Pop.NPop, dtype=int) # index of associated reference direction for each member
    d = zeros((Pop.NPop,nDir)) # distance to the reference directions for each member
    rho = zeros(nDir, dtype=int) # niche count for each reference direction

    for member in arange(Pop.NPop):

        # Calculate distance to every reference direction
        for direct in arange(nDir):
            s = Pop.obj[member,:]
            w = Zr[direct,:]
            vect = s-w.dot(s)*w/(w.dot(w))
            d[member,direct] = sqrt(vect.dot(vect))
        
        ind = argmin(d[member,:])
        pi[member] = ind
        
        # Niche count for members not belonging to the last front
        if(member < nextNPop):
            rho[ind] = rho[ind] + 1
        
    Pop.pi = pi
    Pop.d = d
    Pop.rho = rho



