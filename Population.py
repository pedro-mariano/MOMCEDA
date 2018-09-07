from numpy import *
from objectives import *
from dominate import *
from topsisEval import *
from variation import *
import matplotlib.pyplot as plt

class Population:

##      NPop: Population size
##      Csize: Chromosome size
##      NObj: Number of objectives to optimize
##      Vmin: Limits of chromosome values
##      Vmax
##      minim: minim = 1 if minimizing objectives, minim = 0 otherwise
##      function: type of objective function
##      members: Array of chromosomes
##      obj: Array of objectives values for members of the population
##      rank: Rank of the nondominated front for each member
##      globalRank: Global rank of the members according to the criteria
##      fronts: List of nondominated fronts and their members
##      crowd: Crowding distance for each member
##      pi: Associated reference point for each member
##      d: Distance to the reference point for each member
##      domCount: number of members that dominate each member
   
        # Create a new population
    def __init__(self, NPop, Csize, Vmin, Vmax, NObj,minim, function):
        self.NPop = NPop
        self.Csize = Csize
        self.NObj = NObj
        self.Vmin = Vmin
        self.Vmax = Vmax
        self.function = function
        self.members = random.rand(NPop,Csize)*(Vmax - Vmin) + Vmin
        self.obj = objectives(function,self.members,NPop,NObj)
        self.minim = minim
        self.crowd = zeros(NPop)

    # Union between the current population and 'pop'
    def addPopulation(self, pop):
    
        npop = Population(0,pop.Csize,pop.Vmin,pop.Vmax,pop.NObj,pop.minim,pop.function)
        npop.NPop = self.NPop + pop.NPop
        npop.members = vstack((self.members,pop.members))
        npop.obj = vstack((self.obj,pop.obj))
        npop.crowd = hstack((self.crowd,pop.crowd))
##        npop.rank = vstack((self.rank,pop.rank))
    
        return npop

    # Add members from calling population to another one.
    # indList: list of the members' indexes
    def addMembers(self,indList):
    
        newPop = Population(0,self.Csize,self.Vmin,self.Vmax,self.NObj,self.minim,self.function)
        newPop.NPop = len(indList)
        newPop.members = self.members[indList,:]
        newPop.obj = self.obj[indList,:]
        newPop.crowd = self.crowd[indList]
        newPop.rank = self.rank[indList]
        newPop.domCount = self.domCount[indList]
##        newPop.fronts = []
##        for i in arange(len(sizeEvol)-1):
##            newPop.fronts.append(arange(sizeEvol[i],sizeEvol[i+1]))
##        newPop.globalRank = self.globalRank[indList]

        return newPop

    # Remove listed members from a population
    # indList : list of the members' indexes
    # obj : original objective values
    def removeMembers(self,indList,obj):

        newPop = Population(0,self.Csize,self.Vmin,self.Vmax,self.NObj,self.minim,self.function)
        newPop.NPop = self.NPop - len(indList)
        newPop.members = delete(self.members,indList,0)
        newPop.obj = delete(obj,indList,0)
        newPop.crowd = delete(self.crowd,indList,0)
        newPop.pi = delete(self.pi,indList,0)
        newPop.d = delete(self.d,indList,0)
        newPop.rho = self.rho
        newPop.rank = delete(self.rank,indList,0)
        newPop.domCount = delete(self.domCount,indList,0)
##        newPop.distPareto = delete(self.distPareto,indList,0)
        newPop.hvCont = delete(self.hvCont,indList,0)
##        newPop.fronts = []
##        newPop.fronts.append(self.fronts[:-1])
##        newPop.fronts.append(self.fronts[-1][:-len(indList)])

        return newPop

    # Sort the population into nondominated fronts until Nsort
    # members are found
    def fastNonDominatedSort(self,Nsort=inf):
     
        # Matrix of domination : 
        # M(i,j) = 1 if member i dominates member j
        # M(i,j) = 0 otherwise
        domMatrix = zeros((self.NPop,self.NPop),dtype=int)
        rank = zeros(self.NPop,dtype = int)
        fronts = [] # Array of fronts
        currentFront = array([],dtype=int)
        sortCount = 0 # Counter of sorted members
        frontCount = 0 # Front counter
        domCount = zeros(self.NPop) # Domination counter
        obj = self.obj
        
        for p in arange(self.NPop):
            
            for q in arange(p+1,self.NPop):
                
                domMatrix[p,q] = dominate(obj[p,:],obj[q,:],self.minim)
                
                if(domMatrix[p,q]): # If p dominates q
                    domMatrix[q,p] = 0 
                else:
                    domMatrix[q,p] = dominate(obj[q,:],obj[p,:],self.minim) 
                
            domCount[p] = sum(domMatrix[:,p])
            
            if(domCount[p] == 0):
                
                rank[p] = 1 # If p is nondominated, add p to the first front
                currentFront = append(currentFront,p)
                sortCount = sortCount + 1


        fronts.append(currentFront) # First front #
        frontCount = 1 # Front counter
        
        while(fronts[frontCount-1].size > 0 and sortCount < Nsort):
            currentFront = array([],dtype=int)
            for p in fronts[frontCount-1]:
                q = arange(self.NPop)
                test = (domMatrix[p,:] == 1)
                for q in q[test]:
                    # Update domination counter for members dominated
                    # by the ones on the previous front
                    domCount[q] = domCount[q] - 1 
                    if(domCount[q] == 0):
                        # If q is nondominated, add q to the current front
                        rank[q] = frontCount + 1
                        currentFront = append(currentFront,q)
                        sortCount = sortCount + 1
                    
            frontCount = frontCount + 1
            fronts.append(currentFront) 
        
        if(Nsort == inf):
            fronts = fronts[:-1]
        self.fronts =  fronts
        self.rank = rank
        
        self.domCount = domMatrix.sum(axis=0)


    # Calculate the crowding distance for individuals in a nondominated
    # front
    # minObj,maxObj: array of minimum and maximum values for the objectives
    def crowdingDistanceAssignment(self,front,minObj,maxObj):

        Nfront = len(front) # Size of the front
        crowd = zeros(Nfront) # Array of crowding distances
        obj = self.obj[front,:]

        for i in arange(self.NObj): # For each objective

            ind = argsort(obj[:,i]) # Sort individuals according to objectives

            # Boundaries
            crowd[ind[0]] = Inf 
            crowd[ind[Nfront-1]] = Inf

            norm = maxObj[i] - minObj[i] # Normalization
            for j in arange(1,Nfront - 1):
                crowd[ind[j]] = crowd[ind[j]] + (obj[ind[j+1],i] - obj[ind[j-1],i])/norm

        self.crowd = crowd

    # Generate offspring population
    # If probability pMut and standard deviation sigma are specified,
    # perform recombination and mutation
    # Otherwise, estimate the distribution of the population
    def offspringPop(self,coefGau,sigma,pMut,spread,pSwitch=0,nc=None,nm=None):
        
        newPop = Population(self.NPop,self.Csize,self.Vmin,self.Vmax,self.NObj,self.minim,self.function)
        
        # Estimate distribution
        if(nc == None):
            
            ind = argsort(self.globalRank)
            newMembers = zeros((2,self.Csize))
            
            for i in arange(self.NPop/2):

                pos = zeros(2,dtype=int)
                roulet = random.rand(2)

                sumCoef = coefGau[0]
                # Find gaussian to be sampled             
                while(roulet[0] >= sumCoef):
                    pos[0] = pos[0] + 1
                    sumCoef = sumCoef + coefGau[pos[0]]

                sumCoef = coefGau[0]
                while(roulet[1] >= sumCoef):
                    pos[1] = pos[1] + 1
                    sumCoef = sumCoef + coefGau[pos[1]]

                if(pos[0] == pos[1]):
                    if(pos[0] == 0):
                        pos[1] = pos[0] + 1
                    elif(pos[0] == self.NPop-1):
                        pos[1] = pos[0] - 1
                    else:
                        if(random.rand() <= 0.5):
                            pos[1] = pos[0] - 1
                        else:
                            pos[1] = pos[0] + 1
                
                # Sampling of the gaussians
                ind1 = self.members[ind[pos[0]],:]
                ind2 = self.members[ind[pos[1]],:]
                x1 = vstack((ind1,ind2)).min(axis=0)
                x2 = vstack((ind1,ind2)).max(axis=0)
                deltaMember = x2 - x1

##                stdev = (spread/4.0)*deltaMember
##
##                newMembers[0] = stdev*random.randn(self.Csize) + ind1
##                newMembers[1] = stdev*random.randn(self.Csize) + ind2

                spread1 = (0.1-0.5)*coefGau[pos[0]]/coefGau[0] + 0.5
                spread2 = (0.1-0.5)*coefGau[pos[1]]/coefGau[0] + 0.5
                stdev1 = (spread1/4.0)*deltaMember
                stdev2 = (spread2/4.0)*deltaMember
                
                newMembers[0] = stdev1*random.randn(self.Csize) + ind1
                newMembers[1] = stdev2*random.randn(self.Csize) + ind2

##                pSwitch = -0.5*(coefGau[max(pos)] - coefGau[0])/coefGau[0]

##                pos = zeros(2,dtype=int)
##                roulet = random.rand(2)
##
##                # Find gaussian to be sampled
##
##                sumCoef = coefGau[0]
##                while(roulet[0] >= sumCoef):
##                    pos[0] = pos[0] + 1
##                    sumCoef = sumCoef + coefGau[pos[0]]
##
##                sumCoef = coefGau[0]
##                while(roulet[1] >= sumCoef):
##                    pos[1] = pos[1] + 1
##                    sumCoef = sumCoef + coefGau[pos[1]]
##                    
####                # Sampling of the gaussians
####                stdev = sigma*(self.Vmax - self.Vmin)
####                nMember1 = stdev*random.randn(self.Csize) + self.members[ind[pos[0]],:]
####                nMember2 = stdev*random.randn(self.Csize) + self.members[ind[pos[1]],:]
####                newMember = (coefGau[pos[0]]*nMember1 + coefGau[pos[1]]*nMember2)/(coefGau[pos[0]]+coefGau[pos[1]])
##
##                # Sampling of the gaussians
##                stdev = sigma*(self.Vmax - self.Vmin)
##                meanMember = 0.5*(self.members[ind[pos[0]],:]+self.members[ind[pos[1]],:])
##                newMember = stdev*random.randn(self.Csize) + meanMember
                
                # Fix variables outside their limits
##                indPlus = where(newMember > self.Vmax)[0]
##                indMinus = where(newMember < self.Vmin)[0]
##                newMember[indPlus] = self.Vmax[indPlus] - remainder(newMember[indPlus],self.Vmax[indPlus]-self.Vmin[indPlus])    
##                newMember[indMinus] = self.Vmin[indMinus] - remainder(newMember[indMinus],self.Vmin[indMinus]-self.Vmax[indMinus])

##                newPop.members[i,:] = newMember

                newMembers = maximum(minimum(newMembers,self.Vmax),self.Vmin)

                u = random.rand(self.Csize)
                indLess = where(u <= pSwitch)[0]
                newPop.members[2*i:2*i+2] = copy(newMembers)
                newPop.members[2*i,indLess] = newMembers[1,indLess]
                newPop.members[2*i+1,indLess] = newMembers[0,indLess]

            
##            # Mutation
            stdev = sigma*(self.Vmax - self.Vmin)
            for i in arange(self.NPop):
                ind = newPop.members[i,:]
##                newValue = mutPolynomialBounded(ind, nc-10, self.Vmin.tolist(), self.Vmax.tolist(), 1.0/self.Csize)
##                newPop.members[i,:] = newValue
                for j in arange(self.Csize):
                    if(random.rand() <= pMut):
                        newValue = stdev[j]*random.randn() + ind[j]
                        newValue = maximum(minimum(newValue,self.Vmax[j]),self.Vmin[j])
                        newPop.members[i,j] = newValue


        # Recombination and Mutation
        else:
            # Binary tournament
            for i in arange(self.NPop):
                members = random.randint(self.NPop,size=2)
                rank1 = self.globalRank[members[0]]
                #crowd1 = self.crowd[members[0]]
                rank2 = self.globalRank[members[1]]
                #crowd2 = self.crowd[members[1]]
                if(rank1 < rank2):
                    newPop.members[i,:] = self.members[members[0],:]
                else:
                    newPop.members[i,:] = self.members[members[1],:]
                
            # Recombination
            for i in arange(self.NPop/2):
                coef = random.rand()
                child1 = coef*newPop.members[2*i,:] + (1-coef)*newPop.members[2*i+1,:]
                child2 = (1-coef)*newPop.members[2*i,:] + coef*newPop.members[2*i+1,:]
                newPop.members[2*i,:] = child1
                newPop.members[2*i+1,:] = child2
            
            # Mutation
            for i in arange(self.NPop*self.Csize):
                if(random.rand() <= pMut):
                    ind = i/self.Csize
                    pos = i % self.Csize                   
                    newValue = newPop.members[ind,pos] + sigma[pos]*random.randn()
                    while ((newValue > self.Vmax[pos]) | (newValue < self.Vmin[pos])):
                        newValue = newPop.members[ind,pos] + sigma[pos]*random.randn()
                    
                    newPop.members[ind,pos] = newValue
                
        newPop.obj = objectives(newPop.function,newPop.members,newPop.NPop,newPop.NObj)
        return newPop

    # Apply the TOPSIS sorting to each front of the population
    # rank: 'default' for TOPSIS rank
    #       'refDir' for TOPSIS rank using reference directions
    #        'all' for using ref. directions and front number
    def topsisPop(self,weight=array([]),rank='default'):

        rankFront = zeros(self.NPop,dtype=int) # rank within each front
        globalRank = zeros(self.NPop,dtype=int) # Global rank
        
        if(rank == 'refDir'):
            if(len(weight) != 2):
                weight = ones(2)
            i = arange(self.NPop)
            j = self.pi[i]
            matrixEval = c_[self.d[i,j],self.rho[j]]
            
        if(rank == 'default'):
            if(len(weight) != self.NObj):
                weight = ones(self.NObj)
            matrixEval = self.obj
            

        if(rank == 'hv'):
            if(len(weight) != 5):
                weight = ones(5)
            i = arange(self.NPop)
            j = self.pi[i]
            matrixEval = c_[self.domCount,self.d[i,j],self.rho[j],self.rank,-self.hvCont]
            relSim = topsisEval(matrixEval,self.minim,weight)
            ind = relSim.argsort()[::-1]
            globalRank[ind] = arange(1,self.NPop+1)
            
        elif(rank == 'all'):
            if(len(weight) != 4):
                weight = ones(4)
            i = arange(self.NPop)
            j = self.pi[i]
            matrixEval = c_[self.domCount,self.d[i,j],self.rho[j],self.rank]
            relSim = topsisEval(matrixEval,self.minim,weight)
            ind = relSim.argsort()[::-1]
            globalRank[ind] = arange(1,self.NPop+1)
            
        else:
            relSim = topsisEval(matrixEval,self.minim,weight)

            # For each front i
            for i in arange(max(self.rank)):
                indFront = where(self.rank == i+1)[0]
                relFront = relSim[indFront]
                ind = relFront.argsort()[::-1]
                rankFront[indFront[ind]] = arange(1,len(indFront)+1)

            accFrontSize = zeros(max(self.rank),dtype=int)
            accFrontSize[0] = 0
            for i in arange(1,max(self.rank)):
                accFrontSize[i] = accFrontSize[i-1] + sum(self.rank == i)

            globalRank = accFrontSize[self.rank-1] + rankFront
            
        self.globalRank = globalRank

    # Plot nondominated members of a generation
    def plot(self,color,scale,center,names, countFig):

        obj = self.obj[self.rank == 1]
        if (self.NObj == 1):
            plt.plot(obj[:,0],'.',color=color)
            plt.xlabel('Generations')
            plt.ylabel(names[0])
            plt.draw()
        else:
            pair_n = 0
            for i in arange(self.NObj):
                 for j in arange(i+1,self.NObj):

                     plt.figure(pair_n + countFig)
                     plt.plot(obj[:,i]*scale[i]+center[i],obj[:,j]*scale[j]+center[j],'.',color=color)
                     plt.xlabel(names[i])
                     plt.ylabel(names[j])
                     plt.draw()
                     pair_n = pair_n + 1

        plt.pause(0.05)
