from Population import *
from coefGau import *
from objectiveRecords import *
from normalize import *
from generateRefPoints import *
from associate import *
from niching import *
import matplotlib.pyplot as plt
import cPickle as pickle
from hv import *
from avgDist import *

# NSGAII main loop 

## Parameters ##

NExec = 10
NPop = 100 # Population size
NEval = 20000 # Number of function evaluations
NGer = NEval/NPop - 1 # Number of generations
Csize = 30 # Chromosome size
NObj = 2 # Number of objectives to optimize
minim = 1 # minim = 1 if minimizing objectives, minim = 0 otherwise
function = 'ZDT1' # Type of problem

Vmin = 0.0*ones(Csize) # Limits of chromosome values
Vmax = 1.0*ones(Csize)

minObj = array([0.0,0.0])# Limits of objective values
maxObj = array([1.0,1.0])

refPoint = (1.1,1.1) # Reference point for hypervolume indicator

objRec = objectiveRecords(NObj,minim) # Records for objective values

p = 19 # Number of objective axes divisions to generate structured ref. points
Z = generateRefPoints(NObj, p) # Generate structured reference points
rankType = 'all' # Type of rank used for TOPSIS
weight = array([0.0,1.0,3.0,5.0]) # Array of weights used for TOPSIS
# Weight order: domCount, pi, rho, rank

hvValues = zeros(NGer)

# Plot parameters
color = plt.cm.get_cmap('Reds') 
ObjNames = ['f1', 'f2','f3']
scale = 1.0/(maxObj - minObj)
center = array([0.0,0.0,0.0])

## Mutation parameters

pMut = 0.2
sigma0 = r_[array([1.0]),1.0/10*ones(Csize-1)]
sigmaf = r_[array([1.0/10]),1.0/500*ones(Csize-1)]
#deltaSigma = (sigma0 - sigmaf) / NGer
qSigma = (sigmaf/sigma0)**(1.0/NGer)
sigma = sigma0
dist = zeros((NGer,Csize))

## Distribution variables

# Coeffients of the gaussians of the mixture

##choiceOK = False
##while(not choiceOK):
##    dec = input('Choose decay type:\n(1) Linear\n(2) Exponential\n(3) Logarithmic\n')
##    if(dec in [1,2,3]):
##        choiceOK = True

dec = 1 # decay type:1 - Linear, 2 - Exponential, 3 - Logarithmic
coefGau = coefGau(NPop,dec)

# figure()
# plot(coefGau)

for nexec in arange(NExec):

    ## Initialization ##

    t = 0 # Counter of generations
    Pt = Population(NPop, Csize, Vmin, Vmax, NObj, minim, function) # Initial Population
    Pt.fastNonDominatedSort() # Nondominated sorting
    Pt.topsisPop() # Rank within each front
    Qt = Pt.offspringPop(coefGau,sigma) # Offspring population

    plt.ion()

    ## Main loop ##

    while(t < NGer):
        
        Rt = Pt.addPopulation(Qt) # Combined population
        Rt.fastNonDominatedSort() # Nondominated sorting
    ##    Rt.topsisPop(rank=rankType) # Rank within each front
    ##    Rt.globalRankEval()
        
        indList = zeros(NPop, dtype = int) # List of indexes of members to the next population
        i = 0 # Counter of fronts
        sizeEvol = array([0, len(Rt.fronts[i])]) # Evolution of population's size by adding the fronts
        
        # Fill population with the first fronts
        while(sizeEvol[i+1] < NPop):
    #         Rt.crowd(Rt.fronts[i]) = Rt.crowdingDistanceAssignment(Rt.fronts[i],minObj,maxObj)
            indList[sizeEvol[i]:sizeEvol[i+1]] = Rt.fronts[i] # Add members to the list
            i = i + 1
            sizeEvol = append(sizeEvol, sizeEvol[i]+len(Rt.fronts[i]))
        
        # Sort members of the last front according to
        # crowding distance
    #     Rt.crowd[Rt.fronts[i]] = Rt.crowdingDistanceAssignment(Rt.fronts[i],minObj,maxObj)
    #     ind = argsort(Rt.crowd[Rt.fronts[i]])[::-1]

        # Sort members of the last front according to the rank
    #     ind = argsort(Rt.rank[Rt.fronts[i],1])

        listSt = r_[indList[:sizeEvol[i]], Rt.fronts[i]]
        St = Rt.addMembers(listSt) 

        K = NPop - sizeEvol[i]
        StObj = St.obj # Save original objective values of St
        Zr = normalize(St, objRec, minim, Z, p)
        associate(St, Zr, sizeEvol[i])
        indRemove = niching(K, St, sizeEvol)
        Pt = St.removeMembers(indRemove, StObj) # Next generation's population
        
        Pt.plot(color(float(t)/NGer), scale, center, ObjNames)
        Pt.topsisPop(weight,rank=rankType) # Rank within each front
        Qt = Pt.offspringPop(coefGau,sigma) # Selection, recombination and mutation

        if(remainder(t,1) == 0):
            dist[t] = avgDist(Pt.members)/(2*(Vmax - Vmin))
    ##        if((any(dist[t] > dist[t-1])) & (t > 0)):
    ##            ind = where(dist[t] > dist[t-1])[0]
    ##            dist[t,ind] = dist[t-1,ind]
        else:
            dist[t] = dist[t-1]
    ##        
    ##    dist[t,:] = nicheDist(Pt)/(Vmax - Vmin)
    ##    dist[t,:] = maxDist(Pt.members)/(10*(Vmax - Vmin))
        sigma = dist[t,:]
        
    ##    MEv = sqrt(sum(Pt.members[:,1:]**2)/(NPop*(Csize-1)))
    ##    print 'Generation = ', t, 'Mean Square Value =', MEv

        hv = HyperVolume(refPoint)
        hvValues[t] = hv.compute(Pt.obj)
        
        t = t + 1
    ##    sigma = sigma - deltaSigma
    ##    sigma = sigma*qSigma
        
        
    MEv = sqrt(sum(Pt.members[:,1:]**2)/(NPop*(Csize-1)))
    print 'Execution : ',nexec
    print 'Generation = ', t, 'Mean Square Value =', MEv

    print 'rho =', Pt.rho
    print 'minObj=', objRec.objIdeal

    hv = HyperVolume(refPoint)
    print 'hypervolume=', hv.compute(Pt.obj)

    ##ind = where((St.obj[:,0] != 0) & (St.obj[:,1]!=0))[0]
    ##j = random.choice(ind)
    ##a = (StObj[j] - objRec.objIdeal) / St.obj[j]
    ##ZPlot = Z*a + objRec.objIdeal
    ##
    ##for i in arange(0,len(Z),len(Z)/10):
    ##    plt.plot(vstack((objRec.objIdeal[0],ZPlot[i,0])),vstack((objRec.objIdeal[1],ZPlot[i,1])),'-',color='k')
    ##plt.plot(vstack((objRec.objIdeal[0],ZPlot[-1,0])),vstack((objRec.objIdeal[1],ZPlot[-1,1])),'-',color='k')
    ##inter = diag(a) + objRec.objIdeal
    ##plt.plot(inter[:,0],inter[:,1],'-',color='b')

    ##inter = diag(a) + objRec.objIdeal
    ##plt.plot(inter[:,0],inter[:,1],'-',color='b')

    with open(''.join(['Prt_',function,'.pk1']), 'r') as filename:
        f = pickle.load(filename)
    ##step = len(f)/50
    ##plt.scatter(f[::step,0],f[::step,1],s=1,color='b')
    plt.plot(f[:,0],f[:,1],color='b')

    plt.ioff()
    ##plt.show()

    plt.savefig(''.join([function,'_NSGA3_',str(nexec),'.png']), bbox_inches='tight')

##    plt.figure(2)
##    plt.plot(hvValues)
##    plt.savefig(''.join(['HV_',function,'_NSGA3_',nexec,'.png']), bbox_inches='tight')

    # Save Population
    with open(''.join(['Pop',function,'_NSGA3_',str(nexec),'.pk1']), 'wb') as output:
        pickle.dump(Pt, output, pickle.HIGHEST_PROTOCOL)

    # Save Hypervolume
    with open(''.join(['HV',function,'_NSGA3_',str(nexec),'.pk1']), 'wb') as output:
        pickle.dump(hvValues, output, pickle.HIGHEST_PROTOCOL)
