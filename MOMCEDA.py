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
from metrics import *
import time

# MOMCEDA main loop 

## Parameters ##

NPop = 100 # Population size
NEval = 20000 # Number of function evaluations
NGer = NEval/NPop - 1 # Number of generations
Csize = 30 # Chromosome size
NObj = 2 # Number of objectives to optimize
minim = 1 # minim = 1 if minimizing objectives, minim = 0 otherwise
function = 'ZDT1' # Type of problem

Vmin = 0.0*ones(Csize) # Limits of chromosome values
Vmax = 1.0*ones(Csize)
##Vmin = append(0,-5*ones(Csize-1))
##Vmax = append(1,5*ones(Csize-1))

minObj = array([0.0,0.0])# Limits of objective values
maxObj = array([1.0,1.0])

refPoint = (1.1,1.1) # Reference point for hypervolume indicator

objRec = objectiveRecords(NObj,minim) # Records for objective values

p = 99 # Number of objective axes divisions to generate structured ref. points
Z = generateRefPoints(NObj, p) # Generate structured reference points
rankType = 'hv' # Type of rank used for TOPSIS
weight = array([0.0,3.0,5.0,10.0,1.0]) # Array of weights used for TOPSIS
multiple = True #    multiple: use multiple criteria to select solutions from the last front
sampleAll = True #sampleAll: if true, all members from parent population are sampled with the TOPSIS rank
RTPlot = True # Plot the evolution of the population 
nReps = 1 #number of repetitions
 
# Weight order: domCount, distRefPoint, rho, rank, distPareto

hvValues = zeros((nReps,NGer))
finalPop = []
extime = []

# Plot parameters
color = plt.cm.get_cmap('Reds') 
ObjNames = ['f1', 'f2','f3']
scale = 1.0/(maxObj - minObj)
center = array([0.0,0.0,0.0])
countFig = 0 # counter of figures

## Offspring parameters

pMut = 1.0/Csize # Mutation probability
pSwitch = 0.5 # Probability to switch variables between members
spread = 0.5 # Parameter to control the spread of generated members
nc = 30
nm = 20
sigma = 1.0*ones(Csize)/2.0
##sigma = append(1.0,0.1*ones(Csize-1))
##sigma0 = r_[array([1.0]),1.0/10*ones(Csize-1)]
##sigmaf = r_[array([1.0/10]),1.0/500*ones(Csize-1)]
###deltaSigma = (sigma0 - sigmaf) / NGer
##qSigma = (sigmaf/sigma)**(1.0/NGer)
##sigma = sigma0

spreadf = 0.05
qSpread = (spreadf/spread)**(1.0/NGer)
deltaSpread = (spread - spreadf) / NGer

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


for nExec in arange(nReps,dtype=int):

    start = time.time()

    ## Initialization ##

    t = 0 # Counter of generations
    Pt = Population(NPop, Csize, Vmin, Vmax, NObj, minim, function) # Initial Population
    Pt.fastNonDominatedSort() # Nondominated sorting
    Pt.topsisPop() # Rank within each front
    Qt = Population(NPop, Csize, Vmin, Vmax, NObj, minim, function) # Offspring population

    plt.ion()

    ## Main loop ##

    while(t < NGer):

##        if(t == 150):
##            weight = array([0.0,1.0,3.0,10.0,5.0]) # Array of weights used for TOPSIS
        
        Rt = Pt.addPopulation(Qt) # Combined population
        Rt.fastNonDominatedSort() # Nondominated sorting
        ##    Rt.topsisPop(rank=rankType) # Rank within each front
        ##    Rt.globalRankEval()

        if(sampleAll):

            RtObj = Rt.obj # Save original objective values of St
            Zr,a = normalize(Rt, objRec, minim, Z, p)
            ZRef = Zr*a + objRec.objIdeal
            associate(Rt, Zr,0)
##            distPareto(Rt,ZRef,objRec.objIdeal,a)
            hvContribution(Rt,objRec.objIdeal,a)
            indRemove = niching(NPop,Rt,array([0,2*NPop]),weight,multiple)
            Pt = Rt.removeMembers(indRemove, RtObj)

        else:
         
            indList = zeros(NPop, dtype = int) # List of indexes of members to the next population
            i = 0 # Counter of fronts
            sizeEvol = array([0, len(Rt.fronts[i])]) # Evolution of population's size by adding the fronts
            
            # Fill population with the first fronts
            while(sizeEvol[i+1] <= NPop):
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
            Zr,a = normalize(St, objRec, minim, Z, p)
            ZRef = Zr*a + objRec.objIdeal
            associate(St, Zr, sizeEvol[i])
            ZRef = Zr*a + objRec.objIdeal
        ##    distPareto(St,ZRef,objRec.objIdeal,a)
            hvContribution(St,objRec.objIdeal,a)
            indRemove = niching(K, St, sizeEvol, weight,multiple)
            Pt = St.removeMembers(indRemove, StObj) # Next generation's population
        
        if(RTPlot):
            Pt.plot(color(float(t)/NGer), scale, center, ObjNames, countFig)
            axes = plt.gca()
            axes.set_ylim([0,1])
        Pt.topsisPop(weight,rank=rankType) # Rank within each front
        Qt = Pt.offspringPop(coefGau,sigma,pMut,spread,pSwitch) # Selection, recombination and mutation

        ##    dist[t] = maxDist(Pt.members)/(1*(Vmax - Vmin))
        ##    sigma = dist[t]
##        sigma = random.rand()*append(10.0,1.0*ones(Csize-1))
        
        ##    MEv = sqrt(sum(Pt.members[:,1:]**2)/(NPop*(Csize-1)))
        ##    print 'Generation = ', t, 'Mean Square Value =', MEv

        hv = HyperVolume(refPoint)
        hvValues[nExec,t] = hv.compute(Pt.obj)
        
        t = t + 1
        ##    sigma = sigma - deltaSigma
##        sigma = sigma*qSigma

##        spread = spread*qSpread
        spread = spread - deltaSpread

    end = time.time()
    extime.append(end - start)
        
    print 'Execution ',nExec+1,' - hypervolume = ',hvValues[nExec,-1]
    print 'Elapsed time: ',extime[-1]

    import json
    with open(''.join(['pareto_front/zdt',function[3],'_front.json'])) as optimal_front_data:
        optimal_front = json.load(optimal_front_data)

    print 'Convergence: ',convergence(Pt.obj.tolist(), optimal_front)

##    normIgd,Zint = normIGDmetric(ZRef,objRec.objIdeal,a,Pt.obj,function)
##    igd = IGDmetric(ZRef,objRec.objIdeal,Pt.obj,function)[0]

##    normIgd2,Zint = normIGDmetric2(ZRef,objRec.objIdeal,a,Pt.obj,function)
##    igd2 = IGDmetric2(ZRef,objRec.objIdeal,Pt.obj,function)[0]

##    print 'NormIGD = {0:e}'.format(normIgd)
##    print 'IGD = {0:e}'.format(igd)
##    print 'NormIGD2 = {0:e}'.format(normIgd2)
##    print 'IGD2 = {0:e}'.format(igd2)
       
    countFig = countFig + NObj*(NObj - 1)/2
    with open(''.join(['Pareto/Prt_',function,'.pk1']), 'r') as filename:
        f = pickle.load(filename)
    ##step = len(f)/50
    ##plt.scatter(f[::step,0],f[::step,1],s=1,color='b')
    Pt.plot(color(float(t)/NGer), scale, center, ObjNames, countFig)
    plt.plot(f[:,0],f[:,1],color='b')
    plt.ioff()

    finalPop.append(Pt)

if(nReps == 1):     
    MEv = sqrt(sum(Pt.members[:,1:]**2)/(NPop*(Csize-1)))
    print 'Generation = ', t, 'Mean Square Value =', MEv

    print 'rho =', Pt.rho
    print 'minObj=', objRec.objIdeal

    ##refPoint = objRec.extPoints.max(axis=0)*1.1
    ##print 'refPoint: ',refPoint
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

    ##plt.show()
    axes = plt.gca()
    axes.set_ylim([0,5])

    plt.savefig(''.join([function,'.png']), bbox_inches='tight')

    plt.figure(2)
    plt.plot(hvValues)
    plt.savefig(''.join(['HV_',function,'.png']), bbox_inches='tight')

    # Save Population
    with open(''.join(['Pop_',function,'_NSGA3','.pk1']), 'wb') as output:
        pickle.dump(finalPop, output, pickle.HIGHEST_PROTOCOL)

    # Save Hypervolume
    with open(''.join(['HV_',function,'_NSGA3','.pk1']), 'wb') as output:
        pickle.dump(hvValues, output, pickle.HIGHEST_PROTOCOL)

    # Save Elapsed time
    with open(''.join(['time_',function,'_NSGA3','.pk1']), 'wb') as output:
        pickle.dump(extime, output, pickle.HIGHEST_PROTOCOL)

else:

    print 'average hypervolume=', mean(hvValues[:,-1])
    print 'best hypervolume=', max(hvValues[:,-1])

    plt.figure(countFig)
    plt.plot(hvValues.mean(axis=0))
    plt.savefig(''.join(['meanHV_',function,'.png']), bbox_inches='tight')

    # Save Population
    with open(''.join(['Pop_',function,'_NSGA3','.pk1']), 'wb') as output:
        pickle.dump(finalPop, output, pickle.HIGHEST_PROTOCOL)

    # Save Hypervolume
    with open(''.join(['HV_',function,'_NSGA3','.pk1']), 'wb') as output:
        pickle.dump(hvValues, output, pickle.HIGHEST_PROTOCOL)

    # Save Elapsed time
    with open(''.join(['time_',function,'_NSGA3','.pk1']), 'wb') as output:
        pickle.dump(extime, output, pickle.HIGHEST_PROTOCOL)

