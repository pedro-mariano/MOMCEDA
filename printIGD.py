import cPickle as pickle
from Population import *
import numpy as np
from normalize import *
from objectiveRecords import *
from generateRefPoints import *
from metrics import *
import json

NObj = 2
minim = 1
p = 99
Z = generateRefPoints(NObj, p) # Generate structured reference points
objRec = objectiveRecords(NObj,minim) # Records for objective values


for i in ['1','2','4','6']:
    
    algo = 'NSGA3'
    
    with open(''.join(['results/Pop_ZDT',i,'_',algo,'.pk1']), 'r') as filename:
        Pop = pickle.load(filename)

    nreps = len(Pop)

    igdValues = np.zeros(nreps)
    normIgdValues = np.zeros(nreps)
    function = Pop[0].function

    for run in np.arange(nreps):

        Pt = Pop[run]
        objRec.objIdeal = Pt.obj.min(axis=0)
        Zr,a = normalize(Pt, objRec, minim, Z, p)
        ZRef = Zr*a + objRec.objIdeal

    ##    print ZRef

        normIgd,Zint = normIGDmetric(ZRef,objRec.objIdeal,a,Pt.obj,function)
        igd = IGDmetric(ZRef,objRec.objIdeal,Pt.obj,function)[0]

        normIgdValues[run] = normIgd
        igdValues[run] = igd
##
##        print 'Execution',run+1
##        print 'NormIGD = {0:e}'.format(normIgd)
##        print 'IGD = {0:e}'.format(igd)

    print "ZDT",i
    print algo
    print "Avg. Norm IGD: ",np.mean(normIgdValues),np.std(normIgdValues)
    print "Avg. IGD: ",np.mean(igdValues),np.std(igdValues)

    for algo in ['NSGA2','SPEA2','SMSEMOA']:
        
        with open(''.join(['results/',algo,'_population_',function,'.json']),'r') as filename:
            Pop = json.load(filename)

        nreps = len(Pop)

        igdValues = np.zeros(nreps)
        normIgdValues = np.zeros(nreps)

        for run in np.arange(nreps):

            Pt = Pop[run]
            obj = np.array(Pt)
            objRec.objIdeal = obj.min(axis=0)

            normIgd,Zint = normIGDmetric(ZRef,objRec.objIdeal,a,obj,function)
            igd = IGDmetric(ZRef,objRec.objIdeal,obj,function)[0]

            normIgdValues[run] = normIgd
            igdValues[run] = igd

        print algo
        print "Avg. Norm IGD: ",np.mean(normIgdValues),np.std(normIgdValues)
        print "Avg. IGD: ",np.mean(igdValues),np.std(igdValues)

    
