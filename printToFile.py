import cPickle as pickle
import Population
from metrics import diversity, convergence
from hv import *

algo = 'NSGA3'
function = 'ZDT1'
refPoint = (1.1,1.1)

with open(''.join(['Pop_',function,'_NSGA3','.pk1']),'r') as filename:
    Pop = pickle.load(filename)

with open(''.join(['Pareto/Prt_',function,'.pk1']), 'r') as filename:
    optimal_front = pickle.load(filename)

nreps = len(Pop)

nind = Pop[0].NPop
nobj = Pop[0].NObj

f = open(''.join(['Population_',algo,'_',function,'.txt']), 'w')
for nexec in xrange(nreps):
    
    for ind in xrange(nind):
        for obj in xrange(nobj):
            f.write(str(Pop[nexec].obj[ind][obj]))
            f.write('    ')
        f.write('\n')
    f.write('\n')
f.close()
