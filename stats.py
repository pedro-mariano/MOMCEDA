import cPickle as pickle
import Population
from metrics import diversity, convergence
from hv import *
import json
import numpy as np

algo = 'NSGA3'
function = '2'
refPoint = (1.1,1.1)

with open(''.join(['results/Pop_ZDT',function,'_NSGA3','.pk1']),'r') as filename:
    Pop = pickle.load(filename)

with open(''.join(['pareto_front/zdt',function,'_front.json'])) as optimal_front_data:
        optimal_front = json.load(optimal_front_data)

nreps = len(Pop)
hvValues = np.zeros(nreps)
conv = np.zeros(nreps)
diver = np.zeros(nreps)

for nexec in xrange(nreps):
      
    print 'Execution ',nexec+1
    conv[nexec] = convergence(Pop[nexec].obj, optimal_front)
    print 'Convergence: ',conv[nexec]
    diver[nexec] = diversity(Pop[nexec].obj, optimal_front[0], optimal_front[-1])
    print 'Diversity: ', diver[nexec]
    hv = HyperVolume(refPoint)
    hvValues[nexec] = hv.compute(Pop[nexec].obj)
    print 'Hypervolume :', hvValues[nexec]

f = open(''.join(['results/Stats_',algo,'_ZDT',function,'.txt']), 'w')
f.write('Hypervolume: ');
f.write(''.join([str(np.mean(hvValues)),' ',str(np.std(hvValues)),'\n']))
f.write('Convergence: ');
f.write(''.join([str(np.mean(conv)),' ',str(np.std(hvValues)),'\n']))
f.write('Diversity: ');
f.write(''.join([str(np.mean(diver)),' ',str(np.std(hvValues)),'\n']))

f.close()
