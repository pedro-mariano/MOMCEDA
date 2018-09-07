import cPickle as pickle
import numpy as np
import json

algo = 'NSGA3'
function = 'ZDT6'
with open(''.join(['results/HV_',function,'_',algo,'.pk1']),'r') as filename:
    hv_NSGA3 = pickle.load(filename)

algo = 'NSGA2'
with open(''.join(['results/',algo,'_hvValues_',function,'.json']),'r') as filename:
    hv_NSGA2 = np.asarray(json.load(filename))

algo = 'SPEA2'
with open(''.join(['results/',algo,'_hvValues_',function,'.json']),'r') as filename:
    hv_SPEA2 = np.asarray(json.load(filename))

algo = 'SMSEMOA'
with open(''.join(['results/',algo,'_hvValues_',function,'.json']),'r') as filename:
    hv_SMSEMOA = np.asarray(json.load(filename))

algo = 'proposal'
f = open(''.join(['HV_',algo,'_',function,'.txt']), 'w')

nreps = len(hv_NSGA3)
nger = len(hv_NSGA3[0])

for nexec in np.arange(nreps):
    for ind in np.arange(nger):
        f.write(str(hv_NSGA3[nexec][ind]))
        f.write('\n')
    f.write('\n')
f.close()

algo = 'NSGA2'
f = open(''.join(['HV_',algo,'_',function,'.txt']), 'w')

nreps = len(hv_NSGA2)
nger = len(hv_NSGA2[0])

for nexec in np.arange(nreps):
    for ind in np.arange(nger):
        f.write(str(hv_NSGA2[nexec][ind]))
        f.write('\n')
    f.write('\n')
f.close()

algo = 'SPEA2'
f = open(''.join(['HV_',algo,'_',function,'.txt']), 'w')

nreps = len(hv_SPEA2)
nger = len(hv_SPEA2[0])

for nexec in np.arange(nreps):
    for ind in np.arange(nger):
        f.write(str(hv_SPEA2[nexec][ind]))
        f.write('\n')
    f.write('\n')
f.close()

algo = 'SMSEMOA'
f = open(''.join(['HV_',algo,'_',function,'.txt']), 'w')

nreps = len(hv_SMSEMOA)
nger = len(hv_SMSEMOA[0])

for nexec in np.arange(nreps):
    for ind in np.arange(nger):
        f.write(str(hv_SMSEMOA[nexec][ind]))
        f.write('\n')
    f.write('\n')
f.close()
