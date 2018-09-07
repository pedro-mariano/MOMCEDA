import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt

function = 'ZDT3'
nExec = 10

with open(''.join(['HV_',function,'_NSGA3','.pk1']), 'r') as filename:
    hv = pickle.load(filename)


meanHV = hv.mean(axis=0)
NGer = len(meanHV)
stdHV = hv.std(axis=0)
pos = np.append(np.arange(0,NGer,25),NGer-1)
plt.boxplot(hv[:,pos],positions=pos,widths=5)
plt.plot(meanHV)
axes = plt.gca()
axes.set_xlim([-10,NGer+10])
plt.show()
print function, '- mean hypervolume: ', np.mean(hv[:,-1]),'+/-',np.std(hv[:,-1])
