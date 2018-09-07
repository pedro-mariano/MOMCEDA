import cPickle as pickle
import matplotlib.pyplot as plt

function = 'ZDT3'

with open(''.join(['HV',function,'_NSGA2','.pk1']), 'r') as filename:
    hvNSGA2 = pickle.load(filename)

with open(''.join(['HV',function,'_NSGA3','.pk1']), 'r') as filename:
    hvNSGA3 = pickle.load(filename)
    
plt.plot(hvNSGA2,label='NSGA2')
plt.plot(hvNSGA3,label='NSGA3')
plt.legend(loc='lower right')
plt.xlabel('Generations')
plt.ylabel('Hypervolume')

plt.savefig(''.join(['HV_',function,'.png']), bbox_inches='tight')
