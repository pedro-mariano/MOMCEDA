from Population import *
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import cPickle as pickle

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

NPop = 1024 # Population size
NGer = 1000 # Number of generations
NObj = 3 # Number of objectives to optimize
k = 10
Csize = NObj + k - 1 # Chromosome size
minim = 1 # minim = 1 if minimizing objectives, minim = 0 otherwise
function = 'DTLZ2' # Type of problem
Vmin = 0.0
Vmax = 1.0
Pop = Population(NPop, Csize, k, Vmin, Vmax, NObj, minim, function)

Pop.members[:,-k:] = 0.5
div = NPop**(1.0/(NObj-1))
firstMembers = linspace(0,1.0,div)
f0,f1 = meshgrid(firstMembers,firstMembers)
f2 = np.zeros(f0.shape)
for i in np.arange(f0.shape[0]):
    for j in np.arange(f0.shape[1]):
        if(f0[i,j]**2 + f1[i,j]**2 <= 1.0):
            f2[i,j] = np.sqrt(1.0 - f0[i,j]**2 - f1[i,j]**2)
        else:
            f2[i,j] = 0.0
##Pop.obj = objectives(function,Pop.members,NPop,NObj,k)
##Pop.fastNonDominatedSort()
##
##ndom = where(Pop.rank == 1)[0]
##f = c_[Pop.obj[ndom,0],Pop.obj[ndom,1],Pop.obj[ndom,2]]
ax.plot_surface(f0,f1,f2)
f = array([f0,f1,f2])
plt.show()

with open(''.join(['Prt_',function,'.pk1']), 'wb') as output:
    pickle.dump(f, output, pickle.HIGHEST_PROTOCOL)
