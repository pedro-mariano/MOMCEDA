import numpy as np
import matplotlib.pyplot as plt
import json

function = 'ZDT1'

NPoints = 1000

optimal_front = np.zeros((NPoints,2));

f1 = np.linspace(0,1.0,NPoints)

if(function == 'ZDT1' or function == 'ZDT4'):

    f2 = 1 - np.sqrt(f1)

elif(function == 'ZDT2' or function == 'ZDT6'):

    if(function == 'ZDT6'):
        
        f1 = np.linspace(0.2807753191,1.0,NPoints)

    f2 = 1 - f1**2

elif(function == 'ZDT3'):
    f1[:NPoints/5] = np.linspace(0,0.0830015349,NPoints/5)
    f1[NPoints/5:2*NPoints/5] = np.linspace(0.1822287280,0.2577623634,NPoints/5)
    f1[2*NPoints/5:3*NPoints/5] = np.linspace(0.4093136748,0.4538821041,NPoints/5)
    f1[3*NPoints/5:4*NPoints/5] = np.linspace(0.6183967944,0.6525117038,NPoints/5)
    f1[4*NPoints/5:] = np.linspace(0.8233317983,0.8518328654,NPoints/5)

    f2 = 1 - np.sqrt(f1) -f1*np.sin(10*np.pi*f1)

    
optimal_front[:,0] = f1
optimal_front[:,1] = f2

optimal_front = optimal_front.tolist()

with open(''.join(['Pareto/Prt_',function,'.json']),'w') as outfile:
    json.dump(optimal_front,outfile)
