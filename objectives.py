import numpy as np

def objectives(function,members,NPop,NObj):
    # Evaluate the objectives for each member of the population

    objValues = np.zeros((NPop,NObj))
    n = members.shape[1]

    if(function == 'ZDT1'):
        f1 = members[:,0]
        g = 1 + 9*np.sum(members[:,1:n],axis=1)/(n-1)
        h = 1 - np.sqrt(f1/g)
    elif(function == 'ZDT2'):
        f1 = members[:,0]
        g = 1 + 9*np.sum(members[:,1:n],axis=1)/(n-1)
        h = 1 - (f1/g)**2
    elif(function == 'ZDT3'):
        f1 = members[:,0]
        g = 1 + 9*np.sum(members[:,1:n],axis=1)/(n-1)
        h = 1 - np.sqrt(f1/g)-(f1/g)*np.sin(10*np.pi*f1)
    elif(function == 'ZDT4'):
        f1 = members[:,0]
        g = 1 + 10*(n-1)+ np.sum(members[:,1:n]**2 -10*np.cos(4*np.pi*members[:,1:n]),axis=1)
        h = 1 - np.sqrt(f1/g)
    elif(function == 'ZDT6'):
        f1 = 1 - np.exp(-4*members[:,0])*np.sin(6*np.pi*members[:,0])**6
        g = 1 + 9*(np.sum(members[:,1:n],axis=1)/9)**0.25
        h = 1 - (f1/g)**2
    elif(function == 'DebTiwari'):
        objValues[:,0] = np.sum(np.sin(np.pi*members),axis=1)
        objValues[:,1] = np.sum(np.cos(np.pi*members),axis=1)
    elif(function == 'EBN'):
        objValues[:,0] = np.sum(abs(members),axis=1) / n
        objValues[:,1] = np.sum(abs(members - 1),axis=1) / n
    elif(function == 'LameSuperSphere'):
        eps = np.sum(members[:,1:n],axis=1)/(n-1)
        r = np.sin(np.pi*eps)**2
        objValue[:,0] = (1+r)*np.cos(members[:,0])
        objValues[:,1] = (1+r)*np.sin(members[:,0])
    elif(function == 'TwoOnOne'):
        objValues[:,0] = members[:,0]**4 + members[:,1]**4 - members[:,0]**2 \
                         + members[:,1]**2 - 10* members[:,0] * members[:,1] \
                         + 0.25*members[:,0] + 20
        objValues[:,1] = (members[:,0] - 1)**2 + members[:,1]**2
    elif(function == '3Obj'):
        f1 = members[:,0]
        g = 1 + 9*np.sum(members[:,1:n],axis=1)/(n-1)
        h = 1 - np.sqrt(f1/g)
        objValues[:,0] = f1
        objValues[:,1] = g*h
        objValues[:,2] = np.random.rand(NPop)
    else:
        objValues[:,0] = np.random.rand(NPop)
        objValues[:,1] = np.random.rand(NPop)             

    if(function[:3] == 'ZDT'):
        objValues[:,0] = f1
        objValues[:,1] = g*h

    return objValues
