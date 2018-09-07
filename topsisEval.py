from numpy import *

def topsisEval(matrixEval, minim, weight):
# Apply the TOPSIS algorithm for the entire population

#Input:
# matrixEval: objective values of members of the front
# weight: weight vector for objectives
# minim = 1 if minimizing objectives, minim = 0 otherwise

# Output:
# relSim: relative similarity for each individual

    # Normalization
    scale = sqrt((matrixEval**2).sum(axis=0))
    ind = where(scale != 0)[0]
    normMatrix = matrixEval[:,ind]/scale[ind]

    # Weighted normalized matrix
    sumWeight = sum(weight[ind])
    if(sumWeight > 0):
        weight = weight[ind]/sumWeight
    else:
        weight = weight[ind]
    normMatrix = normMatrix * weight

    # Best and worst cases
    if(minim == 1):
        mBest = normMatrix.min(0)
        mWorst = normMatrix.max(0)
    else:
        mBest = normMatrix.max(0)
        mWorst = normMatrix.min(0)

    # Calculate distances
    dBest = (normMatrix - mBest)**2
    dBest = sqrt(dBest.sum(axis=1))
    dWorst = (normMatrix - mWorst)**2
    dWorst = sqrt(dWorst.sum(axis=1))

    # Relative similarity
    ind = where((dBest == 0) & (dWorst == 0))[0]
    dBest[ind] = 1.0
    relSim = dWorst / (dBest + dWorst)

    return relSim
