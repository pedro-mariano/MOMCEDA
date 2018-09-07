from numpy import *
from generateRefPoints import *
import matplotlib.pyplot as plt
from collections import Counter
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter

def normalize(Pop, objRec, minim, Z, p=0):
# normalize 
#    Normalization procedure of the objective values

#    Input:
#    Pop: population
#    objRec: records for objective values
#    minim = 1 if minimizing objectives, minim = 0 otherwise
#    p = 0 if there are supplied reference points
#    p: number of divisions for each objective axis to generate structured
#       reference points
#    Z : supplied reference points (if given)

#    Output:
#    Zr: reference points on normalized hyper-plane


## Update ideal point

     #   fn: normalized objectives
     if(minim == 1):
         objRec.objIdeal = vstack((Pop.obj,objRec.objIdeal)).min(axis=0)
         fn = Pop.obj - objRec.objIdeal
         extPoints = objRec.extPoints - objRec.objIdeal
     else:
         objRec.objIdeal = vstack((Pop.obj,objRec.objIdeal)).max(axis=0)
         fn = objRec.objIdeal - Pop.obj
         extPoints = objRec.objIdeal - objRec.extPoints
     
## Update extreme points

     epsilon = 1e-6
     w = epsilon*ones((Pop.NObj,Pop.NObj)) #weight matrix
     axisWeight = 1.0
     fill_diagonal(w,axisWeight)

     # Achievement scalarizing function
     ASF = zeros((Pop.NPop,Pop.NObj))
     for i in arange(Pop.NObj):
          objRec.minASF[i] = max(extPoints[i]/w[i]) # Update min ASF
          s = fn/w[i]
          ASF[:,i] = s.max(axis=1)

     nASF = ASF.min(axis=0) # New ASF values
     ind = ASF.argmin(axis=0)
     upObj = where(nASF < objRec.minASF)[0] # Keep the new values that are better
     objRec.minASF[upObj] = nASF[upObj]
     objRec.extPoints[upObj] = Pop.obj[ind[upObj]]
          
     # Verification for repeated extreme points
     rep = [item for item, count in Counter(ind[upObj]).iteritems() if count > 1]
     if(len(rep) > 0):
          for i in rep:
               upPos = where(ind[upObj] == i)[0]
               upPos = upObj[upPos]
               objRec.extPoints[upPos,upPos] = objRec.extPoints[upPos,upPos]*1.1

##     test = Pop.obj[ind]
     
## Find intercepts and normalize

     if(minim == 1):
         extPoints = objRec.extPoints - objRec.objIdeal
     else:
         extPoints = objRec.objIdeal - objRec.extPoints

     # Verification for null extreme points
     indNull = where(extPoints.sum(axis=1) == 0)[0]
     if(len(indNull) > 0):
          for i in indNull:
               nObj = delete(Pop.obj[:,i],where(Pop.obj[:,i] == min(Pop.obj[:,i]))[0])
               extPoints[i,i] = min(nObj)
          if(minim == 1):
               objRec.extPoints = extPoints + objRec.objIdeal
          else:
               objRec.extPoints = objRec.objIdeal - extPoints
     
     b = ones(Pop.NObj)
     a = 1.0/linalg.solve(extPoints,b)

     if(any(a < 0)):
          print 'Negative Intercept! a = ',a
     
     fn = fn/a

     Pop.obj = fn
     
## Generate structured reference points or map given points on the hyperplane
     
     if(p > 0):
         Zr = Z
     else:
         if(minim == 1):
              Z = Z - objRec.objIdeal
         else:
              Z = objRec.objIdeal - Z
        
         Zr = Z/a
    
## Plot

##     fig = plt.figure()
##     ax = fig.add_subplot(111, projection='3d')
##
##     inter = eye(Pop.NObj)
##     inter = diag(a+objRec.objIdeal)
##     ZrPlot = Zr*a + objRec.objIdeal
     
##     plt.plot(fn[:,0],fn[:,1],'.')
##     plt.plot(objRec.extPoints[:,0]/a[0],objRec.extPoints[:,1]/a[1],'o')
##     plt.plot(objRec.extPoints[:,0],objRec.extPoints[:,1],'x',color='b')
##     plt.plot(test[:,0],test[:,1],'o',color='k')
##     plt.plot(inter[:,0],inter[:,1],'x',color='k')
##     for i in arange(len(Zr)):
##          plt.plot(vstack((objRec.objIdeal[0],ZrPlot[i,0])),vstack((objRec.objIdeal[1],ZrPlot[i,1])),'-',color='k')
##     plt.show()

##     plt.plot(fn[:,0],fn[:,1],fn[:,2],'.')
##     plt.plot(objRec.extPoints[:,0]/a[0],objRec.extPoints[:,1]/a[1],objRec.extPoints[:,2]/a[2],'o')
####     plt.plot(objRec.extPoints[:,0],objRec.extPoints[:,1],objRec.extPoints[:,2],'o')
##     plt.plot(inter[:,0],inter[:,1],inter[:,2],'x')
##     plt.show()
     
#      [X1,Y1] = meshgrid(0:0.01:1)
#      Z1 = 1 - X1 - Y1
#      plt.surf(X1,Y1,Z1)
#      plt.show()
     

     return Zr,a

