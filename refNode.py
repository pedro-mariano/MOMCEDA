from numpy import *

class refNode:
    #refNode 
    #   Tree node for generation of structured reference points
    
    #values : objective values of the branch
    #obj : height of the tree node
    #sumValues : sumValues of the branch until the node
        
    def __init__(self,values,obj,sumValues):
        
        self.values = values
        self.obj = obj
        self.sumValues = sumValues
        self.children = []
        
    def DFS(self, nObj,p,Z):
        # Build the tree with p divisions on each objective axe
        
        if(self.obj == (nObj - 1)):
            nRef = refNode(append(self.values,p - self.sumValues),self.obj+1,p)
            self.children.append(nRef)
            nValues = append(self.values[1:],p - self.sumValues)
            if(size(Z) > 0):
                Z = vstack((Z,nValues))
            else:
                Z = nValues
        else:
            values = linspace(0.0,p-self.sumValues,p+1-self.sumValues)
            for i in arange(size(values)):
                nRef = refNode(append(self.values,values[i]),self.obj+1,self.sumValues+values[i])
                self.children.append(nRef)
                Z = nRef.DFS(nObj,p,Z)

        return Z
                
            
  
        
        
    
    


