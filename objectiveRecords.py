import numpy as np

class objectiveRecords:
    # objectiveRecords
    #   Keeps track of objective values ever found since the start of
    #   the simulation
    
    #objIdeal:  Ideal values
    #minASF:  ASF minimum values
    #extPoints:  Extreme points coordinates
    
    def __init__(self,NObj,minim):
        
        # Standard case: Minimizing objectives
        self.objIdeal = np.inf*np.ones(NObj)
        self.minASF = np.inf*np.ones(NObj)
        self.extPoints = np.inf*np.ones((NObj,NObj))
            
        # Maximizing objectives 
        if(minim == 0):
            self.objIdeal = -self.objIdeal
        
            
        
        
        
    
    


