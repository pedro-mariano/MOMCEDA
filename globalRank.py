def = globalRank(self):

    # Calculate the global rank for each member

    frontRank = zeros(size(self.fronts))
    frontRank[0] = 0
    for i in arange(1,size(self.fronts)):
        frontRank[i] = frontRank[i-1] + size(self.fronts[i-1])
    
    self.rank[:,2] = frontRank[self.rank[:,0]-1] + self.rank[:,1] 


