from numpy import *

def coefGau(NPop,x=1):
# Coefficients of the gaussians of the mixture

    # Linear decay
    linear = 2.0*(NPop-arange(1,NPop+1))/(NPop*(NPop-1)) 

    # Exponential decay
    # desired value for k = coefGau(n)/coefGau(1)
    k = 10**(-3) 
    expo = ((1.0-k**(1.0/(NPop-1)))*k**(arange(1,NPop+1,dtype='float')/(NPop-1)))/(k**(1.0/(NPop-1))*(1.0-k**(float(NPop)/(NPop-1))))

    # Logarithmic decay
    loga = log(NPop-arange(1,NPop+1)+1)/math.log(math.factorial(NPop))
    
    return {
        1: linear,
        2: expo,
        3: loga,
    }[x]
