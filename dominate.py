def dominate(x,y,minim):
# x,y: array of objective values
# minim = 1 if minimizing objectives, minim = 0 otherwise
# dom = 1 if x dominates y, dom = 0 otherwise

    dom = 0
    arraySize = len(x)

    if(minim):
        p = (sum(x <= y) == arraySize) # p = prod(x <= y)
        s = sum(x < y)
        if( p and s ):
            dom = 1
    
    else:
        p = (sum(x >= y) == arraySize) # p = prod(x >= y)
        s = sum(x > y)
        if( p and s ):
            dom = 1
    return dom
