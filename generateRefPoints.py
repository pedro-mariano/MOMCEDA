from refNode import *

def generateRefPoints(nObj, p):
#generateRefPoints
#   Generate structured reference points on the normalized hyperplane with
#   p divisions in each objective axe
#   Using Das and Dennis's systematic approach

    Z = []

    tree = refNode(array([0]),0,0)
    Z = tree.DFS(nObj,p,Z)/p

    return Z
