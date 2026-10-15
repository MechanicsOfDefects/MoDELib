import numpy as np
from readEVL import readEVLtxt

# Center and mean radius of each loop of a configuration file, from its loop nodes
class loop:
    def __init__(self,b,n):
        self.b = np.array(b)
        self.n = np.array(n)
        self._c = np.zeros(3)  # center
        self._r = 0.0          # mean distance of the loop nodes from the center

def getLoops(filename):
    evl=readEVLtxt(filename)
    llist={}
    for k in range(len(evl.loopsID)):
        loopID=int(evl.loopsID[k])
        l=loop(evl.burgers[k],evl.normal[k])
        P=evl.loopNodesPos[evl.loopNodesLoopID==loopID]
        if len(P)>0:
            l._c=P.mean(axis=0)
            l._r=np.sqrt(((P-l._c)**2).sum(axis=1)).mean()
        llist[loopID]=l
    return llist
