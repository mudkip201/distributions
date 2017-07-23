'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r

class zerotruncpoisson(Distribution): #zero-truncated poisson
    @staticmethod
    def random(lmbda):
        k=1
        t=math.exp(-lmbda)/((1-math.exp(-lmbda))*lmbda)
        s=t
        u=r.random()
        while(s<u):
            k+=1
            t=t*lmbda/k
            s=s+t
        return k
