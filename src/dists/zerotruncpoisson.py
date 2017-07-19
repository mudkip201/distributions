'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import numpy as np
from numpy import random as r
import scipy.stats as st
import scipy.special as sp
import scipy.optimize as op

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
