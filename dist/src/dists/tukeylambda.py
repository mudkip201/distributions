'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r

class tukeylambda(Distribution):
    @staticmethod
    def random(lmbda):
        p=r.random()
        if(lmbda==0):
            while(p==1):
                p=r.random()
            return math.log(p/(1-p))
        return 1/lmbda*(math.pow(p,lmbda)-math.pow(1-p,lmbda))
    @staticmethod
    def median(lmbda):
        return 0
    @staticmethod
    def ppf(lmbda,q):
        if(lmbda==0):
            return math.log(q/(1-q))
        return 1/lmbda*(math.pow(q,lmbda)-math.pow(1-q,lmbda))