'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r

class geometricextremeexp(Distribution): #geometric extreme-exponential
    @staticmethod
    def random(gmma):
        n=r.random()
        while n==1:
            n=r.random()
        return math.log((gmma)/(1-n)+1-gmma)