'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r
import dists.exponential.exponential as exponential
import dists.gamma.gamma as gamma

class gammashiftedgompertz(Distribution):
    @staticmethod
    def random(b,aa,bb):
        return max(exponential.random(b),b-gamma.random(aa,bb)*math.log(-math.log(r.random())))
