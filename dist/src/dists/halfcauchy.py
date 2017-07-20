'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r

class halfcauchy(Distribution):
    @staticmethod
    def random(x,gmma):
        return x+gmma*math.tan(math.pi*(r.random()/2.0))
    @staticmethod
    def median(x,gmma):
        return x+gmma*math.tan(math.pi*(1/4))
    @staticmethod
    def ppf(x,gmma,q):
        return x+gmma*math.tan(math.pi*(q/2.0))