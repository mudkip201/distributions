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