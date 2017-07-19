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

class geometricextremeexp(Distribution): #geometric extreme-exponential
    @staticmethod
    def random(gmma):
        n=r.random()
        while n==1:
            n=r.random()
        return math.log((gmma)/(1-n)+1-gmma)