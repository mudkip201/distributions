'''
Created on Jul 15, 2017

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
import dists.beta.beta as beta

class betaexp(Distribution): #beta-exponential
    @staticmethod
    def random(a,g):
        return math.log(beta.random(a,g))