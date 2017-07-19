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
import dists.invgamma.invgamma as invgamma
import dists.normal.normal as normal

class norminvgamma(Distribution): #normal-inverse gamma
    @staticmethod
    def random(aa,bb,mu,lmbda):
        n=invgamma.random(aa,bb)
        return normal.random(mu,math.sqrt(n/lmbda))
