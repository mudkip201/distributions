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
import dists.gammaexp.gammaexp as gammaexp

class gengumbel(Distribution): #generalized gumbel
    @staticmethod
    def random(u,lmbda,n):
        return gammaexp.random(u-lmbda*math.log(n),lmbda,n)