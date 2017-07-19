'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Dist
import dists.Distribution as ds
import math
import numpy as np
from numpy import random as r
import scipy.stats as st
import scipy.special as sp
import scipy.optimize as op
import dists.poisson.poisson as poisson
import dists.geometric.geometric as geometric

class bartlett(Dist):
    @staticmethod
    def random(a,q):
        return poisson.random(a)*geometric.random(q)