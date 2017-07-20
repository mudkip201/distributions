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
import dists.gamma.gamma as gamma

class prentice(Distribution):
    @staticmethod
    def random(z,l,a,g):
        return z+l*math.log(gamma.random(g,1)/gamma.random(a,1))