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
import dists.normal.normal as normal
import dists.gamma.gamma as gamma

class pearson7(Distribution):
    @staticmethod
    def random(s,m):
        return normal.random(0,s)/math.sqrt(gamma.random(1/2,m-1/2))
