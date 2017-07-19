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

class pickandsgenpareto(Distribution):
    @staticmethod
    def pdf(sigma,k,x):
        if(x>0 and (k*x)/sigma<1 and k!=0):
            return 1/sigma*math.pow(1-k*x/sigma,(1-k)/k)
        return None