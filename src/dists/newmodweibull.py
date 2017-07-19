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

class newmodweibull(Distribution):
    @staticmethod
    def pdf(a,b,l,x):
        return a*(b+l*x)*math.pow(x,b-1)*math.exp(l*x)*math.exp(-a*math.pow(x,b)*math.exp(x*l))
    @staticmethod
    def cdf(a,b,l,x):
        return 1-math.exp(-a*math.pow(x,b)*math.exp(x*l))