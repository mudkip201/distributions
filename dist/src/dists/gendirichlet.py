'''
Created on Jul 23, 2017

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

class gendirichlet(Distribution):
    @staticmethod
    def pdf(a,b,K,x):
        x0=1-np.sum(x)
        ff=math.pow(x0,b[-1][0]-1)
        for i in range(K):
            ff*=math.pow(x[i][0],a[i+1][0]-1)*math.pow(np.sum(x),b[i][0]-(a[i+1][0]+b[i+1][0]))/sp.beta(a[i+1][0],b[i+1][0])
        return ff
        pass
    @staticmethod
    def cdf():
        pass
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean():
        pass
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance():
        pass
    @staticmethod
    def stddev():
        pass
    @staticmethod
    def kurtosis():
        pass
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness():
        pass
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle():
        pass