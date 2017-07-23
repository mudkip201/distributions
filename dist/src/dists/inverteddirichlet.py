'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import numpy as np
import scipy.special as sp

class inverteddirichlet(Distribution):
    @staticmethod
    def pdf(nu,K,x):
        ff=sp.gamma(np.sum(nu))
        for i in range(K+1):
            ff/=sp.gamma(nu[i][0])
        for i in range(K):
            ff*=math.pow(x[i][0],nu[i][0]-1)
        ff*=math.pow(1+np.sum(x),-np.sum(nu))
        return ff
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