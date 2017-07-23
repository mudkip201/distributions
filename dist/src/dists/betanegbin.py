'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import numpy as np
import scipy.special as sp

class betanegbin(Distribution):
    @staticmethod
    def pdf(a,b,n,k):
        return sp.poch(n,k)*sp.poch(a,n)*sp.poch(b,k)/(math.factorial(k)*sp.poch(a+b,n)*sp.poch(n+a+b,k))
    @staticmethod
    def cdf():
        pass
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(a,b,n):
        if(a>1):
            return n*b/(a-1)
        return None
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(a,b,n):
        if(a>2):
            return n*(n+a-1)*b*(a+b-1)/((a-2)*(a-1)**2)
        return np.Infinity
    @staticmethod
    def stddev(a,b,n):
        return math.sqrt(betanegbin.variance(a,b,n))
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