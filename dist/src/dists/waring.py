'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp

class waring(Distribution):
    @staticmethod
    def pdf(a,b,k):
        return a*sp.poch(b,k)/sp.poch(a+b,1+k)
    @staticmethod
    def cdf(a,b,k):
        return 1-sp.poch(b,1+math.floor(k))/sp.poch(a+b,1+math.floor(k))
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(a,b):
        if(a>1):
            return b/(-1+a)
        return None
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