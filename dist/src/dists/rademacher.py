'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r

class rademacher(Distribution):
    @staticmethod
    def pdf(k):
        if(k!=-1 or k!=1):
            raise ValueError("k must be -1 or 1")
        return 1/2
    @staticmethod
    def cdf(k):
        if(k!=-1 or k!=1):
            raise ValueError("k must be -1 or 1")
        return 1/2+(1+k)/4
    @staticmethod
    def random():
        n=r.random()
        if(n<=0.5):
            return -1
        return 1
    @staticmethod
    def kurtosis():
        return -2
    @staticmethod
    def mean():
        return 0
    @staticmethod
    def median():
        return 0
    @staticmethod
    def mode():
        return None
    @staticmethod
    def variance():
        return 1
    @staticmethod
    def stddev():
        return 1
    @staticmethod
    def entropy():
        return math.log(2)
    @staticmethod
    def skewness():
        return 0
    @staticmethod
    def ppf(q):
        if(q<=0.5):
            return -1
        return 1
