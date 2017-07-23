'''
Created on Jul 18, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.normal.normal as normal

class modburr3normal(Distribution):
    @staticmethod
    def pdf(a,b,g,m,s,x):
        G=normal.cdf(m,s,x)
        gg=normal.pdf(m,s,x)
        return a*b*math.pow(G/(1-G),-b+1)*math.pow(1+g*math.pow(G/(1-G),-b),-a/g-1)*gg/G**2
    @staticmethod
    def cdf(a,b,g,m,s,x):
        return math.pow(1+g*math.pow(normal.cdf(m,s,x)/(1-normal.cdf(m,s,x)),-b),-a/g)
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