'''
Created on Jul 18, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math

class modweibullvar(Distribution):
    @staticmethod
    def pdf(a,b,g,x):
        return (a+b*g*math.pow(x,g-1))*math.exp(-a*x-b*math.pow(x,g))
    @staticmethod
    def cdf(a,b,g,x):
        return 1-math.exp(-a*x-b*math.pow(x,g))
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