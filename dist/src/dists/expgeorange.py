'''
Created on Jul 24, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math

class expgeorange(Distribution):
    @staticmethod
    def pdf(b,p,x):
        return (1-p)*b*math.exp(-b*x)*math.pow(1-p*(1-math.exp(-b*x)),-2)
    @staticmethod
    def cdf(b,p,x):
        return (1-p)/p*(1/(1-p*(1-math.exp(-b*x)))-1)
    @staticmethod
    def random(b,p):
        u=ds.rg0()
        return 1/b*math.log((1-p+u*p)/(1+math.pow(p,u)-p-u))
    @staticmethod
    def mean(b,p):
        return -math.log(1-p)/(b*p)
    @staticmethod
    def median(b,p):
        return 1/b*math.log((2-p)/(1-p))
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