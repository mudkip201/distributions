'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp

class polyaaeppli(Distribution):
    @staticmethod
    def pdf(t,p,n):
        if n==0:
            return math.exp(-t)
        return math.exp(-t)*(1-p)*math.pow(p,n-1)*t*sp.hyp1f1(1-n,2,-(1-p)*t/p)
    @staticmethod
    def cdf():
        pass
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(t,p):
        return t/(1-p)
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(t,p):
        return (1+p)*t/(1-p)**2
    @staticmethod
    def stddev(t,p):
        return math.sqrt((1+p)*t)/(1-p)
    @staticmethod
    def kurtosis(t,p):
        return (p**2+10*p+1)/(t*(p+1))+3
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness(t,p):
        return (p**2+4*p+1)/((p+1)*math.sqrt(t*(p+1)))
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle():
        pass