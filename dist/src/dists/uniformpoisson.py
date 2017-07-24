'''
Created on Jul 24, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp

class uniformpoisson(Distribution):
    @staticmethod
    def pdf(l,x):
        return math.pow(l,x)*math.exp(-l)/math.factorial(x+1)*sp.hyp1f1(1,x+2,l)
    @staticmethod
    def cdf(l,x):
        return 1-math.pow(l,x+1)*math.exp(-l)/math.factorial(x+2)*sp.hyp1f1(2,x+3,l)
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(l):
        return l/2
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(l):
        return l*(l+6)/12
    @staticmethod
    def stddev(l):
        return math.sqrt(l*(l+6)/12)
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