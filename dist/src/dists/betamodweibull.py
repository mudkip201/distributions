'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp

class betamodweibull(Distribution):
    @staticmethod
    def pdf(a,b,aa,g,l,x):
        return (aa*math.pow(x,g-1)*(g+l*x)*math.exp(l*x))/sp.beta(a,b)*math.pow(1-math.exp(-aa*math.pow(x,g)*math.exp(l*x)),a-1)*math.exp(-b*aa*math.pow(x,g)*math.exp(l*x))
    @staticmethod
    def cdf(a,b,aa,g,l,x):
        return sp.betainc(a,b,1-math.exp(-aa*math.pow(x,g)*math.exp(l*x)))
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