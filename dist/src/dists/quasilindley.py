'''
Created on Jul 17, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math

class quasilindley(Distribution):
    @staticmethod
    def pdf(a,t,x):
        return t/(a+1)*(a+t*x)*math.exp(-t*x)
    @staticmethod
    def cdf(a,t,x):
        return 1-(1+a+t*x)/(a+1)*math.exp(-t*x)
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
    def mode(a,t):
        if(abs(a)<1):
            return (1-a)/t
        return 0
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