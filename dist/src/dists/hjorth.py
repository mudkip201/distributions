'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math

class hjorth(Distribution):
    @staticmethod
    def pdf(m,s,f,x):
        return math.exp(-x**2/(2*m**2)*math.pow(1+s*x,-1-f/s))*(f+x*(1+s*x)/m**2)
    @staticmethod
    def cdf(m,s,f,x):
        return 1-math.exp(-x**2/(2*m**2))*math.pow(1+s*x,-f/s)
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