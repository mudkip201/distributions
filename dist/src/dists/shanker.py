'''
Created on Jul 18, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math

class shanker(Distribution):
    @staticmethod
    def pdf(t,x):
        return t**2/(t**2+1)*(t+x)*math.exp(-t*x)
    @staticmethod
    def cdf(t,x):
        return 1-((t**2+1)+t*x)/(t**2+1)*math.exp(-t*x)
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(t):
        return (t**2+2)/(t*(t**2+1))
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(t):
        return (t**4+4*t**2+2)/(t**2*(t**2+1)**2)
    @staticmethod
    def stddev(t):
        return math.sqrt(shanker.variance(t))
    @staticmethod
    def kurtosis(t):
        return 3*(3*t**8+24*t**6+44*t**4+32*t**2+8)/(t**4*(t**2+1)**4)
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness(t):
        return 2*(t**6+6*t**4+6*t**2+2)/(t**3*(t**2+1)**3)
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle():
        pass