'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math

class relativisticbreitwigner(Distribution):
    @staticmethod
    def pdf(M,G,E):
        g=math.sqrt(M**2*(M**2+G**2))
        k=2*math.sqrt(2)*M*G*g/(math.pi*math.sqrt(M**2+g))
        return k/((E**2-M**2)**2+M**2*G**2)
        pass
    @staticmethod
    def cdf():
        pass
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