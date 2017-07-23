'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp

class variancegamma(Distribution):
    @staticmethod
    def pdf(a,b,l,m,x):
        g=math.sqrt(a**2-b**2)
        return math.pow(g,2*l)*math.pow(abs(x-m),l-1/2)*sp.kv(l-1/2,a*abs(x-m))/(math.sqrt(math.pi)*sp.gamma(l)*math.pow(2*a,l-1/2))*math.exp(b*(x-m))
    @staticmethod
    def cdf():
        pass
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(a,b,l,m):
        g=math.sqrt(a**2-b**2)
        return m+2*b*l/g**2
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(a,b,l,m):
        g=math.sqrt(a**2-b**2)
        return 2*l*(1+2*b**2/g**2)/g**2
    @staticmethod
    def stddev(a,b,l,m):
        return math.sqrt(variancegamma.variance(a,b,l,m))
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