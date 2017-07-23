'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp

class davis(Distribution):
    @staticmethod
    def pdf(b,n,m,x):
        return math.pow(b,n)*math.pow(x-m,-n-1)/((math.exp(b/(x-m))-1)*sp.gamma(n)*sp.zetac(n))
    @staticmethod
    def cdf():
        pass
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(b,n,m):
        if(n>2):
            return m+b*sp.zetac(n-1)/((n-1)*sp.zetac(n))
        return None
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(b,n,m):
        if(n>3):
            return b**2*(-(n-2)*sp.zetac(n-1)**2+(n-1)*sp.zetac(n-2)*sp.zetac(n))/((n-2)*(n-1)**2*sp.zetac(n)**2)
        return None
    @staticmethod
    def stddev(b,n,m):
        if(n>3):
            return math.sqrt(b**2*(-(n-2)*sp.zetac(n-1)**2+(n-1)*sp.zetac(n-2)*sp.zetac(n))/((n-2)*(n-1)**2*sp.zetac(n)**2))
        return None
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