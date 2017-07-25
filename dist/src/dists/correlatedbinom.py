'''
Created on Jul 25, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp

class correlatedbinom(Distribution):
    @staticmethod
    def pdf(n,p,phi,x):
        return sp.binom(n,x)*math.pow(p,x)*math.pow(1-p,n-x)*(1+phi/(2*p**2*(1-p)**2)*((x-n*p)**2+x*(2*p-1)-n*p**2))
    @staticmethod
    def cdf():
        pass
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(n,p,phi):
        return n*p
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(n,p,phi):
        return n*p*(1-p)+n*(n-1)*phi
    @staticmethod
    def stddev(n,p,phi):
        return math.sqrt(n*p*(1-p)+n*(n-1)*phi)
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