'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp

class studentsz(Distribution):
    @staticmethod
    def pdf(n,z):
        return sp.gamma(n/2)/(math.sqrt(math.pi)*sp.gamma((n-1)/2))*math.pow(1+z**2,-n/2)
    @staticmethod
    def cdf(n,z):
        d=math.pow(abs(z),1-n)*sp.gamma(n/2)*sp.hyp2f1((n-1)/2,n/2,(n+1)/2,-1/z**2)/(2*math.sqrt(math.pi)*sp.gamma((n+1)/2))
        if z>0:
            return 1-d
        return d
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(n):
        return 0
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(n):
        return 1/(n-3)
    @staticmethod
    def stddev(n):
        return 1/math.sqrt(n-3)
    @staticmethod
    def kurtosis(n):
        return 6/(n-5)
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness(n):
        return 0
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle():
        pass