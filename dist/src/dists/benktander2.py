'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp

class benktander2(Distribution):
    @staticmethod
    def pdf(a,b,x):
        return math.exp(a/b*(1-math.pow(x,b)))*math.pow(x,b-2)*(a*math.pow(x,b)-b+1)
    @staticmethod
    def cdf(a,b,x):
        return 1-math.pow(x,b-1)*math.exp(a/b*(1-math.pow(x,b)))
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(a,b):
        return 1+1/a
    @staticmethod
    def median(a,b):
        if(b==1):
            return math.log(2)/a+1
        return math.pow((1-b)/a*sp.lambertw(math.pow(2,b/(1-b))*a*math.exp(a/(1-b))/(1-b)),1/b)
    @staticmethod
    def mode(a,b):
        return 1
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