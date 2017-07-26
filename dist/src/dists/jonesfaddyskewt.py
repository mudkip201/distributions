'''
Created on Jul 26, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp

class jonesfaddyskewt(Distribution): #Jones and Faddy's skewed t
    @staticmethod
    def pdf(m,s,a,b,x):
        z=(x-m)/s
        v=a+b
        C=math.pow(2,x-1)*sp.beta(a,b)*math.sqrt(v)
        return math.pow(1+z/math.sqrt(v+z**2),a+1/2)/C*math.pow(1-z/math.sqrt(v+z**2),b+1/2)
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