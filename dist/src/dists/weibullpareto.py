'''
Created on Jul 17, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import numpy as np
from numpy import random as r
import scipy.stats as st
import scipy.special as sp
import scipy.optimize as op
import dists.weibull.weibull as weibull

class weibullpareto(Distribution):
    @staticmethod
    def pdf(c,g,k,t,x):
        if(x>t):
            return k*c/(g*x)*math.pow(k/g*math.log(x/t),c-1)*math.exp(math.pow(k/g*math.log(x/t),c))
        return None
    @staticmethod
    def cdf(c,g,k,t,x):
        return 1-math.exp(-math.pow(k/g*math.log(x/t)),c)
    @staticmethod
    def random(c,g,k,t):
        return t*math.exp(math.pow(-math.log(1-ds.rg0()),1/c)/(k/g))
    @staticmethod
    def mean():
        pass
    @staticmethod
    def median(c,g,k,t):
        return t*math.exp(math.pow(-math.log(1/2),1/c)/(k/g))
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
    def ppf(c,g,k,t,q):
        return t*math.exp(math.pow(-math.log(1-q),1/c)/(k/g))
        pass
    @staticmethod
    def mle():
        pass