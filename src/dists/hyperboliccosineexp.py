'''
Created on Jul 18, 2017

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

class hyperboliccosineexp(Distribution):
    @staticmethod
    def pdf(a,l,x):
        return 2*a*math.exp(a)/(math.exp(2*a)-1)*l*math.exp(-l*x)*math.cosh(a*(1-math.exp(-l*x)))
    @staticmethod
    def cdf(a,l,x):
        return 2*math.exp(a)/(math.exp(2*a)-1)*math.sinh(a*(1-math.exp(-l*x)))
    @staticmethod
    def random(a,l):
        u=ds.rg0()
        return -1/l*math.log(1-math.asinh((math.exp(2*a)-1)/(2*math.exp(a))*u)/a)
    @staticmethod
    def mean():
        pass
    @staticmethod
    def median(a,l):
        return -1/l*math.log(1-math.asinh((math.exp(2*a)-1)/(2*math.exp(a))/2)/a)
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
    def ppf(a,l,q):
        return -1/l*math.log(1-math.asinh((math.exp(2*a)-1)/(2*math.exp(a))*q)/a)
    @staticmethod
    def mle():
        pass