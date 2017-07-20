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

class modweibull(Distribution):
    @staticmethod
    def pdf(a,g,l,x):
        return a*math.pow(x,g-1)*(g+l*x)*math.exp(l*x-a*math.pow(x,g)*math.exp(l*x))
    @staticmethod
    def cdf(a,g,l,x):
        return 1-math.exp(-a*math.pow(x,g)*math.exp(l*x))
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