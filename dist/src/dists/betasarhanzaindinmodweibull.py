'''
Created on Jul 26, 2017

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

class betasarhanzaindinmodweibull(Distribution):
    @staticmethod
    def pdf(a,b,k,l,r,x):
        return (l+b*k*math.pow(x,k-1)*math.exp(-r*l*x-r*b*math.pow(x,k)))*math.pow(1-math.exp(-l*x-b*math.pow(x,k)),a-1)/sp.beta(a,r)
    @staticmethod
    def cdf(a,b,k,l,r,x):
        return sp.betainc(1-math.exp(-l*x-b*math.pow(x,k)),a,r)
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