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
import dists.normal.normal as normal

class betaweibullnormal(Distribution):
    @staticmethod
    def pdf(a,b,l,m,s2,x):
        G=normal.cdf(m,s2,x)
        g=normal.pdf(m,s2,x)
        return l*g*math.pow(G,l-1)/(sp.beta(a,b)*math.pow(1-G,l+1))*math.exp(-b*math.pow(G/(1-G),l))*math.pow(1-math.exp(-math.pow(G/(1-G),l)),a-1)
    @staticmethod
    def cdf(a,b,l,m,s2,x):
        G=normal.cdf(m,s2,x)
        return sp.betainc(a,b,1-math.exp(-math.pow(G/(1-G),l)))
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