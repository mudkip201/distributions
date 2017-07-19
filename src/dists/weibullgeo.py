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

class weibullgeo(Distribution):
    @staticmethod
    def pdf(a,b,p,x):
        return a*math.pow(b,a)*(1-p)*math.pow(x,a-1)*math.exp(-math.pow(b*x,a))/math.pow(1-p*math.exp(-math.pow(b*x,a)),2)
    @staticmethod
    def cdf(a,b,p,x):
        return (1-math.exp(-math.pow(b*x,a)))/(1-p*math.exp(-math.pow(b*x,a)))
    @staticmethod
    def random(a,b,p):
        u=ds.rg0()
        return math.pow(math.log((1-p*u)/(1-u)),1/a)/b
    @staticmethod
    def mean():
        pass
    @staticmethod
    def median(a,b,p):
        return math.pow(math.log(1-p),1/a)/b
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
    def ppf(a,b,p,q):
        return math.pow(math.log((1-p*q)/(1-q)),1/a)/b
    @staticmethod
    def mle():
        pass