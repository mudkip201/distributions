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

class dagumpoisson(Distribution):
    @staticmethod
    def pdf(b,d,l,t,x):
        return b*l*d*t*math.pow(x,-d-1)*math.pow(1+l*math.pow(x,-d),-b-1)*math.exp(t*math.pow(1+l*math.pow(x,-d),-b))/(math.exp(t)-1)
    @staticmethod
    def cdf(b,d,l,t,x):
        return (math.exp(t*math.pow(1+l*math.pow(x,-d),-b))-1)/(math.exp(t)-1)
    @staticmethod
    def random(b,d,l,t):
        u=ds.rg0()
        return math.pow((math.pow(math.log(u*(math.exp(t)-1)+1)/t,-1/b)-1)/l,-1/d)
    @staticmethod
    def mean():
        pass
    @staticmethod
    def median(b,d,l,t):
        return math.pow((math.pow(math.log((math.exp(t)-1)+1)/(2*t),-1/b)-1)/l,-1/d)
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
    def ppf(b,d,l,t,q):
        return math.pow((math.pow(math.log(q*(math.exp(t)-1)+1)/t,-1/b)-1)/l,-1/d)
    @staticmethod
    def mle():
        pass