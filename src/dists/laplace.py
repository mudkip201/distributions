'''
Created on Jul 16, 2017

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

class laplace(Distribution):
    @staticmethod
    def pdf(mu,b,x):
        if(b<=0):
            raise ValueError("b must be positive")
        return 1/(2*b)*math.exp(-abs(x-mu)/b)
    @staticmethod
    def cdf(mu,b,x):
        if(b<=0):
            raise ValueError("b must be positive")
        if(x<mu):
            return 1/2*math.exp((x-mu)/b)
        return 1-1/2*math.exp(-(x-mu)/b)
    @staticmethod
    def random(mu,b):
        if(b<=0):
            raise ValueError("b must be positive")
        u=r.random()
        return mu-b*(abs(u-0.5)/(u-0.5))*math.log(1-2*abs(u-0.5))
    @staticmethod
    def kurtosis(mu,b):
        if(b<=0):
            raise ValueError("b must be positive")
        return 3
    @staticmethod
    def mean(mu,b):
        if(b<=0):
            raise ValueError("b must be positive")
        return mu
    @staticmethod
    def median(mu,b):
        if(b<=0):
            raise ValueError("b must be positive")
        return mu
    @staticmethod
    def mode(mu,b):
        if(b<=0):
            raise ValueError("b must be positive")
        return mu
    @staticmethod
    def variance(mu,b):
        if(b<=0):
            raise ValueError("b must be positive")
        return 2*(b**2)
    @staticmethod
    def stddev(mu,b):
        if(b<=0):
            raise ValueError("b must be positive")
        return b*math.sqrt(2)
    @staticmethod
    def entropy(mu,b):
        if(b<=0):
            raise ValueError("b must be positive")
        return math.log(2*b*math.exp(1))
    @staticmethod
    def skewness(mu,b):
        if(b<=0):
            raise ValueError("b must be positive")
        return 0
    @staticmethod
    def ppf(mu,b,q):
        if(b<=0):
            raise ValueError("b must be positive")
        return mu-b*(abs(q-0.5)/(q-0.5))*math.log(1-2*abs(q-0.5))
    @staticmethod
    def mle(x):
        mu=np.median(x)
        b=np.average(np.abs(np.array(x)-mu))/len(x)
        return {'mu':mu,'b':b}