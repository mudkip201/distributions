'''
Created on Jul 15, 2017

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

class exponential(Distribution):
    @staticmethod
    def pdf(lmbda,x):
        if(x<0 or lmbda<=0):
            raise ValueError("lmbda must be bgger than 0 and x must be non-negative")
        return lmbda*math.exp(-lmbda*x)
    @staticmethod
    def cdf(lmbda,x):
        if(x<0 or lmbda<=0):
            raise ValueError("lmbda must be bgger than 0 and x must be non-negative")
        return 1-math.exp(-lmbda*x)
    @staticmethod
    def random(lmbda):
        if(lmbda<=0):
            raise ValueError("lmbda must be bgger than 0")
        return -math.log(ds.rg0())/lmbda
    @staticmethod
    def kurtosis(lmbda):
        if(lmbda<=0):
            raise ValueError("lmbda must be bgger than 0")
        return 6
    @staticmethod
    def mean(lmbda):
        if(lmbda<=0):
            raise ValueError("lmbda must be bgger than 0")
        return 1/lmbda
    @staticmethod
    def median(lmbda):
        if(lmbda<=0):
            raise ValueError("lmbda must be bgger than 0")
        return math.log(2)/lmbda
    @staticmethod
    def mode(lmbda):
        if(lmbda<=0):
            raise ValueError("lmbda must be bgger than 0")
        return 0
    @staticmethod
    def variance(lmbda):
        if(lmbda<=0):
            raise ValueError("lmbda must be bgger than 0")
        return 1/(lmbda**2)
    @staticmethod
    def stddev(lmbda):
        if(lmbda<=0):
            raise ValueError("lmbda must be bgger than 0")
        return 1/lmbda
    @staticmethod
    def entropy(lmbda):
        if(lmbda<=0):
            raise ValueError("lmbda must be bgger than 0")
        return 1-math.log(lmbda)
    @staticmethod
    def skewness(lmbda):
        if(lmbda<=0):
            raise ValueError("lmbda must be bgger than 0")
        return 2
    @staticmethod
    def ppf(lmbda,q):
        if(lmbda<=0):
            raise ValueError("lmbda must be bgger than 0")
        return -math.log(q)/lmbda
    @staticmethod
    def mle(x):
        return {'lambda':len(x)/sum(x)}