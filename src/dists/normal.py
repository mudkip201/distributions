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

class normal(Distribution):
    @staticmethod
    def pdf(mu,sigma2,x):
        if(sigma2<=0):
            raise ValueError("sigma2 must be positive")
        return 1/math.sqrt(2*math.pi*sigma2)*math.exp(-(x-mu)**2/(2*sigma2))
    @staticmethod
    def cdf(mu,sigma2,x):
        if(sigma2<=0):
            raise ValueError("sigma2 must be positive")
        return 1/2*(1+math.erf((x-mu)/math.sqrt(2*sigma2)))
    @staticmethod
    def random(mu,sigma2):
        if(sigma2<=0):
            raise ValueError("sigma2 must be positive")
        a=ds.rg0()
        b=ds.rg0()
        while(a==1):
            a=ds.rg0()
        while(b==1):
            b=ds.rg0()
        return math.sqrt(-2*math.log(a))*math.cos(2*math.pi*b)*math.sqrt(sigma2)+mu
    @staticmethod
    def kurtosis(mu,sigma2):
        if(sigma2<=0):
            raise ValueError("sigma2 must be positive")
        return 0
    @staticmethod
    def mean(mu,sigma2):
        if(sigma2<=0):
            raise ValueError("sigma2 must be positive")
        return mu
    @staticmethod
    def median(mu,sigma2):
        if(sigma2<=0):
            raise ValueError("sigma2 must be positive")
        return mu
    @staticmethod
    def mode(mu,sigma2):
        if(sigma2<=0):
            raise ValueError("sigma2 must be positive")
        return mu
    @staticmethod
    def variance(mu,sigma2):
        if(sigma2<=0):
            raise ValueError("sigma2 must be positive")
        return sigma2
    @staticmethod
    def stddev(mu,sigma2):
        if(sigma2<=0):
            raise ValueError("sigma2 must be positive")
        return math.sqrt(sigma2)
    @staticmethod
    def entropy(mu,sigma2):
        if(sigma2<=0):
            raise ValueError("sigma2 must be positive")
        return 1/2*math.log(2*math.pi*math.exp(1)*sigma2)
    @staticmethod
    def skewness(mu,sigma2):
        if(sigma2<=0):
            raise ValueError("sigma2 must be positive")
        return 0
    @staticmethod
    def mle(x):
        return {'mu':sum(x)/len(x),'sigma2':np.var(x)*(len(x)-1)/len(x)}
