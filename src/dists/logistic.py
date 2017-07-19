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

class logistic(Distribution):
    @staticmethod
    def pdf(mu,s,x):
        if(s<=0):
            raise ValueError("s must be positive")
        z=-(x-mu)/s
        return math.exp(z)/(s*(1+math.exp(z))**2)
    @staticmethod
    def cdf(mu,s,x):
        if(s<=0):
            raise ValueError("s must be positive")
        z=-(x-mu)/s
        return 1/(1+math.exp(z))
    @staticmethod
    def random(mu,s):
        if(s<=0):
            raise ValueError("s must be positive")
        n=ds.rg0()
        while(n==1):
            n=ds.rg0()
        return mu+s*math.log(n/(1-n))
    @staticmethod
    def kurtosis(mu,s):
        if(s<=0):
            raise ValueError("s must be positive")
        return 1.2
    @staticmethod
    def mean(mu,s):
        if(s<=0):
            raise ValueError("s must be positive")
        return mu
    @staticmethod
    def median(mu,s):
        if(s<=0):
            raise ValueError("s must be positive")
        return mu
    @staticmethod
    def mode(mu,s):
        if(s<=0):
            raise ValueError("s must be positive")
        return mu
    @staticmethod
    def variance(mu,s):
        if(s<=0):
            raise ValueError("s must be positive")
        return (s*math.pi)**2/3
    @staticmethod
    def stddev(mu,s):
        if(s<=0):
            raise ValueError("s must be positive")
        return s*math.pi/math.sqrt(3)
    @staticmethod
    def entropy(mu,s):
        if(s<=0):
            raise ValueError("s must be positive")
        return math.log(s)+2
    @staticmethod
    def skewness(mu,s):
        if(s<=0):
            raise ValueError("s must be positive")
        return 0
    @staticmethod
    def ppf(mu,s,q):
        if(s<=0):
            raise ValueError("s must be positive")
        return mu+s*math.log(q/(1-q))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=logistic.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'mu':ret[0],'s':ret[1]}