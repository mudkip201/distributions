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

class extremevalue(Distribution):
    @staticmethod
    def pdf(mu,sigma,x):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return 1/sigma*math.exp(-(x-mu)/sigma)*math.exp(-math.exp(-(x-mu)/sigma))
    @staticmethod
    def cdf(mu,sigma,x):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return math.exp(-(x-mu)/sigma)
    @staticmethod
    def random(mu,sigma):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return mu-sigma*math.log(-math.log(ds.rg0()))
    @staticmethod
    def kurtosis(mu,sigma):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return 2.4
    @staticmethod
    def mean(mu,sigma):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return mu+sigma*ds.euler_gamma
    @staticmethod
    def median(mu,sigma):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return mu-sigma*math.log(math.log(2))
    @staticmethod
    def mode(mu,sigma):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return mu
    @staticmethod
    def variance(mu,sigma):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return sigma**2*math.pi**2/6
    @staticmethod
    def stddev(mu,sigma):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return sigma*math.pi/math.sqrt(6)
    @staticmethod
    def entropy(mu,sigma):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return math.log(sigma)+ds.euler_gamma+1
    @staticmethod
    def skewness(mu,sigma):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return 12*math.sqrt(6)*sp.zeta(3,1)/(math.pi**3)
    @staticmethod
    def ppf(mu,sigma,q):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return mu-sigma*math.log(-math.log(q))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=extremevalue.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'mu':ret[0],'sigma':ret[1]}
