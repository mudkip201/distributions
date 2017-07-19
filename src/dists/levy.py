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

class levy(Distribution):
    @staticmethod
    def pdf(mu,c,x):
        if(c<=0):
            raise ValueError("c must be positive")
        if(x<mu):
            raise ValueError("x must be greater than or equal to mu")
        return math.sqrt(c/(2*math.pi))*math.exp((-c)/(2*(x-mu)))/math.pow(x-mu,3/2)
    @staticmethod
    def cdf(mu,c,x):
        if(c<=0):
            raise ValueError("c must be positive")
        if(x<mu):
            raise ValueError("x must be greater than or equal to mu")
        return math.erfc(math.sqrt(c/(2*(x-mu))))
    @staticmethod
    def random(mu,c):
        if(c<=0):
            raise ValueError("c must be positive")
        return c/math.pow(abs(st.norm.ppf(1-ds.rg0()/2)),2)+mu
    @staticmethod
    def kurtosis(mu,c):
        if(c<=0):
            raise ValueError("c must be positive")
        return None
    @staticmethod
    def mean(mu,c):
        if(c<=0):
            raise ValueError("c must be positive")
        return float("infinity")
    @staticmethod
    def variance(mu,c):
        if(c<=0):
            raise ValueError("c must be positive")
        return float("infinity")
    @staticmethod
    def stddev(mu,c):
        if(c<=0):
            raise ValueError("c must be positive")
        return float("infinity")
    @staticmethod
    def entropy(mu,c):
        if(c<=0):
            raise ValueError("c must be positive")
        return (1+3*ds.euler_gamma+math.log(16*math.pi*c*c))/2
    @staticmethod
    def skewness(mu,c):
        if(c<=0):
            raise ValueError("c must be positive")
        return None
    @staticmethod
    def ppf(mu,c,q):
        if(c<=0):
            raise ValueError("c must be positive")
        return c/math.pow(abs(st.norm.ppf(1-q/2)),2)+mu
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=levy.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'mu':ret[0],'c':ret[1]}
