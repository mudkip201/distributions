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

class gamma(Distribution):
    @staticmethod
    def random(k,t):
        u=r.random()
        x=-2*math.log(1-math.pow(u,1/k))
        v=r.random()
        while(v>(math.pow(x,k/2.0)*math.exp(-x/2.0))/(math.pow(2,k-1)*math.pow((1-math.exp(-x/2.0)),k-1))):
            u=r.random()
            x=-2*math.log(1-math.pow(u,1/k))
            v=r.random()
        return t*v
    @staticmethod
    def pdf(k,t,x):
        return 1/(math.gamma(k)*math.pow(t,k)*math.pow(x,k-1)*math.exp(-x/t))
    @staticmethod
    def kurtosis(k,t):
        return 6/k
    @staticmethod
    def mean(k,t):
        return k*t
    @staticmethod
    def mode(k,t):
        if(k>=1):
            return (k-1)*t
        return None
    @staticmethod
    def variance(k,t):
        return k*t**2
    @staticmethod
    def stddev(k,t):
        return math.sqrt(k)*t
    @staticmethod
    def skewness(k,t):
        return 2/math.sqrt(k)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=gamma.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'k':ret[0],'t':ret[1]}