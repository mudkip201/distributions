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
import dists.exponential.exponential as exponential
import dists.gamma.gamma as gamma

class xgamma(Distribution):
    @staticmethod
    def random(t):
        u=ds.rg0()
        v=exponential.random(t)
        w=gamma.random(3,t)
        if(u<=t/(1+t)):
            return v
        return w
    @staticmethod
    def pdf(t,x):
        return t**2/(1+t)*(1+t/2*x**2)*math.exp(-t*x)
    @staticmethod
    def cdf(t,x):
        return 1-(1+t+t*x+(t*x)**2/2)/(1+t)*math.exp(-t*x)
    @staticmethod
    def mode(t):
        if(0<t and t<1/2):
            return (1+math.sqrt(1-2*t))/t
        return 0
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=xgamma.pdf(args_[0],i)
            return -tomin
        ret=op.minimize(mlefunc,(1),method='Nelder-Mead').x.tolist()
        return {'t':ret[0]}