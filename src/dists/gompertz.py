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

class gompertz(Distribution):
    @staticmethod
    def pdf(eta,b,x):
        if(eta<=0 or b<=0):
            raise ValueError("eta and b must be positive")
        if(x<0):
            raise ValueError("x must be non-negative")
        return b*eta*math.exp(b*x)*math.exp(eta)*math.exp(-eta*math.exp(b*x))
    @staticmethod
    def cdf(eta,b,x):
        if(eta<=0 or b<=0):
            raise ValueError("eta and b must be positive")
        if(x<0):
            raise ValueError("x must be non-negative")
        return 1-math.exp(-eta*(math.exp(b*x)-1))
    @staticmethod
    def random(eta,b):
        if(eta<=0 or b<=0):
            raise ValueError("eta and b must be positive")
        return math.log((-math.log(1-r.random())/eta)+1)/b
    @staticmethod
    def median(eta,b):
        if(eta<=0 or b<=0):
            raise ValueError("eta and b must be positive")
        return (1/b)*math.log(math.log(1/2)/eta+1)
    @staticmethod
    def ppf(eta,b,q):
        if(eta<=0 or b<=0):
            raise ValueError("eta and b must be positive")
        return math.log((-math.log(1-q)/eta)+1)/b
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=gompertz.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'eta':ret[0],'b':ret[1]}