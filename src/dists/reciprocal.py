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

class reciprocal(Distribution):
    @staticmethod
    def pdf(a,b,x):
        if(a<=0 or a>=b):
            raise ValueError("a must be positive and less than b")
        if(x<a or x>b):
            raise ValueError("x must be between a and b inclusive")
        return 1/(x*(math.log(b)-math.log(a)))
    @staticmethod
    def cdf(a,b,x):
        if(a<=0 or a>=b):
            raise ValueError("a must be positive and less than b")
        if(x<a or x>b):
            raise ValueError("x must be between a and b inclusive")
        return (math.log(x)-math.log(a))/(math.log(b)-math.log(a))
    @staticmethod
    def random(a,b):
        if(a<=0 or a>=b):
            raise ValueError("a must be positive and less than b")
        return math.exp(r.random()*(math.log(b)-math.log(a))+math.log(a))
    @staticmethod
    def mean(a,b):
        if(a<=0 or a>=b):
            raise ValueError("a must be positive and less than b")
        return (b-a)/(math.log(b)-math.log(a))
    @staticmethod
    def ppf(a,b,q):
        if(a<=0 or a>=b):
            raise ValueError("a must be positive and less than b")
        return math.exp(q*(math.log(b)-math.log(a))+math.log(a))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=reciprocal.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1]}