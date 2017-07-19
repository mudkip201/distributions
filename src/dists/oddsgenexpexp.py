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

class oddsgenexpexp(Distribution): #odds generalized exponential-exponential
    @staticmethod
    def pdf(l,t,x):
        return l*t*math.exp(x*t)*math.exp(-l*(math.exp(x*t)-1))
    @staticmethod
    def cdf(l,t,x):
        return 1-math.exp(-l*(math.exp(t*x)-1))
    @staticmethod
    def random(l,t):
        return math.log(1-math.log(1-ds.rg0())/l)/t
    @staticmethod
    def median(l,t):
        return math.log(1-math.log(1/2)/l)/t
    @staticmethod
    def ppf(l,t,q):
        return math.log(1-math.log(1-q)/l)/t
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=oddsgenexpexp.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'l':ret[0],'t':ret[1]}