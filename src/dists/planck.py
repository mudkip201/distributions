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

class planck(Distribution):
    @staticmethod
    def random(l):
        return math.ceil(-math.log(1-ds.rg0())/l-1)
    @staticmethod
    def pdf(l,x):
        if(x*l>=0):
            return (1-math.exp(-l))*math.exp(-l*x)
    @staticmethod
    def cdf(l,x):
        if(x*l>=0):
            return 1-math.exp(-l*(math.floor(x)+1))
    @staticmethod
    def median(l):
        return math.ceil(-math.log(1/2)/l-1)
    @staticmethod
    def ppf(l,q):
        return math.ceil(-math.log(1-q)/l-1)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=planck.pdf(args_[0],i)
            return -tomin
        ret=op.minimize(mlefunc,(1),method='Nelder-Mead').x.tolist()
        return {'l':ret[0]}