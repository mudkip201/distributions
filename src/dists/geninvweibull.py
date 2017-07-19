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

class geninvweibull(Distribution):
    @staticmethod
    def random(a,b,g):
        return a*math.pow(-math.log(ds.rg0())/g,-1/b)
    @staticmethod
    def pdf(a,b,g,x):
        return g*b*math.pow(a,b)*math.pow(x,-b-1)*math.exp(-g*math.pow(a/x,b))
    @staticmethod
    def cdf(a,b,g,x):
        return math.exp(-g*math.pow(a/x,b))
    @staticmethod
    def median(a,b,g):
        return a*math.pow(-math.log(1/2)/g,-1/b)
    @staticmethod
    def ppf(a,b,g,q):
        return a*math.pow(-math.log(q)/g,-1/b)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=geninvweibull.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'g':ret[2]}