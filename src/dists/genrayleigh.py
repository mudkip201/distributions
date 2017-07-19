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

class genrayleigh(Distribution):
    @staticmethod
    def random(a,l):
        return math.sqrt(-math.log(1-math.pow(ds.rg0(),1/a)))/l
    @staticmethod
    def pdf(a,l,x):
        return 2*a*l**2*x*math.exp(-(l*x)**2)*math.pow(1-math.exp(-(l*x)**2),a-1)
    @staticmethod
    def cdf(a,l,x):
        return math.pow(1-math.exp(-(l*x)**2),a)
    @staticmethod
    def median(a,l):
        return math.sqrt(-math.log(1-math.pow(1/2,1/a)))/l
    @staticmethod
    def ppf(a,l,q):
        return math.sqrt(-math.log(1-math.pow(q,1/a)))/l
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=genrayleigh.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'aa':ret[0],'mu':ret[1]}