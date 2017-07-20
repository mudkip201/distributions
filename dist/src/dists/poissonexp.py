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

class poissonexp(Distribution): #poisson exponential
    @staticmethod
    def random(l,t):
        u=ds.rg0()
        return -math.log(-math.log(1-(1-math.exp(-t))*u)/t)/l
    @staticmethod
    def pdf(l,t,x):
        return t*l*math.exp(-l*x-t*math.exp(-l*x))/(1-math.exp(-t))
    @staticmethod
    def cdf(l,t,x):
        return 1-(1-math.exp(-t*math.exp(-l*x)))/(1-math.exp(-t))
    @staticmethod
    def median(l,t):
        return -math.log(-math.log(1-(1-math.exp(-t))/2)/t)/l
    @staticmethod
    def ppf(l,t,q):
        return -math.log(-math.log(1-(1-math.exp(-t))*q)/t)/l
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=poissonexp.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'l':ret[0],'t':ret[1]}