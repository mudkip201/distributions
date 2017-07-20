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

class weibullfrechet(Distribution):
    @staticmethod
    def random(a,b,aa,bb):
        u=ds.rg0()
        return aa*math.pow(math.log(1+math.pow(-math.log(1-u)/a,-1/b)),-1/bb)
    @staticmethod
    def pdf(a,b,aa,bb,x):
        return a*b*bb*math.pow(aa,bb)*math.pow(x,-bb-1)*math.exp(-b*math.pow(aa/x,bb))*math.pow(1-math.exp(-math.pow(aa/x,bb))-b-1)*math.exp(-a*math.pow(math.exp(math.pow(aa/x,bb)-1),-b))
    @staticmethod
    def cdf(a,b,aa,bb,x):
        return 1-math.exp(-a*math.pow(math.exp(math.pow(aa/x,bb)-1),-b))
    @staticmethod
    def median(a,b,aa,bb):
        return aa*math.pow(math.log(1+math.pow(-math.log(1/2)/a,-1/b)),-1/bb)
    @staticmethod
    def ppf(a,b,aa,bb,q):
        return aa*math.pow(math.log(1+math.pow(-math.log(1-q)/a,-1/b)),-1/bb)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=weibullfrechet.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'aa':ret[2],'bb':ret[3]}