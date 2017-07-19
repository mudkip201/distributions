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

class burr8(Distribution):
    @staticmethod
    def random(r,l,s):
        return s*math.log(math.tan(math.pi/2*math.pow(ds.rg0(),1/r)))+l
    @staticmethod
    def pdf(r,l,s,x):
        return 1/s*r*math.exp((x-l)/s)*math.pow(2/math.pi,r)*math.pow(math.atan(math.exp((x-l)/s)),r-1)/(1+math.exp(2*(x-l)/s))
    @staticmethod
    def cdf(r,l,s,x):
        return math.pow(2/math.pi*math.atan(math.exp((x-l)/s)),r)
    @staticmethod
    def median(r,l,s):
        return s*math.log(math.tan(math.pi/2*math.pow(1/2,1/r)))+l
    @staticmethod
    def ppf(r,l,s,q):
        return s*math.log(math.tan(math.pi/2*math.pow(q,1/r)))+l
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=burr8.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,0,1),method='Nelder-Mead').x.tolist()
        #ret=op.differential_evolution(mlefunc,[(0.05,100),(-25,25),(0.01,100)]).x.tolist()
        return {'r':ret[0],'l':ret[1],'s':ret[2]}