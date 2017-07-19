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

class gengompertz(Distribution): #generalized gompertz
    @staticmethod
    def random(c,l,t):
        return math.log(1-c/l*math.log(1-math.pow(ds.rg0(),t)))/c
    @staticmethod
    def pdf(c,l,t,x):
        return t*l*math.exp(c*x)*math.exp(-l/c*(math.exp(c*x)-1))*math.pow(1-math.exp(-l/c*(math.exp(c*x)-1)),t-1)
    @staticmethod
    def cdf(c,l,t,x):
        return math.pow(1-math.exp(-l/c*(math.exp(c*x)-1)),t)
    @staticmethod
    def median(c,l,t):
        return math.log(1-c/l*math.log(1-math.pow(1/2,t)))/c
    @staticmethod
    def ppf(c,l,t,q):
        return math.log(1-c/l*math.log(1-math.pow(q,t)))/c
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=gengompertz.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'c':ret[0],'l':ret[1],'t':ret[2]}
