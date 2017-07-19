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

class geninvgenexp(Distribution): #generalized inverse generalized exponential
    @staticmethod
    def random(a,g,l):
        return l*g/(math.log(1-math.pow(1-ds.rg0(),1/a)))
    @staticmethod
    def pdf(a,g,l,x):
        return a*l*g*math.pow(x,-2)*math.exp(-g*l/x)*math.pow(1-math.exp(-g*l/x),a-1)
    @staticmethod
    def cdf(a,g,l,x):
        return 1-math.pow(1-math.exp(-g*l/x),a)
    @staticmethod
    def median(a,g,l):
        return l*g/(math.log(1-math.pow(1/2,1/a)))
    @staticmethod
    def ppf(a,g,l,q):
        return l*g/(math.log(1-math.pow(1-q,1/a)))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=geninvgenexp.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'g':ret[1],'l':ret[2]}
