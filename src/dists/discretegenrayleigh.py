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

class discretegenrayleigh(Distribution):
    @staticmethod
    def random(a,l):
        p=math.exp(-(l**2))
        return 1-math.sqrt(math.log(1-math.pow(ds.rg0(),1/a))/math.log(p))
    @staticmethod
    def pdf(a,l,x):
        p=math.exp(-(l**2))
        return math.pow(1-math.pow(p,(x+1)**2),a)-math.pow(1-math.pow(p,x**2),a)
    @staticmethod
    def cdf(a,l,x):
        p=math.exp(-(l**2))
        return math.pow(1-math.pow(p,(x+1)**2),a)
    @staticmethod
    def median(a,l):
        p=math.exp(-(l**2))
        return 1-math.sqrt(math.log(1-math.pow(1/2,1/a))/math.log(p))
    @staticmethod
    def ppf(a,l,q):
        p=math.exp(-(l**2))
        return 1-math.sqrt(math.log(1-math.pow(q,1/a))/math.log(p))
    @staticmethod
    def mle(x): #not working
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=discretegenrayleigh.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        #ret=op.differential_evolution(mlefunc,[(0.01,10),(-5,5)]).x.tolist()
        return {'a':ret[0],'l':ret[1]}