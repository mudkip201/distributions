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
import dists.normal.normal as normal

class foldednormal(Distribution):
    @staticmethod
    def random(m,s):
        return abs(normal.random(m,s))
    @staticmethod
    def pdf(m,s,x):
        c=abs(m)/s
        return math.sqrt(2/math.pi)*math.cosh(c*x)*math.exp(-(x**2+c**2)/2)
    @staticmethod
    def cdf(m,s,x):
        c=abs(m)/s
        return st.norm.cdf(x-c)-st.norm.cdf(-x-c)
    @staticmethod
    def mean(m,s):
        c=abs(m)/s
        return math.sqrt(2/math.pi)*math.exp(-c**2/2)+c*math.erf(c/math.sqrt(2))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=foldednormal.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'m':ret[0],'s':ret[1]}
