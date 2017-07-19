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

class toppleone(Distribution):
    @staticmethod
    def random(bb):
        return 1-math.sqrt(1-math.pow(r.random(),1/bb))
    @staticmethod
    def pdf(bb,x):
        return bb*(2-2*x)*math.pow(2*x-x**2,bb-1)
    @staticmethod
    def cdf(bb,x):
        return math.pow(2*x-x**2,bb)
    @staticmethod
    def median(bb):
        return 1-math.sqrt(1-math.pow(1/2,1/bb))
    @staticmethod
    def ppf(bb,q):
        return 1-math.sqrt(1-math.pow(q,1/bb))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=toppleone.pdf(args_[0],i)
            return -tomin
        ret=op.minimize(mlefunc,(1),method='Nelder-Mead').x.tolist()
        return {'bb':ret[0]}