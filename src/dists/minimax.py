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

class minimax(Distribution):
    @staticmethod
    def random(aa,bb):
        return math.pow(1-math.pow(1-r.random(),1/aa),1/bb)
    @staticmethod
    def pdf(aa,bb,x):
        return aa*bb*math.pow(x,aa-1)*math.pow(1-math.pow(x,aa),bb-1)
    @staticmethod
    def cdf(aa,bb,x):
        return 1-math.pow(1-math.pow(x,aa),bb)
    @staticmethod
    def median(aa,bb):
        return math.pow(1-math.pow(1/2,1/aa),1/bb)
    @staticmethod
    def ppf(aa,bb,q):
        return math.pow(1-math.pow(1-q,1/aa),1/bb)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=minimax.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'aa':ret[0],'bb':ret[1]}