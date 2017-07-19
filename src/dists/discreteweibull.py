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

class discreteweibull(Distribution):
    @staticmethod
    def random(aa,bb):
        n=r.random()
        while n==1:
            n=r.random()
        return math.pow(-math.log(1-n),1/bb)*aa-1
    @staticmethod
    def pdf(aa,bb,x):
        return math.exp(-math.pow(x/aa,bb))-math.exp(-math.pow((x+1)/aa,bb))
    @staticmethod
    def cdf(aa,bb,x):
        return 1-math.exp(-math.pow((x+1)/aa,bb))
    @staticmethod
    def median(aa,bb):
        return math.pow(-math.log(1/2),1/bb)*aa-1
    @staticmethod
    def ppf(aa,bb,q):
        return math.pow(-math.log(1-q),1/bb)*aa-1
    @staticmethod
    def mle(x): #not working
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=discreteweibull.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(2,2),method='Nelder-Mead').x.tolist()
        #ret=op.differential_evolution(mlefunc,[(0.01,10),(0.01,10)]).x.tolist()
        print(ret[0])
        return {'aa':ret[0],'bb':ret[1]}