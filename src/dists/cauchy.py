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

class cauchy(Distribution):
    @staticmethod
    def pdf(gmma,x0,x):
        if(gmma<=0):
            raise ValueError("gamma must be bigger than 0")
        return 1/(math.pi*gmma*(1+((x-x0)/gmma)**2))
    @staticmethod
    def cdf(gmma,x0,x):
        if(gmma<=0):
            raise ValueError("gamma must be bigger than 0")
        return 1/math.pi*math.atan((x-x0)/gmma)+1/2
    @staticmethod
    def random(gmma,x0):
        if(gmma<=0):
            raise ValueError("gamma must be bigger than 0")
        return x0+gmma*math.tan(math.pi*(r.random()-1/2))
    @staticmethod
    def kurtosis(gmma,x0):
        return None
    @staticmethod
    def mean(gmma,x0):
        return None
    @staticmethod
    def median(gmma,x0):
        if(gmma<=0):
            raise ValueError("gamma must be bigger than 0")
        return x0+gmma*math.tan(0)
    @staticmethod
    def mode(gmma,x0):
        if(gmma<=0):
            raise ValueError("gamma must be bigger than 0")
    @staticmethod
    def variance(gmma,x0):
        return None
    @staticmethod
    def stddev(gmma,x0):
        return None
    @staticmethod
    def entropy(gmma,x0):
        if(gmma<=0):
            raise ValueError("gamma must be bigger than 0")
        return math.log(gmma)+math.log(4*math.pi)
    @staticmethod
    def skewness(gmma,x0):
        return None
    @staticmethod
    def ppf(gmma,x0,q):
        if(gmma<=0):
            raise ValueError("gamma must be bigger than 0")
        return x0+gmma*math.tan(math.pi*(q-1/2))
    @staticmethod
    def mle(x):#not in eclipse
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=cauchy.pdf(args_[0],args_[1],i)
            return -tomin
        '''ret=op.differential_evolution(mlefunc,[(0.01,100),(-25,25)]).x.tolist()
        return {'gamma':ret[0],'x0':ret[1]}
        '''