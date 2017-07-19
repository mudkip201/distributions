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

class exponentiatedweibullexp(Distribution): #exponentiated weibull exponential
    @staticmethod
    def random(a,c,g):
        return -math.log(1-math.pow(1-math.exp(-g*math.pow(-math.log(1-ds.rg0()),1/a)),1/c))
    @staticmethod
    def pdf(a,c,g,x):
        return c*a/g*math.exp(-x)*math.pow(1-math.exp(-x),c-1)/(1-math.pow(1-math.exp(-x),c))*math.pow(-math.log(1-math.pow(1-math.exp(-x),c))/g,a-1)*math.exp(-math.pow(-math.log(1-math.pow(1-math.exp(-x),c))/g),a)
    @staticmethod
    def cdf(a,c,g,x):
        return 1-math.exp(-math.pow(-math.log(1-math.pow(1-math.exp(-x),c))/g,a))
    @staticmethod
    def median(a,c,g):
        return -math.log(1-math.pow(1-math.exp(-g*math.pow(-math.log(1/2),1/a)),1/c))
    @staticmethod
    def ppf(a,c,g,q):
        return -math.log(1-math.pow(1-math.exp(-g*math.pow(-math.log(1-q),1/a)),1/c))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=exponentiatedweibullexp.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'c':ret[1],'g':ret[2]}