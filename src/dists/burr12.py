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

class burr12(Distribution):
    @staticmethod
    def pdf(c,k,x):
        if(c<=0 or k<=0 or x<=0):
            raise ValueError("c, k, and x must be bigger than 0")
        return c*k*math.pow(x,c-1)/math.pow(1+math.pow(x,c),k+1)
    @staticmethod
    def cdf(k,c,x):
        if(c<=0 or k<=0 or x<=0):
            raise ValueError("c, k, and x must be bigger than 0")
        return 1-math.pow(1+math.pow(x,c),-k)
    @staticmethod
    def random(c,k):
        if(c<=0 or k<=0):
            raise ValueError("c and k must be bigger than 0")
        return math.pow(math.pow(1-ds.rg0(),-1/k)-1,1/c)
    @staticmethod
    def median(c,k):
        if(c<=0 or k<=0):
            raise ValueError("c and k must be bigger than 0")
        math.pow(math.pow(2,1/k)-1,1/c)
    @staticmethod
    def ppf(c,k,q):
        if(c<=0 or k<=0):
            raise ValueError("c and k must be bigger than 0")
        math.pow(math.pow(1-q,-1/k)-1,1/c)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=burr12.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        #ret=op.differential_evolution(mlefunc,[(0.01,10),(0.01,10)]).x.tolist()
        return {'c':ret[0],'k':ret[1]}