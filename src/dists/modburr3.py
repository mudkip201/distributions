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

class modburr3(Distribution):
    @staticmethod
    def random(a,b,g):
        return math.pow((math.pow(ds.rg0(),-g/a)-1)/g,-1/b)
    @staticmethod
    def pdf(a,b,g,x):
        return a*b*math.pow(x,-b-1)*math.pow(1+g*math.pow(x,-b),-a/g-1)
    @staticmethod
    def cdf(a,b,g,x):
        return math.pow(1+g*math.pow(x,-b),-a/g)
    @staticmethod
    def median(a,b,g):
        return math.pow((math.pow(1/2,-g/a)-1)/g,-1/b)
    @staticmethod
    def ppf(a,b,g,q):
        return math.pow((math.pow(q,-g/a)-1)/g,-1/b)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=modburr3.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'g':ret[2]}