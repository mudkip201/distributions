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

class burr7(Distribution):
    @staticmethod
    def random(r,l,s):
        return s*math.atanh(math.pow(ds.rg0()*math.pow(2,r),1/r)-1)+l
    @staticmethod
    def pdf(r,l,s,x):
            return 1/s*r*(x**2)*math.pow(1+math.tanh((x-l)/s),r-1)/math.pow(2,r)
    @staticmethod
    def cdf(r,l,s,x):
            return math.pow(2,-r)*math.pow(1+math.tanh((x-l)/s),r)
    @staticmethod
    def median(r,l,s):
        return s*math.atanh(math.pow(math.pow(2,r)/2,1/r)-1)+l
    @staticmethod
    def ppf(r,l,s,q):
        return s*math.atanh(math.pow(q*math.pow(2,r),1/r)-1)+l
    @staticmethod
    def mle(x): #not working
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=burr7.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,0,1),method='Nelder-Mead').x.tolist()
        #ret=op.differential_evolution(mlefunc,[(0.01,100),(-25,25),(0.01,100)]).x.tolist()
        return {'r':ret[0],'l':ret[1],'s':ret[2]}