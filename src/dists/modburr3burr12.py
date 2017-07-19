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

class modburr3burr12(Distribution):
    @staticmethod
    def random(a,b,c,g,k):
        u=ds.rg0()
        return math.pow(math.pow(math.pow((math.pow(u,-g/a)-1)/g,-1/b)+1,1/k)-1,1/c)
    @staticmethod
    def pdf(a,b,c,g,k,x):
        return a*b*c*k*math.pow(x,c-1)*math.pow(1+math.pow(x,c),k-1)*math.pow(math.pow(1+math.pow(x,c),k)-1,-b-1)*math.pow(1+g*math.pow(math.pow(1+math.pow(x,c),k)-1,-b),-a/g-1)
    @staticmethod
    def cdf(a,b,c,g,k,x):
        return math.pow(1+g*math.pow(math.pow(1+math.pow(x,c),k)-1,-b),-a/g)
    @staticmethod
    def median(a,b,c,g,k):
        return math.pow(math.pow(math.pow((math.pow(1/2,-g/a)-1)/g,-1/b)+1,1/k)-1,1/c)
    @staticmethod
    def ppf(a,b,c,g,k,q):
        return math.pow(math.pow(math.pow((math.pow(q,-g/a)-1)/g,-1/b)+1,1/k)-1,1/c)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=modburr3burr12.pdf(args_[0],args_[1],args_[2],args_[3],args_[4],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'c':ret[2],'g':ret[3],'k':ret[4]}