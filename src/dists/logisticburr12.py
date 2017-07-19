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

class logisticburr12(Distribution):
    @staticmethod
    def random(c,k,l,s):
        return s*math.pow(math.exp(math.pow(math.pow(1/ds.rg0()-1,-1/l),1/k))-1,1/c)
    @staticmethod
    def pdf(c,k,l,s,x):
        return l*c*k*math.pow(x,c-1)/(math.pow(s,c)*(1+math.pow(x/s,c)))*math.pow(math.log(math.pow(1+math.pow(x/s,c),k)),-l-1)*1/math.pow(1+math.pow(math.log(math.pow(1+math.pow(x/s,c),k)),-l),2)
    @staticmethod
    def cdf(c,k,l,s,x):
        return 1/(1+math.pow(math.log(math.pow(1+math.pow(x/s,c),k)),-l))
    @staticmethod
    def median(c,k,l,s):
        return s*math.pow(math.exp(math.pow(math.pow(1,-1/l),1/k))-1,1/c)
    @staticmethod
    def ppf(c,k,l,s,q):
        return s*math.pow(math.exp(math.pow(math.pow(1/q-1,-1/l),1/k))-1,1/c)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=logisticburr12.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'c':ret[0],'k':ret[1],'l':ret[2],'s':ret[3]}