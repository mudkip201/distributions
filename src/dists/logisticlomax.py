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
import dists.logisticburr12.logisticburr12 as logisticburr12

class logisticlomax(Distribution):
    @staticmethod
    def random(k,l,s):
        return logisticburr12.random(1,k,l,s)
    @staticmethod
    def pdf(k,l,s,x):
        return logisticburr12.pdf(1,k,l,s,x)
    @staticmethod
    def cdf(k,l,s,x):
        return logisticburr12.cdf(1,k,l,s,x)
    @staticmethod
    def median(k,l,s):
        return logisticburr12.median(1,k,l,s)
    @staticmethod
    def ppf(k,l,s,q):
        return logisticburr12.ppf(1,k,l,s,q)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=logisticlomax.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'k':ret[0],'l':ret[1],'s':ret[2]}
