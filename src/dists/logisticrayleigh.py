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
import dists.logisticweibull.logisticweibull as logisticweibull


class logisticrayleigh(Distribution):
    @staticmethod
    def random(b,l):
        return logisticweibull.random(2,b,l)
    @staticmethod
    def pdf(b,l,x):
        return logisticweibull.pdf(2,b,l)
    @staticmethod
    def cdf(b,l,x):
        return logisticweibull.cdf(2,b,l)
    @staticmethod
    def median(b,l):
        return logisticweibull.median(2,b,l)
    @staticmethod
    def ppf(b,l,q):
        return logisticweibull.ppf(2,b,l,q)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=logisticrayleigh.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'b':ret[0],'l':ret[1]}