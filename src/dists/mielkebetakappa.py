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

class mielkebetakappa(Distribution):
    @staticmethod
    def random(k,theta):
        n=ds.rg0()
        while n==1:
            n=ds.rg0()
        return math.pow((math.pow(n,theta/k))/(1-math.pow(n,theta/k)),1/theta)
    @staticmethod
    def pdf(k,theta,x):
        return k*math.pow(x,k-1)/math.pow(1+math.pow(x,theta),1+k/theta)
    @staticmethod
    def cdf(k,theta,x):
        return math.pow(x,k)/math.pow(1+math.pow(x,theta),k/theta)
    @staticmethod
    def median(k,theta):
        return math.pow((math.pow(1/2,theta/k))/(1-math.pow(1/2,theta/k)),1/theta)
    @staticmethod
    def ppf(k,theta,q):
        return math.pow((math.pow(q,theta/k))/(1-math.pow(q,theta/k)),1/theta)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=mielkebetakappa.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'k':ret[0],'theta':ret[1]}