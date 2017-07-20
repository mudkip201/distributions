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
import dists.exponential.exponential as exponential
import dists.gumbel.gumbel as gumbel

class shiftedgompertz(Distribution):
    @staticmethod
    def random(b,eta):
        return max(exponential.random(b),gumbel.random(b,eta))
    @staticmethod
    def pdf(b,eta,x):
        return b*math.exp(-b*x)*math.exp(-eta*math.exp(-b*x))*(1+eta*(1-math.exp(-b*x)))
    @staticmethod
    def cdf(b,eta,x):
        return (1-math.exp(-b*x))*math.exp(-eta*math.exp(-b*x))
    @staticmethod
    def mode(b,eta):
        if(0<eta and eta<=0.5):
            return 0
        z=(3+eta-math.sqrt(eta**2+2*eta+5))/(2*eta)
        return (-1/b)*math.log(z)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=shiftedgompertz.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'b':ret[0],'eta':ret[1]}