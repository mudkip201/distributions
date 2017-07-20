'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r
import scipy.optimize as op

class logcauchy(Distribution):
    @staticmethod
    def random(mu,sigma):
        return math.exp(mu+sigma*math.tan(math.pi*(r.random()-1/2)))
    @staticmethod
    def pdf(mu,sigma,x):
        return 1/(x*math.pi)*sigma/((math.log(x)-mu)**2+sigma**2)
    @staticmethod
    def cdf(mu,sigma,x):
        return 1/math.pi*math.atan((math.log(x)-mu)/sigma)+1/2
    @staticmethod
    def median(mu,sigma):
        return math.exp(mu)
    @staticmethod
    def variance(mu,sigma):
        return float("infinity")
    @staticmethod
    def stddev(mu,sigma):
        return float("infinity")
    @staticmethod
    def ppf(mu,sigma,q):
        return math.exp(mu+sigma*math.tan(math.pi*(q-1/2)))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=logcauchy.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'mu':ret[0],'sigma':ret[1]}