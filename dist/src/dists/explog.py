'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r
import scipy.optimize as op

class explog(Distribution): #exponentiated logistic
    @staticmethod
    def pdf(p,bb,x):
        if(p<=0 or p>=1):
            raise ValueError("p must be between 0 and 1 exclusive")
        if(bb<=0):
            raise ValueError("bb must be positive")
        if(x<0):
            raise ValueError("x must be non-negative")
        return 1/(-math.log(p))*(bb*(1-p)*math.exp(-bb*x))/(1-(1-p)*math.exp(-bb*x))
    @staticmethod
    def cdf(p,bb,x):
        if(p<=0 or p>=1):
            raise ValueError("p must be between 0 and 1 exclusive")
        if(bb<=0):
            raise ValueError("bb must be positive")
        if(x<0):
            raise ValueError("x must be non-negative")
        return 1-math.log(1-(1-p)*math.exp(-bb*x))/math.log(p)
    @staticmethod
    def random(p,bb):
        if(p<=0 or p>=1):
            raise ValueError("p must be between 0 and 1 exclusive")
        if(bb<=0):
            raise ValueError("bb must be positive")
        return 1/bb * math.log((1-p)/(1-math.pow(p,r.random())))
    @staticmethod
    def median(p,bb):
        if(p<=0 or p>=1):
            raise ValueError("p must be between 0 and 1 exclusive")
        if(bb<=0):
            raise ValueError("bb must be positive")
        return math.log(1+math.sqrt(p))/bb
    @staticmethod
    def mode(p,bb):
        if(p<=0 or p>=1):
            raise ValueError("p must be between 0 and 1 exclusive")
        if(bb<=0):
            raise ValueError("bb must be positive")
        return 0
    @staticmethod
    def ppf(p,bb,q):
        if(p<=0 or p>=1):
            raise ValueError("p must be between 0 and 1 exclusive")
        if(bb<=0):
            raise ValueError("bb must be positive")
        return 1/bb * math.log((1-p)/(1-math.pow(p,q)))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=explog.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(0.5,1),method='Nelder-Mead').x.tolist()
        return {'p':ret[0],'bb':ret[1]}