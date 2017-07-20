'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.special as sp
import scipy.optimize as op

class gumbel(Distribution):
    @staticmethod
    def pdf(mu,bb,x):
        if(bb<=0):
            raise ValueError("bb must be positive")
        z=(x-mu)/bb
        return 1/bb*math.exp(-(z+math.exp(-z)))
    @staticmethod
    def cdf(mu,bb,x):
        if(bb<=0):
            raise ValueError("bb must be positive")
        return math.exp(-math.exp(-(x-mu)/bb))
    @staticmethod
    def random(mu,bb):
        if(bb<=0):
            raise ValueError("bb must be positive")
        return mu-bb*math.log(-math.log(ds.rg0()))
    @staticmethod
    def kurtosis(mu,bb):
        if(bb<=0):
            raise ValueError("bb must be positive")
        return 12/5
    @staticmethod
    def mean(mu,bb):
        if(bb<=0):
            raise ValueError("bb must be positive")
        return mu+bb*ds.euler_gamma
    @staticmethod
    def median(mu,bb):
        if(bb<=0):
            raise ValueError("bb must be positive")
        return mu-bb*math.log(math.log(2))
    @staticmethod
    def mode(mu,bb):
        if(bb<=0):
            raise ValueError("bb must be positive")
        return mu
    @staticmethod
    def variance(mu,bb):
        if(bb<=0):
            raise ValueError("bb must be positive")
        return math.pi**2/6*(bb**2)
    @staticmethod
    def stddev(mu,bb):
        if(bb<=0):
            raise ValueError("bb must be positive")
        return math.pi*bb/math.sqrt(6)
    @staticmethod
    def entropy(mu,bb):
        if(bb<=0):
            raise ValueError("bb must be positive")
        return math.log(bb)+ds.euler_gamma+1
    @staticmethod
    def skewness(mu,bb):
        if(bb<=0):
            raise ValueError("bb must be positive")
        return 12*math.sqrt(6)*sp.zeta(3,1)/(math.pi**3)
    @staticmethod
    def ppf(mu,bb,q):
        if(bb<=0):
            raise ValueError("bb must be positive")
        return mu-bb*math.log(-math.log(q))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=gumbel.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'mu':ret[0],'bb':ret[1]}