'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.optimize as op
import dists.normal.normal as normal

class halfnormal(Distribution):
    @staticmethod
    def pdf(sigma2,x):
        if(sigma2<=0):
            raise ValueError("sigma2 must be positive")
        if(x<0):
            raise ValueError("x must be non-negative")
        return math.sqrt(2/(math.pi*sigma2))*math.exp(-x**2/(2*sigma2))
    @staticmethod
    def cdf(sigma2,x):
        if(sigma2<=0):
            raise ValueError("sigma2 must be positive")
        if(x<0):
            raise ValueError("x must be non-negative")
        return math.erf(math.sqrt(x**2/(2*sigma2)))
    @staticmethod
    def random(sigma2):
        if(sigma2<=0):
            raise ValueError("sigma2 must be positive")
        n=normal.random(0,sigma2)
        while(n<0):
            n=normal.random(0,sigma2)
        return n
    @staticmethod
    def mean(sigma2):
        if(sigma2<=0):
            raise ValueError("sigma2 must be positive")
        return math.sqrt(2*sigma2/math.pi)
    @staticmethod
    def variance(sigma2):
        if(sigma2<=0):
            raise ValueError("sigma2 must be positive")
        return sigma2*(1-2/math.pi)
    @staticmethod
    def stddev(sigma2):
        if(sigma2<=0):
            raise ValueError("sigma2 must be positive")
        return math.sqrt(sigma2*(1-2/math.pi))
    @staticmethod
    def entropy(sigma2):
        if(sigma2<=0):
            raise ValueError("sigma2 must be positive")
        return 1/2*math.log(math.pi*sigma2/2)+1/2
    @staticmethod
    def skewness(sigma2):
        if(sigma2<=0):
            raise ValueError("sigma2 must be positive")
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=halfnormal.pdf(args_[0],i)
            return -tomin
        ret=op.minimize(mlefunc,(1),method='Nelder-Mead').x.tolist()
        return {'sigma2':ret[0]}
