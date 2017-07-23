'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.optimize as op

class raisedcosine(Distribution):
    @staticmethod
    def pdf(mu,s,x):
        if(s<=0):
            raise ValueError("s must be positive")
        if(abs(x-mu)>s):
            raise ValueError("x must be between mu-s and mu+s")
        return 1/(2*s)*(1+math.cos((x-mu)/s*math.pi))
    @staticmethod
    def cdf(mu,s,x):
        if(s<=0):
            raise ValueError("s must be positive")
        if(abs(x-mu)>s):
            raise ValueError("x must be between mu-s and mu+s")
        return 1/2*(1+(x-mu)/s+1/math.pi*math.sin((x-mu)/s*math.pi))
    @staticmethod
    def kurtosis(mu,s):
        if(s<=0):
            raise ValueError("s must be positive")
        return 6/5*(90-math.pi**4)/((math.pi**2-6)**2)
    @staticmethod
    def mean(mu,s):
        if(s<=0):
            raise ValueError("s must be positive")
        return mu
    @staticmethod
    def median(mu,s):
        if(s<=0):
            raise ValueError("s must be positive")
        return mu
    @staticmethod
    def mode(mu,s):
        if(s<=0):
            raise ValueError("s must be positive")
        return mu
    @staticmethod
    def variance(mu,s):
        if(s<=0):
            raise ValueError("s must be positive")
        return s**2*(1/3-2/math.pi**2)
    @staticmethod
    def stddev(mu,s):
        if(s<=0):
            raise ValueError("s must be positive")
        return s*math.sqrt(1/3-2/math.pi**2)
    @staticmethod
    def skewness(mu,s):
        if(s<=0):
            raise ValueError("s must be positive")
        return 0
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=raisedcosine.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'mu':ret[0],'s':ret[1]}
