'''
Created on Jul 15, 2017

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
import dists.extremevalue.extremevalue as extremevalue
import dists.explog.explog as explog

class genextremevalue(Distribution):
    @staticmethod
    def pdf(mu,sigma,xi,x):
        if(xi==0):
            return extremevalue.pdf(mu,sigma,x)
        if(xi>0 and x<mu-sigma/xi):
            raise ValueError("x must be bigger than or equal to mu-sigma/xi")
        if(xi<0 and x>mu-sigma/xi):
            raise ValueError("x must be less than or equal to mu-sigma/xi")
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        t=math.pow(1+xi*(x-mu)/sigma,-1/xi)
        return 1/sigma*math.pow(t,xi+1)*math.exp(t)
    @staticmethod
    def cdf(mu,sigma,xi,x):
        if(xi==0):
            return extremevalue.cdf(mu,sigma,x)
        if(xi>0 and x<mu-sigma/xi):
            raise ValueError("x must be bigger than or equal to mu-sigma/xi")
        if(xi<0 and x>mu-sigma/xi):
            raise ValueError("x must be less than or equal to mu-sigma/xi")
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        t=math.pow(1+xi*(x-mu)/sigma,-1/xi)
        return math.exp(-t)
    @staticmethod
    def random(mu,sigma,xi):
        if(xi==0):
            return extremevalue.random(mu,sigma)
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return mu+sigma*(math.pow(math.log(1/ds.rg0()),-xi)-1)/xi
    @staticmethod
    def kurtosis(mu,sigma,xi):
        if(xi==0):
            return extremevalue.random(mu,sigma)
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        if(xi>=1/4):
            return float("infinity")
        return (math.gamma(1-4*xi)-4*math.gamma(1-xi)*math.gamma(1-3*xi)+6*math.gamma(1-2*xi)*(math.gamma(1-xi)**2)-3*(math.gamma(1-xi)**4))/math.pow(math.gamma(1-2*xi)-math.gamma(1-xi)**2,2)-3
    @staticmethod
    def mean(mu,sigma,xi):
        if(xi==0):
            return extremevalue.random(mu,sigma)
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        if(xi>=1):
            return float("infinity")
        return mu+sigma*(math.gamma(1-xi)-1)/xi
    @staticmethod
    def median(mu,sigma,xi):
        if(xi==0):
            return extremevalue.random(mu,sigma)
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return mu+sigma*(math.pow(math.log(2),-xi)-1)/xi
    @staticmethod
    def mode(mu,sigma,xi):
        if(xi==0):
            return extremevalue.random(mu,sigma)
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return mu+sigma*(math.pow(1+xi,-xi)-1)/xi
    @staticmethod
    def variance(mu,sigma,xi):
        if(xi==0):
            return extremevalue.random(mu,sigma)
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        if(xi>=1/2):
            return float("infinity")
        return sigma**2*(math.gamma(1-2*xi)-math.gamma(1-xi)**2)/xi**2
    @staticmethod
    def stddev(mu,sigma,xi):
        if(xi==0):
            return extremevalue.random(mu,sigma)
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        if(xi>=1/2):
            return float("infinity")
        return sigma/xi*math.sqrt(math.gamma(1-2*xi)-math.gamma(1-xi)**2)
    @staticmethod
    def entropy(mu,sigma,xi):
        if(xi==0):
            return extremevalue.random(mu,sigma)
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return math.log(sigma)+ds.euler_gamma*(xi+1)+1
    @staticmethod
    def skewness(mu,sigma,xi):
        if(xi==0):
            return extremevalue.random(mu,sigma)
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        if(xi>=1/3):
            return float("infinity")
        return (abs(xi)/xi)*(math.gamma(1-3*xi)-3*math.gamma(1-xi)*math.gamma(1-2*xi)+2*math.gamma(1-xi)**3)/math.pow(math.gamma(1-2*xi)-math.gamma(1-xi)**2,3/2)
    @staticmethod
    def ppf(mu,sigma,xi,q):
        if(xi==0):
            return extremevalue.random(mu,sigma)
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return mu+sigma*(math.pow(math.log(1/q),-xi)-1)/xi
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=explog.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'mu':ret[0],'sigma':ret[1],'xi':ret[2]}
