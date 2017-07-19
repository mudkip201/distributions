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

class genpareto(Distribution):
    @staticmethod
    def pdf(mu,sigma,xi,x):
        if(sigma<=0):
            raise ValueError("sigma must be greater than 0")
        if(x<mu and xi>=0):
            raise ValueError("x must be bigger than or equal to mu")
        if(xi<0 and (x<mu or x>(mu-sigma/xi))):
            raise ValueError("x must be between mu and mu-sigma/xi")
        return 1/sigma*math.pow(1+xi*(x-mu)/sigma,-1/xi+1)
    @staticmethod
    def cdf(mu,sigma,xi,x):
        if(sigma<=0):
            raise ValueError("sigma must be greater than 0")
        if(x<mu and xi>=0):
            raise ValueError("x must be bigger than or equal to mu")
        if(xi<0 and (x<mu or x>(mu-sigma/xi))):
            raise ValueError("x must be between mu and mu-sigma/xi")
        return 1-math.pow(1+xi*(x-mu)/sigma,-1/xi)
    @staticmethod
    def random(mu,sigma,xi):
        if(sigma<=0):
            raise ValueError("sigma must be greater than 0")
        if(xi==0):
            return mu-sigma*math.log(ds.rg0())
        return mu+(sigma*(math.pow(ds.rg0(),-xi)-1)/xi)
    @staticmethod
    def kurtosis(mu,sigma,xi):
        if(sigma<=0):
            raise ValueError("sigma must be greater than 0")
        if(xi<1/4):
            return 3*(1-2*xi)*(2*(xi**2)+xi+3)/((1-3*xi)*(1-4*xi))-3
        return None
    @staticmethod
    def mean(mu,sigma,xi):
        if(sigma<=0):
            raise ValueError("sigma must be greater than 0")
        if(xi<1):
            return mu+sigma/(1+xi)
        return None
    @staticmethod
    def median(mu,sigma,xi):
        if(sigma<=0):
            raise ValueError("sigma must be greater than 0")
        return mu+sigma*(math.pow(2,xi)-1)/xi
    @staticmethod
    def variance(mu,sigma,xi):
        if(sigma<=0):
            raise ValueError("sigma must be greater than 0")
        if(xi<1/2):
            return sigma**2/((1-xi)**2*(1-2*xi))
        return None
    @staticmethod
    def stddev(mu,sigma,xi):
        if(sigma<=0):
            raise ValueError("sigma must be greater than 0")
        if(xi<1/2):
            return sigma/math.sqrt((1-xi)**2*(1-2*xi))
        return None
    @staticmethod
    def entropy(mu,sigma,xi):
        if(sigma<=0):
            raise ValueError("sigma must be greater than 0")
        return math.log(sigma)+xi+1
    @staticmethod
    def skewness(mu,sigma,xi):
        if(sigma<=0):
            raise ValueError("sigma must be greater than 0")
        if(xi<1/3):
            return 2*(1+xi)*math.sqrt(1-2*xi)/(1-3*xi)
        return None
    @staticmethod
    def ppf(mu,sigma,xi,q):
        if(sigma<=0):
            raise ValueError("sigma must be greater than 0")
        if(xi==0):
            return mu-sigma*math.log(q)
        return mu+(sigma*(math.pow(q,-xi)-1)/xi)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=genpareto.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'mu':ret[0],'sigma':ret[1],'xi':ret[2]}