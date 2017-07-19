'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''

import math
import numpy as np
from numpy import random as r
import scipy.stats as st
import scipy.special as sp
import scipy.optimize as op


class Distribution:
    @staticmethod
    def pdf():
        raise NotImplementedError("Should have implemented pdf function.")
    @staticmethod
    def cdf():
        raise NotImplementedError("Should have implemented cdf function.")
    @staticmethod
    def random():
        raise NotImplementedError("Should have implemented random number generator function.")
    @staticmethod
    def mean():
        raise NotImplementedError("Should have implemented mean function.")
    @staticmethod
    def median():
        raise NotImplementedError("Should have implemented median function.")
    @staticmethod
    def mode():
        raise NotImplementedError("Should have implemented mode function.")
    @staticmethod
    def variance():
        raise NotImplementedError("Should have implemented variance function.")
    @staticmethod
    def stddev():
        raise NotImplementedError("Should have implemented standard deviation function.")
    @staticmethod
    def kurtosis():
        raise NotImplementedError("Should have implemented kurtosis function.")
    @staticmethod
    def entropy():
        raise NotImplementedError("Should have implemented entropy function.")
    @staticmethod
    def skewness():
        raise NotImplementedError("Should have implemented skewness function.")
    @staticmethod
    def ppf():
        raise NotImplementedError("Should have implemented ppf function.")
    @staticmethod
    def mle():
        raise NotImplementedError("Should have implemented mle function.")
    
    
euler_gamma=2.71828182845904523536028747135266249

def kexp(kppa,x):
    return math.pow(math.sqrt(1+kppa**2*x**2)+kppa*x,1/kppa)

def qexp(x,q): #q-exponential
    if(q!=0):
        return math.pow(1+(1-q)*x,1/(1-q))
    return math.exp(x)

def qlog(x,q): #q-logarithm
    if(x>=0 and q==1):
        return math.log(x)
    if(x>=0 and q!=1):
        return (math.pow(x,1-q)-1)/(1-q)

def rg0():
    n=r.random()
    while(n==0):
        n=r.random()
    return n

class lefttrunc(Distribution):
    @staticmethod
    def random(dist,maxval,*args,**kwargs):
        u=dist.random(*args,**kwargs)
        while(u>maxval):
            u=dist.random(*args,**kwargs)
        return u
    @staticmethod
    def pdf(dist,maxval,*args,**kwargs):
        if(args[-1]>maxval):
            return 0
        aargs=list(args)[:-1]
        return dist.pdf(*args,**kwargs)/dist.cdf(tuple(aargs),x=maxval,**kwargs)
    @staticmethod
    def cdf(dist,maxval,*args,**kwargs):
        if(args[-1]>maxval):
            return 1
        aargs=list(args)[:-1]
        return dist.cdf(*args,**kwargs)/dist.cdf(tuple(aargs),x=maxval,**kwargs)
    @staticmethod
    def mean(dist,maxval,*args,**kwargs):
        return dist.ppf(dist.cdf(*args,x=maxval,**kwargs)/2)

class righttrunc(Distribution):
    @staticmethod
    def random(dist,minval,*args,**kwargs):
        u=dist.random(*args,**kwargs)
        while(u<minval):
            u=dist.random(*args,**kwargs)
        return u
    @staticmethod
    def pdf(dist,minval,*args,**kwargs):
        if(args[-1]<minval):
            return 0
        aargs=list(args)[:-1]
        return dist.pdf(*args,**kwargs)/(1-dist.cdf(tuple(aargs),x=minval,**kwargs))
    @staticmethod
    def cdf(dist,minval,*args,**kwargs):
        if(args[-1]<minval):
            return 0
        aargs=list(args)[:-1]
        return dist.cdf(*args,**kwargs)/(1-dist.cdf(tuple(aargs),x=minval,**kwargs))
    @staticmethod
    def mean(dist,minval,*args,**kwargs):
        return dist.ppf(1-dist.cdf(x=minval,*args,**kwargs)/2)

class lrtrunc(Distribution):
    @staticmethod
    def random(dist,minval,maxval,*args,**kwargs):
        u=dist.random(*args,**kwargs)
        while(u<minval or u>maxval):
            u=dist.random(*args,**kwargs)
        return u
    @staticmethod
    def pdf(dist,minval,maxval,*args,**kwargs):
        if(args[-1]>maxval or args[-1]<minval):
            return 0
        aargs=list(args)[:-1]
        return dist.pdf(*args,**kwargs)/(dist.cdf(tuple(aargs),x=maxval,**kwargs)-dist.cdf(tuple(aargs),x=minval,**kwargs))
    @staticmethod
    def cdf(dist,minval,maxval,*args,**kwargs):
        if(args[-1]>maxval):
            return 1
        if(args[-1]<minval):
            return 0
        aargs=list(args)[:-1]
        return (dist.cdf(*args,**kwargs)-dist.cdf(tuple(aargs),x=minval,**kwargs))/(dist.cdf(tuple(aargs),x=maxval,**kwargs)-dist.cdf(tuple(aargs),x=minval,**kwargs))
    @staticmethod
    def mean(dist,minval,maxval,*args,**kwargs):
        return dist.ppf((dist.cdf(x=maxval,*args,**kwargs)-dist.cdf(x=minval,*args,**kwargs))/2)