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

class erlang(Distribution):
    @staticmethod
    def pdf(k,lmbda,x):
        if(k%1!=0 or k<=0):
            raise ValueError("k must be a positive integer")
        if(lmbda<=0 or x<=0):
            raise ValueError("lmbda and x must be bigger than 0")
        return math.pow(lmbda,k)*math.pow(x,k-1)*math.exp(-lmbda*x)/math.factorial(k-1)
    @staticmethod
    def random(k,lmbda):
        if(k%1!=0 or k<=0):
            raise ValueError("k must be a positive integer")
        if(lmbda<=0):
            raise ValueError("lmbda must be bigger than 0")
        total=1
        for _ in range(k):
            total*=ds.rg0()
        #total=1
        #for _ in range(k):
        #total+=exponential.random(lmbda/(float(k)))
        return -math.log(total)/lmbda
    @staticmethod
    def kurtosis(k,lmbda):
        if(k%1!=0 or k<=0):
            raise ValueError("k must be a positive integer")
        if(lmbda<=0):
            raise ValueError("lmbda must be bigger than 0")
        return 6/k
    @staticmethod
    def mean(k,lmbda):
        if(k%1!=0 or k<=0):
            raise ValueError("k must be a positive integer")
        if(lmbda<=0):
            raise ValueError("lmbda must be bigger than 0")
        return k/lmbda
    @staticmethod
    def mode(k,lmbda):
        if(k%1!=0 or k<=0):
            raise ValueError("k must be a positive integer")
        if(lmbda<=0):
            raise ValueError("lmbda must be bigger than 0")
        return 1/lmbda*(k-1)
    @staticmethod
    def variance(k,lmbda):
        if(k%1!=0 or k<=0):
            raise ValueError("k must be a positive integer")
        if(lmbda<=0):
            raise ValueError("lmbda must be bigger than 0")
        return k/lmbda**2
    @staticmethod
    def stddev(k,lmbda):
        if(k%1!=0 or k<=0):
            raise ValueError("k must be a positive integer")
        if(lmbda<=0):
            raise ValueError("lmbda must be bigger than 0")
        return math.sqrt(k)/lmbda
    @staticmethod
    def skewness(k,lmbda):
        if(k%1!=0 or k<=0):
            raise ValueError("k must be a positive integer")
        if(lmbda<=0):
            raise ValueError("lmbda must be bigger than 0")
        return 2/math.sqrt(k)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=erlang.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'k':ret[0],'lmbda':ret[1]}