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

class yulesimon(Distribution):
    @staticmethod
    def random(rho):
        a=r.random()
        i=1
        while True:
            b=1-i*sp.beta(i,a+1)
            if(b>a):
                return i
            i+=1
    @staticmethod
    def pdf(rho,x):
        return rho*sp.beta(x,rho+1)
    @staticmethod
    def cdf(rho,x):
        return 1-x*sp.beta(x,rho+1)
    @staticmethod
    def kurtosis(rho):
        if(rho>4):
            return rho+3+(11*rho**3-49*rho-22)/((rho-4)*(rho-3)*rho)
    @staticmethod
    def mean(rho):
        if(rho>1):
            return rho/(rho-1)
    @staticmethod
    def mode(rho):
        return 1
    @staticmethod
    def variance(rho):
        if(rho>2):
            return rho**2/((rho-1)**2*(rho-2))
    @staticmethod
    def stddev(rho):
        if(rho>2):
            return rho/math.sqrt((rho-1)**2*(rho-2))
    @staticmethod
    def skewness(rho):
        if(rho>3):
            return (rho+1)**2*math.sqrt(rho-2)/((rho-3)*rho)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=yulesimon.pdf(args_[0],i)
            return -tomin
        ret=op.minimize(mlefunc,(1),method='Nelder-Mead').x.tolist()
        return {'rho':ret[0]}