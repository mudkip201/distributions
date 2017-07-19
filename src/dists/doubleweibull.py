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

class doubleweibull(Distribution):
    @staticmethod
    def random(a,b,c):
        y=r.random()
        if y<=0.5:
            return -b*math.pow(math.log(1/(2*y)),1/c)+a
        return b*math.pow(math.log(1/(2*y-1)),1/c)+a
    @staticmethod
    def pdf(a,b,c,x):
        z=(x-a)/b
        return c/2*math.pow(abs(z),c-1)*math.exp(-math.pow(abs(z),c))
    @staticmethod
    def cdf(a,b,c,x):
        z=(x-a)/b
        if(z<=0):
            return 1/2*math.exp(-math.pow(abs(z),c))
        return 1-1/2*math.exp(-math.pow(abs(z),c))
    @staticmethod
    def median(a,b,c):
        return a
    @staticmethod
    def mean(a,b,c):
        return a
    @staticmethod
    def variance(a,b,c):
        return math.gamma((c+2)/c*b*b)
    @staticmethod
    def stddev(a,b,c):
        return math.sqrt(math.gamma((c+2)/c*b*b))
    @staticmethod
    def skewness(a,b,c):
        return 0
    @staticmethod
    def ppf(a,b,c,q):
        if q<=0.5:
            return -b*math.pow(math.log(1/(2*q)),1/c)+a
        return b*math.pow(math.log(1/(2*q-1)),1/c)+a
    @staticmethod
    def mle(x): #not working
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=doubleweibull.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(0,5,5),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'c':ret[2]}