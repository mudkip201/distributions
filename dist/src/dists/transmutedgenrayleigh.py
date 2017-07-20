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

class transmutedgenrayleigh(Distribution): #transmuted generalized rayleigh
    @staticmethod
    def random(a,b,l):
        i=((1+l)-math.sqrt((1+l)**2-4*l*ds.rg0()))/(2*l)
        return math.sqrt(-math.log(1-math.pow(i,1/a))/b)
    @staticmethod
    def pdf(a,b,l,x):
        return 2*a*b**2*x*math.exp(-(b*x)**2)*math.pow(1-math.exp(-(b*x)**2),a-1)*(1+l-2*l*math.pow(1-math.exp(-(b*x)**2),a))
    @staticmethod
    def cdf(a,b,l,x):
        return math.pow(1-math.exp(-(b*x)**2),a)*(1+l-2*l*math.pow(1-math.exp(-(b*x)**2),a))
    @staticmethod
    def median(a,b,l):
        i=((1+l)-math.sqrt(1+l**2)/(2*l))
        return math.sqrt(-math.log(1-math.pow(i,1/a))/b)
    @staticmethod
    def ppf(a,b,l,q):
        i=((1+l)-math.sqrt((1+l)**2-4*l*q))/(2*l)
        return math.sqrt(-math.log(1-math.pow(i,1/a))/b)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=transmutedgenrayleigh.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'l':ret[2]}