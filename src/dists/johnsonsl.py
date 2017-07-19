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
import dists.normal.normal as normal

class johnsonsl(Distribution):
    @staticmethod
    def random(delta,gmma,xi,lmbda):
        v=(normal.random(0,1)-gmma)/delta
        return xi+lmbda*math.exp(v)
    @staticmethod
    def pdf(delta,gmma,xi,lmbda,x):
        return delta*math.exp(-1/2*(gmma+delta*math.log((x-xi)/lmbda)))/(math.sqrt(2*math.pi)*(x-xi))
    @staticmethod
    def cdf(delta,gmma,xi,lmbda,x):
        if(xi<x and x<=xi+lmbda):
            return 1/2*math.erfc(-(gmma+delta*math.log((x-xi)/lmbda))/math.sqrt(2))
        if(x>xi+lmbda):
            return 1/2*(1+math.erf((gmma+delta*math.log((x-xi)/lmbda))/math.sqrt(2)))
    @staticmethod
    def kurtosis(delta,gmma,xi,lmbda):
        return -3+math.exp(2/delta**2)*(3+math.exp(1/delta**2)*(2+math.exp(1/delta**2)))
    @staticmethod
    def mean(delta,gmma,xi,lmbda):
        return xi+lmbda*math.exp((1-2*gmma*delta)/(2*delta**2))
    @staticmethod
    def median(delta,gmma,xi,lmbda):
        return xi+lmbda*math.exp(-gmma/delta)
    @staticmethod
    def variance(delta,gmma,xi,lmbda):
        return math.exp((1-2*gmma*delta)/(delta**2))*lmbda**2*(-1+math.exp(1/delta**2))
    @staticmethod
    def stddev(delta,gmma,xi,lmbda):
        return math.sqrt(math.exp((1-2*gmma*delta)/(delta**2))*lmbda**2*(-1+math.exp(1/delta**2)))
    @staticmethod
    def skewness(delta,gmma,xi,lmbda):
        return (2+math.exp(1/delta**2))*math.sqrt(-1+math.exp(1/delta**2))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=johnsonsl.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'delta':ret[0],'gamma':ret[1],'xi':ret[2],'lmbda':ret[3]}