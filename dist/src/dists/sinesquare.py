'''
Created on Nov 20, 2017

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
import math.pi as pi
import dists.uniform.uniform as uniform

class sinesquare(Distribution):
    @staticmethod
    def pdf(lmbda,x):
        if(0<x and x<pi*lmbda and lmbda>0):
            return 2/(lmbda*pi)*math.sin(x/(2*lmbda))**2
        return 0
    @staticmethod
    def cdf(lmbda,x):
        if(x<=0):
            return 0
        elif(0<x and x<lmbda*pi):
            return (x/(lmbda*pi)-1/pi*math.sin(x/lmbda))
        else:
            return 1
    @staticmethod
    def random(lmbda):
        while(True):
            u1=uniform.random(0,1)
            u2=uniform.random(0,lmbda*pi)
            y=lmbda*pi*u2
            if(u1<=math.sin(pi*u2/2)**2):
                return y
        pass
    @staticmethod
    def mean(lmbda):
        return lmbda*(pi**2+4)/(2*pi)
    @staticmethod
    def median(lmbda):
        return 1.659019676*lmbda
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(lmbda):
        return (lmbda**2*pi**4-48)/(12*pi**2)
    @staticmethod
    def stddev(lmbda):
        return math.sqrt((lmbda**2*pi**4-48)/(12*pi**2))
    @staticmethod
    def kurtosis(lmbda):
        return lmbda**4*(pi**4/80+24/pi**2+128/(5*pi)-1)
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness(lmbda):
        return lmbda**3*(4*pi/3+64/(3*pi**3)+96/pi)
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle(x):
        return {'lmbda':(2*pi/(pi**2+4))*np.average(x)}