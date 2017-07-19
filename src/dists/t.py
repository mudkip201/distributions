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
import dists.f.f as f

class t(Distribution):
    @staticmethod
    def random(nu):
        return math.sqrt(f.random(1,nu))
    @staticmethod
    def pdf(nu,x):
        return math.gamma((nu+1)/2)/(math.sqrt(nu*math.pi)*math.gamma(nu/2))*math.pow(1+x**2/nu,-(nu+1)/2)
    @staticmethod
    def kurtosis(nu):
        if(nu>4):
            return 6/(nu-4)
        if(2<nu and nu<=4):
            return float("infinity")
        return None
    @staticmethod
    def mean(nu):
        if(nu>1):
            return 0
        return None
    @staticmethod
    def median(nu):
        return 0
    @staticmethod
    def mode(nu):
        return 0
    @staticmethod
    def variance(nu):
        if(nu>2):
            return nu/(nu-2)
        if(1<nu and nu<=2):
            return float("infinity")
        return None
    @staticmethod
    def stddev(nu):
        if(nu>2):
            return math.sqrt(nu/(nu-2))
        if(1<nu and nu<=2):
            return float("infinity")
        return None
    @staticmethod
    def skewness(nu):
        if(nu>3):
            return 0
        return None
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=t.pdf(args_[0],i)
            return -tomin
        ret=op.minimize(mlefunc,(1),method='Nelder-Mead').x.tolist()
        return {'nu':ret[0]}