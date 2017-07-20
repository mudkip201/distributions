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
import dists.chi.chi as chi

class nakagami(Distribution):
    @staticmethod
    def random(m,omega):
        return math.sqrt(omega/(2*m))*chi.random(2*m)
    @staticmethod
    def pdf(m,omega,x):
        return 2*math.pow(m,m)/(math.gamma(m)*math.pow(omega,m))*math.pow(x,2*m-1)*math.exp(-m/omega*x**2)
    @staticmethod
    def mean(m,omega):
        return math.gamma(m+1/2)/math.gamma(m)*math.sqrt(omega/m)
    @staticmethod
    def mode(m,omega):
        return math.sqrt(1/2)*math.sqrt((2*m-1)*omega/m)
    @staticmethod
    def variance(m,omega):
        return omega*(1-1/m*(math.gamma(m+1/2)/math.gamma(m))**2)
    @staticmethod
    def stddev(m,omega):
        return math.sqrt(omega*(1-1/m*(math.gamma(m+1/2)/math.gamma(m))**2))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=nakagami.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'m':ret[0],'omega':ret[1]}