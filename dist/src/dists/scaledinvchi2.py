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
import dists.invchi2.invchi2 as invchi2

class scaledinvchi2(Distribution): #scaled inverse chi-squared
    @staticmethod
    def random(nu,t2):
        return invchi2(nu)*t2*nu
    @staticmethod
    def pdf(nu,t2,x):
        return math.pow(t2*nu/2,nu/2)/math.gamma(nu/2)*math.exp((-nu*t2/(2*x)))/math.pow(x,1+nu/2)
    @staticmethod
    def kurtosis(nu,t2):
        if(nu>8):
            return 12*(5*nu-22)/((nu-6)*(nu-8))
    @staticmethod
    def mean(nu,t2):
        if(nu>2):
            return nu*t2/(nu-2)
    @staticmethod
    def mode(nu,t2):
        return nu*t2/(nu+2)
    @staticmethod
    def variance(nu,t2):
        if(nu>4):
            return 2*nu**2*t2**2/((nu-2)**2*(nu-4))
    @staticmethod
    def stddev(nu,t2):
        if(nu>4):
            return 2*nu*t2/math.sqrt((nu-2)**2*(nu-4))
    @staticmethod
    def skewness(nu,t2):
        if(nu>6):
            return 4/(nu-6)*math.sqrt(2*(nu-4))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=scaledinvchi2.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'nu':ret[0],'t2':ret[1]}