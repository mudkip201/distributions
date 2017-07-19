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
import dists.chi2.chi2 as chi2

class invchi2(Distribution): #inverse chi-squared
    @staticmethod
    def random(nu):
        k=chi2.random(nu)
        while(k==0):
            k=chi2.random(nu)
        return 1/k
    @staticmethod
    def pdf(nu,x):
        return math.pow(2,-nu/2)/math.gamma(nu/2)*math.pow(x,-nu/2-1)*math.exp(-1/(2*x))
    @staticmethod
    def kurtosis(nu):
        return 12*(5*nu-22)/((nu-6)*(nu-8))
    @staticmethod
    def mean(nu):
        if(nu>2):
            return 1/(nu-2)
    @staticmethod
    def mode(nu):
        return 1/(nu+2)
    @staticmethod
    def variance(nu):
        return 2/((nu-2)**2*(nu-4))
    @staticmethod
    def stddev(nu):
        return math.sqrt(2/((nu-2)**2*(nu-4)))
    @staticmethod
    def skewness(nu):
        return 4/(nu-6)*math.sqrt(2*(nu-4))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=invchi2.pdf(args_[0],i)
            return -tomin
        ret=op.minimize(mlefunc,(1),method='Nelder-Mead').x.tolist()
        return {'nu':ret[0]}