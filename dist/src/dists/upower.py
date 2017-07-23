'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r
import scipy.optimize as op

class upower(Distribution):
    @staticmethod
    def random(k):
        return math.pow((2*r.random()-1),1/(2*k+1))
    @staticmethod
    def pdf(k,x):
        return (2*k+1)/2*math.pow(x,2*k)
    @staticmethod
    def cdf(k,x):
        return 1/2*(1+math.pow(x,2*k+1))
    @staticmethod
    def kurtosis(k):
        return (2*k+3)**2/((2*k+5)*(2*k+1))
    @staticmethod
    def mean(k):
        return 0
    @staticmethod
    def median(k):
        return 0
    @staticmethod
    def variance(k):
        return (2*k+1)/(2*k+3)
    @staticmethod
    def stddev(k):
        return math.sqrt((2*k+1)/(2*k+3))
    @staticmethod
    def skewness(k):
        return 0
    @staticmethod
    def ppf(k,q):
        return math.pow((2*q-1),1/(2*k+1))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=upower.pdf(args_[0],i)
            return -tomin
        ret=op.minimize(mlefunc,(1),method='Nelder-Mead').x.tolist()
        return {'k':ret[0]}
