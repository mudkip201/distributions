'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r
import scipy.stats as st
import scipy.optimize as op
import dists.normal.normal as normal

class birnbaumsaunders(Distribution): #birnbaum-saunders
    '''Birnbaum-Saunders'''
    @staticmethod
    def random(aa,bb):
        n=normal.random(0,1)
        return bb*(aa*n/2+math.sqrt((aa*n/2)**2+1))**2
    @staticmethod
    def pdf(aa,bb,x):
        return (math.sqrt(x/bb)+math.sqrt(bb/x))/(2*aa*x)*normal.pdf((math.sqrt(x/bb)-math.sqrt(bb/x))/aa)
    @staticmethod
    def cdf(aa,bb,x):
        return normal.cdf(1/aa*(math.sqrt(x/bb)-math.sqrt(bb/x)))
    @staticmethod
    def mean(aa,bb):
        return bb*(1+aa**2/2)
    @staticmethod
    def median(aa,bb):
        return bb
    @staticmethod
    def variance(aa,bb):
        return (aa*bb)**2*(1+5*aa**2/2)
    @staticmethod
    def stddev(aa,bb):
        return aa*bb*math.sqrt(1+5*aa**2/2)
    @staticmethod
    def skewness(aa,bb):
        return (4*aa)*(11*aa**2+6)/math.pow(5*aa**2+4,3/2)
    @staticmethod
    def kurtosis(aa,bb):
        return 6*aa**2*(93*aa**2+41)/(5*aa**2+4)**2
    @staticmethod
    @staticmethod
    def ppf(aa,bb,q):
        return 1/4*((aa*normal.ppf(0,1,q)+math.sqrt(4+(aa*normal.ppf(0,1,q))**2))/bb)**2
    @staticmethod
    def mle(x):#not working
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=birnbaumsaunders.pdf(args_[0],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        #ret=op.differential_evolution(mlefunc,[(-10,10)]).x.tolist()
        return {'gamma':ret[0]}