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

class birnbaumsaunders(Distribution): #standard fatigue-life
    '''Fatigue-life'''
    @staticmethod
    def random(gmma):
        n=r.random()
        return 1/4.0*math.pow(gmma*st.norm.ppf(n)+math.sqrt(4+math.pow(gmma*st.norm.ppf(n),2)),2)
    @staticmethod
    def pdf(gmma,x):
        return (math.sqrt(x)+math.sqrt(1/x))/(2*gmma*x)*st.norm.pdf((math.sqrt(x)-math.sqrt(1/x))/gmma)
    @staticmethod
    def cdf(gmma,x):
        return st.norm.cdf((math.sqrt(x)-math.sqrt(1/x))/gmma)
    @staticmethod
    def mean(gmma):
        return 0
    @staticmethod
    def median(gmma):
        return 1/4.0*math.pow(gmma*st.norm.ppf(1/2)+math.sqrt(4+math.pow(gmma*st.norm.ppf(1/2),2)),2)
    @staticmethod
    def ppf(gmma,q):
        return 1/4.0*math.pow(gmma*st.norm.ppf(q)+math.sqrt(4+math.pow(gmma*st.norm.ppf(q),2)),2)
    @staticmethod
    def mle(x):#not working
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=birnbaumsaunders.pdf(args_[0],i)
            return -tomin
        ret=op.minimize(mlefunc,(2),method='Nelder-Mead').x.tolist()
        #ret=op.differential_evolution(mlefunc,[(-10,10)]).x.tolist()
        return {'gamma':ret[0]}