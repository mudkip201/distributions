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

class gentoppleone(Distribution): #generalized topp-leone
    @staticmethod
    def random(aa,bb):
        if(aa==1):
            return math.pow(r.random(),1/bb)
        return (-aa+math.sqrt(math.pow(aa,2)-4*(1-aa)*(-math.pow(r.random(),1/bb))))/(2-2*aa)
    @staticmethod
    def pdf(aa,bb,x):
        return bb*math.pow(aa*x-(aa-1)*x**2,bb-1)*(aa-2*(aa-1)*x)
    @staticmethod
    def cdf(aa,bb,x):
        return math.pow(x,bb)*math.pow(aa-(aa-1)*x,bb)
    @staticmethod
    def median(aa,bb):
        if(aa==1):
            return math.pow(1/2,1/bb)
        return (-aa+math.sqrt(math.pow(aa,2)-4*(1-aa)*(-math.pow(1/2,1/bb))))/(2-2*aa)
    @staticmethod
    def ppf(aa,bb,q):
        if(aa==1):
            return math.pow(q,1/bb)
        return (-aa+math.sqrt(math.pow(aa,2)-4*(1-aa)*(-math.pow(q,1/bb))))/(2-2*aa)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=gentoppleone.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'aa':ret[0],'bb':ret[1]}
