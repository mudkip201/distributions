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

class exponentiatedweibulllog(Distribution): #exponentiated weibull logarithmic
    @staticmethod
    def random(a,b,g,t):
        u=ds.rg0()
        return math.pow(-math.log(1-math.pow((1-math.pow(1-t,u))/t,1/a)),1/g)/b
    @staticmethod
    def pdf(a,b,g,t,x):
        return a*t*g*math.pow(b,g)*math.pow(x,g-1)*math.exp(-math.pow(b*x,g))*math.pow(1,math.exp(-math.pow(b*x,g)),a-1)/(math.log(1-t)*(t*math.pow(1-math.exp(-math.pow(b*x,g)),a)-1))
    @staticmethod
    def cdf(a,b,g,t,x):
        return math.log(1-t*math.pow(1-math.exp(-math.pow(b*x,g)),a))/math.log(1-t)
    @staticmethod
    def median(a,b,g,t):
        return math.pow(-math.log(1-math.pow((1-math.pow(1-t,1/2))/t,1/a)),1/g)/b
    @staticmethod
    def ppf(a,b,g,t,q):
        return math.pow(-math.log(1-math.pow((1-math.pow(1-t,q))/t,1/a)),1/g)/b
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=exponentiatedweibulllog.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'g':ret[2],'t':ret[3]}
