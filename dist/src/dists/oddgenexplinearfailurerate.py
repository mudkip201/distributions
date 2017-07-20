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

class oddgenexplinearfailurerate(Distribution): #odd generalized exponentiated linear failure rate
    @staticmethod
    def random(a,b,aa,bb):
        i=1+math.log(1/math.pow(1-math.pow(ds.rg0(),1/bb),1/aa))
        return (-a+math.sqrt(a*a+2*b*math.log(i)))/b
    @staticmethod
    def pdf(a,b,aa,bb,x):
        aa*bb*(a+b*x)*math.exp(a*x+b/2*x**2)*math.exp(-aa*(math.exp(a*x+b/2*x**2-1)))*math.pow(1-math.exp(-aa*(math.exp(a*x+b/2*x**2-1))),bb-1)
    @staticmethod
    def cdf(a,b,aa,bb,x):
        return math.pow(1-math.exp(-aa*(math.exp(a*x+b/2*x**2-1))),bb)
    @staticmethod
    def median(a,b,aa,bb):
        i=1+math.log(1/math.pow(1-math.pow(1/2,1/bb),1/aa))
        return (-a+math.sqrt(a*a+2*b*math.log(i)))/b
    @staticmethod
    def ppf(a,b,aa,bb,q):
        i=1+math.log(1/math.pow(1-math.pow(q,1/bb),1/aa))
        return (-a+math.sqrt(a*a+2*b*math.log(i)))/b
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=oddgenexplinearfailurerate.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'aa':ret[2],'bb':ret[3]}
