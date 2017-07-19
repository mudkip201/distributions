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

class exponentiatedexpgeo(Distribution): #exponentiated exponential-geometric
    @staticmethod
    def random(a,l,t):
        u=ds.rg0()
        return -math.log(1-math.pow((u*t)/(1-u+u*t),1/a))/l
    @staticmethod
    def pdf(a,l,t,x):
        return a*l*t*math.exp(-l*x)*math.pow(1-math.exp(-l*x),a-1)/math.pow(1-(1-t)*math.pow(1-math.exp(-l*x),a),2)
    @staticmethod
    def cdf(a,l,t,x):
        return 1-(1-math.pow(1-math.exp(-l*x),a))/(1-(1-t)*math.pow(1-math.exp(-l*x),a))
    @staticmethod
    def median(a,l,t):
        return -math.log(1-math.pow(t/(1+t),1/a))/l
    @staticmethod
    def ppf(a,l,t,q):
        return -math.log(1-math.pow((q*t)/(1-q+q*t),1/a))/l
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=exponentiatedexpgeo.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'l':ret[1],'t':ret[2]}