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

class expkumaraswamyinvweibull(Distribution):
    @staticmethod
    def random(a,b,e,l,t):
        return a/math.pow(-math.log(math.pow(1-math.pow(1-math.pow(ds.rg0(),1/t),1/e),1/l)),1/b)
    @staticmethod
    def pdf(a,b,e,l,t,x):
        return b*l*t*e*math.pow(a,b)*math.pow(x,-b-1)*math.exp(-l*math.pow(a/x,b))*math.pow(1-math.exp(-l*math.pow(a/x,b)),e-1)
    @staticmethod
    def cdf(a,b,e,l,t,x):
        return math.pow(1-math.pow(1-math.exp(-l*math.pow(a/x,b)),e),t)
    @staticmethod
    def median(a,b,e,l,t):
        return a/math.pow(-math.log(math.pow(1-math.pow(1-math.pow(1/2,1/t),1/e),1/l)),1/b)
    @staticmethod
    def ppf(a,b,e,l,t,q):
        return a/math.pow(-math.log(math.pow(1-math.pow(1-math.pow(q,1/t),1/e),1/l)),1/b)
    @staticmethod
    def mle(x): #not working
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=expkumaraswamyinvweibull.pdf(args_[0],args_[1],args_[2],args_[3],args_[4],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'e':ret[2],'l':ret[3],'t':ret[4]}