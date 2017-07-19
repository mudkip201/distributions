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

class kumaraswamyexponentiatedpareto(Distribution):
    @staticmethod
    def random(a,b,l,m,t):
        u=ds.rg0()
        return l/math.pow(1-math.pow(1-math.pow(1-u,1/b),1/(t*a)),1/m)
    @staticmethod
    def pdf(a,b,l,m,t,x):
        return a*b*t*m*math.pow(l,m)/math.pow(x,m+1)*math.pow(1-math.pow(l/x,m),t*a-1)*math.pow(1-math.pow(1-math.pow(l/x,m),t*a),b-1)
    @staticmethod
    def cdf(a,b,l,m,t,x):
        return 1-math.pow(1-math.pow(1-math.pow(l/x,m),t*a),b)
    @staticmethod
    def median(a,b,l,m,t):
        return l/math.pow(1-math.pow(1-math.pow(1/2,1/b),1/(t*a)),1/m)
    @staticmethod
    def ppf(a,b,l,m,t,q):
        return l/math.pow(1-math.pow(1-math.pow(1-q,1/b),1/(t*a)),1/m)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=kumaraswamyexponentiatedpareto.pdf(args_[0],args_[1],args_[2],args_[3],args_[4],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'l':ret[2],'m':ret[3],'t':ret[4]}