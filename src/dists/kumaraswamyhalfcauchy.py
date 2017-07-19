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

class kumaraswamyhalfcauchy(Distribution):
    @staticmethod
    def random(a,b,d):
        return d*math.tan(math.pi/2*math.pow(1-math.pow(1-ds.rg0(),1/b),1/a))
    @staticmethod
    def pdf(a,b,d,x):
        return a*b*math.pow(2/math.pi,a)/d*1/(1+(x/d)**2)*math.pow(math.atan(x/d),a-1)*math.pow(1-math.pow(2/math.pi*math.atan(x/d),a),b-1)
    @staticmethod
    def cdf(a,b,d,x):
        return 1-math.pow(1-math.pow(2/math.pi*math.atan(x/d),a),b)
    @staticmethod
    def median(a,b,d):
        return d*math.tan(math.pi/2*math.pow(1-math.pow(1/2,1/b),1/a))
    @staticmethod
    def ppf(a,b,d,q):
        return d*math.tan(math.pi/2*math.pow(1-math.pow(1-q,1/b),1/a))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=kumaraswamyhalfcauchy.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'d':ret[2]}