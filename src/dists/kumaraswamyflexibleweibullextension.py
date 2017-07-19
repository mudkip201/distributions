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

class kumaraswamyflexibleweibullextension(Distribution):
    @staticmethod
    def random(a,b,aa,bb):
        u=ds.rg0()
        i=math.log(-math.log(1-math.pow(1-math.pow(1-u,1/b),1/a)))
        return 1/(2*aa)*(i+math.sqrt(i**2+4*aa*bb))
    @staticmethod
    def pdf(a,b,aa,bb,x):
        return a*b*(aa+bb/x**2)*math.exp(aa*x-bb/x)*math.exp(-math.exp(aa*x-bb/x))*math.pow(1-math.exp(-math.exp(aa*x-bb/x)),a-1)*math.pow(1-math.pow(1-math.exp(-math.exp(aa*x-bb/x)),a),b-1)
    @staticmethod
    def cdf(a,b,aa,bb,x):
        return 1-math.pow(1-math.pow(1-math.exp(-math.exp(aa*x-bb/x)),a),b)
    @staticmethod
    def median(a,b,aa,bb):
        i=math.log(-math.log(1-math.pow(1-math.pow(1/2,1/b),1/a)))
        return 1/(2*aa)*(i+math.sqrt(i**2+4*aa*bb))
    @staticmethod
    def ppf(a,b,aa,bb,q):
        i=math.log(-math.log(1-math.pow(1-math.pow(1-q,1/b),1/a)))
        return 1/(2*aa)*(i+math.sqrt(i**2+4*aa*bb))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=kumaraswamyflexibleweibullextension.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'aa':ret[2],'bb':ret[3]}
