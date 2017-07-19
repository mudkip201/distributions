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

class kumaraswamyweibullpoisson(Distribution):
    @staticmethod
    def random(a,b,bb,c,l):
        u=ds.rg0()
        return math.pow(-math.log(1+math.log(1-u*(1-math.exp(-l)))/l),1/c)/bb
    @staticmethod
    def pdf(a,b,bb,c,l,x):
        return 1-math.exp(-l*(1-math.exp(-math.pow(bb*x,c))))/(1-math.exp(-l))
    @staticmethod
    def cdf(a,b,bb,c,l,x):
        z=1-math.exp(-math.pow(bb*x,c))
        return l*a*b*c*math.pow(bb,c)/(math.exp(l)-1)*math.pow(x,c-1)*math.pow(z,a-1)*math.pow(1-math.pow(z,a),b-1)*math.exp(l*math.pow(1-math.pow(z,a),b)-math.pow(bb*x,c))
    @staticmethod
    def median(a,b,bb,c,l):
        return math.pow(-math.log(1+math.log(1-(1-math.exp(-l))/2)/l),1/c)/bb
    @staticmethod
    def ppf(a,b,bb,c,l,q):
        return math.pow(-math.log(1+math.log(1-q*(1-math.exp(-l)))/l),1/c)/bb
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=kumaraswamyweibullpoisson.pdf(args_[0],args_[1],args_[2],args_[3],args_[4],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'bb':ret[2],'c':ret[3],'l':ret[4]}