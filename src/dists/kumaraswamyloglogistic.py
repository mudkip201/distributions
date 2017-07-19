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

class kumaraswamyloglogistic(Distribution):
    @staticmethod
    def random(a,b,aa,g):
        return aa*(math.pow(1-math.pow(1-math.pow(1-ds.rg0(),1/b),1/a),-1/g)-1)
    @staticmethod
    def pdf(a,b,aa,g,x):
        return a*b*g/math.pow(aa,a*g)*math.pow(x,a*g-1)*math.pow(1+math.pow(x/aa,g),-a-1)*math.pow(1-math.pow(1-1/(1+(x/aa)*g),a),b-1)
    @staticmethod
    def cdf(a,b,aa,g,x):
        return 1-math.pow(1-math.pow(1-1/(1+(x/aa)*g),a),b)
    @staticmethod
    def median(a,b,aa,g):
        return aa*(math.pow(1-math.pow(1-math.pow(1/2,1/b),1/a),-1/g)-1)
    @staticmethod
    def ppf(a,b,aa,g,q):
        return aa*(math.pow(1-math.pow(1-math.pow(1-q,1/b),1/a),-1/g)-1)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=kumaraswamyloglogistic.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'aa':ret[2],'g':ret[3]}
