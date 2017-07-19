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

class kumaraswamydagum(Distribution):
    @staticmethod
    def random(a,b,bb,d,l):
        u=ds.rg0()
        return math.pow((math.pow(1-math.pow(1-u,1/b),1/(-bb*a))-1)/l,-1/d)
    @staticmethod
    def pdf(a,b,bb,d,l,x):
        return a*b*bb*l*d*math.pow(x,-d-1)*math.pow(1+l*math.pow(x,-d),-bb*a-1)*math.pow(1-math.pow(1+l*math.pow(x,-d),-bb*a),-b-1)
    @staticmethod
    def cdf(a,b,bb,d,l,x):
        return 1-math.pow(1-math.pow(1+l*math.pow(x,-d),-bb*a),b)
    @staticmethod
    def median(a,b,bb,d,l):
        return math.pow((math.pow(1-math.pow(1/2,1/b),1/(-bb*a))-1)/l,-1/d)
    @staticmethod
    def ppf(a,b,bb,d,l,q):
        return math.pow((math.pow(1-math.pow(1-q,1/b),1/(-bb*a))-1)/l,-1/d)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=kumaraswamydagum.pdf(args_[0],args_[1],args_[2],args_[3],args_[4],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'bb':ret[2],'d':ret[3],'l':ret[4]}