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

class expkumaraswamydagum(Distribution): #exponentiated kumaraswamy-dagum
    @staticmethod
    def random(a,d,l,p,t):
        return math.pow(l,1/d)*math.pow(math.pow(1-math.pow((1-math.pow(ds.rg0(),1/t)),1/p),-1/a)-1,-1/d)
    @staticmethod
    def pdf(a,d,l,p,t,x):
        return math.pow(1-math.pow(1-math.pow(1+l*math.pow(x,-d),-a),p),t)
    @staticmethod
    def cdf(a,d,l,p,t,x):
        return a*l*d*p*t*math.pow(x,-d-1)*math.pow(1+l*math.pow(x,-d),-a-1)*math.pow(1-math.pow(1+l*math.pow(x,-d),-a),p-1)*math.pow(1-math.pow(1-math.pow(1+l*math.pow(x,-d),-a),p),t-1)
    @staticmethod
    def median(a,d,l,p,t):
        return math.pow(l,1/d)*math.pow(math.pow(1-math.pow((1-math.pow(1/2,1/t)),1/p),-1/a)-1,-1/d)
    @staticmethod
    def ppf(a,d,l,p,t,q):
        return math.pow(l,1/d)*math.pow(math.pow(1-math.pow((1-math.pow(q,1/t)),1/p),-1/a)-1,-1/d)
    @staticmethod
    def mle(x): #not working
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=expkumaraswamydagum.pdf(args_[0],args_[1],args_[2],args_[3],args_[4],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'d':ret[1],'l':ret[2],'p':ret[3],'t':ret[4]}