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

class kumaraswamygenhalfnormal(Distribution): #kumaraswamy generalized half-normal
    @staticmethod
    def random(a,b,aa,t):
        return t*math.pow(st.norm.ppf((1+math.pow(1-math.pow(1-ds.rg0(),1/b),1/a))/2),1/aa)
    @staticmethod
    def pdf(a,b,aa,t,x):
        pp=math.erf(math.pow(x/t,aa)/math.sqrt(2))
        return a*b*math.sqrt(2/math.pi)*(aa/x)*math.pow(x/t,aa)*math.exp(-1/2*math.pow(x/t,2*aa))*math.pow(pp-1,a-1)*math.pow(1-math.pow(pp-1,a),b-1)
    @staticmethod
    def cdf(a,b,aa,t,x):
        return 1-math.pow(1-math.pow(math.erf(math.pow(x/t,aa)/math.sqrt(2)),a),b)
    @staticmethod
    def median(a,b,aa,t):
        return t*math.pow(st.norm.ppf((1+math.pow(1-math.pow(1/2,1/b),1/a))/2),1/aa)
    @staticmethod
    def ppf(a,b,aa,t,q):
        return t*math.pow(st.norm.ppf((1+math.pow(1-math.pow(1-q,1/b),1/a))/2),1/aa)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=kumaraswamygenhalfnormal.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'aa':ret[2],'t':ret[3]}
