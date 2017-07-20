'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r
import scipy.optimize as op

class expweibull(Distribution): #exponentiated weibull
    @staticmethod
    def cdf(aa,k,lmbda,x):
        if(aa<=0 or k<=0 or lmbda<=0 or x<=0):
            raise ValueError("All inputs must be positive")
        return math.pow(1-math.exp(-math.pow(x/lmbda,k)),aa)
    @staticmethod
    def pdf(aa,k,lmbda,x):
        if(aa<=0 or k<=0 or lmbda<=0 or x<=0):
            raise ValueError("All inputs must be positive")
        return aa*(k/lmbda)*math.pow(x/lmbda,k-1)*math.pow(1-math.exp(-math.pow(x/lmbda,k)),aa-1)*math.exp(-math.pow(x/lmbda,k))
    @staticmethod
    def random(aa,k,lmbda):
        if(aa<=0 or k<=0 or lmbda<=0):
            raise ValueError("All inputs must be positive")
        return aa*math.pow(-math.log(1-math.pow(r.random(),1/k)),1/lmbda)
    @staticmethod
    def median(aa,k,lmbda):
        if(aa<=0 or k<=0 or lmbda<=0):
            raise ValueError("All inputs must be positive")
        return aa*math.pow(-math.log(1-math.pow(1/2,1/k)),1/lmbda)
    @staticmethod
    def ppf(aa,k,lmbda,q):
        if(aa<=0 or k<=0 or lmbda<=0):
            raise ValueError("All inputs must be positive")
        return aa*math.pow(-math.log(1-math.pow(q,1/k)),1/lmbda)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=expweibull.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'aa':ret[0],'k':ret[1],'lmbda':ret[2]}