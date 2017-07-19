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

class lomax(Distribution):
    @staticmethod
    def pdf(lmbda,aa,x):
        if(lmbda<=0 or aa<=0):
            raise ValueError("lambda and aa must be positive")
        if(x<0):
            raise ValueError("x must be non-negative")
        return aa/lmbda*math.pow(1+x/lmbda,-aa-1)
    @staticmethod
    def cdf(lmbda,aa,x):
        if(lmbda<=0 or aa<=0):
            raise ValueError("lambda and aa must be positive")
        if(x<0):
            raise ValueError("x must be non-negative")
        return 1-math.pow(1+x/lmbda,-aa)
    @staticmethod
    def random(lmbda,aa):
        if(lmbda<=0 or aa<=0):
            raise ValueError("lambda and aa must be positive")
        n=r.random()
        while(n==1):
            n=r.random()
        return lmbda*(math.pow(1/(1-n),1/aa)-1)
    @staticmethod
    def kurtosis(lmbda,aa):
        if(lmbda<=0 or aa<=0):
            raise ValueError("lambda and aa must be positive")
        if(aa>4):
            return 6*(aa**3+aa**2-6*aa-2)/(aa*(aa-3)*(aa-4))
        return None
    @staticmethod
    def mean(lmbda,aa):
        if(lmbda<=0 or aa<=0):
            raise ValueError("lambda and aa must be positive")
        if(aa>1):
            return lmbda/(aa-1)
        return None
    @staticmethod
    def median(lmbda,aa):
        if(lmbda<=0 or aa<=0):
            raise ValueError("lambda and aa must be positive")
        return lmbda*(math.pow(2,1/aa)-1)
    @staticmethod
    def mode(lmbda,aa):
        if(lmbda<=0 or aa<=0):
            raise ValueError("lambda and aa must be positive")
        return 0
    @staticmethod
    def variance(lmbda,aa):
        if(lmbda<=0 or aa<=0):
            raise ValueError("lambda and aa must be positive")
        if(aa>2):
            return lmbda**2*aa/((aa-1)**2*(aa-2))
        if(aa>1):
            return float("infinity")
        return None
    @staticmethod
    def stddev(lmbda,aa):
        if(lmbda<=0 or aa<=0):
            raise ValueError("lambda and aa must be positive")
        if(aa>2):
            return lmbda*math.sqrt(aa/(aa-2))/(aa-1)
        if(aa>1):
            return float("infinity")
        return None
    @staticmethod
    def skewness(lmbda,aa):
        if(lmbda<=0 or aa<=0):
            raise ValueError("lambda and aa must be positive")
        if(aa>3):
            return 2*(1+aa)/(aa-3)*math.sqrt((aa-2)/aa)
        return None
    @staticmethod
    def ppf(lmbda,aa,q):
        if(lmbda<=0 or aa<=0):
            raise ValueError("lambda and aa must be positive")
        return lmbda*(math.pow(1/(1-q),1/aa)-1)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=lomax.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'lmbda':ret[0],'aa':ret[1]}