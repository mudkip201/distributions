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

class pareto(Distribution):
    @staticmethod
    def pdf(aa,xm,x):
        if(aa<=0 or xm<=0):
            raise ValueError("alpha and xm must be positive")
        if(x<xm):
            raise ValueError("x must be greater than or equal to xm")
        return aa*math.pow(xm,aa)/math.pow(x,aa+1)
    @staticmethod
    def cdf(aa,xm,x):
        if(aa<=0 or xm<=0):
            raise ValueError("alpha and xm must be positive")
        if(x<xm):
            raise ValueError("x must be greater than or equal to xm")
        return 1-math.pow(xm/x,aa)
    @staticmethod
    def random(aa,xm):
        if(aa<=0 or xm<=0):
            raise ValueError("alpha and xm must be positive")
        return xm/math.pow(ds.rg0(),1/aa)
    @staticmethod
    def kurtosis(aa,xm):
        if(aa<=0 or xm<=0):
            raise ValueError("alpha and xm must be positive")
        if(aa>4):
            return 6*(aa**3+aa**2-6*aa-2)/(aa*(aa-3)*(aa-4))
        return None
    @staticmethod
    def mean(aa,xm):
        if(aa<=0 or xm<=0):
            raise ValueError("alpha and xm must be positive")
        if(aa<=1):
            return float("infinity")
        return aa*xm/(aa-1)
    @staticmethod
    def median(aa,xm):
        if(aa<=0 or xm<=0):
            raise ValueError("alpha and xm must be positive")
        return xm*math.pow(2,1/aa)
    @staticmethod
    def mode(aa,xm):
        if(aa<=0 or xm<=0):
            raise ValueError("alpha and xm must be positive")
        return xm
    @staticmethod
    def variance(aa,xm):
        if(aa<=0 or xm<=0):
            raise ValueError("alpha and xm must be positive")
        if(aa<=2):
            return float("infinity")
        return (xm**2)*aa/((aa-1)**2*(aa-2))
    @staticmethod
    def stddev(aa,xm):
        if(aa<=0 or xm<=0):
            raise ValueError("alpha and xm must be positive")
        if(aa<=2):
            return float("infinity")
        return xm*math.sqrt(aa/(aa-2))/(aa-1)
    @staticmethod
    def entropy(aa,xm):
        if(aa<=0 or xm<=0):
            raise ValueError("alpha and xm must be positive")
        return math.log((xm/aa)*math.exp(1+1/aa))
    @staticmethod
    def skewness(aa,xm):
        if(aa<=0 or xm<=0):
            raise ValueError("alpha and xm must be positive")
        if(aa>3):
            return 2*(1+aa)/(aa-3)*math.sqrt((aa-2)/aa)
        return None
    @staticmethod
    def ppf(aa,xm,q):
        if(aa<=0 or xm<=0):
            raise ValueError("alpha and xm must be positive")
        return xm/math.pow(q,1/aa)
    @staticmethod
    def mle(x):
        params={}
        params['xm']=min(x)
        params['aa']=len(x)/np.sum(np.log(x)-math.log(params['xm']))
        return params