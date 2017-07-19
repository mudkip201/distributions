'''
Created on Jul 17, 2017

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

class gammalomax(Distribution):
    @staticmethod
    def pdf(aa,b,a,x):
        return aa*math.pow(b,aa)/math.gamma(a)*math.pow(b+x,-aa-1)*math.pow(-aa*math.log(b/(b+x)),a-1)
    @staticmethod
    def cdf(aa,b,a,x):
        return sp.gammainc(a,-aa*math.log(b/(b+x)))
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(aa,b,a):
        return b/(aa-1)
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode(aa,b,a):
        return b*(math.exp((a-1)/(aa-1))-1)
        pass
    @staticmethod
    def variance(aa,b,a):
        return aa*b**2/((aa-1)**2*(aa-2))
    @staticmethod
    def stddev(aa,b,a):
        return math.sqrt(gammalomax.variance(aa,b,a))
    @staticmethod
    def kurtosis():
        pass
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness():
        pass
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle():
        pass