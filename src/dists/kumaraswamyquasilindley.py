'''
Created on Jul 18, 2017

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

class kumaraswamyquasilindley(Distribution):
    @staticmethod
    def pdf(a,b,aa,t,x):
        return a*b*t/(aa+1)*(aa+t*x)*math.exp(-t*x)*math.pow(1-math.exp(-t*x)*(1+t*x/(aa+1)),a-1)*math.pow(1-math.pow(1-math.exp(-t*x)*(1+t*x/(aa+1)),a),b-1)
    @staticmethod
    def cdf(a,b,aa,t,x):
        return 1-math.pow(1-math.pow(1-math.exp(-t*x)*(1+t*x/(aa+1)),a),b)
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean():
        pass
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance():
        pass
    @staticmethod
    def stddev():
        pass
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