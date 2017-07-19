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

class hyperbolicsecant(Distribution):
    @staticmethod
    def pdf(x):
        return 1/(2*np.cosh(math.pi*x/2))
    @staticmethod
    def cdf(x):
        return 2/math.pi*math.atan(math.exp(math.pi*x/2))
    @staticmethod
    def random():
        return 2/math.pi*math.log(math.tan(math.pi*r.random()/2))
    @staticmethod
    def kurtosis():
        return 2
    @staticmethod
    def mean():
        return 0
    @staticmethod
    def median():
        return 0
    @staticmethod
    def mode():
        return 0
    @staticmethod
    def variance():
        return 1
    @staticmethod
    def stddev():
        return 1
    @staticmethod
    def entropy():
        return 4/math.pi*0.915965594177219015054603
    @staticmethod
    def skewness():
        return 0
    @staticmethod
    def ppf(q):
        return 1/(2*math.cosh(math.pi*q/2))