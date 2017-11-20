'''
Created on Nov 20, 2017

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

class correlatedstduniformgeomin(Distribution):
    @staticmethod
    def pdf(theta,x):
        if(theta<=x and x<=1):
            return theta/((1-theta)*x**2)
        return 0
    @staticmethod
    def cdf():
        pass
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(theta):
        return (-theta*math.log(theta))/(1-theta)
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(theta):
        return (theta**3-2*theta**2-theta**2*math.log(theta)**2+theta)/(1-theta)**2
    @staticmethod
    def stddev(theta):
        return math.sqrt((theta**3-2*theta**2-theta**2*math.log(theta)**2+theta)/(1-theta)**2)
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
    def mle(x):
        return {'theta':max(x)}