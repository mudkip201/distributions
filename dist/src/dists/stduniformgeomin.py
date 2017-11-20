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

class MyClass(Distribution):
    @staticmethod
    def pdf(theta,x):
        if(0<x and x<1):
            return theta/(theta+(1-theta)*x)**2
        return 0
    @staticmethod
    def cdf(theta,x):
        pass
    @staticmethod
    def random(theta):
        pass
    @staticmethod
    def mean(theta):
        return theta*(theta-1-math.log(theta))/(1-theta)**2
    @staticmethod
    def median(theta):
        pass
    @staticmethod
    def mode(theta):
        pass
    @staticmethod
    def variance(theta):
        return (theta**3-2*theta**2-theta**2*math.log(theta)**2+theta)/(1-theta)**4
    @staticmethod
    def stddev(theta):
        return math.sqrt((theta**3-2*theta**2-theta**2*math.log(theta)**2+theta)/(1-theta)**4)
    @staticmethod
    def kurtosis(theta):
        pass
    @staticmethod
    def entropy(theta):
        pass
    @staticmethod
    def skewness(theta):
        pass
    @staticmethod
    def ppf(theta,q):
        pass
    @staticmethod
    def mle(x):
        pass