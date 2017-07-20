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

class uniform(Distribution):
    @staticmethod
    def pdf(a,b,x):
        if(b<=a):
            raise ValueError("b must be greater than a")
        if(x<a or x>b):
            raise ValueError("x must be between a and b inclusive")
        return 1/(b-a)
    @staticmethod
    def cdf(a,b,x):
        if(b<=a):
            raise ValueError("b must be greater than a")
        if(x<a or x>b):
            raise ValueError("x must be between a and b inclusive")
        return (x-a)/(b-a)
    @staticmethod
    def random(a,b):
        if(b<=a):
            raise ValueError("b must be greater than a")
        return ((b-a)*r.random()+a)
    @staticmethod
    def kurtosis(a,b):
        if(b<=a):
            raise ValueError("b must be greater than a")
        return -6/5
    @staticmethod
    def mean(a,b):
        if(b<=a):
            raise ValueError("b must be greater than a")
        return 1/2*(a+b)
    @staticmethod
    def median(a,b):
        if(b<=a):
            raise ValueError("b must be greater than a")
        return 1/2*(a+b)
    @staticmethod
    def variance(a,b):
        if(b<=a):
            raise ValueError("b must be greater than a")
        return 1/12*((b-a)**2)
    @staticmethod
    def stddev(a,b):
        if(b<=a):
            raise ValueError("b must be greater than a")
        return 1/math.sqrt(12)*(b-a)
    @staticmethod
    def entropy(a,b):
        if(b<=a):
            raise ValueError("b must be greater than a")
        return math.log(b-a)
    @staticmethod
    def skewness(a,b):
        if(b<=a):
            raise ValueError("b must be greater than a")
        return 0
    @staticmethod
    def ppf(a,b,q):
        if(b<=a):
            raise ValueError("b must be greater than a")
        return ((b-a)*q+a)
    @staticmethod
    def mle(x):
        return {'a':min(x),'b':max(x)}