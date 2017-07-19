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

class geometric(Distribution):
    @staticmethod
    def pdf(p,k):
        if(p<=0 or p>1):
            raise ValueError("p must be greater than zero and less than or equal to 1")
        if(k%1!=0 or k<1):
            raise ValueError("k must be a positive integer")
        return math.pow(1-p,k-1)*p
    @staticmethod
    def cdf(p,k):
        if(p<=0 or p>1):
            raise ValueError("p must be greater than zero and less than or equal to 1")
        if(k%1!=0 or k<1):
            raise ValueError("k must be a positive integer")
        return 1-math.pow(1-p,k)
    @staticmethod
    def random(p):
        if(p<=0 or p>1):
            raise ValueError("p must be greater than zero and less than or equal to 1")
        return math.floor(math.log(ds.rg0())/math.log(1-p))
    @staticmethod
    def kurtosis(p):
        if(p<=0 or p>1):
            raise ValueError("p must be greater than zero and less than or equal to 1")
        return 6+(p**2)/(1-p)
    @staticmethod
    def mean(p):
        if(p<=0 or p>1):
            raise ValueError("p must be greater than zero and less than or equal to 1")
        return 1/p
    @staticmethod
    def median(p):
        if(p<=0 or p>1):
            raise ValueError("p must be greater than zero and less than or equal to 1")
        return math.ceil(-1/math.log(1-p,2))
    @staticmethod
    def mode(p):
        if(p<=0 or p>1):
            raise ValueError("p must be greater than zero and less than or equal to 1")
        return 1
    @staticmethod
    def variance(p):
        if(p<=0 or p>1):
            raise ValueError("p must be greater than zero and less than or equal to 1")
        return (1-p)/(p**2)
    @staticmethod
    def stddev(p):
        if(p<=0 or p>1):
            raise ValueError("p must be greater than zero and less than or equal to 1")
        return math.sqrt(1-p)/p
    @staticmethod
    def entropy(p):
        if(p<=0 or p>1):
            raise ValueError("p must be greater than zero and less than or equal to 1")
        return (-(1-p)*math.log(1-p,2)-p*math.log(p,2))/p
    @staticmethod
    def skewness(p):
        if(p<=0 or p>1):
            raise ValueError("p must be greater than zero and less than or equal to 1")
        return (2-p)/math.sqrt(1-p)
    @staticmethod
    def ppf(p,q):
        if(p<=0 or p>1):
            raise ValueError("p must be greater than zero and less than or equal to 1")
        return math.floor(math.log(q)/math.log(1-p))
    @staticmethod
    def mle(x):
        return {'p':len(x)/sum(x)}