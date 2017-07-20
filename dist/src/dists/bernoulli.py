'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r

class bernoulli(Distribution):
    @staticmethod
    def pdf(p,k):
        if((k!=0 and k!=1) or p<0 or p>1):
            raise ValueError("k must be 0 or 1, and p must be between 0 and 1 inclusive")
        if(k==0):
            return 1-p
        return p
    @staticmethod
    def cdf(p,k):
        if((k!=0 and k!=1) or p<0 or p>1):
            raise ValueError("k must be 0 or 1, and p must be between 0 and 1 inclusive")
        if(k==0):
            return 1-p
        return 1
    @staticmethod
    def random(p):
        if(p<0 or p>1):
            raise ValueError("p must be between 0 and 1 inclusive")
        if(r.random()>=p):
            return 1
        return 0
    @staticmethod
    def kurtosis(p):
        if(p<0 or p>1):
            raise ValueError("p must be between 0 and 1 inclusive")
        return (1-6*p*(1-p))/(p*(1-p))
    @staticmethod
    def mean(p):
        if(p<0 or p>1):
            raise ValueError("k must be 0 or 1, and p must be between 0 and 1 inclusive")
        return p
    @staticmethod
    def median(p):
        if(p<0 or p>1):
            raise ValueError("p must be between 0 and 1 inclusive")
        if(p<1/2):
            return 0
        if(p==1/2):
            return 0.5
        if(p==1):
            return 1
    @staticmethod
    def variance(p):
        if(p<0 or p>1):
            raise ValueError("p must be between 0 and 1 inclusive")
        return p*(1-p)
    @staticmethod
    def stddev(p):
        if(p<0 or p>1):
            raise ValueError("p must be between 0 and 1 inclusive")
        return math.sqrt(p*(1-p))
    @staticmethod
    def entropy(p):
        if(p<0 or p>1):
            raise ValueError("p must be between 0 and 1 inclusive")
        return 1/(p*(1-p))
    @staticmethod
    def skewness(p):
        if(p<0 or p>1):
            raise ValueError("p must be between 0 and 1 inclusive")
        return (1-2*p)/math.sqrt(p*(1-p))
    @staticmethod
    def ppf(p,q):
        if(p<0 or p>1):
            raise ValueError("p must be between 0 and 1 inclusive")
        if(q>=p):
            return 1
        return 0
    @staticmethod
    def mle(x):
        return {'p':sum(x)/len(x)}