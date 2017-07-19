'''
Created on Jul 15, 2017

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

class binomial(Distribution):
    @staticmethod
    def random(n,p):
        if(n%1!=0 or n<0 or p<0 or p>1):
            raise ValueError("n must be a non-negative integer and p must be between 0 and 1 inclusive")
        count=0
        for _ in range(n):
            if r.random()<p:
                count+=1
        return count
    @staticmethod
    def kurtosis(n,p):
        if(n%1!=0 or n<0 or p<0 or p>1):
            raise ValueError("n must be a non-negative integer and p must be between 0 and 1 inclusive")
        return (1-6*p*(1-p))/(n*p*(1-p))
    @staticmethod
    def mean(n,p):
        if(n%1!=0 or n<0 or p<0 or p>1):
            raise ValueError("n must be a non-negative integer and p must be between 0 and 1 inclusive")
        return n*p
    @staticmethod
    def median(n,p):
        if(n%1!=0 or n<0 or p<0 or p>1):
            raise ValueError("n must be a non-negative integer and p must be between 0 and 1 inclusive")
        return int(round(n*p))
    @staticmethod
    def variance(n,p):
        if(n%1!=0 or n<0 or p<0 or p>1):
            raise ValueError("n must be a non-negative integer and p must be between 0 and 1 inclusive")
        return n*p*(1-p)
    @staticmethod
    def stddev(n,p):
        if(n%1!=0 or n<0 or p<0 or p>1):
            raise ValueError("n must be a non-negative integer and p must be between 0 and 1 inclusive")
        return math.sqrt(n*p*(1-p))
    @staticmethod
    def skewness(n,p):
        if(n%1!=0 or n<0 or p<0 or p>1):
            raise ValueError("n must be a non-negative integer and p must be between 0 and 1 inclusive")
        return (1-2*p)/math.sqrt(n*p*(1-p))